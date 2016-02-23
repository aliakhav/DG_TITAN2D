c------------------------------------------------------------------
c     routine name melem      
c     lastest revision: April, 1993
c     arguments:
c       in:    nel,kind - number and kind of the element
c                  ifg - 1 not do element stffness matrix
c                      - otherwise: do element stiffness matrix
c               Norder - order of element      
c               nedofr - dof number of element
c                   Nc - length of the column of ek      
c                 nequ -       
c        out:       ek - the element stiffness matrix
c            ef - element load vector      
c    comment: this routine establish the element stiffness matrix
c             using Hierarchical shape functions. The element geometry
c             is described by the normal quadratic shape functions.
c------------------------------------------------------------------
c      
      subroutine pelem(ifg,Norder,nedofr,nequ,Nc,xnod,ek,ef,matid)
c      
      include '../com_fem/syscom.blk'
      include '../com_fem/paramf.blk'
      include '../com_fem/cint.blk'
c
c      include '../com_fem/celemf.blk'
c


      dimension ek(Nc,*),ef(*)
      dimension xi(2),xnod(2,9),Norder(*),xy(2)
      dimension gpsi(9),dgpsi(2,9),dxds(2,2),dsdx(2,2)
      dimension dpsiy(121),dpsix(121),psi(121),dpsi(2,121)
c
      dimension a11(MAXEQS,MAXEQS),a12(MAXEQS,MAXEQS),
     .          a21(MAXEQS,MAXEQS),a22(MAXEQS,MAXEQS),
     .           b1(MAXEQS,MAXEQS), b2(MAXEQS,MAXEQS),
     .           d1(MAXEQS,MAXEQS), d2(MAXEQS,MAXEQS), 
     .            c(MAXEQS,MAXEQS),
     .            f(MAXEQS),fx(MAXEQS),fy(MAXEQS)
c      
c  ...local geometry function
      x(i,eta1,eta2) = Xnod(i,1)*(eta1**2-eta1)*(eta2**2-eta2)/4.d0
     .               + Xnod(i,2)*(eta1**2+eta1)*(eta2**2-eta2)/4.d0
     .               + Xnod(i,3)*(eta1**2+eta1)*(eta2**2+eta2)/4.d0
     .               + Xnod(i,4)*(eta1**2-eta1)*(eta2**2+eta2)/4.d0
     .               + Xnod(i,5)*(1.d0-eta1**2)*(eta2**2-eta2)/2.d0
     .               + Xnod(i,6)*(eta1**2+eta1)*(1.d0-eta2**2)/2.d0
     .               + Xnod(i,7)*(1.d0-eta1**2)*(eta2**2+eta2)/2.d0
     .               + Xnod(i,8)*(eta1**2-eta1)*(1.d0-eta2**2)/2.d0
     .               + Xnod(i,9)*(1.d0-eta1**2)*(1.d0-eta2**2)
c
      zero = 0.d0
      small = 1.e-20
      one = 1.d0      
c
c  ...initialize element arrays
      if (ifg.ne.1) then
        do 10 i=1,nedofr
          call setz(ek(1,i),nedofr)
   10   continue
      endif
      call setz(ef,nedofr)
c
c  ...find the order of polynormial in two directions
      maxdg1 = max0(Norder(1),Norder(3),Norder(5))
      maxdg2 = max0(Norder(2),Norder(4),Norder(5))
c
      if (((maxdg1.lt.1).or.(maxdg2.lt.1)).or.
     .  ((maxdg1.gt.9).or.(maxdg2.gt.9)))  then
        write(*,*) 'WRONG INPUT IN ELEM, MAXDG1,MAXDG2 = ',
     .              maxdg1,maxdg2
        stop
      endif
c
      nedof = nedofr/nequ
c
c  ...determine the element geometry
c
c....define the number of integration point
      nl1 = maxdg1+1
      nl2 = maxdg2+1
c
c ...initiate dpsix and dpsiy
      do i=1,121
        dpsix(i) = 0.d0
        dpsiy(i) = 0.d0
      enddo
c

c ...loop through integration points
      do 30 l1 = 1, nl1
        do 40 l2= 1, nl2
          xi(1) = XIGAUS(l1,nl1)
          xi(2) = XIGAUS(l2,nl2)
          call gshape(xi, gpsi,dgpsi)
          xy(1) = x(1,xi(1),xi(2))
          xy(2) = x(2,xi(1),xi(2))
c
c      ...calculate Jacobi matrix dx/ds(i,j)
          do 50 i=1,2
            do 50 j=1,2
              dxds(i,j)=0.d0
              do 50 k=1,9
                dxds(i,j)=dxds(i,j)+dgpsi(j,k)*xnod(i,k)
   50     continue
c          
c      ...inverse jacobi matrix dx/ds(i,j) to get ds/dx(i,j)
          detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
          if (detj.le.0.0) go to 99
c
          dsdx(1,1)=dxds(2,2)/detj
          dsdx(2,2)=dxds(1,1)/detj
          dsdx(1,2)=-dxds(1,2)/detj
          dsdx(2,1)=-dxds(2,1)/detj          
c
c      ...accumulate integration point value of integrals
          fac=detj*WAGAUS(l1,nl1)*WAGAUS(l2,nl2)          
c
c      ...evaluate values and derivatives of the shape functions
          call shape2(Norder,xi(1),xi(2), psi)
          call dshap2(Norder,xi(1),xi(2), dpsi)
          do i=1,nedof
            dpsix(i)=dpsi(1,i)*dsdx(1,1)+dpsi(2,i)*dsdx(2,1)
            dpsiy(i)=dpsi(1,i)*dsdx(1,2)+dpsi(2,i)*dsdx(2,2)
          enddo          
c
c      ...call for coefficients defining bilinear and linear form
          call mcoeff(nequ,MAXEQS,
     .      xy(1),xy(2),a11,a12,a21,a22,b1,b2,d1,d2,c,f,fx,fy,matid)
c
c      ...first loop through components
          do 60 i=1,nequ
c
c        ...if this is a resolution skip the stiffness calculations
            if (ifg.eq.1) goto 80            
            do 70 j=1,nequ
c
c          ...coefficient a11(i,j)
              if (a11(i,j).ne.zero) then
                temp = a11(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * dpsix(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * dpsix(l)
                  enddo
                enddo
              endif
c
c          ...coefficient a12(i,j)
              if (a12(i,j).ne.zero) then
                temp = a12(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * dpsix(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * dpsiy(l)
                  enddo
                enddo
              endif
c
c          ...coefficient a21(i,j)
              if (a21(i,j).ne.zero) then
                temp = a21(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * dpsiy(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * dpsix(l)
                  enddo
                enddo
              endif
c
c          ...coefficient a22(i,j)
              if (a22(i,j).ne.zero) then
                temp = a22(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * dpsiy(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * dpsiy(l)
                  enddo
                enddo
              endif
c
c          ...coefficient d1(i,j)
              if (d1(i,j).ne.zero) then
                temp = d1(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * dpsix(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * psi(l)
                  enddo
                enddo
              endif              
c
c          ...coefficient d2(i,j)
              if (d2(i,j).ne.zero) then
                temp = d2(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * dpsiy(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * psi(l)
                  enddo
                enddo
              endif
c
c          ...coefficient b1(i,j)
              if (b1(i,j).ne.zero) then
                temp = b1(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * psi(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * dpsix(l)
                  enddo
                enddo
              endif
c
c          ...coefficient b2(i,j)
              if (b2(i,j).ne.zero) then
                temp = b2(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * psi(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * dpsiy(l)
                  enddo
                enddo
              endif
c
c          ...coefficient c(i,j)
              if (c(i,j).ne.zero) then
                temp = c(i,j) * fac
                do k = 1, Nedof
                  temp1 = temp * psi(k)
                  do l = 1, Nedof
                    ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .              ek(nequ*(k-1)+i,nequ*(l-1)+j) + temp1 * psi(l)
                  enddo
                enddo
              endif                
   70       continue
c
c        ...skip to here if doing a resolution           
   80       continue
c
c        ...accumulate for the load vector...
c        ...coefficient f(i) 
            if (f(i).ne.zero) then
              temp = f(i) * fac
              do k = 1, nedof
                ef(nequ*(k-1)+i) = ef(nequ*(k-1)+i) + temp*psi(k)
              enddo
            endif
c
c        ...coefficient f(i) 
            if (fx(i).ne.zero) then
              temp = fx(i) * fac
              do k = 1, nedof
                ef(nequ*(k-1)+i) = ef(nequ*(k-1)+i) + temp*dpsix(k)
              enddo
            endif
c            
c        ...coefficient f(i) 
            if (fy(i).ne.zero) then
              temp = fy(i) * fac
              do k = 1, nedof
                ef(nequ*(k-1)+i) = ef(nequ*(k-1)+i) + temp*dpsiy(k)
              enddo
            endif
   60   continue
c
c.......end of first loop through components        
   40   continue
   30 continue
c
c.....end of loop through gaussian points

      return
c      
99    write(*,100) detj,xi
100   format(' in elem: bad jacobian = ',e10.3/2(2x,e10.3))
      do i=1,9
        write(*,*) xnod(1,i),xnod(2,i)
      enddo
c      
      stop
      end
c
