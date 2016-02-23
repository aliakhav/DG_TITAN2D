c------------------------------------------------------------------
c     routine name calcvol      
c     lastest revision: April, 1993
c     arguments:
c       in:     Norder - order of element (for DG all values are the same)     
c               nedofr - number of dof of element
c                   Nc - length of the column of ek      
c                 nequ - number of equations
c                 xnod - the nodal coordinates
c                  sol - the solution vector for the mass shape functions
c                              
c        out:      ave - the average of each component of the solution
c                        values for an element
c
c    comment: this routines calculates the total volume in the element
c------------------------------------------------------------------
c      
c      
      subroutine calcave(Norder,nedofr,nequ,Nc,
     1     xnod,sol, ave)

      include '../com_fem/syscom.blk'
      include '../com_fem/paramf.blk'
      include '../com_fem/cint.blk'
c
c      include '../com_fem/celemf.blk'
c


      dimension sol(Nc)
      dimension xi(2),xnod(2,9),xy(2)
      dimension gpsi(9),dgpsi(2,9),dxds(2,2)
      dimension psi(121), Norder(*),ave(3)
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
      small = 1.e-5
      one = 1.d0 
      area = 0.d0
      do i = 1, 3
         ave(i) = 0.d0
      enddo

c
c  ...find the order of polynormial in two directions
c
      if (Norder(1) .lt.  1 .or. Norder(1) .gt. 9)  then
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
      nl1 = norder(1)+1
c
c ...loop through integration points
      do 30 l1 = 1, nl1
         do 40 l2= 1, nl1
            xi(1) = XIGAUS(l1,nl1)
            xi(2) = XIGAUS(l2,nl1)
            call gshape(xi, gpsi,dgpsi)
            xy(1) = x(1,xi(1),xi(2))
            xy(2) = x(2,xi(1),xi(2))
c     
c     ...calculate Jacobi matrix dx/ds(i,j)
            do 50 i=1,2
               do 50 j=1,2
                  dxds(i,j)=0.d0
                  do 50 k=1,9
                     dxds(i,j)=dxds(i,j)+dgpsi(j,k)*xnod(i,k)
 50         continue
c     
c     ...inverse jacobi matrix dx/ds(i,j) to get ds/dx(i,j)
            detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
            if (detj.le.0.0) go to 99
c     
c
c      ...accumulate integration point value of integrals
            fac=detj*WAGAUS(l1,nl1)*WAGAUS(l2,nl1)          
c
c      ...evaluate values and derivatives of the shape functions
            call shape2(Norder,xi(1),xi(2), psi)

c     ...calculate the volume at this point
            do i=1,3
               sum2 = 0.d0
               do k=1, nedof
                  sum2 = sum2 + psi(k) * sol(k+(i-1)*nequ)
               enddo
               ave(i) = ave(i) + fac*sum2
            enddo
            area = area + fac
            
c     
c.......end of first loop through components        
 40      continue
 30   continue
c
c.....end of loop through gaussian points
      do i=1,3
         ave(i) = ave(i)/fac
      enddo

      return
c      
 99   write(*,100) detj,xi
 100  format(' in elem: bad jacobian = ',e10.3/2(2x,e10.3))
      do i=1,9
         write(*,*) xnod(1,i),xnod(2,i)
      enddo
c     
      stop
      end


      
