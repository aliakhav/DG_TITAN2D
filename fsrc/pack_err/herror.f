c--------------------------------------------------------------------
c   routine name        - herror
c   writer              - Kristofer Max
c   latest revision     - modified AP Mar '99
c   purpose             - calculate the numerical error with the 
c                         Hierarchical error estimator
c---------------------------------------------------------------------
c
      subroutine herror(nequ,norder,xnod,sol,err,solnrm,matid)
c
c nequ == no. of equations/sh.fn.
c norder== polynom. ord. for elt.
c xnod == node coord
c sol == solution for elt.
c
      include '../com_fem/syscom.blk'
      include '../com_fem/cint.blk'
c
      dimension norder(5),morder(5),ek(242,242),ef(242),
     .     fb(242),sol(2,121),ae(242),ab(242),ff(242),abb(2,121)
c
      double precision kbb(242,242),keb(242,242),kebae(242)
c
      dimension vshape(121),xi(2),uf(2),er(2),erd(4),ufd(4)
      dimension gpsi(9),dgpsi(2,9),dxds(2,2),dsdx(2,2),xnod(2,9)
      dimension dpsiy(121),dpsix(121),dpsi(2,121),ip(242)
c


c


c         
c         write(6,*)'In herror, nel = ', nel
c
c     initialize arrays
c
         solnrm=0.
         err=0.d0
         do i=1,242
            kebae(i)=0.
            enddo
c increase local polynom order for err. calc.
         morder(1)=norder(1)
         morder(2)=norder(2)
         morder(3)=norder(3)
         morder(4)=norder(4)
         morder(5)=norder(5)+2
c
         nrdofr=norder(1)+norder(2)+norder(3)+norder(4)+(norder(5)-1)**2
         mrdofr=morder(1)+morder(2)+morder(3)+morder(4)+(morder(5)-1)**2
         nrdof=nrdofr*nequ
         mrdof=mrdofr*nequ
c
c     element stiffness matrix
c
c         call elem(nel,kind,2,morder,mrdof,nequ,242,  ek,ef)
         ifg=2
        
          call pelem(ifg,morder,mrdof,nequ,242,xnod,ek,ef,matid)
c
c     determine submatrices
c     
         nsize=mrdof-nrdof
         do  i=1,nsize
            fb(i)=ef(nrdof+i)
            do  j=1,nsize
               kbb(i,j)=ek(nrdof+i,nrdof+j)
               enddo
               enddo
         do  i=1,nsize
            do  j=1,nrdof
               keb(i,j)=ek(nrdof+i,j)
                enddo
               enddo
c
c     find solution coefficients
c
         j=1
         do 40 i=1,nrdof
            ae(j)=sol(1,i)
            j=j+1
            ae(j)=sol(2,i)
            j=j+1
 40      continue
c
c     set up equation 
c
         do 50 i=1,242
            do 50 j=1,242
               kebae(i)=kebae(i)+keb(i,j)*ae(j)
 50      continue
c
         do 60 i=1,242
            ff(i)=fb(i)-kebae(i)
 60      continue
c
c     solve for ab
c

         call gausse(kbb,242,ff,nsize,ab)
c         call gauss2(nsize,242,kbb,ip,ff,ab)
c
         i=1
         do 70 j=1,121
            abb(1,j)=ab(i)
            i=i+1
            abb(2,j)=ab(i)
            i=i+1
 70      continue
c
c     calculate error norm
c     
c     define the number of integration point
c
         maxdg1 = max0(morder(1),morder(3),morder(5))
         maxdg2 = max0(morder(2),morder(4),morder(5))
c
         if (((maxdg1.lt.1).or.(maxdg2.lt.1)).or.
     .        ((maxdg1.gt.9).or.(maxdg2.gt.9)))  then
            write(*,*) 'WRONG INPUT IN ELEM, MAXDG1,MAXDG2 = ',
     .           maxdg1,maxdg2
            stop
         endif
c
         nl1 = maxdg1+1
         nl2 = maxdg2+1
c
c     begin integration point loop
c
         do i=1,121
            dpsix(i) = 0.d0
            dpsiy(i) = 0.d0
         enddo
c     
         do 230 int1 = 1, nl1
            do 240 int2= 1, nl2
c
               xi(1) = XIGAUS(int1,nl1)
               xi(2) = XIGAUS(int2,nl2)
               call gshape(xi, gpsi,dgpsi)
c     
c     to calculate shape func
c
               call shape2(morder,xi(1),xi(2),vshape)
c
c     calculate fem  solution @ gauss points
c    
               uf(1) = 0.
               uf(2) = 0.
               do jj =1,nequ
                  do j =1,nrdofr
                     uf(jj) = uf(jj) + vshape(j)*sol(jj,j)
                  enddo
               end do
c                
c     calculate error at Gauss points
c     
               er(1) = 0.
               er(2) = 0.
               do j =1,nsize/nequ
                  er(1) = er(1) + vshape(j+nrdofr)*abb(1,j)
                  er(2) = er(2) + vshape(j+nrdofr)*abb(2,j)
               enddo
c
c     calculate Jacobi matrix dx/ds(i,j)                                   
c     
               do 250 i=1,2
                  do 250 j=1,2
                     dxds(i,j)=0.d0
                     do 250 k=1,9
                        dxds(i,j)=dxds(i,j)+dgpsi(j,k)*xnod(i,k)
 250           continue
c     
c     inverse jacobi matrix dx/ds(i,j) to get ds/dx(i,j)
c
               detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
               if (detj .le. 0.0) stop
c
c     accumulate integration point value of integrals
c     
               fac=detj*WAGAUS(int1,nl1)*WAGAUS(int2,nl2)
c 
               dsdx(1,1)=dxds(2,2)/detj
               dsdx(2,2)=dxds(1,1)/detj
               dsdx(1,2)=-dxds(1,2)/detj
               dsdx(2,1)=-dxds(2,1)/detj
c     
c     to calculate dpsi                                                   
c     
               call dshap2(morder,xi(1),xi(2), dpsi)
               
c     calculate d(psi)/dx equation                                        
c     
               do 270 i=1,mrdofr
                  dpsix(i)=dpsi(1,i)*dsdx(1,1)+dpsi(2,i)*dsdx(2,1)
                  dpsiy(i)=dpsi(1,i)*dsdx(1,2)+dpsi(2,i)*dsdx(2,2)
 270           continue
c               
               ufd(1) = 0.
               ufd(2) = 0.
               ufd(3) = 0.
               ufd(4) = 0.
               do 280 jj = 1,nequ
                  do 280 i =1,nrdofr
                     ufd(2*(jj-1)+1)=ufd(2*(jj-1)+1)+sol(jj,i)*dpsix(i)
                     ufd(2*(jj-1)+2)=ufd(2*(jj-1)+2)+sol(jj,i)*dpsiy(i)
 280           continue
c
               erd(1) = 0.
               erd(2) = 0.
               erd(3) = 0.
               erd(4) = 0.
               do 290 i =1,nsize/nequ
                  erd(1) = erd(1) + abb(1,i)*dpsix(i+nrdofr)
                  erd(2) = erd(2) + abb(2,i)*dpsiy(i+nrdofr)
                  erd(3) = erd(3) + abb(2,i)*dpsix(i+nrdofr)
                  erd(4) = erd(4) + abb(1,i)*dpsiy(i+nrdofr)
 290           continue
c     
c     h1err norms
c         
               err=err+(er(1)**2+er(2)**2+erd(1)**2+erd(2)**2+erd(3)**2
     &              +erd(4)**2)*fac

               solnrm = solnrm + (uf(1)**2+uf(2)**2+ufd(1)**2+ufd(2)**2
     &              + ufd(3)**2 + ufd(4)**2)*fac
c
 240        continue
 230     continue
c
c
      return
      end
c


