c---------------------------------------------------------------------
c   routine name        -  err
c
c   computer            - machine independent
c
c   latest revision     - July,1992
c
c   purpose             - compute and plot exact error
c
c   usage               - call err
c
c   arguments
c       in
c         nequ          - Number of equations
c         Norder(5)     - Order of edges and bubble node
c         Xnod(2,9)     - Coordinates of nodes
c         Utemp(2,*)    - Solution values: vertices, edges, and bubble
c       out
c         errorsq(2)    - Square of element's error
c         solsq(2)      - Square of element's solution
c
c   required   routines - exsol
c---------------------------------------------------------------------
c........0.........0.........0.........0.........0.........0.........0..
c
      subroutine elemerrp(nequ,Norder,Xnod,Utemp,errorsq,solsq)
c
      include '../com_fem/syscom.blk'
      include '../com_fem/cint.blk'
c
      dimension Xnod(2,9),Utemp(2,*),Norder(5)
c
C      dimension xy(2,5)
      dimension Vshape(121),xi(2),xya(2),uf(2),uex(2),uexd(4),ufd(4)
c
      dimension gpsi(9),dgpsi(2,9),dxds(2,2),dsdx(2,2)
      dimension dpsiy(121),dpsix(121),dpsi(2,121)
C     psi(121)
c
c------------------------------------------------------------------
c
c  define local functions defining geometry of the element
c
c...transformation of element geometry
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
c
      ntdof =  Norder(1)+Norder(2)+Norder(3)+Norder(4)
     .      +  (Norder(5)-1)*(Norder(5)-1)
c
c....find the order of polynormial in two directions
c
      maxdg1 = max0(Norder(1),Norder(3),Norder(5))
      maxdg2 = max0(Norder(2),Norder(4),Norder(5))
c
      if (((maxdg1.lt.1).or.(maxdg2.lt.1)).or.
     .  ((maxdg1.gt.9).or.(maxdg2.gt.9)))  then
        write(*,*) 'WRONG INPUT IN ELEM, MAXDG1,MAXDG2 = ',
     .              maxdg1,maxdg2
        stop
      endif

c....define the number of integration point
c
      nl1 = maxdg1+1
C       nl1 = 10
      nl2 = maxdg2+1
C       nl2 = 10
c....begin integration point loop
c
      solsq = 0.0d0
      errorsq = 0.0d0

c

c  debugging
      do i=1,5
         xi(1) = -1.4d0+(i*.4d0)
         do jjj=1,5
            xi(2) =  -1.4d0+(jjj*.4d0)
c  global position
            xya(1) = x(1,xi(1),xi(2))
            xya(2) = x(2,xi(1),xi(2))            

            call shape2(Norder,xi(1),xi(2),vshape)

            uf(1) = 0.
            uf(2) = 0.
            do jj =1,nequ
               do j =1,ntdof
                  uf(jj) = uf(jj) + vshape(j)*Utemp(jj,j)
               enddo
            end do

            call exsol(xya,uex,uexd)

            uexd(1) = dabs(uex(1) - uf(1))
            uexd(2) = dabs(uex(2) - uf(2))
         enddo
      enddo



c  end debugging
      do 230 l1 = 1, nl1
        do 240 l2= 1, nl2

          xi(1) = XIGAUS(l1,nl1)
          xi(2) = XIGAUS(l2,nl2)
          call gshape(xi, gpsi,dgpsi)
          xya(1) = x(1,xi(1),xi(2))
          xya(2) = x(2,xi(1),xi(2))
c
c........to calculate shape func

          call shape2(Norder,xi(1),xi(2),vshape)

c
c calculate fem  solution @ gauss points
c
          uf(1) = 0.
          uf(2) = 0.
           do jj =1,nequ
              do j =1,ntdof
                 uf(jj) = uf(jj) + vshape(j)*Utemp(jj,j)
             enddo
          end do
c
c compute exact solution @gauss pt

             call exsol(xya,uex,uexd)
c
c         calculate Jacobi matrix dx/ds(i,j)
c
          do 250 i=1,2
            do 250 j=1,2
              dxds(i,j)=0.d0
              do 250 k=1,9
                dxds(i,j)=dxds(i,j)+dgpsi(j,k)*xnod(i,k)
  250        continue
c
c........inverse jacobi matrix dx/ds(i,j) to get ds/dx(i,j)
c
          detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
          if (detj .le. 0.0) stop
c
c........accumulate integration point value of integrals
c

           fac=detj*WAGAUS(l1,nl1)*WAGAUS(l2,nl2)
c
            dsdx(1,1)=dxds(2,2)/detj
            dsdx(2,2)=dxds(1,1)/detj
            dsdx(1,2)=-dxds(1,2)/detj
            dsdx(2,1)=-dxds(2,1)/detj
c
c..........to calculate dpsi
c
            call dshap2(Norder,xi(1),xi(2), dpsi)

c..........calculate d(psi)/dx equation
c
            do 270 i=1,ntdof
              dpsix(i)=dpsi(1,i)*dsdx(1,1)+dpsi(2,i)*dsdx(2,1)
              dpsiy(i)=dpsi(1,i)*dsdx(1,2)+dpsi(2,i)*dsdx(2,2)
  270      continue

           ufd(1) = 0.
           ufd(2) = 0.
           ufd(3) = 0.
           ufd(4) = 0.

           do 280 jj = 1,nequ

           do 280 i =1,ntdof

         ufd(2*(jj-1)+1) = ufd(2*(jj-1)+1) + Utemp(jj,i)*dpsix(i)
         ufd(2*(jj-1)+2) = ufd(2*(jj-1)+2) + Utemp(jj,i)*dpsiy(i)

  280    continue

c
c h1err norms
c
         errorsq = errorsq + ((uf(1)-uex(1))**2
     &        + (uexd(1)-ufd(1))**2 + (uexd(2)-ufd(2))**2)*fac
     &        + ((uf(2)-uex(2))**2
     &        + (uexd(3)-ufd(3))**2 + (uexd(4)-ufd(4))**2)*fac
c l2 err norms
c
c         errorsq = errorsq + (uf(1)-uex(1))**2+(uf(2)-uex(2))**2
c
         solsq = solsq + (uf(1)**2 + ufd(1)**2 + ufd(2)**2)*fac
     &        + (uf(2)**2 + ufd(3)**2 + ufd(4)**2)*fac
c
  240   continue
  230 continue

 999  return
      end
c
c........0.........0.........0.........0.........0.........0.........0..

