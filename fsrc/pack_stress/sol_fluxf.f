c------------------------------------------------------------------
c         routine name: sol_fluxf.f
c     lastest revision: April, 1993
c              purpose: give the derivative interpolation by least
c                       square method
c     arguments:
c                   in: nel,kind - number and kind of the element
c                       xi1,xi2 - two master coordinates
c                  out: du(1,*),du(2,*)
c------------------------------------------------------------------
c      
      subroutine sol_fluxf(lcol,nequ, np1,du,xnod,sol,norder)
c      
      include '../com_fem/syscom.blk'
c      include '../com_fem/paramf.blk'
      include '../com_fem/cint.blk'      
c
      dimension du(2,lcol,*)
      dimension Xnod(2,9),xi(2),gpsi(9),dgpsi(2,9),
     .          dxds(2,2),dsdx(2,2),psi(121),dpsi(2,121)
      dimension Norder(5),Nporder(5),sol(2,121)
c
      dimension A(100,100),B1(100),B2(100),C(100)
c
      zero = 0.d0
      small = 1.e-20
      one = 1.d0      
c
c  ...order of the element
c      call findap(nel,kind, Norder)
      nedof =  Norder(1)+Norder(2)+Norder(3)+Norder(4)
     .      + (Norder(5)-1)**2
c
c  ...number of superconvergent points in one direction      
      np1 = min0(Norder(1),Norder(2),Norder(3),Norder(4))
c
c  ...give
      np2 = np1*np1
      Nporder(1) = np1-1
      Nporder(2) = np1-1
      Nporder(3) = np1-1
      Nporder(4) = np1-1
      Nporder(5) = np1-1
c
c  ...initiating A     
      do 10  i=1,np2
        do 10 j=1,np2
   10     A(i,j) = 0.d0
        
      do 20 i1=1,np1
      do 20 i2=1,np1      
        xi(1) = xigaus(i1,np1)
        xi(2) = xigaus(i2,np1)
        call shape2(Nporder,xi(1),xi(2), psi)
c
        do 30 i=1,np2
          do 30 j=1,np2
   30       A(i,j) = A(i,j) + psi(i)*psi(j)
   20 continue
c
      call tri(A,100,np2)
c
c  ...determine the coordinates of element nodes
c      call nodcor(nel,kind, xnod)
c
c  ...element solution 
c      call solelm(nel,kind, sol)
c
c  ...loop through superconvergent point
      do 100 i1=1,nequ
        do i=1,np2
          B1(i) = 0.d0
          B2(i) = 0.d0
        enddo
        do 110 l1 = 1, np1
        do 110 l2 = 1, np1
c
c      ...coordinates of superconvergent points      
          xi(1) = xigaus(l1,np1)
          xi(2) = xigaus(l2,np1)
          call gshape(xi, gpsi,dgpsi)          
c
c      ...calculate Jacobi matrix dx/ds(i,j)
          do 120 i=1,2
            do 120 j=1,2
              dxds(i,j)=0.d0
              do 120 k=1,9
                dxds(i,j)=dxds(i,j)+dgpsi(j,k)*xnod(i,k)
  120     continue
c          
c      ...inverse jacobi matrix dx/ds(i,j) to get ds/dx(i,j)
          detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
          if (detj.le.0.0) then
            write(*,*) 'In sol_fluxf: detj = ',detj
            stop
          endif
c
          dsdx(1,1) =  dxds(2,2)/detj
          dsdx(2,2) =  dxds(1,1)/detj
          dsdx(1,2) = -dxds(1,2)/detj
          dsdx(2,1) = -dxds(2,1)/detj
c
c      ...derivatives of shape functions        
          call dshap2(Norder,xi(1),xi(2), dpsi)
          call shape2(Nporder,xi(1),xi(2), psi)
c
          dudxi1 = 0.d0
          dudxi2 = 0.d0        
          do j=1, nedof
            dudxi1 = dudxi1 + sol(i1,j)*dpsi(1,j)
            dudxi2 = dudxi2 + sol(i1,j)*dpsi(2,j)
          enddo
          do i=1,np2
            B1(i) = B1(i)+psi(i)
     .            * (dudxi1*dsdx(1,1)+dudxi2*dsdx(2,1))
            B2(i) = B2(i)+psi(i)
     .            * (dudxi1*dsdx(1,2)+dudxi2*dsdx(2,2))
          enddo
  110   continue
c
        call rhsub(A,B1,100,np2, C)                
        do j=1,np2
          du(1,i1,j) = C(j)
        enddo
c        
        call rhsub(A,B2,100,np2, C)
        do j=1,np2
          du(2,i1,j) = C(j)
        enddo        
c        
  100 continue
c
      return
      end
c


