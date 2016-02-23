C file:  ./elem/dshap2.f
C=========================
      subroutine dshap2 (Idegr, X, Y, Dgxy)

       IMPLICIT REAL*8 (A-H,O-Z)

c
c  double precision declaration
c
c-----------------------------------------------------------------------
c    routine name   -  dshap2
c
c-----------------------------------------------------------------------
c
c   computer            - ibm/pc
c
c   latest revision      - october 29, 1990
c
c   purpose            - routine computes derivatives of 2-d shape
c                    functions ( on master element )
c
c   usage            - call dshap2 (Idegr, X, Y, Dgx, Dgy )
c
c
c   arguments           in Idegr - degree of approximation
c                    X,Y - coordinates (between -1. and 1.)
c
c              out Dgx,Dgy - values of derivatives of shape
c                        functions
c
c   NOTE:  routine modified by tim to allow for unstructured
c          meshes.  common block must be filled for element
c          by a call to infelm
c
c   required   routines - shape1,dshap1
c
c-----------------------------------------------------------------------
      common /cunstruc/ NODNIC(13)
      dimension Dgxy(2,*),Idegr(*),dfx(11),dfy(11),
     *                                  fx(11), fy(11)
c
c....................................................................
c     derivatives of shape functions at x,y on master element
c     2-d case
c     Idegr(i)- degree of polynomials for nodal points
c....................................................................
c
      maxdgx = max(Idegr(1),Idegr(3),Idegr(5))
      maxdgy = max(Idegr(2),Idegr(4),Idegr(5))     
      call dshap1(maxdgx, X, dfx)
      call dshap1(maxdgy, Y, dfy)
      call shape1(maxdgx, X, fx)
      call shape1(maxdgy, Y, fy)

c      write(*,*) fx(3),fy(3),dfx(3),dfy(3)
c
      Dgxy(1,1) = dfx(1)*fy(1)
      Dgxy(1,2) = dfx(2)*fy(1)
      Dgxy(1,3) = dfx(2)*fy(2)
      Dgxy(1,4) = dfx(1)*fy(2)
      Dgxy(2,1) = fx(1)*dfy(1)
      Dgxy(2,2) = fx(2)*dfy(1)
      Dgxy(2,3) = fx(2)*dfy(2)
      Dgxy(2,4) = fx(1)*dfy(2)
c
      idg = 4
c
      idgst = idg
      do 10 i = 3,Idegr(1)+1
         idg = idg + 1
         Dgxy(1,idg) = dfx(i)* fy(1)
 10      Dgxy(2,idg) =  fx(i)*dfy(1)
c
coh      iskip = 0
coh      if (NODNIC(5) .lt. 0) iskip = 9
coh      if (NODNIC(1 + iskip) .gt. NODNIC(2 + iskip)) then
coh         do 20 i = 4,Idegr(5) + 1,2
coh            idgst = idgst + 2
coh            Dgx(idgst) = -Dgx(idgst)
coh 20         Dgy(idgst) = -Dgy(idgst)
coh      endif
c
      idgst = idg
      do 30 i = 3,Idegr(2)+1
         idg = idg + 1
         Dgxy(1,idg) = dfx(2)* fy(i)
 30      Dgxy(2,idg) = fx(2)*dfy(i)
c
coh      iskip = 0
coh      if (NODNIC(6) .lt. 0) iskip = 9
coh      if (NODNIC(2 + iskip) .gt. NODNIC(3 + iskip)) then
coh         do 40 i = 4,Idegr(6) + 1,2
coh            idgst = idgst + 2
coh            Dgx(idgst) = -Dgx(idgst)
coh 40         Dgy(idgst) = -Dgy(idgst)
coh      endif
c
      idgst = idg
      do 50 i = 3,Idegr(3)+1
         idg = idg + 1
         Dgxy(1,idg) = dfx(i)* fy(2)
 50      Dgxy(2,idg) = fx(i)*dfy(2)
c
coh      iskip = 0
coh      if (NODNIC(7) .lt. 0) iskip = 9
coh      if (NODNIC(4 + iskip) .gt. NODNIC(3 + iskip)) then
coh         do 60 i = 4,Idegr(7) + 1,2
coh            idgst = idgst + 2
coh            Dgx(idgst) = -Dgx(idgst)
coh 60         Dgy(idgst) = -Dgy(idgst)
coh      endif
c
      idgst = idg
      do 70 i = 3,Idegr(4)+1
         idg = idg + 1
         Dgxy(1,idg) = dfx(1)* fy(i)
 70      Dgxy(2,idg) = fx(1)*dfy(i)
c
coh      iskip = 0
coh      if (NODNIC(8) .lt. 0) iskip = 9
coh      if (NODNIC(1 + iskip) .gt. NODNIC(4 + iskip)) then
coh         do 80 i = 4,Idegr(8) + 1,2
coh            idgst = idgst + 2
coh            Dgx(idgst) = -Dgx(idgst)
coh 80         Dgy(idgst) = -Dgy(idgst)
coh      endif
c
      do 90 j = 3,Idegr(5)+1
         do 90 i = 3,Idegr(5)+1
            idg = idg + 1
            Dgxy(1,idg) = dfx(i)* fy(j)
 90         Dgxy(2,idg) = fx(i)*dfy(j)
c
      ndof = idg
c
      return
      end

C =bottom========= ./elem/dshap2.f
