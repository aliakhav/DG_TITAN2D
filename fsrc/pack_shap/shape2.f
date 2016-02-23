C file:  ./elem/shape2.f
C=========================
      subroutine shape2 (Idegr, X, Y, F )

       IMPLICIT REAL*8 (A-H,O-Z)
c
c  double precision declaration
c
c-----------------------------------------------------------------------
c    routine name   -  shape2
c
c-----------------------------------------------------------------------
c
c   computer            - ibm/pc
c
c   latest revision      - october 29, 1990
c
c   purpose            - routine computes values  of 2-d shape
c                    functions ( on master element )
c
c   usage            - call      shape2 (Idegr, X, Y, F )
c
c
c   arguments           in   Idegr - degrees of approximation for all
c                          nodal points
c                      X,Y - coordinates (between -1. and 1.)
c
c                out       F - values of shape functions
c
c
c   NOTE: routine modified by tim to account for unstrucured meshes
c         common block cunstruc must be filled for the element
c         by a call to infelm
c
c   required   routines - shape1
c
c-----------------------------------------------------------------------
      common /cunstruc/ NODNIC(13)
      dimension F(*),Idegr(*),fx(11),fy(11)
C
c....................................................................
c     values        of shape functions at x,y on master element
c     2-d case
c....................................................................
c
      maxdgx = max(Idegr(1),Idegr(3),Idegr(2),Idegr(4),Idegr(5))
      maxdgy = max(Idegr(2),Idegr(4),Idegr(1),Idegr(3),Idegr(5))
c
      call shape1(maxdgx, X, fx)
      call shape1(maxdgy, Y, fy)
c
      F(1) = fx(1)*fy(1)
      F(2) = fx(2)*fy(1)
      F(3) = fx(2)*fy(2)
      F(4) = fx(1)*fy(2)
      idg = 4
c
      idgst = idg
      do 10 i = 3,Idegr(1)+1
         idg = idg + 1
 10      F(idg) = fx(i)*fy(1)
c
coh      iskip = 0
coh      if (NODNIC(5) .lt. 0) iskip = 9
coh      if (NODNIC(1 + iskip) .gt. NODNIC(2 + iskip)) then
coh         do 20 i = 4,Idegr(5) + 1,2
coh            idgst = idgst + 2
coh 20         F(idgst) = -F(idgst)
coh      endif
c
      idgst = idg
      do 30 i = 3,Idegr(2)+1
         idg = idg + 1
 30      F(idg) = fx(2)*fy(i)
c
coh      iskip = 0
coh      if (NODNIC(6) .lt. 0) iskip = 9
coh      if (NODNIC(2 + iskip) .gt. NODNIC(3 + iskip)) then
coh         do 40 i = 4,Idegr(6) + 1,2
coh           idgst = idgst + 2
coh 40         F(idgst) = -F(idgst)
coh      endif
c
      idgst = idg
      do 50 i = 3,Idegr(3)+1
         idg = idg + 1
 50      F(idg) = fx(i)*fy(2)
c
coh      iskip = 0
coh      if (NODNIC(7) .lt. 0) iskip = 9
coh      if (NODNIC(4 + iskip) .gt. NODNIC(3 + iskip)) then
coh         do 60 i = 4,Idegr(7) + 1,2
coh            idgst = idgst + 2
coh 60         F(idgst) = -F(idgst)
coh      endif
c
      idgst = idg
      do 70 i = 3,Idegr(4)+1
         idg = idg + 1
 70      F(idg) = fx(1)*fy(i)
c
coh      iskip = 0
coh      if (NODNIC(8) .lt. 0) iskip = 9
coh      if (NODNIC(1 + iskip) .gt. NODNIC(4 + iskip)) then
coh         do 80 i = 4,Idegr(8) + 1,2
coh            idgst = idgst + 2
coh 80         F(idgst) = -F(idgst)
coh      endif
c
      do 90 j = 3,Idegr(5)+1
        do 90 i = 3,Idegr(5)+1
            idg = idg + 1
 90         F(idg) = fx(i)*fy(j)

      ndof=idg
c
      return
      end

C =bottom========= ./elem/shape2.f
