C file:  ./pack_shap/dshap2.f
C=========================
      subroutine dshap2dg (Idegr, X, Y, Dgxy)
       IMPLICIT REAL*8 (A-H,O-Z)

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
c
c   required   routines - shape1,dshap1
c
c-----------------------------------------------------------------------
      dimension Dgxy(2,*)
c
c....................................................................
c     derivatives of shape functions at x,y on master element
c     2-d case
c     Idegr(i)- degree of polynomials for nodal points
c....................................................................
c
      Dgxy(1,1) = 0.d0
      Dgxy(2,1) = 0.d0
      if(idegr .gt. 2) goto 99
      do i=1,6
         do j=1,2
            Dgxy(j,i) = 0.d0
         enddo
      enddo
      if(idegr .ge. 1) then
         Dgxy(1,2) = 1.d0
         Dgxy(2,2) = 0.d0
         Dgxy(1,3) = 0.d0
         Dgxy(2,3) = 1.d0
         Dgxy(1,4) = y
         Dgxy(2,4) = x
      endif
      if(idegr .ge. 2) then
         Dgxy(1,5) = 2.d0*x
         Dgxy(2,5) = 0.d0
         Dgxy(1,6) = 0.d0
         Dgxy(2,6) = 2.d0*y
         Dgxy(1,7) = -2.d0*x*y
         Dgxy(2,7) = 1-x*x
         Dgxy(1,8) = 1-y*y
         Dgxy(2,8) = -2.d0*x*y
         Dgxy(1,9) = -2.d0*x*(1-y*y)
         Dgxy(2,9) = -2.d0*y*(1-x*x)
      endif

      return
99    write (*,100)i
100   format(/' maximum of 2nd order ',i5)
      end

C =bottom========= ./elem/dshap2.f
