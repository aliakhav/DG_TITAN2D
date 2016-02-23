C file:  ./pack_shape/shape2dg.f
C=========================
      subroutine shape2dg (Idegr, x, y, F )

       IMPLICIT double precision (a-h,o-z)
       dimension F(*)
c
c  double precision declaration
c
c-----------------------------------------------------------------------
c    routine name   -  shape2dg
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
c
c-----------------------------------------------------------------------
C
c....................................................................
c     values        of shape functions at x,y on master element
c     2-d case
c....................................................................
c
      if(idegr .gt. 2) goto 99

      F(1) = 1.d0

      if(idegr .ge. 1) then
         F(2) = x
         F(3) = y
         F(4) = x*y
      endif

      if(idegr .ge. 2) then
         F(5) = x*x-1.0/3.0 
         F(6) = y*y-1.0/3.0
         F(7) = (1-x*x)*y
         F(8) = (1-y*y)*x
         f(9) = (1-x*x)*(1-y*y)-4.0/9.0
      endif

      return
99    write (*,100)i
100   format(/' maximum of 2nd order ',i5)
      end

C =bottom========= ./elem/shape2.f
