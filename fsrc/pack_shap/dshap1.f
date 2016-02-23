C file:  ./elem/dshap1.f
C=========================
      subroutine dshap1(Idg, X, F)
c---------------------------------------------------------------------
c    routine name   - dshap1
c
c---------------------------------------------------------------------
c
c   computer            - independent
c
c   latest revision     - 
c
c   purpose             - derivatives of shape functions
c                         
c                         
c
c
c   usage               - 
c
c
c   required   routines -
c
c---------------------------------------------------------------------
c   default header supplied by (TL) 10/14/91
c---------------------------------------------------------------------


       IMPLICIT REAL*8 (A-H,O-Z)


c
c  double precision declaration
c
c
      common /dlagr/ DLAG(9,8)
      dimension F(*)
c...................................................................
c     derivatives of shape functions at X on master element
c     1-d case
c     Idg max degree of polynomials
c...................................................................
      F(1) = -0.5
      F(2) = 0.5
c
      x2 = X*X
c
c  ...x**3, x**5, x**7
      do 20 nrf = 4,Idg+1,2
c
         F(nrf) = DLAG(nrf,nrf-1)*x2 + DLAG(nrf,nrf-3)
         do 10 ip = nrf-5,1,-2
            F(nrf) = F(nrf)*x2 + DLAG(nrf,ip)
 10         continue
c
 20      continue
c
c  ...x**2, x**4, x**6, x**8
c
      F(3) = DLAG(3,2)*X
c
      do 40 nrf = 5,Idg+1,2
c
         F(nrf) = DLAG(nrf,nrf-1)*x2 + DLAG(nrf,nrf-3)
         do 30 ip = nrf-5,2,-2
            F(nrf) = F(nrf)*x2 + DLAG(nrf,ip)
 30         continue
c
         F(nrf) = F(nrf)*X
c
 40      continue
c
      return
      end

C =bottom========= ./elem/dshap1.f
