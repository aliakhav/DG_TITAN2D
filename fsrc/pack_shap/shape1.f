C file:  ./elem/shape1.f
C=========================
      subroutine shape1(Idg, X, F )
c-----------------------------------------------------------------------
c    routine name   - 
c
c-----------------------------------------------------------------------
c
c   computer            - machine independent
c
c   latest revision     - 
c
c   purpose             - shape functions of Legendre type
c
c   arguments in:  
c
c   arguments out: 
c
c   required   routines -
c
c-----------------------------------------------------------------------
c - default header supplied by TL 10/11/91
c-----------------------------------------------------------------------
       IMPLICIT REAL*8 (A-H,O-Z)
c
c  double precision declaration
c
c
      dimension F(*)
      common /clagr/ CLAG(9,9)
c
      F(2) = (1 + X)*0.5
      F(1) = (1 - X)*0.5
c
      x2 = X*X
c
c  ...x**3, x**5, x**7
      do 20 nrf = 4,Idg+1,2
         F(nrf) = CLAG(nrf,nrf)*x2 + CLAG(nrf,nrf-2)
         do 10 ip = nrf-4,2,-2
            F(nrf) = F(nrf)*x2 + CLAG(nrf,ip)
 10         continue
c
         F(nrf) = F(nrf)*X
 20      continue
c
c  ...x**2, x**4, x**6, x**8
      do 40 nrf = 3,Idg+1,2
         F(nrf) = CLAG(nrf,nrf)*x2 + CLAG(nrf,nrf-2)
         do 30 ip = nrf-4,1,-2
            F(nrf) = F(nrf)*x2 + CLAG(nrf,ip)
 30         continue
 40      continue
c
      return
      end

C =bottom========= ./elem/shape1.f
