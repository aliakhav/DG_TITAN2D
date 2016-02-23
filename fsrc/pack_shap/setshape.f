C file:  ./elem/setshape.f
C=========================
      SUBROUTINE setshape
c-----------------------------------------------------------------------
c    routine name   - 
c
c-----------------------------------------------------------------------
c
c   computer            - machine independent
c
c   latest revision     - 
c
c   purpose             - define coefficients for legendre shape functions
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


      COMMON /CRRR  / RRR(10,10,2)

CWR+
ctim+
C     MATRIX RRR CONTAINS THE constraint COEFFICIENTS for ensuring
c     continuty of the f.e. approximation  across elements of
c     differing size (different h-refinement levels).  The shape
c     functions from the smaller element are multiplied by these
c     coefficients.
c
c              -------------------------------
c                                |
c                                |  Small
c                                |  Element 2
c                  BIG           |
c                 ELEMENT        |------------
c                                |
c                                |  Small
c                                |  Element 1
c                                |
c              -------------------------------
c
c
c     RRR(I,J,K) : Coefficient for Small Element K's
c                  shape function J in relation to the
c                  BIG ELEMENT's shape function I.
c
c     Note:   The constraint process is done only on the
c             constrained side and is a 1-D process and so
c             it only involves the shape functions on the side.
c             Therefore the shape functions are the two appropriate
c             corner node shape functions plus the shape functions
c             associated with the midside node (if it exists).
c
C     NOTE:   THIS MATRIX DEPENDS SOLELY UPON THE CHOICE OF 1-D
c             HIERARCHICAL SHAPE FUNCTIONS.
c
c     Note:   This array is filled by a call to routine CNST.
C

      common /clagr/ CLAG(9,9)
      common /dlagr/ DLAG(9,8)
      common /elagr/ ELAG(9,7)
      common /flagr/ FLAG(9,6)

      dimension aul(9,8),a(9,8),b(9,8),ip(8),rhs(8)

      dimension FACT(10)
      dimension RR(10,10,2) ,cappsi(20),psi(20)


      call setz(clag,81)

c.....set coefficients defining shape functions:
c     meaning of clag( , ):
c     clag( fun#,degree+1 ) i.e.:
c                 shape(fun#,x) =  clag(fun#,1)*x     +
c                                  clag(fun#,2)*x**2  +
c                                  clag(fun#,3)*x**3  +
c                                  .................
c                                  clag(fun#,9)*x**8
c------------------------------------------------------


c----------------------------------------------------
C      WRITE(*,*)'***********************************************'
C      WRITE(*,*)'***                                         ***'
C      write(*,*)'***   NEW  SHAPE  FUNCTIONS  ARE  SET:      ***'
C      WRITE(*,*)'***        LEGENDRE POLYNOMIALS             ***'
C      WRITE(*,*)'***                                         ***'
C      WRITE(*,*)'***********************************************'
c.....L E G E N D R E  P O L Y N O M I A L S:
c      clag(3,1)=  -1.5 d0
c      clag(3,3)=   1.5 d0
c
c      clag(4,2)=  -2.5 d0
c      clag(4,4)=   2.5 d0
c
c      clag(5,1)=  -0.625 d0
c      clag(5,3)=  -3.75  d0
c      clag(5,5)=   4.375 d0
c
c      clag(6,2)=   0.875 d0
c      clag(6,4)=  -8.75  d0
c      clag(6,6)=   7.875 d0
c
c      clag(7,1)=  -1.3125  d0
c      clag(7,3)=   6.5625  d0
c      clag(7,5)= -19.6875  d0
c      clag(7,7)=  14.4375  d0
c
c      clag(8,2)=  -3.1875  d0
c      clag(8,4)=  19.6875  d0
c      clag(8,6)= -43.3125  d0
c      clag(8,8)=  26.8125  d0
c
c      clag(9,1)=  -0.7265625 d0
c      clag(9,3)=  -9.84375   d0
c      clag(9,5)=  54.140625  d0
c      clag(9,7)= -93.84375   d0
c      clag(9,9)=  50.2734375 d0
c------------------------------------------------------------------
c.....I N T E G R A T E D   L E G E N D R E   P O L Y N O M I A L S

c     WRITE(*,*)'***********************************************'
c     WRITE(*,*)'***                                         ***'
c     write(*,*)'***   NEW  SHAPE  FUNCTIONS  ARE  SET:      ***'
c     WRITE(*,*)'***   PRIMITIVES OF LEGENDRE POLYNOMIALS    ***'
c     WRITE(*,*)'***                                         ***'
c     WRITE(*,*)'***********************************************'

      clag(3,1)=   -0.5  d0
      clag(3,3)=    0.5  d0

      clag(4,2)=   -0.5  d0
      clag(4,4)=    0.5  d0

      clag(5,1)=    0.125 d0
      clag(5,3)=   -0.75  d0
      clag(5,5)=    0.625 d0

      clag(6,2)=    0.375  d0
      clag(6,4)=   -1.25   d0
      clag(6,6)=    0.8750 d0

      clag(7,1)=   -0.0625 d0
      clag(7,3)=    0.9375 d0
      clag(7,5)=   -2.1875 d0
      clag(7,7)=    1.3125 d0

      clag(8,2)=   -0.3125  d0
      clag(8,4)=    2.1875  d0
      clag(8,6)=   -3.9375  d0
      clag(8,8)=    2.0625  d0

      clag(9,1)=    0.0390625 d0
      clag(9,3)=   -1.09375   d0
      clag(9,5)=    4.921875  d0
      clag(9,7)=   -7.21875   d0
      clag(9,9)=    3.3515625 d0
c---------------------------------------------------------------


c  ...ceefficients for derivatives:
      call setz(DLAG,72)
      do 10 nrf = 3,9
         do 10 kp = 1,nrf-1
            DLAG(nrf,kp) = CLAG(nrf,kp+1) * kp
 10         continue
c
c  ...coefficients for second derivatives:
      call setz(ELAG,63)
      do 15 nrf = 3,9
         do 15 kp = 1,nrf-2
            ELAG(nrf,kp) = DLAG(nrf,kp+1) * kp
 15   continue
c
c  ---coefficients for third derivatives
      call setz(FLAG,54)
      do 17 nrf = 4,9
         do 17 kp = 1,nrf-3
            FLAG(nrf,kp) = ELAG(nrf,kp+1)*kp
 17         continue
c
c****************************************************************
c****************************************************************
c.....establishing matrix r:

c.....matrix a:
      call setz(a  ,72)
      call setz(aul,72)
      do 20 i=2,8
            do 20 j=2,8
                  a  (i,j) = clag(i+1,j+1)
                  aul(i,j) = clag(i+1,j+1)
 20               continue
      a  (1,1)=1.d0
      aul(1,1)=1.d0

c  ...b = invers of a
      N = 8
      NA= 9
      CALL DECOMP(N,NA,aul,IP,IFLAG)
      IF(IFLAG.EQ.1) WRITE(*,*)'FLAG=',IFLAG
      DO 40 I=1,N
         CALL SETZ(RHS,8)
         RHS(I)=1.D0
         CALL GAUSS2(N,NA,aul,IP,RHS,B(1,I) )
 40   CONTINUE


c.....computing rrr itself:

      FACT(1)=1
      do 60 i=2,10
            FACT(i)=FACT(i-1)/i
 60         continue

      call setz(rr ,200)
      call setz(rrr,200)

c.....old version of rrr:
      na=10
      rr(1,1,1)=.5
      rr(1,1,2)=.5
      do 70 i=2,na,2
          rr(i,1,1)= -1.
 70          rr(i,1,2)= -1.
      xmlt= .5
      do 80 i=2,na
          xmlt=xmlt*      .5
          rr(i,i,1)=xmlt
          rr(i,i,2)=xmlt
          one=1
          do 80 j=i-1,2,-1
              one=-one
              rr(i,j,2)=xmlt    *FACT(i-j)
 80              rr(i,j,1)=xmlt*one*FACT(i-j)
      do 90 i=2,na
          do 90 k=2,i
              rr(i,k,1)=rr(i,k,1)*FACT(k)/FACT(i)
 90              rr(i,k,2)=rr(i,k,2)*FACT(k)/FACT(i)


c.....new version expressed by the old one:
      do 100 i=1,8
            do 101 j=1,8
                  do 102 k=1,8
                        do 103 l=1,8
                              rrr(i,j,1) = rrr(i,j,1) +
     .                                     a(i,k)*rr(k,l,1)*b(l,j)
                              rrr(i,j,2) = rrr(i,j,2) +
     .                                     a(i,k)*rr(k,l,2)*b(l,j)
 103                          continue
 102                    continue
 101              continue
 100        continue
c4999 write(*,*) 'enter n,x'
c      read(*,*) n,x
c      call shape1(n,x,psi)
c      do 5000 k=1,2
c        capxi=-0.5*(1.0-x)
c        if (k.eq.2) capxi=0.5*(1.0+x)
c write(*,*)'------------------------------------'
c 	 write(*,*) 'k,x,capxi',k,x,capxi
c call shape1(n,capxi,cappsi)
c sum=0.
c do 5001 j=0,n
c           sum=sum+rrr(n,j+1,k)*psi(j+1)
c     write(*,*)'rrr,psi',rrr(n,j+1,k),psi(j+1),rrr(n,j+1,k)*psi(j+1)
c5001    continue
c        write(*,*)'k,n,err',k,n,sum,cappsi(n+1),sum-cappsi(n+1)
c5000 continue
c     goto 4999

      RETURN
      END

C =bottom========= ./elem/setshape.f
