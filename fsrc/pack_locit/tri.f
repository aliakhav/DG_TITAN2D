c-----------------------------------------------------------------------
c   routine name       - tri
c
c   computer           - machine independent
c
c   latest revision    - april 2, 1988
c
c   purpose            - forward gauss elimination
c
c   usage              - call tri(a,n,m)
c
c   arguments in: a    - m by m matrix
c                 n    - row dimension of matrix a exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c             out:a    - a after gauss factorization
c
c   required routines  - none
c-----------------------------------------------------------------------
c
      subroutine tri(a,n,m)
c
      include '../com_fem/syscom.blk'
c
      dimension a(n,n)
c
cc      do i=1,m
cc         do j=1,m
cc            write(*,*) a(j,i)
cc         enddo
cc      enddo

      m1=m-1
      tiny=1.e-30
      zero= 0.d0
c
c.....eliminate degree of freedom i
      do 30  i=1,m1
c.....check for excessively small pivot
      if(dabs(a(i,i)).lt.tiny) go to 99
      j1=i+1
c.....modify row j
      do 20  j=j1,m
      if(a(j,i).eq.zero) go to 20
      fac=a(j,i)/a(i,i)
      do 10 k=j1,m
10    a(j,k)=a(j,k)-a(i,k)*fac
20    continue
30    continue

c      do j=1,m
c         write(*,*) a(j,j)
c      enddo
      
      return
99    write (*,100)i,dabs(a(i,i))
100   format(/' reduction failed due to small pivot'/
     .   ' equation no.',i5,' pivot ',e10.3)
      stop
      end
c





