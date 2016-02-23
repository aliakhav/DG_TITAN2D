c-----------------------------------------------------------------------
c
c   routine name       - rhsub
c
c   computer           - machine independent
c
c   latest revision    - april 2, 1988
c
c   purpose            - backward gauss substitution
c
c   usage              - call rhsub(a,b,n,m,x)
c
c   arguments in: a    - m by m matrix
c                 b    - vector of length m containing the right
c                        hand side
c                 n    - row dimension of matrix a exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 m    - number of equations
c             out:x    - vector of length m containing the solution
c
c   required routines  - none
c
c-----------------------------------------------------------------------
c
      subroutine rhsub(a,b,n,m, x)
c
      include '../com_fem/syscom.blk'
c
      dimension a(n,*),x(*),b(*)

cc      do j=1,m
cc         write(*,*) "predictor,el_rhs", x(j),b(j)
cc      enddo
c
      m1=m-1
c.....begin forward reduction of right hand side
      do 20  i=1,m1
      j1=i+1
      do 10 j=j1,m
10    b(j)=b(j)-b(i)*a(j,i)/a(i,i)
c      if(a(i,i).eq.0) then
c         PAUSE
c      endif
      
20    continue
c.....begin back substitution
      x(m)=b(m)/a(m,m)
      do 40 i=1,m1
      ib=m-i
      j1=ib+1
      do 30  j=j1,m
30    b(ib)=b(ib)-a(ib,j)*x(j)
      x(ib)=b(ib)/a(ib,ib)
c      if(a(ib,ib).eq.0) then
c         PAUSE
c      endif
c      write(*,*) x(ib), a(ib,ib),b(ib)
40    continue
      return
      end
c
