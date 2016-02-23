c-----------------------------------------------------------------------
c   routine name       - gausse
c
c   computer           - machine independent
c
c   latest revision    - April 2, 1988
c
c   purpose            - solution of a full matrix linear complex
c                        system of equations using gauss elimination
c                        without pivoting
c
c   usage              - call gausse(gk,igk,gf,u,n)
c
c   arguments in: gk   - n by n complex matrix
c                 igk  - row dimension of matrix gk exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 gf   - RHS vector of length n 
c                 n    - number of equations
c             out:u    - a vector of length n containing
c                        the solution
c
c-----------------------------------------------------------------------
c
      subroutine gausse(gk,igk,gf,n,  u)
c
      include '../com_fem/syscom.blk'
      dimension gk(igk,*),gf(*),u(*)
c
c ...set debug print flag  
c      
      iprint=0
c
c ...debug print
c      
      if (iprint.eq.1) then
        open(unit=14,file='file2',
     .  form='formatted',access='sequential',status='unknown')
        do 10 i=1,n
c          write(14,*) '--------------------------------------'
c          write(14,*) 'NO OF EQUATION AND RHS = ',
c     .                i,gf(i)
c          write(14,*) 'STIFFNESS MATRIX'
c          write(14,1000) (gk(i,j),j=1,n)
   10   continue
 1000   format(8(e14.8,1x))        
      endif     
c
c ...end of debug print
c      
      call tri(gk,igk,n)
c     
      call rhsub(gk,gf,igk,n,  u)
c
      return
      end
c
