c---------------------------------------------------------------------
c   routine name       - eval2
c   latest revision    - May, 1993
c   purpose            - routine evaluates the projection of a 1-D solution
c                        at a point -0.5 and 0.5
c   arguments :
c     in:     Norder   - order of approximation
c             iposi    - 1 for left
c             iposi    - 2 for right
c             Dof      - dof in 2-D Hierarchical shape functions
c     out:    Val      - values of the solution
c----------------------------------------------------------------------
c
      subroutine eval2(Norder,dof,iposi, val)
c      
      include '../com_fem/syscom.blk'
      include '../com_fem/crrr.blk'
c
      dimension Dof(2,*),Val(2,*)
c
      if ((Norder.lt.1).or.(Norder.gt.10)) then
        write(*,*)'WRONG INPUT TO EVAL2 ! Norder = ',Norder
      endif
c
c  ...accumulate for the final values
      do 10 ivar=1,2
        s = 0.d0
        do 20 i=1, Norder-1
          s = 0.d0
          do 30 k=i,Norder-1
            s = s + Dof(ivar,k)*RRR(k+1,i+1,iposi)
   30     continue
          Val(ivar,i) = s
   20   continue
   10 continue   
c
      return
      end
c
    



