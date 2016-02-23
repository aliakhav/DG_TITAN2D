c----------------------------------------------------------------------
c   routine name       - eval1
c   latest revision    - May, 1993
c   purpose            - routine evaluates value of a 1-D solution
c                        at a point with master coodinate Eta
c   arguments :
c     in:     Norder   - order of approximation
c             Eta      - the master coordinate 
c             Dof      - dof in hierarchical shape function
c     out:    Val      - values of the solution
c----------------------------------------------------------------------
c
      subroutine eval1(Norder,Eta,Dof, Val)
c      
      include '../com_fem/syscom.blk'
c
      dimension Dof(2,*),Val(2)
      dimension vshap(11)
c
c  ...explanation of local variables:
c     vshap   - 1-D hierarchical shape functions
c
      if ((Norder.lt.1).or.(Norder.gt.9)) then
        write(*,*)'WRONG INPUT TO EVAL1 ! Norder = ',Norder
      endif
c
c  ...evaluate shape functions at the point
      call shape1(Norder,Eta, vshap)
c
c  ...accumulate for the final values
      do 20 ivar=1,2
        s = 0.d0
        do 10 k=1,Norder+1
          s = s + Dof(ivar,k)*vshap(k)
   10   continue
        Val(ivar) = s
   20 continue
c
      return
      end
c


