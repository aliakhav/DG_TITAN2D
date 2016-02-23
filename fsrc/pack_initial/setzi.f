c---------------------------------------------------------------------
c   routine name       - setzi
c   latest revision    - jan, 1990
c   purpose            - routine initiates an integer vector
c   arguments :
c     in:        Ia    - an integer vector
c                 N    - number of elements to be zeroed      
c---------------------------------------------------------------------
c
      subroutine setzi(Ia,N)
c      
      include '../com_fem/syscom.blk'
      dimension Ia(*)
c      
      do 10 i=1,N
        Ia(i)=0
   10 continue
c      
      return
      end
c
      
