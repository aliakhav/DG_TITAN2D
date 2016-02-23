c----------------------------------------------------------------------
c   routine name       - setz
c   latest revision    - jan 1990   by  Geng Po
c   purpose            - initiates a real array with zeros
c   arguments :
c     in:         a    - a real array
c                 n    - total length of the array 
c-----------------------------------------------------------------------
c
      subroutine setz(a,n)
c      
      include '../com_fem/syscom.blk'
      dimension a(*)
c      
      do 10 i=1,n
        a(i)=0.d0
   10 continue
c      
      return
      end
c
