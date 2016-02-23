c-----------------------------------------------------------------------
c   routine name       - shape1
c
c   computer           - machine independent
c
c   latest revision    - jan 1990
c
c   purpose            - routine calculates 1-d master element
c                        hierarchical shape functions      
c
c   usage              - call shape1(Idg,X,F)
c
c   arguments :
c     in:
c               Idg    - order of approximation
c               X      - local coordinate (between -1 and 1)      
c     out:
c               F      - a one dimensional array containing
c                        values of the shape functions      
c
c   required  routines - none
c-----------------------------------------------------------------------
c
      subroutine shapeb(xb,yb,idg, p)
c
      include '../com_fem/syscom.blk' 
c
      dimension F(11)
c
      if (dabs(dabs(yb)-1.d0).le.1.d-8) then
        x=xb
      else
        x=yb
      endif  
c  ...linear shape functions      
      f(1)=(1-x)*0.5
      f(2)=(1+x)*0.5
c
c  ...higher order shape functions      
      xp=x
      x1=x
      x2=1
      do 10 i=3,idg+1
        c =x1
        x1=x2
        x2=c
        xp=xp*x
        f(i)=(xp-x1)
   10 continue
c      
      p=F(idg)
c
      return
      end
c

