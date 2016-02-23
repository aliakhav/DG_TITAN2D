c-------------------------------------------------------------------
c     Routine name: elemcom
c     Lastest revision: May, 1993
c     Arguments
c       in:        ifg -  <1> return rhs only
c                         <2> both lhs and rhs
c                 nequ -  number of equations
c                 ndff -  number of degrees of freedom
c                 Nelb -  boundary condition type
c                matid -  material flag 
c              bcvalue -  boundary conditionr array
c      out:         ek -  element coefficient matrix
c                   ef -  element rhs
c   explation: routine elemcom combines three related element routines:
c                   elem
c                 bcelem
c                 elmcon    
c              to give the element stiffness matrix and rhs      
c-----------------------------------------------------------------------
c-----789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
      subroutine elemcom(ifg,nequ,ndff,Nc,Norder,Nelb,Icon,
     .  xnod,ek,ef, matid, bcvalue)
c
      include '../com_fem/syscom.blk'
c      include '../com_fem/paramf.blk'
c      include '../com_fem/configf.blk'

c
      dimension ek(Nc,Nc),ef(*),xnod(2,9),bcvalue(*) 
      dimension Norder(5),Nelb(4),Icon(4)

c  ...element stiffness matrix      
      call pelem(ifg,Norder,ndff,nequ,Nc,xnod,ek,ef,matid)
c
c  ...apply boundary condition 
      call bcelem(ifg,Norder,ndff,nequ,Nc,Nelb,xnod,ek,ef, bcvalue)
c  ...take care of constrained node
      call elmcon(ifg,Norder,ndff,nequ,Nc,Icon,ek,ef)



      return
      end
c





