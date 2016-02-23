c---------------------------------------------------------------------
c     computer: alliant fx/8
c
c     lastest revision: October, 1991
c
c     routine name: gshape
c
c     usage: call gshape (xi,n,psi,dpsi)
c
c     arguments: in - xi:   the master element coordinate
c      
c                out- psi:  the values of the shape functions
c                     dpsi: the v. of the shape functions' derivatives
c
c     dimension: psi(9)
c                dpsi(2,9)
c
c     explanation: this routine returns the values of the shape
c                  functions (psi) and their derivatives (dpsi)
c                  with respect to the element coordinate (xi)
c                  for 9 node quadrilateral elements.
c---------------------------------------------------------------------
c
      subroutine gshape (xi, psi,dpsi)
c
      include '../com_fem/syscom.blk'
      dimension psi(9),dpsi(2,9),xi(2)
c
      psi(1)=(1./4.)*(xi(1)**2-xi(1))*(xi(2)**2-xi(2))
      psi(2)=(1./4.)*(xi(1)**2+xi(1))*(xi(2)**2-xi(2))
      psi(3)=(1./4.)*(xi(1)**2+xi(1))*(xi(2)**2+xi(2))
      psi(4)=(1./4.)*(xi(1)**2-xi(1))*(xi(2)**2+xi(2))
      psi(5)=(1./2.)*(1.-xi(1)**2)*(xi(2)**2-xi(2))
      psi(6)=(1./2.)*(xi(1)**2+xi(1))*(1.-xi(2)**2)
      psi(7)=(1./2.)*(1.-xi(1)**2)*(xi(2)**2+xi(2))
      psi(8)=(1./2.)*(xi(1)**2-xi(1))*(1.-xi(2)**2)
      psi(9)=(1.-xi(1)**2)*(1.-xi(2)**2)
      dpsi(1,1)=(0.5*xi(1)-0.25)*(xi(2)**2-xi(2)) 
      dpsi(2,1)=(0.5*xi(2)-0.25)*(xi(1)**2-xi(1)) 
      dpsi(1,2)=(0.5*xi(1)+0.25)*(xi(2)**2-xi(2)) 
      dpsi(2,2)=(0.5*xi(2)-0.25)*(xi(1)**2+xi(1)) 
      dpsi(1,3)=(0.5*xi(1)+0.25)*(xi(2)**2+xi(2)) 
      dpsi(2,3)=(0.5*xi(2)+0.25)*(xi(1)**2+xi(1)) 
      dpsi(1,4)=(0.5*xi(1)-0.25)*(xi(2)**2+xi(2)) 
      dpsi(2,4)=(0.5*xi(2)+0.25)*(xi(1)**2-xi(1)) 
      dpsi(1,5)=-xi(1)*(xi(2)**2-xi(2)) 
      dpsi(2,5)=(xi(2)-0.5)*(1.-xi(1)**2)
      dpsi(1,6)=(xi(1)+0.5)*(1.-xi(2)**2)
      dpsi(2,6)=-xi(2)*(xi(1)**2+xi(1)) 
      dpsi(1,7)=-xi(1)*(xi(2)**2+xi(2)) 
      dpsi(2,7)=(xi(2)+0.5)*(1.-xi(1)**2)
      dpsi(1,8)=(xi(1)-0.5)*(1.-xi(2)**2)
      dpsi(2,8)=-xi(2)*(xi(1)**2-xi(1)) 
      dpsi(1,9)=-2.*xi(1)*(1.-xi(2)**2) 
      dpsi(2,9)=-2.*xi(2)*(1.-xi(1)**2) 
c
      return
      end 
c






