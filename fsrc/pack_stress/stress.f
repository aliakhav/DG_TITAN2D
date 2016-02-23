c------------------------------------------------------------------
c         routine name: stress.f
c     lastest revision: April, 1993
c              purpose: return stresses at a given point (xi1,xi2)   
c     arguments:
c                   in: nel,kind - number and kind of the element
c                       xi1,xi2 - two master coordinates
c                  out: sigma_x,sigma_y,sigma_xy 
c------------------------------------------------------------------
c      
      subroutine stress(matid,xnod,xi1,xi2,sol,sigma_x,sigma_y,sigma_xy,
     &	sigma_vm, No_order)
c      
      include '../com_fem/syscom.blk'
c      include '../com_fem/paramf.blk'
      include '../com_fem/elatcon.blk'
c      include '../com_fem/cmat.blk'
c
      dimension Sol_du(2,2,121),Norder(5),du(2,2),psi(121),xnod(2,9)
      dimension Sol(2,121),No_order(5)
c
      do i=1,5
      Norder(i)=No_order(i)
      end do
c  set xlam and xmi by matid
      xlam = xlam_array(matid)
      xmi = xmi_array(matid)
c        end if

      call sol_fluxf(2,2, np1,Sol_du,xnod,sol,norder)
c
      Norder(1) =  np1-1
      Norder(2) =  np1-1
      Norder(3) =  np1-1
      Norder(4) =  np1-1
      Norder(5) =  np1-1
c
      du(1,1) = 0.d0
      du(2,1) = 0.d0
      du(1,2) = 0.d0
      du(2,2) = 0.d0
c      
      call shape2(Norder,xi1,xi2, psi)
      np2 = np1*np1
      do i=1,np2
        du(1,1) = du(1,1) + Sol_du(1,1,i)*psi(i)
        du(2,1) = du(2,1) + Sol_du(2,1,i)*psi(i)
        du(1,2) = du(1,2) + Sol_du(1,2,i)*psi(i)
        du(2,2) = du(2,2) + Sol_du(2,2,i)*psi(i)        
      enddo
c

c
c  ...case of plain strain
      sigma_x = xlam*(du(1,1)+du(2,2)) + xmi*(du(1,1)+du(1,1))
      sigma_y = xlam*(du(1,1)+du(2,2)) + xmi*(du(2,2)+du(2,2))
      sigma_xy = xmi*(du(1,2)+du(2,1))
c
      sigma_vm = sqrt((sigma_x - sigma_y)**2 + sigma_xy**2)

      return
      end
c









