c------------------------------------------------------------------
c     routine name mass_elem      
c     lastest revision: April, 1993
c     arguments:
c       in:        ifg - 1 do not do element stffness matrix
c                      - otherwise: do element stiffness matrix
c                 ison - >= 0: do rhs using ison as the son number (used for refinement)
c                      - = -1: do rhs as going from 4 small elements to 1 big element (unrefinement)
c                      - -2:  do rhs as normal element stiffness rhs
c                      - <= -3: don't do rhs calculation
c               Norder - order of element
c               nedofr - number of dof of element
c                   Nc - length of the column of ek      
c                 nequ - number of equations
c                 xnod - the nodal coordinates
c     sol1, sol2, sol3, sol4 - the solution vectors from the original elements 
c                              (sol1 is lower left, sol2 is lower right, sol3 is upper right, ...)
c                   
c        out:       ek - the element stiffness matrix
c            ef - element load vector      
c    comment: this routine establish the mass matrix
c             using Hierarchical shape functions. The element geometry
c             is described by the normal quadratic shape functions.
c             it is assumed that the elements were split 1 to 4 with no 
c             grading (i.e. each element is the same geometry)
c------------------------------------------------------------------
c      
c      
      subroutine masselem(ifg,Norder,nedofr,nequ,Nc,
     1     xnod,ek,ef, sol1, sol2, sol3, sol4, ison)

      include '../com_fem/syscom.blk'
      include '../com_fem/paramf.blk'
      include '../com_fem/cint.blk'
c
c      include '../com_fem/celemf.blk'
c
      dimension ek(Nc,Nc),ef(Nc), sol1(Nc),sol2(Nc),sol3(Nc),sol4(Nc)
      dimension xi(2),xnod(2,9),xy(2), xy_interp(2)
      dimension gpsi(9),dgpsi(2,9),dxds(2,2),dsdx(2,2)
      dimension psi(9)
c      
c  ...local geometry function
      x(i,eta1,eta2) = Xnod(i,1)*(eta1**2-eta1)*(eta2**2-eta2)/4.d0
     .               + Xnod(i,2)*(eta1**2+eta1)*(eta2**2-eta2)/4.d0
     .               + Xnod(i,3)*(eta1**2+eta1)*(eta2**2+eta2)/4.d0
     .               + Xnod(i,4)*(eta1**2-eta1)*(eta2**2+eta2)/4.d0
     .               + Xnod(i,5)*(1.d0-eta1**2)*(eta2**2-eta2)/2.d0
     .               + Xnod(i,6)*(eta1**2+eta1)*(1.d0-eta2**2)/2.d0
     .               + Xnod(i,7)*(1.d0-eta1**2)*(eta2**2+eta2)/2.d0
     .               + Xnod(i,8)*(eta1**2-eta1)*(1.d0-eta2**2)/2.d0
     .               + Xnod(i,9)*(1.d0-eta1**2)*(1.d0-eta2**2)
c
      zero = 0.d0
      small = 1.e-5
      one = 1.d0    
c     
c     ...initialize element arrays
      if (ifg.ne.1) then
         do 10 i=1,nedofr
            call setz(ek(1,i),nedofr)
 10      continue
      endif
      call setz(ef,nedofr)
c     
c     ...find the order of polynormial in two directions
c     
      if (Norder .lt.  0 .or. Norder .gt. 3)  then
         write(*,*) 'WRONG order in masselem.f '
         stop
      endif
c     
      nedof = nedofr/nequ
c     
c     ...determine the element geometry
c     
c     .... define the number of integration point
      nl1 = norder+1
c     
c     ...loop through integration points
      do 30 l1 = 1, nl1
         do 40 l2= 1, nl1
            xi(1) = XIGAUS(l1,nl1)
            xi(2) = XIGAUS(l2,nl1)
            call gshape(xi, gpsi,dgpsi)
            xy(1) = x(1,xi(1),xi(2))
            xy(2) = x(2,xi(1),xi(2))
c     
c     ...calculate Jacobi matrix dx/ds(i,j)
            do 50 i=1,2
               do 50 j=1,2
                  dxds(i,j)=0.d0
                  do 50 k=1,9
                     dxds(i,j)=dxds(i,j)+dgpsi(j,k)*xnod(i,k)
 50         continue
c     
c     ...inverse jacobi matrix dx/ds(i,j) to get ds/dx(i,j)
            detj=dxds(1,1)*dxds(2,2)-dxds(1,2)*dxds(2,1)
            if (detj.le.0.0) go to 99
c     
            dsdx(1,1)=dxds(2,2)/detj
            dsdx(2,2)=dxds(1,1)/detj
            dsdx(1,2)=-dxds(1,2)/detj
            dsdx(2,1)=-dxds(2,1)/detj          
c     
c     ...accumulate integration point value of integrals
            fac=detj*WAGAUS(l1,nl1)*WAGAUS(l2,nl1)          
c     
c     ...evaluate values and derivatives of the shape functions
            call shape2dg(Norder,xi(1),xi(2), psi)
            
c     
c     ...first loop through components
            do 60 i=1,nequ
c     
c     ...if this is a resolution skip the stiffness calculations
               if (ifg.eq.1) goto 80            
               do 70 j=1,nequ
c     
               if(i.eq. j) then	
                  temp = fac
                  do k = 1, Nedof
                     temp1 = temp * psi(k)
                     do l = 1, Nedof
                        ek(nequ*(k-1)+i,nequ*(l-1)+j) =
     .                       ek(nequ*(k-1)+i,nequ*(l-1)+j)+temp1*psi(l)
                     enddo
                  enddo
               endif
 70         continue
c     
c     ...skip to here if doing a resolution           
 80         continue
c
c     ...accumulate for the load vector...
c     ...coefficient f(i) 
               
c     ...is ison >= 0, projecting to a small element            
            if(ison .ge. 0) then
c     ...first figure out coordinates of interpolating element
               if(ison .eq. 0) then
                  xy_interp(1) = xi(1)*.5d0-.5d0
                  xy_interp(2) = xi(2)*.5d0-.5d0
               else if(ison .eq. 1) then
                  xy_interp(1) = xi(1)*.5d0+.5d0
                  xy_interp(2) = xi(2)*.5d0-.5d0
               else if(ison .eq. 2) then
                  xy_interp(1) = xi(1)*.5d0+.5d0
                  xy_interp(2) = xi(2)*.5d0+.5d0
               else if(ison .eq. 3) then
                  xy_interp(1) = xi(1)*.5d0-.5d0
                  xy_interp(2) = xi(2)*.5d0+.5d0
               endif
               call sumshape(Norder,xy_interp,fac,sol1,temp,nedof)
c     ...small to big element
            else if(ison .eq. -1) then
               if(xi(1) .lt. -small .and. xi(2) .lt. -small) then
                  xy_interp(1) = xi(1)*2.d0+1.d0
                  xy_interp(2) = xi(2)*2.d0+1.d0
                  call sumshape(Norder,xy_interp,fac,sol1,temp,nedof)                 
               else if(xi(1) .gt. small .and. xi(2) .lt. -small) then
                  xy_interp(1) = xi(1)*2.d0-1.d0
                  xy_interp(2) = xi(2)*2.d0+1.d0
                  call sumshape(Norder,xy_interp,fac,sol2,temp,nedof)
               else if(xi(1) .gt. small .and. xi(2) .gt. small) then
                  xy_interp(1) = xi(1)*2.d0-1.d0
                  xy_interp(2) = xi(2)*2.d0-1.d0
                  call sumshape(Norder,xy_interp,fac,sol3,temp,nedof)
               else if(xi(1) .lt. -small .and. xi(2) .gt. small) then
                  xy_interp(1) = xi(1)*2.d0+1.d0
                  xy_interp(2) = xi(2)*2.d0-1.d0
                  call sumshape(Norder,xy_interp,fac,sol4,temp,nedof)
               else if(xi(1) .lt. -small) then
c     ...between element 0 and element 4 (average at interface because of DG)
                  xy_interp(1) = xi(1)*2.d0+1.d0
                  xy_interp(2) = 1.d0
                  call sumshape(Norder,xy_interp,fac,sol1,temp2,nedof)
                  temp = .5d0*temp2
                  xy_interp(2) = -1.d0
                  call sumshape(Norder,xy_interp,fac,sol4,temp2,nedof)
                  temp = temp + .5d0*temp2
               else if(xi(1) .gt. small) then
c     ...between element 1 and element 2 (average at interface because of DG)
                  xy_interp(1) = xi(1)*2.d0-1.d0
                  xy_interp(2) = 1.d0
                  call sumshape(Norder,xy_interp,fac,sol2,temp2,nedof)
                  temp = .5d0*temp2
                  xy_interp(2) = -1.d0
                  call sumshape(Norder,xy_interp,fac,sol3,temp2,nedof)
                  temp = temp + .5d0*temp2
               else if(xi(2) .lt. -small) then
c     ...between element 0 and element 1 (average at interface because of DG)
                  xy_interp(1) = 1.d0
                  xy_interp(2) = xi(2)*2.d0+1.d0
                  call sumshape(Norder,xy_interp,fac,sol1,temp2,nedof)
                  temp = .5d0*temp2
                  xy_interp(1) = -1.d0
                  call sumshape(Norder,xy_interp,fac,sol2,temp2,nedof)
                  temp = temp + .5d0*temp2
               else if(xi(2) .gt. small) then
c     ...between element 2 and element 3 (average at interface because of DG)
                  xy_interp(1) = 1.d0
                  xy_interp(2) = xi(2)*2.d0-1.d0
                  call sumshape(Norder,xy_interp,fac,sol2,temp2,nedof)
                  temp = .5d0*temp2
                  xy_interp(1) = -1.d0
                  call sumshape(Norder,xy_interp,fac,sol3,temp2,nedof)
                  temp = temp + .5d0*temp2
               else
c     ... at (0,0) (average at interface because of DG)
                  xy_interp(1) = 1.d0
                  xy_interp(2) = 1.d0
                  call sumshape(Norder,xy_interp,fac,sol1,temp2,nedof)
                  temp = .25d0*temp2
                  xy_interp(1) = -1.d0
                  call sumshape(Norder,xy_interp,fac,sol2,temp2,nedof)
                  temp = temp + .25d0*temp2
                  xy_interp(2) = -1.d0
                  call sumshape(Norder,xy_interp,fac,sol3,temp2,nedof)
                  temp = temp + .25d0*temp2
                  xy_interp(1) = 1.d0
                  call sumshape(Norder,xy_interp,fac,sol4,temp2,nedof)
                  temp = temp + .25d0*temp2
               endif
            else if(ison .eq. -2) then
c     ... normal rhs
               temp2 = 0.d0
               do k=1,nedof
                  temp2 = temp2 + sol1(nequ*(k-1)+i) * psi(k)
               enddo
               temp = temp2 * fac
            endif
            if(ison .ge. -2) then
               do k = 1, nedof
                  ef(nequ*(k-1)+i) = ef(nequ*(k-1)+i) + temp*psi(k)
                enddo
            endif
            
 60      continue
c
c.......end of first loop through components        
 40   continue
 30   continue
c
c.....end of loop through gaussian points
      
      return
c      
 99   write(*,100) detj,xi
 100  format(' in elem: bad jacobian = ',e10.3/2(2x,e10.3))
      do i=1,9
         write(*,*) xnod(1,i),xnod(2,i)
      enddo
c     
      stop
      end
c
c     routine sums up the shape functions multiplied by their solution values to 
c     give the solution at a point
      subroutine sumshape(Norder, xy_interp, fac, sol, sum, nedof)
      implicit none
      double precision sum, xy_interp(*), fac, sol(*), psi_interp(121)
      integer k, nedof, norder
      sum = 0
      call shape2dg(Norder,xy_interp(1),xy_interp(2),psi_interp)
      
      do k = 1, nedof
         sum = sum + psi_interp(k) * sol(k)
      enddo
      sum = sum * fac

      return
      end

      
