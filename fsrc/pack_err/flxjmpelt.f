c--------------------------------------------------------------------
c     Routine Name      - flxjmp
c     latest revision   - FEb 2000
c     purpose           - routine determines the boundary integrals
c                         of flux jumps
c     usage             - call flxjmp(matid,xy,Neig,Nside,Norder,nequ,Rhs)
c     arguments:
c     input             
c                        Norder    - order of element

c      output            Rhs - element stiffness matrix and Rhs vector
c                                 after adding boundary condition       
c-------------------------------------------------------------------
c      
      subroutine flxjmp(matid,xy,Neig_xy,Neig,Nside,Norder,ssol,tsol,
     &                  neq,u,Rhs,Neigh_order)
c     xy -cordinates of element
c     Neig_xy:cordinate of neighboring elements
c     neig :to indicate neighboring generation
c     Nside:side no of neighbor on which the current element is
c     ssol :element solution
c     tsol :solution of neighboring elements
c     nequ : no of equations
c     Rhs returns the value of flx jmp
c     order of neighboring elements

      include '../com_fem/syscom.blk'
      include '../com_fem/paramf.blk'
      include '../com_fem/cint.blk'

c
      dimension Norder(*),Rhs(*),ssol(2,121)
      dimension xy(2,9),Xnod(2,3),Vshap(11),xiloc(10),waloc(10)
      double precision Neig_xy(18,8),N_xy(2,9)
      dimension frh(11*MAXEQS),tsol(242,8),Nnorder(5)
      dimension in(11),sol(2,121), Neigh_order(5,8)

c
c
      dimension Neig(2,4),Nelb(4),Nside(4),g(2)
c
      x(i,eta) = Xnod(i,1)*eta*(eta-1.d0)/2.d0
     .         + Xnod(i,2)*(1.d0-eta*eta)
     .         + Xnod(i,3)*eta*(eta+1.d0)/2.d0
c
      xeta(i,eta) = Xnod(i,1)*(eta-0.5d0)
     .            - Xnod(i,2)*2.d0*eta
     .            + Xnod(i,3)*(eta+0.5d0)
c
      Nelb_con = 100
c

        call setz(Rhs,nequ*(Norder(is)+1))
c
c  ...loop through sides of the element
c
      do 20 is=1,4

        if (is.le.2) then 
          Xnod(1,1) = xy(1,is)
          Xnod(2,1) = xy(2,is)
          Xnod(1,3) = xy(1,is+1)
          Xnod(2,3) = xy(2,is+1)
        elseif (is.eq.3) then
          Xnod(1,1) = xy(1,4)
          Xnod(2,1) = xy(2,4)
          Xnod(1,3) = xy(1,3)
          Xnod(2,3) = xy(2,3)
        else
          Xnod(1,1) = xy(1,1)
          Xnod(2,1) = xy(2,1)
          Xnod(1,3) = xy(1,4)
          Xnod(2,3) = xy(2,4)            
        endif
        Xnod(1,2) = xy(1,is+4)
        Xnod(2,2) = xy(2,is+4)
c
c  ...Initiating f and frh  
c        
        call setz(frh,nequ*(Norder(is)+1))

c
        call bointg(Norder(is)+1,0, xiloc,waloc,igs) 
c
c  ...Specify the integral rule
c
        do 30 i=1,igs
          xglobl = x(1,xiloc(i))
          yglobl = x(2,xiloc(i))
          rtx = xeta(1,xiloc(i))
          rty = xeta(2,xiloc(i))
          ds = sqrt(rtx*rtx+rty*rty)
c
c      ...make sure the tangential direction is always count-clockwise
          if (is.ge.3) then
            rtx = -rtx
            rty = -rty
          endif
c
          rtx = rtx/ds
          rty = rty/ds
          rnx = rty
          rny = -rtx
c          
cc
cc compute flux jump

          if (is .eq. 1) then
             rnita = -1.d0
             rsita = xiloc(i)
             end if
          if (is .eq. 2) then
             rsita = 1.d0
             rnita = xiloc(i)
             end if
          if (is .eq. 3) then
             rnita = 1.d0
             rsita = xiloc(i)
             end if
          if (is .eq. 4) then
             rsita = -1.d0
             rnita = xiloc(i)
             end if
         call stress(matid,xy,rsita,rnita,ssol,sigma_x,sigma_y,sigma_xy,
     &  sigma_vm, Norder)

c
cc normal traction
           trac_nx = sigma_x*rnx + sigma_xy*rny
           trac_ny = sigma_xy*rnx + sigma_y*rny
           trac_n = sqrt (trac_nx*trac_nx+ trac_ny*trac_ny)

           js= Nside(is)
           js=js+1

cc for no neighbor i.e. bdry nside(is) < 0
cc case of boundary
           if (js .lt. 0) then
                 signa_x= -sigma_x
                 signa_xy=-sigma_y
                 signa_y=-sigma_xy
                 end if
                 if (js.lt.0) goto 60

           if ((Neig(2,is) .eq. 0) .and. (Neig(1,is) .gt. 0)) then
cc
cc case 1 : equal sized neighbor
cc              

          if (js .eq. 1) then
             rnita = -1.d0
             rsita = xiloc(i)
             end if
          if (js .eq. 2) then
             rsita = 1.d0
             rnita = xiloc(i)
             end if
          if (js .eq. 3) then
             rnita = 1.d0
             rsita = xiloc(i)
             end if
          if (js .eq. 4) then
             rsita = -1.d0
             rnita = xiloc(i)
             end if

cc get neighbor sol

                      do jjj=1,121
                         sol(1,jjj)=tsol(jjj*2-1,is)
                         sol(2,jjj)=tsol(jjj*2,is)
                         end do

cc get neighbor order
                         do iii=1,5
                            Nnorder(iii) = Neigh_order(iii,is)
                            end do
cc get neighbor cordinates
                          do iii=1,9
                            N_xy(1,iii)   = Neig_xy(iii*2-1,is)
                            N_xy(2,iii) = Neig_xy(iii*2,is)
                            end do

        call stress(matid,N_xy,rsita,rnita,sol,sigma_x,sigma_y,sigma_xy,
     &  sigma_vm, Nnorder)


              
         else
cc case of 2 small neighbors 
cc
                 if (Neig(2,is) .gt. 0) then




            if (xiloc(i) .le. 0) then 
            js= Nside(is)


           if (js .eq. 1) then
             rnita = -1.d0
             rsita = xiloc(i)
             rsita = rsita*2. +1.
             end if
          if (js .eq. 2) then
             rsita = 1.d0
             rnita = xiloc(i)
             rnita = rnita*2. +1.
             end if
          if (js .eq. 3) then
             rnita = 1.d0
             rsita = xiloc(i)
             rsita = rsita*2. +1.
             end if
          if (js .eq. 4) then
             rsita = -1.d0
             rnita = xiloc(i)
             rnita = rnita*2. +1.
             end if                   
c remap the gauss coordinate             

             if ((rsita .gt. 1.) .or. (rsita .lt. -1.)) then
                write(*,*) 'trouble flxjmp ',nel,kind,is,mnel,js,rsita
                end if

 

                      do jjj=1,121
                         sol(1,jjj)=tsol(jjj*2-1,is)
                         sol(2,jjj)=tsol(jjj*2,is)
                         end do
cc get neighbor order
                         do iii=1,5
                            Nnorder(iii) = Neigh_order(iii,is)
                            end do

cc get neighbor cordinates
                          do iii=1,9
                            N_xy(1,iii)   = Neig_xy(iii*2-1,is)
                            N_xy(2,iii) = Neig_xy(iii*2,is)
                            end do
        call stress(matid,N_xy,rsita,rnita,sol,sigma_x,sigma_y,sigma_xy,
     &  sigma_vm, Nnorder)
          end if
c     now for 2nd neighbor
          if (xiloc(i) .gt. 0.) then
          
c          call decode(neig(2,is),mnel,mknd)
           if (js .eq. 1) then
             rnita = -1.d0
             rsita = xiloc(i)
c remap the gauss coordinate             
             rsita = rsita*2. -1.
             end if
          if (js .eq. 2) then
             rsita = 1.d0
             rnita = xiloc(i)
c remap the gauss coordinate             
             rnita = rnita*2. -1.
             end if
          if (js .eq. 3) then
             rnita = 1.d0
             rsita = xiloc(i)
c remap the gauss coordinate             
             rsita = rsita*2. -1.
             end if
          if (js .eq. 4) then
             rsita = -1.d0
             rnita = xiloc(i)
c remap the gauss coordinate             
             rnita = rnita*2. -1
             end if                   

             if ((rsita .gt. 1.) .or. (rsita .lt. -1.)) then
                write(*,*) 'trouble flxjmp ',nel,kind,is,mnel,js,rsita
                end if

                      do jjj=1,121
                         sol(1,jjj)=tsol(jjj*2-1,is+4)
                         sol(2,jjj)=tsol(jjj*2,is+4)
                         end do
cc get neighbor order is+4 for 2nd neighbor
                         do iii=1,5
                            Nnorder(iii) = Neigh_order(iii,is+4)
                            end do

cc get neighbor cordinates
                          do iii=1,9
                            N_xy(1,iii)   = Neig_xy(iii*2-1,is+4)
                            N_xy(2,iii) = Neig_xy(iii*2,is+4)
                            end do
        call stress(matid,N_xy,rsita,rnita,sol,sigma_x,sigma_y,sigma_xy,
     &  sigma_vm, Nnorder)

          
               end if
                    end if

cc
cc case of big neighbor
cc
                    
                    if (Neig(1,is).lt.0) then
                    js= Nside(is)


           if (js .eq. 1) then
             rnita = -1.d0
             rsita = xiloc(i)
             end if
          if (js .eq. 2) then
             rsita = 1.d0
             rnita = xiloc(i)
             end if
          if (js .eq. 3) then
             rnita = 1.d0
             rsita = xiloc(i)
             end if
          if (js .eq. 4) then
             rsita = -1.d0
             rnita = xiloc(i)
             end if                   
c remap the gauss coordinate -- not correct for all meshes             
             rsita = rsita/2. -0.5


                      do jjj=1,121
                         sol(1,jjj)=tsol(jjj*2-1,is)
                         sol(2,jjj)=tsol(jjj*2,is)
                         end do
cc get neighbor order
                         do iii=1,5
                            Nnorder(iii) = Neigh_order(iii,is)
                            end do
cc get neigbor cordinates
                         do iii=1,9
                            N_xy(1,iii)   = Neig_xy(iii*2-1,is+4)
                            N_xy(2,iii) = Neig_xy(iii*2,is+4)
                            end do

        call stress(matid,N_xy,rsita,rnita,sol,sigma_x,sigma_y,sigma_xy,
     &  sigma_vm, Nnorder)

                   end if
cc

                 end if
c
   60     continue
cc normal traction
           tracn_nx = signa_x*rnx + signa_xy*rny
           tracn_ny = signa_xy*rnx + signa_y*rny
           tracn_n = sqrt (trac_nx*trac_nx+ trac_ny*trac_ny)

           g(1) = (trac_nx - tracn_nx)/2.
           g(2) = (trac_ny - tracn_ny)/2.

          call shape1(Norder(is),xiloc(i),vshap)
          do 40 k=1,Norder(is)+1
          do 40 ivar1 = 1, nequ
 40          frh(nequ*(k-1)+ivar1) = frh(nequ*(k-1)+ivar1)
     .                            + vshap(k)*g(ivar1)*ds*waloc(i)

   30   continue
c
        if (is.le.2) then
          in(1) = is-1
          in(2) = is
        elseif (is.eq.3) then
          in(1) = 3
          in(2) = 2
        else
          in(1) = 0
          in(2) = 3
        endif
c        
        iprev = 3
        do i=1,is-1
          iprev = iprev + Norder(i)-1
        enddo
        do i=1,Norder(is)-1
          in(2+i) = iprev+i
        enddo        
c
c  ...add f and frh to 
c          
        do 50 i=1,Norder(is)+1
        do 50 ivar1 = 1, nequ
          Rhs(nequ*in(i)+ivar1) = Rhs(nequ*in(i)+ivar1)
     .                       + frh(nequ*(i-1)+ivar1)
   50   continue  
   20 continue
c
      return
      end
c

      
