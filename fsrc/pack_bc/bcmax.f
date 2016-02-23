c---------------------------------------------------------------------
c     computer          - machine independent                        
c     arguments:                                                    
c     input        nell,kind   - element number and kind                       
c                     kindbc   - flag identifying type of boundary cond
c                        is    - side of element      
c                       x,y    - cartesioan coordinates of a point cons
c                        xi    - master coordinate of a point cons
c                        xval  - x-dir constraint
c                        xval  - y-dir constraint
c                   rnx,rny    - unit normal vector on the boundary    
c                   rtx,rty    - unit tangent vector                    
c    output              c,g   - values of required functions    
c----------------------------------------------------------------------
c
      subroutine bcmax(lcol,kindbc,is,
     .                  x,y,xi,rnx,rny,rtx,rty,xval, yval,  c,g)      
c
      include '../com_fem/syscom.blk'
      include '../com_fem/elatcon.blk'
      include '../com_fem/penalty_bc.blk'
c
      dimension c(lcol,lcol),g(lcol)
c                                     
c---------------------------------------------------------------------      
c    E X A M P L E S
c                                                                
c       - plane linear elasticity ,                                     
c                                                                       
c                         static,                                       
c                         kinematic,                                    
c                         elastic support boundary conditions           
c
c       - the data below corresponds to the exact solution:                 
c
c                     ...................(1,1)                          
c       (kinematic)   .        .        .                               
c       bckind=5      .        .        .                               
c un = 0   .                   .(0,0)   .                               
c ut = 0              ...........................  bckind=1,2 (static)         
c                     .        .        .     sigman = 4*xmi*x + xlam*(2
c                     .        .        .     sigmat = 2*xmi*( x + y )  
c                     .        .        .                               
c                     ...................                               
c                              .                                        
c                           bckind=3,4 ("no penetration")                   
c                           un= 0
c                           sigmat=0
c....................................................................
c

c
c.....force and dispacement on boundary
c      


      sigmat = 0.
      sigman = 0.
c 
c....."no penetration" : un and sigmat are specified                    
c   
      if ((kindbc.eq.3).or.(kindbc.eq. 4)) then
        c(1,1) = 1/eps                         
        c(1,2) = 1/eps            
        c(2,1) = 1/eps
        c(2,2) = 1/eps
        g(1) = sigmat + 1/eps*xval
        g(2) = sigmat + 1/eps*yval
      endif                                                           
c                                                                        
c.....static boundary conditions: sigmat,sigman are specified           
c                                                                       
      if (kindbc.eq.1) then
        g(1) = xval 
        g(2)= yval
      endif
c
      if ((kindbc.eq.2)) then
c         sigman= 0.e6
c         sigmat=0.
        g(1) =  xval                        
        g(2)= yval                           
      endif


c.....elastic foundation -- not correct
c      
      if (kindbc.eq.6) then
         write(*,*) 'warning in bcmax.f'
c          
c   ...set coefficients characterizing the foundation
c          
        rk11=2                                                      
        rk12=1                                                      
        rk21=1                                                      
        rk22=2                                                      
c
c   ...find the inverse of the above                               
c
        rl11=.6666667                                               
        rl12=-.3333333                                              
        rl21=-.3333333                                              
        rl22=.6666667                                               
c                                                                        
c    ...construct un0,ut0 for the assumed exact solution            
c                                             
        un0 = un + rk11*sigman + rk12*sigmat                        
        ut0 = ut + rk21*sigman + rk22*sigmat                 
        c(1,1) = rl11*rnx*rnx + rl12*rtx*rnx + rl21*rnx*rtx        
     .         + rl22*rtx*rtx        
        c(1,2) = rl11*rnx*rny + rl12*rtx*rny + rl21*rnx*rty          
     .         + rl22*rtx*rty          
        c(2,1) = rl11*rny*rnx + rl12*rty*rnx + rl21*rny*rtx         
     .         + rl22*rty*rtx         
        c(2,2) = rl11*rny*rny + rl12*rty*rny + rl21*rny*rty         
     .         + rl22*rty*rty
        g(1) = (rl11*un0+rl12*ut0)*rnx + (rl21*un0+rl22*ut0)*rtx               
        g(2) = (rl11*un0+rl12*ut0)*rny + (rl21*un0+rl22*ut0)*rty   
      endif        
c
c.....kinematic boundary conditions: ut and un are specified  
c            
      if (kindbc.eq.5) then                        
        c(1,1) = 1.d0/eps
        c(2,2) = 1.d0/eps
        g(1) = 1.d0/eps*xval
        g(2) = 1.d0/eps*yval
      endif
c
      return
      end
c
