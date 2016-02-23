c---------------------------------------------------------------------
c     computer          - machine independent                        
c     arguments:                                                    
c     input        nell,kind   - element number and kind                       
c                     kindbc   - flag identifying type of boundary cond
c                        is    - side of element      
c                       x,y    - cartesioan coordinates of a point cons
c                        xi    - master coordinate of a point cons 
c                   rnx,rny    - unit normal vector on the boundary    
c                   rtx,rty    - unit tangent vector                    
c    output              c,g   - values of required functions    
c----------------------------------------------------------------------
c
      subroutine bcexact(lcol,is,
     .     x,y,xi,rnx,rny,rtx,rty,c,g)      
c
      include '../com_fem/syscom.blk'
      include '../com_fem/elatcon.blk'
      include '../com_fem/penalty_bc.blk'
      include '../com_fem/cflag.blk'
c
      dimension c(lcol,lcol),g(lcol)
c                                     
c---------------------------------------------------------------------      
c    E X E M P L E S
c                                                                
c       - plane linear elasticity ,                                     
c                                                                       
c                         static,                                       
c                         kinematic,                                    
c                         elastic support boundary conditions           
c
c       - the data below corresponds to the exact solution:                 
c                           u = (x**2 + y**2)/20                             
c                           v = (x**2 + y**2)/20                             
c                                                              
c                              .   bckind=3 (elastic support)             
c                              . -(un-un0)=(rk11*sigman+rk12*sigmat)    
c                              . -(ut-ut0)=(rk21*sigman+rk22*sigmat)    
c                     ...................(1,1)                          
c       (kinematic)   .        .        .                               
c       bckind=4      .        .        .                               
c un = -(x**2+y**2)   .        .(0,0)   .                               
c ut = -(x**2+y**2)...........................  bckind=2 (static)         
c                     .        .        .     sigman = 4*xmi*x + xlam*(2
c                     .        .        .     sigmat = 2*xmi*( x + y )  
c                     .        .        .                               
c                     ...................                               
c                              .                                        
c                           bckind=1 ("no penetration")                   
c                           un=-(x**2+y**2)                             
c                           sigmat=-2*xmi*(x+y)      
c....................................................................
c
c  set pure dirichlet bc's 
      kindbc = 4

c
c.....force and dispacement on boundary
c      
      xlam = xlam_array(matid)
      xmi = xmi_array(matid)
      if(iflag_bf. eq. 10) then
         sigmax = (2*(x+y)*xlam+4*x*xmi)/20.d0
         sigmay = (2*(x+y)*xlam+4*y*xmi)/20.d0
         sigmaxy = 2*xmi*(x+y)/20.d0
         sigman = sigmax*rnx**2+sigmay*rny**2+2*sigmaxy*rnx*rny
         sigmat = sigmaxy*(rnx**2-rny**2)+(sigmay-sigmax)*rnx*rny
         un = (x**2+y**2)*(rnx+rny)/20.d0
         ut = (x**2+y**2)*(rtx+rty)/20.d0 
      endif
      if(iflag_bf .eq. 20) then
         un = (x**3+y**3)*(rnx+rny)/30.d0
         ut = (x**3+y**3)*(rtx+rty)/30.d0 
      endif
c 
c....."no penetration" : un and sigmat are specified                    
c   
      if (kindbc.eq.1) then
        c(1,1) = 1/eps*rnx*rnx                         
        c(1,2) = 1/eps*rnx*rny               
        c(2,1) = 1/eps*rny*rnx                                    
        c(2,2) = 1/eps*rny*rny
        g(1) = sigmat*rtx + 1/eps*un*rnx
        g(2) = sigmat*rty + 1/eps*un*rny                           
      endif                                                           
c                                                                        
c.....static boundary conditions: sigmat,sigman are specified           
c                                                                       
      if (kindbc.eq.2) then
        g(1) = sigman*rnx + sigmat*rtx                           
        g(2)= sigman*rny + sigmat*rty                           
      endif
c
c.....elastic foundation
c      
      if (kindbc.eq.3) then
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
      if (kindbc.eq.4) then                        
        c(1,1) = 1.d0/eps
        c(2,2) = 1.d0/eps
        g(1) = 1.d0/eps*(un*rnx + ut*rtx)
        g(2) = 1.d0/eps*(un*rny + ut*rty) 
      endif
c
      return
      end
c
