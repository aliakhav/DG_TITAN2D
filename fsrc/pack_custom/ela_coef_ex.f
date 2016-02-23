c
      subroutine ela_coef(lcol,
     .                    x,y, a11,a12,a21,a22,b1,b2,d1,d2,c,f,fx,fy)
c          
      include '../com_fem/syscom.blk'
      include '../com_fem/elatcon.blk'
      include '../com_fem/cflag.blk'
c
      dimension a11(lcol,*),a12(lcol,*),
     .          a21(lcol,*),a22(lcol,*),
     .           b1(lcol,*), b2(lcol,*),
     .           d1(lcol,*), d2(lcol,*),
     .            c(lcol,*),
     .            f(lcol),fx(lcol),fy(lcol)
      dimension xy(2)
c
      do 10 j=1,lcol
        do 9 i=1,lcol
          a11(i,j) = 0.d0       
          a12(i,j) = 0.d0                  
          a21(i,j) = 0.d0         
          a22(i,j) = 0.d0    
           b1(i,j) = 0.d0                   
           b2(i,j) = 0.d0                              
           d1(i,j) = 0.d0
           d2(i,j) = 0.d0
            c(i,j) = 0.d0
    9    continue
         f(j)  = 0.d0                       
        fx(j)  = 0.d0                        
        fy(j)  = 0.d0
   10 continue
c                                                       
c...................................................................
c     example :                               
c     iflag_bg = 10
c     - vector valued equation : plane linear elasticity problem 
c       coefficients correspond to the exact solution:     
c                            u = (x**2 + y**2)/20.d0    
c                            v = (x**2 + y**2)/20.d0  
c       constant lame coefficients xmi, xlam are assumed.            
c
c     iflag_bg = 20
c     - vector valued equation : plane linear elasticity problem 
c       coefficients correspond to the exact solution:     
c                            u = (x**3 + y**3)/30.d0    
c                            v = (x**3 + y**3)/30.d0  
c       constant lame coefficients xmi, xlam are assumed.            
c..................................................................
c                                                                
c     linear elasticity                                          
c                                               
      xlam = xlam_array(matid)
      xmi = xmi_array(matid)

      a11(1,1)=2.*xmi+xlam                                       
      a22(1,1)=xmi
                                                          
      a12(1,2)=xlam                                    
      a21(1,2)=xmi                                         

      a12(2,1)=xmi                                          
      a21(2,1)=xlam                                          

      a11(2,2)=xmi                                   
      a22(2,2)=2.*xmi+xlam                         
c    
c      iflag_bf = 20
      iflag_bf = 0

      if (iflag_bf.eq.0) then
        f(1) = 0.d0
        f(2) = 0.d0        
        goto 1000
      endif
c
      if (iflag_bf.eq.1) then
        f(1) = 0.d0
        f(2) = 0.d0        
        goto 1000
      endif      
c
c  
      if (iflag_bf.eq.10) then
        f(1) = -(6.d0*xmi+2.d0*xlam)/20.d0
        f(2) = -(6.d0*xmi+2.d0*xlam)/20.d0        
        goto 1000
      endif
c
      if (iflag_bf.eq.20) then
        f(1) = -(12.d0*xmi*x+6.d0*xmi*y+6.d0*xlam*x)/30.d0
        f(2) = -(6.d0*xmi*x+12.d0*xmi*y+6.d0*xlam*y)/30.d0        
        goto 1000
      endif
   
c       
 1000 return
      end                                                             
c                              
