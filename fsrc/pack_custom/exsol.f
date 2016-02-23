cccc
ccc
c     xya - coordinates
c     uex - exact solution
c    uexd - derivative of exact solution 
c           uexd(1) = du/dx
c           uexd(2) = du/dy
c           uexd(3) = dv/dx
c           uexd(4) = dv/dy
cccc
      subroutine exsol(xya,uex,uexd)
c
      include '../com_fem/syscom.blk'
      include '../com_fem/cflag.blk'

c

      dimension xya(2),uex(2),uexd(4)

      if(iflag_bf .eq. 10) then
c       - the data below corresponds to the exact solution:                 
c                           u = (x**2 + y**2)/20                             
c                           v = (x**2 + y**2)/20 
      
         uex(1) = (xya(1)*xya(1) + xya(2)*xya(2))/20.
         uex(2) = (xya(1)*xya(1) + xya(2)*xya(2))/20.
         uexd(1) = xya(1)/10.
         uexd(2) = xya(2)/10.
         uexd(3) = xya(1)/10.
         uexd(4) = xya(2)/10.
      endif

      if(iflag_bf .eq. 20) then
c       - the data below corresponds to the exact solution:                 
c                           u = (x**3 + y**3)/30                             
c                           v = (x**3 + y**3)/30 
      
         uex(1) = (xya(1)**3 + xya(2)**3)/30.
         uex(2) = (xya(1)**3 + xya(2)**3)/30.
         uexd(1) = (xya(1)**2)/10.
         uexd(2) = (xya(2)**2)/10.
         uexd(3) = (xya(1)**2)/10.
         uexd(4) = (xya(2)**2)/10.
      endif


      return
      end
