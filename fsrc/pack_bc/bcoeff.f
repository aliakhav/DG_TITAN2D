c---------------------------------------------------------------------
c     latest revision   - September, 1992                                             
c     purpose           - routine specifies functions whose integrals   
c                         over the boundary of the element modify the   
c                         stiffness matrix and the rhs vector according 
c                         to type of boundary conditions required   
c     arguments:         
c     input             lcol   -
c                  nell,kind   - element number and kind           
c                     kindbc   - flag identifying type of boundary cond
c                     x,y      - cartesioan coordinates of a point cons
c                     rnx,rny  - unit normal vector on the boundary    
c                     rtx,rty  - unit tangent vector   
c    output           att,bt,pt, 
c                     c,gt,g   - values of required functions 
c----------------------------------------------------------------------
c
      subroutine bcoeff(nequ,
     .               kindbc,is,x,y,xi,rnx,rny,rtx,rty, 
     .               att,bt,pt,c,gt,g, xval, yval)
c
      include '../com_fem/syscom.blk'
      include '../com_fem/paramf.blk'
c      include '../com_fem/cflag.blk'
c  
      parameter (lcol=MAXEQS)
      dimension  att(lcol,lcol),
     .            bt(lcol,lcol),
     .            pt(lcol,lcol),
     .             c(lcol,lcol),
     .            gt(lcol     ),
     .             g(lcol     )      
c
c     set debug print flag
c      
      iprint = 0
c
      do i=1, lcol
        do j=1, lcol
          att(i,j) = 0.d0
           bt(i,j) = 0.d0
           pt(i,j) = 0.d0
            c(i,j) = 0.d0
        enddo
        gt(i) = 0.d0
         g(i) = 0.d0
      enddo
c

c        call bc_poisson(lcol,kindbc,is,x,y,xi,rnx,rny,rtx,
c     .                  rty,  c,g)

       call bcmax(lcol,kindbc,is,
     .                  x,y,xi,rnx,rny,rtx,rty,xval, yval,c,g)      
c       call bcexact(lcol,is, x,y,xi,rnx,rny,rtx,rty,c,g)      
c
 1000 continue
c      
      if (iprint.eq.1) then
        write(*,*) 'kindbc =',kindbc
        write(*,*) 'x,y = ',x,y
        write(*,*) 'rnx, rny = ',rnx,rny
        write(*,*) 'rtx, rty = ',rtx,rty        
        write(*,*) c(1,1), c(1,2)
        write(*,*) c(2,1),c(2,2)
        write(*,*) g(1),g(2)
      endif
c
      return
      end                  
c




