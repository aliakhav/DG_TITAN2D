      subroutine createfunky(xmax, xmin, ymax, ymin, ny)
      integer nx,ny
c     X will be the down-slope direction -- keep consistent with savage
      double precision tiny
      parameter (tiny=0.0000001d0)
      
c     simple code to generate topographical data
c     aug '00  ---  flat surface
      double precision xa(10000), ya(10000)
      integer i,j, icounter
      double precision xmin, xmax, ymin, ymax, xx, yy, zz
      integer i1, j1

c     topo data
      double precision xt(10000), yt(10000)
      double precision xlength, ylength

      xlength = .98d0*(xmax-xmin)
      ylength = .98d0*(ymax - ymin)
      xmin = xmin+.1d0*xlength
      xmax = xmax - .1d0*xlength
      xlength = xmax - xmin
      ymin = ymin+.1d0*ylength
      ymax = ymax - .1d0*ylength
      ylength = ymax - ymin
 
      nx = ny*(xlength/ylength)
      if(nx .le.  0) then
         nx = 10
         ny = nx*(ylength/xlength)
      endif
      do i=1, nx+1
         xt(i) = xlength*(1.d0*(i-1)/(1.d0*nx))+xmin
c         write(*,*) xt(i)
      enddo
      do i=1, ny+1
         yt(i) = ylength*(1.d0*(i-1)/(1.d0*ny))+ymin
c         write(*,*) yt(i)
      enddo

      do i=1, nx
         xa(2*i-1) = xt(i)
         xa(2*i) = .5*(xt(i) + xt(i+1))
      enddo
      xa(2*nx+1) = xt(nx+1)

      do i=1, ny
         ya(2*i-1) = yt(i)
         ya(2*i) = .5*(yt(i) + yt(i+1))
      enddo
      ya(2*ny+1) = yt(ny+1)         

c*************************************************************
c     afeapi input file
c*************************************************************
      open(unit=2,file='funky.dat')
c     nodes
      write(2,1977) (2*nx+1)*(2*ny+1)-nx*ny
c     elements
      write(2,1977) nx*ny
c     essential bc's
      write(2,1977) 2*(nx+ny)
c     natural bc's
      write(2,1977) 0
c     number of materials
      write(2,1977) 1

c     nodes
      icounter = 1
      do j=1,2*ny+1
         do i=1,2*nx+1
c     format is: node #, x coord, y coord
            if(mod(i,2) .ne. 0 .or. mod(j,2) .ne. 0) then
               write(2,1988) icounter,xa(i),ya(j)
               icounter = icounter + 1
            endif
         enddo
      enddo

c     elements
      icounter = 1
      do j=1,ny
         do i=1,nx
            write(2,1966) icounter, 2*i-1+(j-1)*(3*nx+2), 
     1           2*i+1+(j-1)*(3*nx+2), 2*i+1+j*(3*nx+2),
     2           2*i-1+j*(3*nx+2), 2*i+(j-1)*(3*nx+2), 
     3           i+1+j*(2*nx+1)+(j-1)*(nx+1), 2*i+j*(3*nx+2),
     4           i+j*(2*nx+1)+(j-1)*(nx+1), 1, i-1, j-1        
            icounter = icounter + 1
         enddo
      enddo

c     boundary conditions
      do i=1,nx
         write(2,1955) 2*i, 0.d0, 0.d0
         write(2,1955) 2*i+(2*nx+1)*2*ny-nx*ny, 0.d0, 0.d0
      enddo
      do j=1,ny
         write(2,1955) j*(2*nx+1)+(j-1)*nx+j, 0.d0, 0.d0
         write(2,1955) j*(2*nx+1)+(j-1)*nx+j+nx, 0.d0, 0.d0         
      enddo
      write(2,*) 'stones'


      close(2)
 1977 format('',i8)
 1999 format('',3(f9.3,1x),f9.6)
 1988 format('',i7,2(f19.9,1x))
 1966 format('',12(i8,1x))
 1955 format('',i6, 2(1x,f7.3))

      stop
      end
      

