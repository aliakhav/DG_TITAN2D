c-------------------------------------------------------------------
c-------------------------------------------------------------------
c
	subroutine initial(matcount,xlam_in,xmi_in)
c
	include '../com_fem/syscom.blk'
	include '../com_fem/penalty_bc.blk'
	include '../com_fem/cboundary.blk'
	include '../com_fem/elatcon.blk'

	dimension xlam_in(matcount),xmi_in(matcount)
c
c
c  ...set some constants for boundary condition and
c     control parameter      
c
	eps = 1.d-12
	Nelb_con1 = 200
	Nelb_con2 = 100
c
c  ...set material properties
	if(matcount .gt. MAXMATS) then
	   write(*,*) 'Too many materials, MAXMATS is',MAXMATS
	   stop
	endif
	do i=1,matcount
	   xlam_array(i) = xlam_in(i)
	   xmi_array(i) = xmi_in(i)
	enddo	
c
c  ...Number of equation      
c
c      
        call cnst
	call setshape
c       
	return
	end
c










