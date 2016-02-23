	function minval(arr,n,nx)

c	finds minimum value of array
	implicit real*8(a-h,o-z)
	real*8 minval,arr(nx) 
	minval = 1.d30 
	do i =1,n    
	if (minval .gt. arr(i)) then
	minval = abs(arr(i))
	end if
	end do

	return
	end

