	function maxval(arr,n,nx)
c	finds maximum value of array
	implicit real*8(a-h,o-z)
	real*8 maxval,arr(nx) 
	maxval = 0. 
	do i =1,n    
	if (maxval .lt. arr(i)) then
	maxval = abs(arr(i))
	end if
	end do

	return
	end

