C file:  ./mutil/blas/decomp.f
C=========================
      subroutine decomp(n,na,a,ip,Iflag)
      include '../com_fem/syscom.blk'
      dimension a(na,n),ip(n)

      pivmin=1.d30

      Iflag=0
      inc=1

      do 10 k=1,n
      imax=idamax(n-k+1,a(k,k),inc)

      ip(k)=k+imax-1
      if(  abs( a(ip(k),k)) .le. 1.e-10 )
     .  then
            write(*,*)'PIVOT=', abs(a(ip(k),k))
            pause
        endif
      pivmin=dmin1(pivmin,abs(a(ip(k),k)) )
      if(  abs( a(ip(k),k)) .ne. 0.e0 ) go to 20
      Iflag=1
      return

 20   c=a(ip(k),k)
      a(ip(k),k)=a(k,k)
      a(k,k)=c
      c=1.d0/c
c      call sscal(n-k,c,a(k+1,k),inc)
      ink = n-k
      call dscal(ink,c,a(k+1,k),inc)

      do 30 j=k+1,n
      c=a(ip(k),j)
      a(ip(k),j)=a(k,j)
      a(k,j)=c
      
      call daxpy( n-k ,-a(k,j),  a(k+1,k), inc, a(k+1,j)  ,inc)

 30   continue

 10   continue

c      write(*,*)'PIVMIN=',pivmin

      return
      end

C =bottom========= ./mutil/blas/decomp.f
