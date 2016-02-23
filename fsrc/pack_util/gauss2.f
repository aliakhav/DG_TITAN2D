C file:  ./mutil/blas/gauss2.f
C=========================
      subroutine gauss2(n,na,a,ip,b,x)
      include '../com_fem/syscom.blk'
      dimension a(na,n),b(n),x(n),ip(n)

      inc=1

      do 10 k=1,n-1
      c=b(ip(k))
      b(ip(k))=b(k)
      b(k)=c
 10   call daxpy(n-k,-b(k),a(k+1,k),inc,b(k+1),inc)

      do 20 j=n,1,-1
      x(j)=b(j)/a(j,j)
 20   call daxpy(j-1,-x(j),a(1,j),inc,b,inc)

      return
      end

C =bottom========= ./mutil/blas/gauss2.f
