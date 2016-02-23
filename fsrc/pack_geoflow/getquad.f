      subroutine getquad(nint,xil,wal)
      implicit real*8(a-h,o-z)
      integer nint
      include '../com_fem/cint.blk'
      dimension xil(10),wal(10)

c      call cnst
	do ii=1,nint
	   xil(ii)=xigaus(ii,nint)
	   wal(ii)=wagaus(ii,nint)
	enddo
c
        return
        end

      
