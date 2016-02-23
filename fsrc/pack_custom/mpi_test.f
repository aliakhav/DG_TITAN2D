      subroutine mpi_test(iam)
      integer iam
      include 'mpif.h'
      double precision dd(2), djunk(2)


      dd(1)=5.
      dd(2)=10.
      djunk(1)= 0.
      djunk(2) = 0.
      call MPI_ALLREDUCE(dd,djunk,2,
     1     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(iam .eq. 0) then
         open(5, file='djunk.txt')
         write(5,*) djunk(1), djunk(2)
         close(5)
      endif

      return
      end

      subroutine mat_view(a, rows, columns)
      integer rows, columns, i, j, k
      double precision a(rows,columns), junk, junk2
      
      k=1
      junk2 = 0.00001
      if(k .eq. 1. and. rows.eq.columns) then
         junk = 0.000
         do i = 1, rows
            do j=i+1, columns
               if( dabs(a(i,j) - a(j,i)) .gt. dabs(a(i,j)/100) .and.
     1              dabs(a(i,j)) .gt. junk2) then
                  junk = junk+ dabs(a(i,j) - a(j,i))
c                  write(*,*) "unsymmetric mat, i, j= ",i,j
                  k = k+1
               endif
            enddo
         enddo
         junk = junk/k
         write(*,*) "junk = ",junk
      endif

      return
      end


c  zeros out small values in matrix
      subroutine mat_zero(a, rows, columns)
      integer rows, columns, i, j
      double precision a(rows,columns), junk 

      junk = 0.00001
      
      do i= 1, rows
            do j=i+1, columns
               if((dabs(a(i,j))+dabs(a(j,i))) .lt. junk ) then
                  a(i,j) = 0.
                  a(j,i) = 0.
               endif
            enddo
         enddo


      return
      end


