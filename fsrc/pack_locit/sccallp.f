c---------------------------------------------------------------------
c    routine name     -sccall
c
c     purpose : 
c      
c    original :ys feng  mod: abani patra
c      
c   in: 
c       igk : actual array dimension  c : precond
c         n : actual no. of equs
c     purpose :  call schur complement for sub dom level Xforms
c                and write out sub dom stiif and Xform mat into files
c
c  out    u: solution      
c----------------------------------------------------------------------
      subroutine sccallp(numd,gf,pc,u,pmat,nglob,iprec,gi,gku,gkl,
     .igl_int,gk,gfi,gfk,ifre_int,ix,em,numint,numext,dummy,uu,rhs)



C      include 'mpif.h'
      implicit real*8(a-h,o-z)


      
      dimension gi(numext,numext),gku(numext,numint),gkl(numint,numext),
     . gk(numint,numint),gfi(*),gfk(*),uu(*),rhs(*)
      dimension pc(*),ifre_int(*),ix(*),gf(*),u(*)
      dimension pmat(numext,numext),em(numint,numext)
      double precision dummy(*)
      dimension igl_int(*)


c initialize



         int_lim = ifre_int(1)+ifre_int(2)
         nudof = int_lim + ifre_int(3)+ifre_int(4)+ifre_int(5)
         ned = numext
         nb = nudof - int_lim
         ibb = 1
c         write(*,*)'interface unknowns on proc',numd,int_lim,numext
c         write(*,*)'interior unknowns on proc',numd,nb,numint

		
      call schur5(gi,gku,gkl,gk,int_lim,ned,nb,numint,gfi,gfk,em,uu,rhs) 


c      write(*,*)'back from Schur'

c assemble preconditioners,rhs

c
c      do i =1, nudof
c         write(*,*)'ix(',i,')= ',ix(i)
c         end do
c
             do i =1,int_lim
                pc(ix(i)) = pc(ix(i)) + gi(i,i)
                gf(ix(i)) = gf(ix(i)) + gfi(i)
                end do

             do i =1,nudof - int_lim
              pc(ix(i+int_lim)) = pc(ix(i+int_lim)) + gk(i,i)
              u(ix(i+int_lim)) = u(ix(i+int_lim)) + gfk(i)
                end do

c
c      do i =1, nudof
c         write(*,*)numd,' diag(',ix(i),')=',pc(ix(i)),
c     .             '  gf(',ix(i),')=',gf(ix(i))
c                  
c         end do 
c

                if (iprec .eq. 1) return

                
                do i = 1,ifre_int(1)+ifre_int(2)
                   do j = 1,ifre_int(1)+ifre_int(2)
                 pmat(ix(i),ix(j)) = gi(i,j)+pmat(ix(i),ix(j))
                   end do
                      end do

cccc replace by smart assembly
                      

C234567
                call MPI_ALLREDUCE(pmat,dummy,numext*numext,
     .  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD)
                      

                      
c                call mpiared(pmat,dummy,numext*numext)

                call dcopy(numext*numext,dummy,1,pmat,1)

                   return
                   end










