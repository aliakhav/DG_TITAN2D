      subroutine schur5(gid,gku,gkl,gk,ne,ned,nb,nbd,Fe,Fb,u,rhs)
c***********************************************************
c    The routine performs Schur complement operation
c      on stiffness matrix at sub-dom level.
c  On entry:
c      gid   gku |  Fe
c
c      gkl   gk  |  Fb
c  On return:
c                       -1            -1
c      Fe  = Fe - gku*gk   *Fb  Fb = gk *Fb
c                       -1                 
c     gid  = gid + gku*gk  *gkl   
c***************************************************
      implicit real*8 (a-h,o-z)
      dimension gid(ned,ned),gku(ned,nbd),gkl(nbd,ned),gk(nbd,nbd)
      dimension Fe(*),Fb(*),u(*),rhs(*)

      character A,B
c
c       write(6,*) '----in schur2----',nb,nbd
      A='N'
      B='N'
      alfa=-1.d0
      beta=1.d0


      
c factorize local stiffness
      call tri(gk,nbd,nb)

c multiple rhs solve      

      do i = 1,ne
         do j =1,nb
            rhs(j) = gkl(j,i)
            end do
      call rhsub(gk,rhs,nbd,nb,u)
         do j=1,nb
            gkl(j,i) = u(j)
            end do
      end do

      


c ...gid=gid - gku*gkl~

       call dgemm(A,B,ne,ne,nb,alfa,gku,ned,gkl,nbd,beta,gid,ned)






c-------------------------------------------------------
c solve local problem

              call rhsub(gk,Fb,nbd,nb,u)       

c
c    calculate modification of right hand side Fe~
              
              alfa = -1.d0
              
         
       call dgemv(A,ne,nb,alfa,gku,ned,u,1,beta,Fe,1)

c===============================================================          
       do i = 1,nb
         Fb(i) = u(i)
              end do
              return
              end
















