      subroutine transf(nump,uhat,ifre_int,em,ix,uh,numint,numext)
      implicit real*8 (a-h,o-z)

c      dimension uhat(*),uh(numint),em(numint,numint) // big bug!
      dimension uhat(*),uh(*),em(numint,numint)
      dimension ifre_int(5),ix(numext+numint)
      

           
            int_lim = ifre_int(1)+ifre_int(2)
            nudof = int_lim + ifre_int(3)+ifre_int(4) + ifre_int(5)
            
            
             do ii=1,nudof
                uh(ii) = uhat(ix(ii))
                end do





            do i =1,nudof -int_lim
               do j =1,int_lim
        uhat(ix(i+int_lim)) = uhat(ix(i+int_lim))+em(i,j)*uh(j)
               end do
               end do



         
             return
             end
             
             

      

      







