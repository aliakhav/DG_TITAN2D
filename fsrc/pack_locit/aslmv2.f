c------------------------------------------------------------------
c     Routine Name: aslmv
c     Latest revision: Jun '93
c     Argument
c         in : p
c      
c        out : ap = A*p
c
c    
c    purpose : form block stiffness matrix and form 
c             global matrix vector product
c------------------------------------------------------------------
c
      subroutine  aslmv2(ap,p,iint,ibub,k,gi,ix,ifre_int,numext)
c

      include '../com_fem/syscom.blk'
 

      dimension p(*),ap(*),ix(*),ifre_int(5),gi(numext,numext)
c
	int_lim = ifre_int(1)+ifre_int(2)
               do ii = 1,int_lim
                     do kk=1,int_lim
                        ap(ix(ii)) = ap(ix(ii)) +gi(ii,kk)*p(ix(kk))
                        end do
                        end do
      return
      end





