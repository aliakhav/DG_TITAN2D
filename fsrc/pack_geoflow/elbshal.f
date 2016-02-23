c----------------------------------------------------------------------
c
c       Routine: elbmat
c
c       Purpose: compute element vector associated with the
c		  element integrals of the bilinear form and source terms
c		
c
c	Usage: call elbmat(nel,kind,nord,ndofe,uK,r)
c
c       In:   xnod = node coordinates
c           nord = array of orders of approximation
c          ndofe = total DOFs for the element
c	      uK = vector of element unknowns at t_n
c
c       Out: r = element vector of sources - K*uprev
c           bm = matrix 
c----------------------------------------------------------------------
c
      subroutine elbshal(Nc,xnod,nord,ndofe,uK,bm,r,tiny,gravnd,
     &     curvature,bedfrictang,intfrictang,nequ,dlength,dheight,
     &     phi, dphi, artvisc, src, diffmat,nkstep,ntime_step)
c
      implicit real*8(a-h,o-z)
      integer ndofee
      include '../com_fem/cint.blk'
c
      external sgn
      double precision uK(Nc*4),r(Nc)
c     
      double precision bm(Nc,Nc)
      double precision xnod(2,9), phi_grav(4), dphi_grav(2,4)
      double precision xil(10),wal(10),DUDX(3),DUDY(3)
      double precision src(Nc),src2(Nc),xi(2)

      double precision kactx,kacty,grav(3),dgdx(2),inv_curv(2)
      double precision bedfrictang,intfrictang,gravnd(12),curvature(8)
      double precision phi(10),dphi(2,10),artvisc(Nc)
      double precision diffmat(Nc,Nc)
      double precision rjac(2,2),rjacinv(2,2),gpsi(9),dgpsi(2,9)
      double precision coef(3,3,2),U(3)
c      double precision Ul(3),Ur(3),Ut(3),Ub(3)
      double precision constant,temp
      double precision tp(10)
      integer iord_grav(5)
      constant = 0.d0
c
      cosphi= dcos(intfrictang) 
      tandel= dtan(bedfrictang)
c
cc      write(*,*) "in ELBSHAL"
cc      do i=1,Nc
cc         write(*,*) "elrhs", r(i)
cc      enddo
c
c      write(*,*) "timestep",ntime_step
c     
c     overintegration because of non-linearity
      nint= nord+2
c      nint= nord+4
      do ii=1,nint
         xil(ii)=xigaus(ii,nint)
         wal(ii)=wagaus(ii,nint)
      enddo

      do i = 1,5
         iord_grav(i) = 1
      enddo

c      if(uK(1) .gt. tiny) then
c         do i = 1,9
c            write(*,*)  uK(i)
c     enddo
c      endif

      do i = 1,10
         tp(i) = 0
      enddo
      
c     call nodcor(nel,kind,xnod)
      dx=xnod(1,3)-xnod(1,1)
      dy=xnod(2,3)-xnod(2,1)
      h1=sqrt(dx*dx+dy*dy)
      dx=xnod(1,4)-xnod(1,2)
      dy=xnod(2,4)-xnod(2,2)
      h2=sqrt(dx*dx+dy*dy)
      hk=max(h1,h2)
c     hk=2.d0/sqrt(float(nreles))
      np=nord
      pk2=float(np*np)
      hkopk2=hk/pk2
c     hkopk2=0.d0
      
      beta=0.3
c      epsilon=1.E-01*hk**(2-beta)
c      epsilon=0.001E-03*hk**(2-beta)
      epsilon = 0.d0

      do 10 i=1,Nc
         artvisc(i)=0.d0
         src(i)=0.d0
         src2(i)=0.d0
         
         do 10 j=1,Nc
            bm(i,j)=0.d0
            diffmat(i,j)=0.0d0
 10      continue
c
      do 30 ii=1,nint
         do 30 jj=1,nint

c      ...calculate Jacobi matrix dx/ds(i,j)
            xi(1)=xil(ii)
            xi(2)=xil(jj)
            call gshape(xi, gpsi,dgpsi)
            
            do 50 i=1,2
               do 50 j=1,2
                  rjac(i,j)=0.d0
                  do 50 k=1,9
                     rjac(i,j)=rjac(i,j)+dgpsi(j,k)*xnod(i,k)
 50         continue
c
c     ...inverse jacobi matrix rjac(i,j) to get rjacinv(i,j)
            detj=rjac(1,1)*rjac(2,2)-rjac(1,2)*rjac(2,1)
            if (detj.le.0.0) stop

            call shape2dg(nord,xil(ii),xil(jj),phi)

            call dshap2dg(nord,xil(ii),xil(jj),dphi)
c            write(*,*) "check in elbshal, phi",phi(1),phi(2),phi(3)
c            write(*,*) "check in elbshal, dphi ",dphi(1,1),dphi(1,2),
c     &         dphi(1,3),dphi(2,1),dphi(2,2), dphi(2,3) 
            call shape2(iord_grav, xil(ii),xil(jj),phi_grav)
            call dshap2(iord_grav, xil(ii),xil(jj),dphi_grav)
            rjacinv(1,1)=rjac(2,2)/detj
            rjacinv(2,2)=rjac(1,1)/detj
            rjacinv(1,2)=-rjac(1,2)/detj
            rjacinv(2,1)=-rjac(2,1)/detj 
            do in=1,3
               U(in)=0.
               DUDX(in)=0.
               DUDY(in)=0.
               grav(in)=0.

               do i=1,ndofe
                  U(in)=  U(in) + phi(i)*uK((i-1)*nequ+in)
              
                  DUDX(in)= DUDX(in)+uK((i-1)*nequ+in)
     &                 *(dphi(1,i)*rjacinv(1,1)+dphi(2,i)*rjacinv(1,2))
                  DUDY(in)= DUDY(in)+uK((i-1)*nequ+in)
     &                 *(dphi(1,i)*rjacinv(2,1)+dphi(2,i)*rjacinv(2,2))
               end do

c     gravity -- also goes 1-3 for x,y,z components -- interpolate from nodal values using linear basis
              
               do i=1,4
                  grav(in)=grav(in)+gravnd((i-1)*nequ+in)*phi_grav(i) 
               end do

            end do

            if(U(1) .gt. tiny) then
c               write(*,*) U(1) 
            endif

            DNORM = U(2)*U(2)+U(3)*U(3) 
            DNORM = SQRT(DNORM)+.0000000000001d0
c     derivatives of gravity,curvature

            dgdx(1)=0.
            dgdx(2)=0.
            inv_curv(1)=0.
            inv_curv(2)=0.
            
            do i=1,4
               dgdx(1)=dgdx(1)+gravnd(i*3)*(dphi_grav(1,i)*rjacinv(1,1)
     &              +dphi_grav(2,i)*rjacinv(1,2))
               dgdx(2)=dgdx(2)+gravnd(i*3)*(dphi_grav(1,i)*rjacinv(2,1)
     &              +dphi_grav(2,i)*rjacinv(2,2))
               
               inv_curv(1)=inv_curv(1)+curvature((i-1)*2+1)*phi_grav(i)
               inv_curv(2)=inv_curv(2)+curvature(i*2)*phi_grav(i)
            end do

            if (U(1).gt. tiny) then
               vel=(DUDX(2)/U(1) - U(2)*DUDX(1)/U(1)**2+
     1              DUDY(3)/U(1) - U(3)*DUDY(1)/U(1)**2)
            else
               vel=0.
            endif

c     vel=(DUDX(2)+DUDY(3))
                
            if (U(1) .gt. tiny) then
               kactx=(2.d0/cosphi**2)*(1.d0-sgn(vel,tiny)*
     1              dsqrt(1.d0-(1.d0+tandel**2)*cosphi**2) ) -1.d0
               kacty=(2.d0/cosphi**2)*(1.d0-sgn(vel,tiny)*
     1              dsqrt(1.d0-(1.d0+tandel**2)*cosphi**2) ) -1.d0

c     if there is no yielding...
               if(dabs(U(2)/U(1)) .lt. tiny .and.
     1              dabs(U(3)/U(1)) .lt. tiny) then
                  kactx=1.d0
                  kacty=1.d0
               end if
            else 
               vel = 0.d0
               kactx = 1.d0
               kacty = 1.d0
            endif

            kactx=kactx*dheight/dlength
            kacty=kacty*dheight/dlength
            
            do iii=1,nequ
               do jjj=1,nequ
                  do kkk=1,2
                     coef(iii,jjj,kkk)=0.
                  end do
               end do
            end do

            if (U(1).gt.tiny) then
               coef(1,2,1)=1.
               coef(1,3,2)=1.
               coef(2,1,1)=U(1)*grav(3)*kactx*.5d0
               coef(3,1,2)=U(1)*grav(3)*kacty*.5d0
               coef(2,2,1)=U(2)/U(1)
               coef(2,3,2)=U(2)/U(1)
               coef(3,2,1)=U(3)/U(1)
               coef(3,3,2)=U(3)/U(1)
            end if

            wght=wal(ii)*wal(jj)*detj
            
            beta1 = U(1)
            beta2 = U(2)
            
            if (U(1).gt.tiny) then
               h_inv=1.d0/U(1)
               tmp = h_inv*(DUDY(2)-U(2)*h_inv*DUDY(1))
               sgn_dudy = sgn(tmp, tiny)
               tmp = h_inv*(DUDX(3)-U(3)*h_inv*DUDX(1))
               sgn_dvdx = sgn(tmp,tiny)
               tp(1)=grav(1)*U(1)
               tp(2)=-(U(2)/DNORM)*(grav(3)+(U(2)*h_inv)**2*
     &              inv_curv(1))*U(1)*dtan(bedfrictang)     
               tp(3)=-sgn_dudy*U(1)*kactx*
     &              (grav(3)*DUDY(1)+dgdx(2)*U(1))
     &              *dsin(intfrictang)
               tp(4)=grav(2)*U(1)
               tp(5)=-(U(3)/DNORM)*(grav(3)+(U(3)*h_inv)**2*
     &              inv_curv(2))*U(1)*dtan(bedfrictang)
               tp(6)=-sgn_dvdx*U(1)*kacty*
     &              (grav(3)*DUDX(1)+dgdx(1)*U(1))
     &              *dsin(intfrictang)

            else
               sgn_dudy = 0.d0
               sgn_dvdx = 0.d0
            endif

            if (U(1).gt.tiny) then
               do 21 i=1,ndofe
                  do 20 j=1,ndofe
                     do i1=1,nequ
                        do j1=1,nequ

                       bm((i-1)*nequ+i1,(j-1)*nequ+j1)=
     &                      bm((i-1)*nequ+i1,(j-1)*nequ+j1)
     &                      +wght*phi(j)*(dphi(1,i)*rjacinv(1,1)+
     &                      dphi(2,i)*rjacinv(1,2))*coef(i1,j1,1)
     &                      + wght*phi(j)*(dphi(1,i)*rjacinv(2,1)+
     &                      dphi(2,i)*rjacinv(2,2))*coef(i1,j1,2)

ccc                       if(j.eq.1.and.1.eq.j1.and.i1.eq.1) then
ccc                          constant = constant + wght*phi(i)
ccc                       endif
                          

cc                           diffmat((i-1)*nequ+i1,(j-1)*nequ+j1)=  
cc     &                          diffmat((i-1)*nequ+i1,(j-1)*nequ+j1)
cc     &                          +wght*((dphi(1,i)*rjacinv(1,1) +
cc     &                          dphi(2,i)*rjacinv(1,2))*(dphi(1,j)*
cc     &                          rjacinv(1,1)+dphi(2,j)*rjacinv(1,2)) +
cc     &                          (dphi(1,i)*rjacinv(2,1)+dphi(2,i)*
cc     &                          rjacinv(2,2))*(dphi(1,j)*rjacinv(2,1)+
cc     &                          dphi(2,j)*rjacinv(2,2)))
                        enddo
                     enddo
 20               continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     source terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     x direction
                  i1=2                 
c     
                  src((i-1)*nequ+i1)=src((i-1)*nequ+i1)+
     &                 wght*phi(i)*(tp(1)+tp(2)+tp(3))     
c     
c     y direction
                  i1=3
c     
                  src((i-1)*nequ+i1)=src((i-1)*nequ+i1)+
     &                 wght*phi(i)*(tp(4)+tp(5)+tp(6))
                  
                  
 21            continue
            endif           
 30      continue
c     
      do i=1,ndofe*nequ
         do  j=1,ndofe*nequ
            r(i)=r(i) + bm(i,j)*uK(j)
         enddo
      enddo     
c     
      do i=1,ndofe*nequ
c         artvisc(i)=epsilon*diffmat(i,i)*uK(i)
         r(i)=r(i)+src(i)
c     r(i)=r(i)+src(i) 
c     CN just a check here
c     r(i)=r(i)+artvisc(i)
c     if (src(i).gt. 0) then
c     write(*,*)i,uk(i),r(i),src(i),diffmat(i,i),artvisc(i)
c     endif
      end do
c     write(*,*) "in ELBSHAL(end)"
c      do i=1,Nc
c        write(*,*) "elrhs", r(i)
c      enddo
      return
      end
