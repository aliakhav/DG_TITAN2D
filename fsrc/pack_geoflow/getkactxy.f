cc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getkactxy(ndofe,nord,xi1,eta1,xnod,intfrictang,
     &     bedfrictang,uK,tiny,kactxy)
c
      implicit real*8(a-h,o-z)
      dimension xi(2),gpsi(9),dgpsi(2,9),rjac(2,2),rjacinv(2,2)
      dimension U(3),DUDX(3),DUDY(3),xnod(2,9),uk(ndofe*3),phi(9)
      double precision kactxy(2),intfrictang,bedfrictang,dphi(2,9)

      nequ=3
      xi(1)=xi1
      xi(2)=eta1
      call gshape(xi, gpsi,dgpsi)

      do 50 i=1,2
         do 50 j=1,2
            rjac(i,j)=0.d0
            do 50 k=1,9
               rjac(i,j)=rjac(i,j)+dgpsi(j,k)*xnod(i,k)
 50   continue
c
c     ...inverse jacobi matrix rjac(i,j) to get rjacinv(i,j)
      detj=rjac(1,1)*rjac(2,2)-rjac(1,2)*rjac(2,1)
      if (detj.le.0.0) then
         write(*,*) 'Non positive jacobian in getkactxy'
         stop
      endif

      call shape2dg(nord,xi1,eta1,phi)
      call dshap2dg(nord,xi1,eta1,dphi)
      rjacinv(1,1)=rjac(2,2)/detj
      rjacinv(2,2)=rjac(1,1)/detj
      rjacinv(1,2)=-rjac(1,2)/detj
      rjacinv(2,1)=-rjac(2,1)/detj 
      do in=1,3
         U(in)=0.

         DUDX(in)=0.
         DUDY(in)=0.

         do i=1,ndofe
            U(in)= U(in)+phi(i)*uK((i-1)*nequ+in)
            DUDX(in)=DUDX(in)+uK((i-1)*nequ+in)
     &           *(dphi(1,i)*rjacinv(1,1)+dphi(2,i)*rjacinv(1,2))
            DUDY(in)=DUDY(in)+uK((i-1)*nequ+in)
     &           *(dphi(1,i)*rjacinv(2,1)+dphi(2,i)*rjacinv(2,2))
         end do

c         DUDX(in)=uK(nequ+in)*rjacinv(1,1)
c         DUDY(in)=uK(2*nequ+in)*rjacinv(2,2)

      end do
      cosphi= dcos(intfrictang) 
      tandel=dtan(bedfrictang)
c     vel=(DUDX(2)+DUDY(3))
      if (U(1) .ge. tiny) then
         vel=(DUDX(2)/U(1) - U(2)*DUDX(1)/U(1)**2+
     1        DUDY(3)/U(1) - U(3)*DUDY(1)/U(1)**2)
      else
         vel =0.
      end if
      if (U(1) .gt. tiny) then
         kactxy(1)=(2.d0/cosphi**2)*(1.d0-sgn(vel,tiny)*
     1        dsqrt(1.d0-(1.d0+tandel**2)*cosphi**2) ) -1.d0
         kactxy(2)=(2.d0/cosphi**2)*(1.d0-sgn(vel,tiny)*
     1        dsqrt(1.d0-(1.d0+tandel**2)*cosphi**2) ) -1.d0
c     if there is no yielding...
         if(dabs(U(2)/U(1)) .lt. tiny .and.
     1        dabs(U(3)/U(1)) .lt. tiny) then
            kactxy(1) = 1.d0
            kactxy(2) = 1.d0
         endif
      else
         kactxy(1)=1.
         kactxy(2)=1.
      end if

      kactxy(1)=1.
      kactxy(2)=1.



      return
      end
