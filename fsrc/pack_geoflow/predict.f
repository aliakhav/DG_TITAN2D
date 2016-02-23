C***********************************************************************
      subroutine predict(Uvec,dUdx,dUdy,Uprev,tiny,kactxy,dt2, g, 
     1     inv_curv, bedfrictang, intfrictang,
     2     dgdx, frict_tiny)
C***********************************************************************
c     this routine only calculates Uvec for a half step

***********************************************************************
******NOTE:  d(g(3)*Uvec(1))/dx is approximated by g(3)*dUvec(1)/dx !!!
***********************************************************************
      include 'rnr.h'
      double precision dUdx(3),dUdy(3),kactxy,Uvec(3)
      double precision tiny,dt2, Uprev(3), dgdx(2)
      double precision c_sq, h_inv, g(3),inv_curv(2)
      double precision intfrictang, bedfrictang, sgn, frict_tiny
      double precision sgn_dudy, sgn_dvdx, tmp
      external sgn

c     inv_curv := inverse of curvature
      
      uprev(1) = uvec(1)
      uprev(2) = uvec(2)
      uprev(3) = uvec(3)
c      c_sq = kactxy*g(3)*Uvec(1)
c     h_inv := 1/Uvec(1)                         

c      Uvec(1)=Uvec(1) -dt2*(dUdx(2)+dUdy(3))
c      Uvec(1)=dmax1(Uvec(1),0.d0)
      
c     dF/dU, dG/dU and S terms if Uvec(1) > TINY !
c      if(Uprev(1) .gt. TINY) then
c         h_inv = 1.d0/Uprev(1)
c     ****** X-dir ******
c     dF/dU and dG/dU terms
c         Uvec(2) = Uvec(2) - dt2*((c_sq-(Uprev(2)*h_inv)**2)*dUdx(1) + 
c     1        2.d0*Uprev(2)*h_inv*dUdx(2) -
c     2        Uprev(2)*Uprev(3)*h_inv**2 *dUdy(1)+Uprev(3)*h_inv*dUdy(2)
c     3        + Uprev(2)*h_inv*dUdy(3)) 
c     S terms
c         tmp = h_inv*(dudy(2)-uprev(2)*h_inv*dudy(1))
c         sgn_dudy = sgn(tmp, frict_tiny)
c         Uvec(2) = Uvec(2) + dt2*(g(1)*Uprev(1) - 
c     1        sgn(Uprev(2),frict_tiny)*(g(3)*Uprev(1)+
c     2        Uprev(2)*Uprev(2)*h_inv *
c     2        inv_curv(1))*dtan(bedfrictang) -
c     3        sgn_dudy*Uprev(1)*kactxy*(g(3)*dUdy(1)+dgdx(2)*Uprev(1))*
c     4        dsin(intfrictang))
cc     ****** Y-dir ******
c     dF/dU and dG/dU terms
c         Uvec(3) = Uvec(3) - dt2*((c_sq-(Uprev(3)*h_inv)**2)*dUdy(1) +
c     1        2.d0*Uprev(3)*h_inv*dUdy(3) -
c     2        Uprev(2)*Uprev(3)*h_inv**2 *dUdx(1)+Uprev(3)*h_inv*dUdx(2)
c     3        + Uprev(2)*h_inv*dUdx(3)) 
c     S terms
c         tmp = h_inv*(dudx(3)-uprev(3)*h_inv*dudx(1))
c         sgn_dvdx = sgn(tmp, frict_tiny)
c         Uvec(3) = Uvec(3) + dt2*(g(2)*Uprev(1) - 
c     1        sgn(Uprev(3),frict_tiny)*(g(3)*Uprev(1)+
c     2       Uprev(3)*Uprev(3)*h_inv *
c     2        inv_curv(2))*dtan(bedfrictang) -
c     3        sgn_dvdx*Uprev(1)*kactxy*(g(3)*dUdx(1)+dgdx(1)*Uprev(1))
c     4        *dsin(intfrictang))        
c      endif
            
      return
      end
