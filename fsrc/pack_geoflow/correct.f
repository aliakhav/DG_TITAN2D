C***********************************************************************
      subroutine correct(Uvec,Uprev,fluxxp,fluxyp, fluxxm, fluxym,
     1     tiny,dtdx,dtdy,dt,dUdx, dUdy, 
     2     curv, intfrictang, bedfrictang,g,kactxy, dgdx,
     3     frict_tiny)
C***********************************************************************

      include 'rnr.h'
      double precision fluxxp(3),fluxyp(3),tiny, Uprev(3),Ustore(3)
      double precision fluxxm(3), fluxym(3)
      double precision uvec(3), dUdx(3), dUdy(3)
      double precision sgn, h_inv, curv(2), frict_tiny
      double precision intfrictang, bedfrictang, kactxy, dgdx(2)
      double precision dtdx, dtdy, dt, g(3), sgn_dudy, sgn_dvdx, tmp
c     curv := 1/radius of curvature
      external sgn      
       
      Ustore(1)=Uprev(1)
     1     -dtdx*(fluxxp(1)-fluxxm(1))
     2     -dtdy*(fluxyp(1)-fluxym(1))
      Ustore(1)=dmax1(Ustore(1),0.d0)
      
      Ustore(2)=Uprev(2)
     1     -dtdx*(fluxxp(2)-fluxxm(2))
     2     -dtdy*(fluxyp(2)-fluxym(2))
      
      Ustore(3)=Uprev(3)
     1     -dtdx*(fluxxp(3)-fluxxm(3))
     2     -dtdy*(fluxyp(3)-fluxym(3))

      if(Uvec(1). gt. tiny) then
c     S terms
         h_inv = 1.d0/Uvec(1)
         tmp = h_inv*(dudy(2)-uvec(2)*h_inv*dudy(1))
         sgn_dudy = sgn(tmp, frict_tiny)
         Ustore(2) = Ustore(2) + dt*(g(1)*Uvec(1) - 
     1        sgn(Uvec(2),frict_tiny)*(g(3)+(Uvec(2)*h_inv)**2 *
     2        curv(1))*Uvec(1)*dtan(bedfrictang) -
     3        sgn_dudy*Uvec(1)*kactxy*(g(3)*dUdy(1)+dgdx(2)*Uvec(1))
     4        *dsin(intfrictang))

         tmp = h_inv*(dudx(3)-uvec(3)*h_inv*dudx(1))
         sgn_dvdx = sgn(tmp, frict_tiny)
         Ustore(3) = Ustore(3) + dt*(g(2)*Uvec(1) - 
     1        sgn(Uvec(3),frict_tiny)*(g(3)+(Uvec(3)*h_inv)**2 *
     2        curv(2))*Uvec(1)*dtan(bedfrictang) -
     3        sgn_dvdx*Uvec(1)*kactxy*(g(3)*dudx(1)+dgdx(1)*Uvec(1))
     4        *dsin(intfrictang))


      endif

      Uvec(1) = Ustore(1)
      Uvec(2) = Ustore(2)
      Uvec(3) = Ustore(3)
      
      return
      end
