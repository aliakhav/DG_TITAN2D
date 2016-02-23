C***********************************************************************
      double precision function sgn(zz, tiny)
C***********************************************************************

      double precision zz, tiny
      
      if (dabs(zz) .lt. tiny) then
         sgn=0.d0
      else if (zz .ge. tiny) then
         sgn=1.d0
      else
         sgn=-1.d0
      endif
      
      return
      end
