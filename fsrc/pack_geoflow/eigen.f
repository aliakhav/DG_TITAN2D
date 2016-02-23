C***********************************************************************
      subroutine eigen(Uvec,eigenvxmax,eigenvymax,evalue,tiny,
     1     kactxy,gravity)
C***********************************************************************

      include 'rnr.h'
      double precision eigenvxmax,eigenvymax
      double precision evalue, tiny, Uvec(3), kactxy(2),gravity(3)
         
      if (Uvec(1) .gt. tiny) then
c     iverson and denlinger
         if(kactxy(1) .lt. 0.d0) then
	      write(*,*) 'negative kactxy'
  	      kactxy(1) = -kactxy(1)
	 endif
         eigenvxmax=dabs(Uvec(2)/Uvec(1)) +
     +        dsqrt(kactxy(1)*gravity(3)*Uvec(1))
         eigenvymax=dabs(Uvec(3)/Uvec(1)) +
     +        dsqrt(kactxy(1)*gravity(3)*Uvec(1))

      else
         eigenvxmax=tiny
         eigenvymax=tiny
      endif
         
      evalue=dmax1(eigenvxmax,eigenvymax)
         
      return
      end
