c---------------------------------------------------------------------
c    routine name     -pcg_sol.f
c
c     purpose : block wise cg solver
c      
c      
c   in: gi  u = gfi
c	ix : mapping array 
c       numext:  sub-domain interface unknowns(dimension of gi)
c	igk: array size for globals(amount of vel unknowns)
c	prec:  diagonal preconditioner
c	rm,rmnew,zm,zmnew,pm: vectors used in PCG as in TJR Hughes "FEM"p.485
c       denom:  gi*pm vector in TJR Hughes "FEM"p.485
c       dummy: communication buffer
c       nump:  # of processors
c       iam:  processor number
c       (C^-1 is prec)
c  out    u, pres_sub: solution  
c  note:  bcb, store_p, rm, rmnew, zm, zmnew, pm, prec are 
c            all in global ordering    
c----------------------------------------------------------------------
	subroutine pcg_sol(gi,gfi,u, numext,prec,rm,rmnew,zm,
     1       zmnew,pm,denom,dummy,nump,iam)
c
        include 'mpif.h'	
	implicit double precision (a-h,o-z)	
	dimension gi(numext,numext),gfi(numext),u(numext)
	dimension prec(numext),pm(numext),denom(numext)
	dimension rm(numext),rmnew(numext),zm(numext),zmnew(numext)
	dimension dummy(numext)
	integer n,  nump, iam

c  convergence criteria (epsil and n are convergence criteria)
	epsil=0.5e-8
	resid_min = 1.e+100
	n = numext*4

	write(*,*) 'top of pcg_sol, proc: ',iam
	do i=1, numext
	   if(dabs(gi(i,i)) .gt. .01) then
	      prec(i) = 1.0/gi(i,i)
	   else
	      prec(i) = 1.
	   endif
	   rm(i) = gfi(i)
	   u(i) = 0.
	enddo

cc test stuff!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	do i=1,numext
c	   prec(i) = 1.
c	enddo
c	do i=1, numext
c	   do j=1, numext
c	      gi(i,j) = 0.
c	   enddo
c	   gi(i,i) = 100.
c	enddo
c  end test stuff------------------------------------	

	do i=1, numext
	   zm(i) = rm(i)*prec(i)
	   pm(i) = zm(i)
	enddo
c  k is iteration counter
	k=0
c-------------------------------------start of pcg iteration
 15	alpham=0.0

	k=k+1
c  dot product of rm and zm
	alpham=ddot(numext, rm, 1, zm, 1)
	
c	alpham = 0.
c	do i=1, numext
c	   alpham = alpham + rm(i)*zm(i)
c	enddo

	do i=1, numext
	   denom(i) = 0.
	enddo
	do i=1,numext
	   do j=1, numext
	      denom(i) = denom(i) + gi(i,j)*pm(j)
	   enddo
	enddo
c	call aslmv2(denom,pm,gi,ix,numext)
	

	sum = 0.
c	do i=1,numext
c	   sum = sum + denom(i)*pm(i)
c	enddo
	
	sum = ddot(numext,denom, 1, pm, 1)
	
	alpham = alpham/sum
	if(iam .eq. 0) then
	   write(*,*) 'alpham = ',alpham
	endif
	do i=1,numext
	   u(i)=u(i)+alpham*pm(i)
	   rmnew(i) = rm(i) - alpham*denom(i)
	enddo

c  convergence check using L2 vector norm for rm
	crit = 0.
	do i=1, numext
	   crit = crit + rmnew(i)*rmnew(i)
	enddo
	crit = sqrt(crit)
c  convergence check using L1 vector norm for rm
	crit = 0.
	do i=1, numext
	   crit = crit + dabs(rmnew(i))
	enddo

	if(iam.eq. 0) then    
	   write(*,*) 'crit = ',crit, ' iter=', k
	endif

	if(resid_min .gt.  crit) then
	   resid_min = crit
	   iter_min = k
	endif
	if ((crit.le.epsil) .or. (k.eq.n) .or. (alpham.eq.0.d0) ) then
	   if(iam.eq. 0) then
	      write(*,*) 'solution criteria met, crit = ',crit, ' iter=', k
	      write(*,*) 'minimum residual was: ',resid_min,' iter=',iter_min
	   endif
	   return
	end if


c  find new velocity 
	do i=1, numext
	   zmnew(i) = prec(i)*rmnew(i)
	enddo

	top = ddot(numext, rmnew, 1, zmnew, 1)
	bottom = ddot(numext, rm, 1, zm, 1)

	betam = top/bottom

	do i=1,numext
	   pm(i) = zmnew(i) + betam*pm(i)
	enddo   

	do i=1, numext
	   rm(i) = rmnew(i)
	   zm(i) = zmnew(i)
	enddo

	goto 15
	end











