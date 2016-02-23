
      subroutine triev(d,e,n,np)
c-------------------------------------------------
c   routine to find eigenvalue of a tri-diagnal matrix
c
c    input : d -- diagnoal of the matrix
c            e -- subdiagnal of the matrix
c    output: d -- eigenvalue
c             z -- eigenvector
c--------------------------------------------------  
      implicit real*8 (a-h,o-z)
      dimension d(np),e(np)
c  renumber element of e
      do 11 i=2,n
      e(i-1)=e(i)
 11   continue
      e(n)=0.
c
      do 15 l=1,n
      iter=0
c      
c  look for a single small sub-diag element
c   to split the matrix
c      
 1    do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if(abs(e(m))+dd.eq.dd) goto 2
   12 continue
      m=n
    2 if(m.ne.l) then
          if(iter.eq.300) return
          iter=iter+1
c  form shift
          g=(d(l+1)-d(l))/(2.*e(l))
          r=sqrt(g**2+1.)
c  this is dm-ks
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
c  a plane rotation as in the original QL
c     followed by Givens rotation to restore tri-diag
c
      do 14 i=m-1,l,-1
          f=s*e(i)
          b=c*e(i)
          if(abs(f).ge.abs(g)) then
              c=g/f
              r=sqrt(c**2+1.)
              e(i+1)=f*r
              s=1./r
              c=c*s
          else
              s=f/g
              r=sqrt(s**2+1.)
              e(i+1)=g*r
              c=1./r
              s=s*c
          endif
          g=d(i+1)-p
          r=(d(i)-g)*s+2.*c*b
          p=s*r
          d(i+1)=g+p
          g=c*r-b          
   14   continue
       d(l)=d(l)-p
       e(l)=g
       e(m)=0.
       goto 1
       endif
   15  continue
       return
       end
