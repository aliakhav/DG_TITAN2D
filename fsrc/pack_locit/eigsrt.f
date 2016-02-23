c-------------------------------------------------
      subroutine eigsrt(d,n,np)
      implicit real*8 (a-h,o-z)
      dimension d(np)
      do 13 i=1,n-1
          k=i
          p=d(i)
          do 11 j=i+1,n
              if(d(j).ge.p) then
                  k=j
                  p=d(j)
              endif
   11     continue
          if(k.ne.i) then
              d(k)=d(i)
              d(i)=p
          endif
   13 continue
      return
      end
