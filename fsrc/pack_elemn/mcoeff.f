c
      subroutine mcoeff(nequ,lcol,
     .         x,y, a11,a12,a21,a22,b1,b2,d1,d2,c,f,fx,fy,matid)
c          
      include '../com_fem/syscom.blk'
c      include '../com_fem/cflag.blk'
c
      dimension a11(lcol,*),a12(lcol,*),
     .          a21(lcol,*),a22(lcol,*),
     .           b1(lcol,*), b2(lcol,*),
     .           d1(lcol,*), d2(lcol,*),
     .            c(lcol,*),
     .            f(lcol),fx(lcol),fy(lcol)      
c

c        call poisson_coef(nequ, lcol,
c     .           x,y, a11,a12,a21,a22,b1,b2,d1,d2,c,f,fx,fy )
         call ela_coef(lcol,
     .     x,y, a11,a12,a21,a22,b1,b2,d1,d2,c,f,fx,fy,matid)

c      
      return
      end                                                             
c                              
