c--------------------------------------------------------------------
c     Routine Name      - bcelem
c     latest revision   - March, 1993
c     purpose           - routine determines the boundary integrals
c                         supplementing the stiffness matrix and Rhs vec
c                         according to the type of boundary conditions
c                         required
c     usage             - call bcelem(Nel,Idg,Ndof,xy, nc,xk,Rhs)
c     arguments:
c     input              Nel,kind  - element number and kind
c                        Norder    - order of element
c                        ifg -  no use at present 
c                        nc     - decleared length of the column of xk
c                        xk,Rhs - element stiffness matrix and Rhs vector
c                                 before adding boundary condition
c      output            xk,Rhs - element stiffness matrix and Rhs vector
c                                 after adding boundary condition       
c-------------------------------------------------------------------
c      
      subroutine bcelem(ifg,Norder,nedof,nequ,Nc,Nelb,xy,xk,Rhs,value)
c
      include '../com_fem/syscom.blk'
      include '../com_fem/paramf.blk'
c      include '../com_fem/cdoff.blk'
c      include '../com_fem/celemf.blk'
      include '../com_fem/cint.blk'
c
      dimension Norder(*),xk(Nc,*),Rhs(*), value(*)
      dimension xy(2,9),Xnod(2,3),Vshap(11),xiloc(10),waloc(10)
      dimension f(11*MAXEQS,11*MAXEQS),frh(11*MAXEQS)
      dimension in(11)
c
      dimension att(MAXEQS,MAXEQS),
     .           bt(MAXEQS,MAXEQS),
     .           pt(MAXEQS,MAXEQS),
     .            c(MAXEQS,MAXEQS),
     .           gt(MAXEQS       ),
     .            g(MAXEQS       ) 
c
c      dimension Neig(2,4),Nelb(4),Nside(4)
      dimension Nelb(4),Nside(4)
c
      x(i,eta) = Xnod(i,1)*eta*(eta-1.d0)/2.d0
     .         + Xnod(i,2)*(1.d0-eta*eta)
     .         + Xnod(i,3)*eta*(eta+1.d0)/2.d0
c
      xeta(i,eta) = Xnod(i,1)*(eta-0.5d0)
     .            - Xnod(i,2)*2.d0*eta
     .            + Xnod(i,3)*(eta+0.5d0)
c
      Nelb_con = 100
c
c  ...find types of boundary conditions
c
c      call neigbr(Nel,Kind, Neig,Nside,Nelb)
c      call nodcor(Nel,Kind, xy)
c
c  ...loop through sides of the element
c
      do 20 is=1,4
        if(Nelb(is).eq.0) go to 20
c        call findap(Nel,Kind, Norder)
        if (is.le.2) then 
          Xnod(1,1) = xy(1,is)
          Xnod(2,1) = xy(2,is)
          Xnod(1,3) = xy(1,is+1)
          Xnod(2,3) = xy(2,is+1)
        elseif (is.eq.3) then
          Xnod(1,1) = xy(1,4)
          Xnod(2,1) = xy(2,4)
          Xnod(1,3) = xy(1,3)
          Xnod(2,3) = xy(2,3)
        else
          Xnod(1,1) = xy(1,1)
          Xnod(2,1) = xy(2,1)
          Xnod(1,3) = xy(1,4)
          Xnod(2,3) = xy(2,4)            
        endif
        Xnod(1,2) = xy(1,is+4)
        Xnod(2,2) = xy(2,is+4)
c
c  ...Initiating f and frh  
c        
        do i = 1, nequ*(Norder(is)+1)
          call setz(f(1,i),nequ*(Norder(is)+1))
        enddo
        call setz(frh,nequ*(Norder(is)+1))
c
        call bointg(Norder(is)+1,Nelb(is), xiloc,waloc,igs) 
c
c  ...Natural boundary condition and normal boundary condition      
c  ...Specify the integral rule
c
        do 30 i=1,igs
          xglobl = x(1,xiloc(i))
          yglobl = x(2,xiloc(i))
          rtx = xeta(1,xiloc(i))
          rty = xeta(2,xiloc(i))
          ds = sqrt(rtx*rtx+rty*rty)
c
c      ...make sure the tangential direction is always count-clockwise
          if (is.ge.3) then
            rtx = -rtx
            rty = -rty
          endif
c
          rtx = rtx/ds
          rty = rty/ds
          rnx = rty
          rny = -rtx
c          
          call bcoeff(MAXEQS,nelb(is),is,xglobl,yglobl,xiloc(i),rnx,
     .         rny,rtx,rty,att,bt,pt,c,gt,g,value(2*is-1),value(2*is))          
c
          call shape1(Norder(is),xiloc(i),vshap)
          do 40 k=1,Norder(is)+1
          do 40 ivar1 = 1, nequ
            frh(nequ*(k-1)+ivar1) = frh(nequ*(k-1)+ivar1)
     .                            + vshap(k)*g(ivar1)*ds*waloc(i)
            do 40 l=1,Norder(is)+1
            do 40 ivar2 = 1, nequ
              f(nequ*(k-1)+ivar1,nequ*(l-1)+ivar2) =
     .        f(nequ*(k-1)+ivar1,nequ*(l-1)+ivar2)
     .      + vshap(k)*vshap(l)*c(ivar1,ivar2)*ds*waloc(i)
   40     continue
   30   continue
c
        if (is.le.2) then
          in(1) = is-1
          in(2) = is
        elseif (is.eq.3) then
          in(1) = 3
          in(2) = 2
        else
          in(1) = 0
          in(2) = 3
        endif
c        
        iprev = 3
        do i=1,is-1
          iprev = iprev + Norder(i)-1
        enddo
        do i=1,Norder(is)-1
          in(2+i) = iprev+i
        enddo        
c
c  ...add f and frh to element stiffness matrix
c          
        do 50 i=1,Norder(is)+1
        do 50 ivar1 = 1, nequ
          Rhs(nequ*in(i)+ivar1) = Rhs(nequ*in(i)+ivar1)
     .                       + frh(nequ*(i-1)+ivar1)
          do 60 j=1,Norder(is)+1
          do 60 ivar2 = 1, nequ
            xk(nequ*in(i)+ivar1,nequ*in(j)+ivar2) = 
     .      xk(nequ*in(i)+ivar1,nequ*in(j)+ivar2)
     .        + f(nequ*(i-1)+ivar1,nequ*(j-1)+ivar2)
   60     continue
   50   continue  
   20 continue
c
      return
      end
c

      
