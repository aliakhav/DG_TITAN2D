c-----------------------------------------------------------------------
c     routine name   - elmcon
c   latest revision  - Feb. 14, 1994
c   purpose          - rearrange the stiffness matrix and load vector
c                      of the constraint element
c   arguments
c       in:           ifg -
c                  nedofr -
c                    nequ -
c                      Nc -
c                    Icon - Edge constraint flags
c                      ek -
c                      ef -      
c      out:        ek, ef -
c-----------------------------------------------------------------------
c   Child element numbering scheme matches the parent's
c vertex numbering scheme:
c
c ---------
c | 4 | 3 |
c |-------|
c | 1 | 2 |
c ---------
c
c-----------------------------------------------------------------------
c-----789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c
      subroutine elmcon(ifg,Norder,nedofr,nequ,Nc,Icon,ek,ef)
c
      include '../com_fem/syscom.blk'
      include '../com_fem/crrr.blk'      
c
      dimension ek(Nc,*),ef(*)
      dimension Norder(*),Icon(*)
c
      dimension Nn0(111),Nn1(10),Nn2(10)
      parameter (Nc1 = 121)
      double precision xk(Nc1,Nc1),rh(Nc1)
c
      do i=1,4
        if ( Icon(i) .ne. 0 ) goto 100
      enddo
      return
c
c  ...The case with constrained nodes
  100 nedof = nedofr/nequ
      j_free = 0
      n1 = 0
      n2 = 0
c
c ----------------------------------------------------------------------
      if ( Icon(1) .eq. -1 ) then
c
c       Side #1, Vertex #1  (parent's side #1, Child #2)
c
        if (j_free.eq.0) then              
          iposi = 2
          j_free = 2
          Nn1(1) = 1
          do ii=2,Norder(1)  
            Nn1(ii) = 4+ii-1
          enddo
          n1 = Norder(1)
        else
          iposi1 = 2
          Nn2(1) = 1
          do ii=2,Norder(1)  
            Nn2(ii) = 4+ii-1
          enddo
          n2 = Norder(1)           
        endif
c
      elseif ( Icon(1) .eq. 1 ) then
c
c       Side #1, Vertex #2  (parent's side #1, Child #1)
c
        if (j_free.eq.0) then  
          iposi = 1
          j_free = 1
          Nn1(1) = 2
          do ii=2,Norder(1)  
            Nn1(ii) = 4+ii-1
          enddo
          n1 = Norder(1)
        else
          iposi1 = 1
          Nn2(1) = 2
          do ii=2,Norder(1)  
            Nn2(ii) = 4+ii-1
          enddo
          n2 = Norder(1)
        endif
      endif
c
      if ( Icon(2) .eq. -1 ) then
c
c       Side #2, Vertex #2  (parent's side #2, Child #3)
c
        if (j_free.eq.0) then
          iposi = 2
          j_free = 3
          Nn1(1) = 2
          do ii=2,Norder(2)  
            Nn1(ii) = 3+Norder(1)+ii-1
          enddo
          n1 = Norder(2)
        else
          iposi1 = 2
          Nn2(1) = 2
          do ii=2,Norder(2)  
            Nn2(ii) = 3+Norder(1)+ii-1
          enddo
          n2 = Norder(2)
        endif    
c
      elseif ( Icon(2) .eq. 1 ) then
c
c       Side #2, Vertex #3  (parent's side #2, Child #2)
c
        if (j_free.eq.0) then
          iposi = 1
          j_free = 2
          Nn1(1) = 3          
          do ii=2,Norder(2)  
            Nn1(ii) = 3+Norder(1)+ii-1
          enddo
          n1 = Norder(2)
        else
          iposi1 = 1
          Nn2(1) = 3
          do ii=2,Norder(2)  
            Nn2(ii) = 3+Norder(1)+ii-1
          enddo
          n2 = Norder(2)              
        endif    
      endif
c
      if ( Icon(3) .eq. -1 ) then
c
c       Side #3, Vertex #3  (parent's side #3, Child #4)
c
        if (j_free.eq.0) then
          iposi = 1
          j_free = 4
          Nn1(1) = 3
          do ii=2,Norder(3)  
            Nn1(ii) = 2+Norder(1)+Norder(2)+ii-1
          enddo
          n1 = Norder(3)
        else
          iposi1 = 1
          Nn2(1) = 3
          do ii=2,Norder(3)  
            Nn2(ii) = 2+Norder(1)+Norder(2)+ii-1
          enddo
          n2 = Norder(3)
        endif
c
      elseif ( Icon(3) .eq. 1 ) then
c
c       Side #3, Vertex #4  (parent's side #3, Child #3)
c
        if (j_free.eq.0) then
          iposi = 2
          j_free = 3
          Nn1(1) = 4
          do ii=2,Norder(3)
            Nn1(ii) = 2+Norder(1)+Norder(2)+ii-1
          enddo
          n1 = Norder(3)
        else
          iposi1 = 2
          Nn2(1) = 4
          do ii=2,Norder(3)
            Nn2(ii) = 2+Norder(1)+Norder(2)+ii-1
          enddo
          n2 = Norder(3)              
        endif
      endif
c
      if ( Icon(4) .eq. -1 ) then
c
c       Side #4, Vertex #4  (parent's side #4, Child #1)
c
        if (j_free.eq.0) then
          iposi = 1
          j_free = 1
          Nn1(1) = 4
          do ii=2,Norder(4)  
            Nn1(ii) = 1+Norder(1)+Norder(2)+Norder(3)+ii-1
          enddo
          n1 = Norder(4)
        else
          iposi1 = 1
          Nn2(1) = 4
          do ii=2,Norder(4)  
            Nn2(ii) = 1+Norder(1)+Norder(2)+Norder(3)+ii-1
          enddo
          n2 = Norder(4)
        endif
c
      elseif ( Icon(4) .eq. 1 ) then
c
c       Side #4, Vertex #1  (parent's side #4, Child #4)
c
        if (j_free.eq.0) then
          iposi = 2
          j_free = 4
          Nn1(1) = 1
          do ii=2,Norder(4)  
            Nn1(ii) = 1+Norder(1)+Norder(2)+Norder(3)+ii-1
          enddo
          n1 = Norder(4)
        else
          iposi1 = 2
          Nn2(1) = 1
          do ii=2,Norder(4)  
            Nn2(ii) = 1+Norder(1)+Norder(2)+Norder(3)+ii-1
          enddo
          n2 = Norder(4)
        endif
      endif
c ----------------------------------------------------------------------
c
      n0 = 0
      do 12 ii = 1, nedof
        if (ii.eq.j_free) goto 12
        do jj =1, n1
          if (ii.eq.Nn1(jj)) goto 12
        enddo
        do jj =1, n2
          if (ii.eq.Nn2(jj)) goto 12
        enddo        
        n0 = n0+1
        Nn0(n0) = ii
   12 continue      
c
c  ...end of defining constraining information              
c
      do 20 i_r=1,nequ
      do 20 j_c=1,nequ
        call setz( xk, Nc1*nedof )
        call setz( rh, nedof )          
c
c    ...the case with one hanging node                
        if (i_r.eq.j_c) then
          rh( j_free ) = ef(nequ*(j_free-1)+i_r)
     .                 + 0.5*ef(nequ*(Nn1(1)-1)+i_r)
        endif
        xk( j_free, j_free ) = 
     .    ek(nequ*(j_free-1)+i_r, nequ*(j_free-1)+j_c)
     .    + 0.25 * ek(nequ*(Nn1(1)-1)+i_r,nequ*(Nn1(1)-1)+j_c)
     .    + 0.5 * ( ek(nequ*(j_free-1)+i_r,nequ*(Nn1(1)-1)+j_c)
     .             +ek(nequ*(Nn1(1)-1)+i_r,nequ*(j_free-1)+j_c) )
        do 110 i1 = 1, n0
          if (i_r.eq.j_c) then
            rh( Nn0(i1) ) = ef(nequ*(Nn0(i1)-1)+i_r)
          endif
          do 110 j1 = 1, n0
  110       xk( Nn0(i1), Nn0(j1) )
     .        = ek(nequ*(Nn0(i1)-1)+i_r, nequ*(Nn0(j1)-1)+j_c)          
c
        do 120 j1 = 1, n0
          xk( j_free, Nn0(j1) ) =
     .      ek(nequ*(j_free-1)+i_r, nequ*(Nn0(j1)-1)+j_c)
     .      + 0.5*ek(nequ*(Nn1(1)-1)+i_r, nequ*(Nn0(j1)-1)+j_c)
  120     xk( Nn0(j1), j_free ) =
     .      ek(nequ*(Nn0(j1)-1)+i_r, nequ*(j_free-1)+j_c)
     .      + 0.5*ek(nequ*(Nn0(j1)-1)+i_r, nequ*(Nn1(1)-1)+j_c)
c
        do 130 l1 = 1, n1
          do 130 j1 = l1, n1
            xk( j_free, Nn1(j1) ) = xk( j_free, Nn1(j1) ) +
     .        rrr(j1,l1,iposi) *
     .        (  ek(nequ*(j_free-1)+i_r,nequ*(Nn1(l1)-1)+j_c) 
     .         + 0.5*ek(nequ*(Nn1(1)-1)+i_r,nequ*(Nn1(l1)-1)+j_c))
  130       xk( Nn1(j1), j_free ) = xk( Nn1(j1), j_free ) +
     .        rrr(j1,l1,iposi) *
     .        ( ek(nequ*(Nn1(l1)-1)+i_r, nequ*(j_free-1)+j_c) 
     .         + 0.5*ek(nequ*(Nn1(l1)-1)+i_r, nequ*(Nn1(1)-1)+j_c) )
c
        do 140 l1 = 1, n1
          do 140 j1=l1, n1
            if (i_r.eq.j_c) then
              rh( Nn1(j1) ) = rh( Nn1(j1) )
     .          + rrr(j1,l1,iposi)*ef(nequ*(Nn1(l1)-1)+i_r)
            endif
            do 140 i1 = 1, n0
              xk( Nn0(i1), Nn1(j1) ) = xk( Nn0(i1), Nn1(j1) )            
     .          + rrr(j1,l1,iposi) 
     .          * ek(nequ*(Nn0(i1)-1)+i_r, nequ*(Nn1(l1)-1)+j_c)
  140         xk( Nn1(j1), Nn0(i1) ) = xk( Nn1(j1), Nn0(i1) )            
     .          + rrr(j1,l1,iposi) 
     .          * ek(nequ*(Nn1(l1)-1)+i_r, nequ*(Nn0(i1)-1)+j_c)
c
        do 150 l1 = 1, n1
        do 150  n = 1, n1
        do 150 i1 =l1, n1
        do 150 j1 = n, n1
  150     xk( Nn1(i1), Nn1(j1) ) = xk( Nn1(i1), Nn1(j1) )
     .      + rrr(i1,l1,iposi)*rrr(j1,n,iposi) 
     .      * ek(nequ*(Nn1(l1)-1)+i_r, nequ*(Nn1(n)-1)+j_c)
c
        if (n2.eq.0) goto 300
c
c    ...the case with two hanging nodes
        if (i_r.eq.j_c) then
          rh( j_free ) = rh( j_free )
     .                   + 0.5*ef(nequ*(Nn2(1)-1)+i_r)  
        endif
        xk( j_free, j_free ) = xk( j_free, j_free )
     .    + 0.5*( ek(nequ*(j_free-1)+i_r, nequ*(Nn2(1)-1)+j_c)
     .           +ek(nequ*(Nn2(1)-1)+i_r, nequ*(j_free-1)+j_c) )
     .    +0.25*( ek(nequ*(Nn1(1)-1)+i_r, nequ*(Nn2(1)-1)+j_c)
     .           +ek(nequ*(Nn2(1)-1)+i_r, nequ*(Nn1(1)-1)+j_c) )
     .    +0.25*ek(nequ*(Nn2(1)-1)+i_r, nequ*(Nn2(1)-1)+j_c)
c
        do 160 j1 = 1, n0
          xk( j_free, Nn0(j1) ) = xk( j_free, Nn0(j1) )    
     .    + 0.5*ek(nequ*(Nn2(1)-1)+i_r, nequ*(Nn0(j1)-1)+j_c)
  160     xk( Nn0(j1), j_free ) = xk( Nn0(j1), j_free )      
     .    + 0.5*ek(nequ*(Nn0(j1)-1)+i_r, nequ*(Nn2(1)-1)+j_c)
c
        do 170 l1 = 1, n1
          do 170 j1 = l1, n1
            xk( j_free, Nn1(j1) ) = xk( j_free, Nn1(j1) )
     .        + rrr(j1,l1,iposi)*0.5*
     .        ek(nequ*(Nn2(1)-1)+i_r, nequ*(Nn1(l1)-1)+j_c)
  170       xk( Nn1(j1), j_free ) = xk( Nn1(j1), j_free ) 
     .        + rrr(j1,l1,iposi)*0.5*
     .        ek(nequ*(Nn1(l1)-1)+i_r, nequ*(Nn2(1)-1)+j_c)

        do 180 l1 = 1, n2
          do 180 j1 = l1, n2
            if (i_r.eq.j_c) then
              rh( Nn2(j1) ) = rh( Nn2(j1) )
     .         + rrr(j1,l1,iposi1)*ef(nequ*(Nn2(l1)-1)+i_r)
            endif
            do 180 i1 = 1, n0
              xk( Nn0(i1), Nn2(j1) ) = xk( Nn0(i1), Nn2(j1) )
     .          + rrr(j1,l1,iposi1)
     .          * ek(nequ*(Nn0(i1)-1)+i_r, nequ*(Nn2(l1)-1)+j_c)
  180         xk( Nn2(j1), Nn0(i1)) = xk( Nn2(j1), Nn0(i1) )  
     .          + rrr(j1,l1,iposi1)
     .          * ek(nequ*(Nn2(l1)-1)+i_r, nequ*(Nn0(i1)-1)+j_c)

        do 190 l1 = 1, n2
          do 190 j1 = l1, n2
            xk( j_free, Nn2(j1) ) = xk( j_free, Nn2(j1) )
     .        + rrr(j1,l1,iposi1) *
     .        ( ek(nequ*(j_free-1)+i_r,nequ*(Nn2(l1)-1)+j_c)
     .         +0.5*ek(nequ*(Nn1(1)-1)+i_r,nequ*(Nn2(l1)-1)+j_c)
     .         +0.5*ek(nequ*(Nn2(1)-1)+i_r,nequ*(Nn2(l1)-1)+j_c)  )
  190       xk( Nn2(j1), j_free ) = xk( Nn2(j1), j_free )
     .        + rrr(j1,l1,iposi1)*
     .        ( ek(nequ*(Nn2(l1)-1)+i_r,nequ*(j_free-1)+j_c)
     .         +0.5*ek(nequ*(Nn2(l1)-1)+i_r,nequ*(Nn1(1)-1)+j_c)
     .         +0.5*ek(nequ*(Nn2(l1)-1)+i_r,nequ*(Nn2(1)-1)+j_c)  )

        do 200 l1 = 1, n2
        do 200  n = 1, n2
        do 200 i1 =l1, n2
        do 200 j1 = n, n2
 200      xk( Nn2(i1), Nn2(j1) ) = xk( Nn2(i1), Nn2(j1) )
     .      + rrr(i1,l1,iposi1)*rrr(j1,n,iposi1)
     .      * ek(nequ*(Nn2(l1)-1)+i_r,nequ*(Nn2(n)-1)+j_c)

        do 210 l1 = 1, n1
        do 210  n = 1, n2
        do 210 i1 =l1, n1
        do 210 j1 = n, n2
          xk( Nn1(i1), Nn2(j1) ) = xk( Nn1(i1), Nn2(j1) )
     .      + rrr(i1,l1,iposi)*rrr(j1,n,iposi1)
     .      * ek(nequ*(Nn1(l1)-1)+i_r,nequ*(Nn2(n)-1)+j_c)
 210      xk( Nn2(j1), Nn1(i1) ) = xk( Nn2(j1), Nn1(i1) )
     .      + rrr(j1,n,iposi1)*rrr(i1,l1,iposi)
     .      * ek(nequ*(Nn2(n)-1)+i_r,nequ*(Nn1(l1)-1)+j_c)
c
c  ...copy the constrained stiffness matrix from the workspace xk,
c     rh to ek, ef
c          
  300 do i=1,nedof
        if (i_r.eq.j_c) ef( nequ*(i-1)+i_r ) = rh(i)   
        do j=1, nedof
          ek( nequ*(i-1)+i_r, nequ*(j-1)+j_c ) = xk(i,j)
        enddo
      enddo
c
c  ...end of transformation          
   20 continue

c
      return
      end
c

