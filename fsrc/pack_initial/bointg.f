c----------------------------------------------------------------------
c    routine name   -  bointg
c   latest revision   - jan 1989
c   purpose           - this is a user routine selecting an integra-
c                       tion rule to be used when calculating boun-
c                       dary integrals
c   usage             - call bointg(norder,kind,xiloc,waloc,nint)
c   arguments 
c        in:     norder - order of approximation along a particular
c                         element side
c                kind   - boundary conditions flag
c       out:      xiloc - a vector containing the coordinates of the
c                         integration points
c                 waloc - a vector containing the weights
c                 nint  - number of integration points
c   required   routines - none
c-----------------------------------------------------------------------
c
      subroutine bointg(norder,kind,xiloc,waloc,nint)
c
      include '../com_fem/syscom.blk'
      include '../com_fem/cint.blk'
c
c.....this commons contains a flag controlling a possibility of setting
c     nint=norder e.g. in case of integrating penalty  terms  in error
c     estimation process:
      common /cboint/ iflint
c
      dimension xiloc(*),waloc(*)
c
      iflint = 0
      if ((kind.eq.2).or.(kind.eq.3).or.(kind.eq.1).or.(kind.eq.4))
     .  then
c          
c    ...gauss-lobatto underintegration (for the penalty terms):
        if(iflint.eq.0) then
c      ...while solving the problem:
          nint = norder + 1
        else
c            
c      ...while estimating the error:
          nint = norder
        endif
c
        do 10 i=1,nint
          xiloc(i) = xiloba(i,nint)
          waloc(i) = wlobat(i,nint)
   10   continue
        return
      endif
c
      if ((kind.eq.0).or.(kind.eq.5)) then
c          
c    ...gauss integration rule:
        nint = norder + 1
        do 30 i=1,nint
          xiloc(i) = xigaus(i,nint)
          waloc(i) = wagaus(i,nint)
   30   continue
      endif
c
      return
      end
c
