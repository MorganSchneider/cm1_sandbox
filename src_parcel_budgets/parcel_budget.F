 MODULE parcel_budget_module

  implicit none

  private
  public :: parcel_budgets,setup_parcel_budget_vars,parcel_budget_write,writepbdata_nc,tri_interp2

  CONTAINS

      subroutine parcel_budgets(dt,mtime,prcltime1,xh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                               zh,mh,rmh,zf,mf,   &
                               u3d,v3d,w3d,rho,prs,xv,yv,zv,   &
                               xvt,xvs,xvb,yvt,yvs,yvb,zvt,zvs,zvb,   &
                               pdata,pbdata)
      use input
      use constants
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none


      !-----------------------------------------------------------------------

      real, intent(in) :: dt
      double precision, intent(in) :: mtime
      real, intent(in) :: prcltime1
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u3d
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v3d
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,prs,xv,yv,zv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: xvt,xvs,xvb,yvt,yvs,yvb,zvt,zvs,zvb
      real, intent(in), dimension(nparcels,npvals) :: pdata
      real, intent(inout), dimension(nparcels,npbvals) :: pbdata

      real :: x3d,y3d,z3d
      real :: xv0,yv0,zv0
      real :: drx,dry,drz,dpx,dpy,dpz
      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tem,tem1
      real :: rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: var

      logical, parameter :: debug = .false.




      !----------------------------------------------------------------------
      !  Get xvort, yvort, zvort at appropriate C-grid location:
      !  (assuming no terrain)
      !  cm1r18:  below lowest model level:
      !           Use extrapolated velocities for bbc=1,2
      !           Use log-layer equations for bbc=3 (see below)

      IF( pbxv.ge.1 )THEN
      do i=0,ni+1
      do j=1,nj+1
      do k=1,nk+1
        xv(i,j,k) = (w3d(i,j,k)-w3d(i,j-1,k))*rdy*vf(j)   &
                    -(v3d(i,j,k)-v3d(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
      enddo
      enddo
      enddo
      ENDIF

      IF( pbyv.eq.1 )THEN
      do i=1,ni+1
      do j=0,nj+1
      do k=1,nk+1
        yv(i,j,k) = (u3d(i,j,k)-u3d(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))   &
                    -(w3d(i,j,k)-w3d(i-1,j,k))*rdx*uf(i)
      enddo
      enddo
      enddo
      ENDIF

      IF( pbzv.ge.1)THEN
      do k=0,nk+1
      do j=1,nj+1
      do i=1,ni+1
        zv(i,j,k) = (v3d(i,j,k)-v3d(i-1,j,k))*rdx*uf(i)   &
                   -(u3d(i,j,k)-u3d(i,j-1,k))*rdy*vf(j)
      enddo
      enddo
      enddo
      ENDIF


      !----------------------------------------------------------------------
      !  Get xvort tendency terms at C-grid locations:
      !  (assuming no terrain)
      !  cm1r18:  below lowest model level:
      !           Use extrapolated velocities for bbc=1,2
      !           Use log-layer equations for bbc=3 (see below)


      do i=0,ni+1
      do j=0,nj+1
      do k=0,nk+1
        drx = (rho(i,j,k)-rho(i-1,j,k))*rdx*uf(i)
        dry = (rho(i,j,k)-rho(i,j-1,k))*rdy*vf(j)
        drz = (rho(i,j,k)-rho(i,j,k-1))*rdz*0.25*(mf(i,j,k)+mf(i,j-1,k)+mf(i-1,j,k)+mf(i-1,j-1,k))
        dpx = (prs(i,j,k)-prs(i-1,j,k))*rdx*uf(i)
        dpy = (prs(i,j,k)-prs(i,j-1,k))*rdy*vf(j)
        dpz = (prs(i,j,k)-prs(i,j,k-1))*rdz*0.25*(mf(i,j,k)+mf(i,j-1,k)+mf(i-1,j,k)+mf(i-1,j-1,k))


        IF( prxvten.ge.1 )THEN
        ! Tilting
        xvt(i,j,k) = ( yv(i,j,k)*(u3d(i,j,k)-u3d(i,j-1,k))*rdy*vf(j) )    &
                        + ( zv(i,j,k)*(u3d(i,j,k)-u3d(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k)) )
        ! Stretching
        xvs(i,j,k) = xv(i,j,k)*(u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)
        ! Baroclinic
        xvb(i,j,k) = (1./(rho(i,j,k)**2.))*((dry*dpz)-(drz*dpy))
        ENDIF

        IF( pryvten.ge.1 )THEN
        ! Tilting
        yvt(i,j,k) = ( xv(i,j,k)*(v3d(i,j,k)-v3d(i-1,j,k))*rdx*uf(i) )    &
                        + ( zv(i,j,k)*(v3d(i,j,k)-v3d(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k)) )
        ! Stretching
        yvs(i,j,k) = yv(i,j,k)*(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)
        ! Baroclinic
        yvb(i,j,k) = (1./(rho(i,j,k)**2.))*((drz*dpx)-(drx*dpz))
        ENDIF

        IF( przvten.ge.1 )THEN
        ! Tilting
        zvt(i,j,k) = ( xv(i,j,k)*(w3d(i,j,k)-w3d(i-1,j,k))*rdx*uf(i) )    &
                        + ( yv(i,j,k)*(w3d(i,j,k)-w3d(i,j-1,k))*rdy*vf(j) )
        ! Stretching
        zvs(i,j,k) = zv(i,j,k)*(w3d(i,j,k+1)-w3d(i,j,k))*rdz*mh(1,1,k)
        ! Baroclinic
        zvb(i,j,k) = (1./(rho(i,j,k)**2.))*((drx*dpy)-(dry*dpx))
        ENDIF
      enddo
      enddo
      enddo




        !----------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, get interpolated info:

    nploop:  &
    DO np=1,nparcels

      pbdata(np,pbtime) = mtime

      x3d = pdata(np,prx)
      y3d = pdata(np,pry)
      z3d = pdata(np,prz)


      iflag = -100
      jflag = -100
      kflag = 0

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  haveit2:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  haveit2

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. pdata(np,prx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. pdata(np,pry).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      myprcl:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                   (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        kflag = 1
        do while( pdata(np,prz).ge.zf(iflag,jflag,kflag+1) )
          kflag = kflag+1
        enddo

        x3d = pdata(np,pbx)
        y3d = pdata(np,pby)
        z3d = pdata(np,pbz)

        pbdata(np,pbx) = x3d
        pbdata(np,pby) = y3d
        pbdata(np,pbz) = z3d
        
        if( mtime.gt.prcltime1 )then
          xv0 = pbdata(np,pbxv)
          yv0 = pbdata(np,pbyv)
          zv0 = pbdata(np,pbzv)
        else
          xv0 = 0.0
          yv0 = 0.0
          zv0 = 0.0
        endif



      !----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        if( z3d.lt.zh(iflag,jflag,k) )then
          k=k-1
        endif
        rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )

        ! saveit:
        rxu = rx
        rys = ry
        rzs = rz

        !----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( z3d.lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
        rz = rzs

        ! saveit:
        rxs = rx
        ryv = ry

        !----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        rx = rxs
        ry = rys
        rz = ( z3d-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )

        ! saveit:
        rzw = rz



        !----------------------------------------------------------------------
      !  Data on xvort points

      IF( pbxv.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( pbdata(np,pbx).lt.xh(i) )then
          i=i-1
        endif

        rx = rxs
        ry = ryv
        rz = rzw

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,xv,pbdata(np,pbxv))

      ENDIF


      IF( prxvten.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( pbdata(np,pbx).lt.xh(i) )then
          i=i-1
        endif

        rx = rxs
        ry = ryv
        rz = rzw

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,xvt,pbdata(np,prxvt))
        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,xvs,pbdata(np,prxvs))
        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,xvb,pbdata(np,prxvb))

        if( mtime.gt.prcltime1 )then
          pbdata(np,prxvten) = (pbdata(np,pbxv)-xv0)*(1./prclbfrq)
          pbdata(np,prxvf) = pbdata(np,prxvten)-pbdata(np,prxvt)-pbdata(np,prxvs)-pbdata(np,prxvb)
        else
          pbdata(np,prxvten) = 0.0
          pbdata(np,prxvf) = 0.0
        endif

      ENDIF


      !----------------------------------------------------------------------
!  Data on yvort points

      IF( pbyv.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( pbdata(np,pby).lt.yh(j) )then
          j=j-1
        endif

        rx = rxu
        ry = rys
        rz = rzw

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,yv,pbdata(np,pbyv))

      ENDIF


      IF( pryvten.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( pbdata(np,pby).lt.yh(j) )then
          j=j-1
        endif

        rx = rxu
        ry = rys
        rz = rzw

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,yvt,pbdata(np,pryvt))
        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,yvs,pbdata(np,pryvs))
        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,yvb,pbdata(np,pryvb))
        if( mtime.gt.prcltime1 )then
          pbdata(np,pryvten) = (pbdata(np,pbyv)-yv0)*(1./prclbfrq)
          pbdata(np,pryvf) = pbdata(np,pryvten)-pbdata(np,pryvt)-pbdata(np,pryvs)-pbdata(np,pryvb)
        else
          pbdata(np,pryvten) = 0.0
          pbdata(np,pryvf) = 0.0
        endif

      ENDIF




      !----------------------------------------------------------------------
!  Data on zvort points

      IF( pbzv.ge.1 )THEN

        i=iflag
        j=jflag
        k=kflag

        if( pbdata(np,pbz).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = rxu
        ry = ryv
        rz = rzs

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.1 .or. i.gt.(ni+1)   .or.        &
            j.lt.1 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15559: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zvt,pbdata(np,przvt))
        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zvs,pbdata(np,przvs))
        call tri_interp2(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,zvb,pbdata(np,przvb))
        if( mtime.gt.prcltime1 )then
          pbdata(np,przvten) = (pbdata(np,pbzv)-zv0)*(1./prclbfrq)
          pbdata(np,przvf) = pbdata(np,przvten)-pbdata(np,przvt)-pbdata(np,przvs)-pbdata(np,przvb)
        else
          pbdata(np,przvten) = 0.0
          pbdata(np,przvf) = 0.0
        endif


      ENDIF



#ifdef MPI
      ELSE

        ! set to really small number (so we can use the allreduce command below)
        do n=4,npbvals
          pbdata(np,n) = -1.0e30
        enddo
#endif

      ENDIF  myprcl

    ENDDO  nploop


    !----------------------------------------------------------------------
!  communicate data
#ifdef MPI
  do np=4,npbvals
    if( myid.eq.0 )then
      call MPI_REDUCE(MPI_IN_PLACE,pbdata(1,np),nparcels,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    else
      call MPI_REDUCE( pbdata(1,np),pbdata(1,np),nparcels,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    endif
  enddo
#endif
!----------------------------------------------------------------------

      end subroutine parcel_budgets







    






        
        













!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setup_parcel_budget_vars(name_prclb,desc_prclb,unit_prclb)
      use input
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_prclb,desc_prclb,unit_prclb

      integer :: n,n2
      character(len=8) :: text1
      character(len=30) :: text2

      prclb_out = 0 !this comes from input

      prclb_out = prclb_out+1
      name_prclb(prclb_out) = 'x'
      desc_prclb(prclb_out) = 'x position'
      unit_prclb(prclb_out) = 'm'

      prclb_out = prclb_out+1
      name_prclb(prclb_out) = 'y'
      desc_prclb(prclb_out) = 'y position'
      unit_prclb(prclb_out) = 'm'

      prclb_out = prclb_out+1
      name_prclb(prclb_out) = 'z'
      desc_prclb(prclb_out) = 'z position (above sea level)'
      unit_prclb(prclb_out) = 'm'

      prclb_out = prclb_out+1
      name_prclb(prclb_out) = 'mtime'
      desc_prclb(prclb_out) = 'model time (seconds since beginning of simulation)'
      unit_prclb(prclb_out) = 's'

      if( pbxv.ge.1 )then
        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'xvort'
        desc_prclb(prclb_out) = 'zonal vorticity'
        unit_prclb(prclb_out) = '1/s'
      endif

      if( pbyv.ge.1 )then
        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'yvort'
        desc_prclb(prclb_out) = 'meridional vorticity'
        unit_prclb(prclb_out) = '1/s'
      endif

      if( pbzv.ge.1 )then
        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'zvort'
        desc_prclb(prclb_out) = 'vertical vorticity'
        unit_prclb(prclb_out) = '1/s'
      endif


      if(prxvten.ge.1)then
        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'xvten'
        desc_prclb(prclb_out) = 'total x-vorticity tendency'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'xvtilt'
        desc_prclb(prclb_out) = 'x-vorticity tendency from tilting'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'xvstretch'
        desc_prclb(prclb_out) = 'x-vorticity tendency from stretching'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'xvbcl'
        desc_prclb(prclb_out) = 'x-vorticity tendency from baroclinity'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'xvfric'
        desc_prclb(prclb_out) = 'x-vorticity tendency from viscosity'
        unit_prclb(prclb_out) = '1/s/s'
      endif

      if( pryvten.ge.1)then
        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'yvten'
        desc_prclb(prclb_out) = 'total y-vorticity tendency'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'yvtilt'
        desc_prclb(prclb_out) = 'y-vorticity tendency from tilting'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'yvstretch'
        desc_prclb(prclb_out) = 'y-vorticity tendency from stretching'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'yvbcl'
        desc_prclb(prclb_out) = 'y-vorticity tendency from baroclinity'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'yvfric'
        desc_prclb(prclb_out) = 'y-vorticity tendency from viscosity'
        unit_prclb(prclb_out) = '1/s/s'
      endif

      if(przvten.ge.1)then
        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'zvten'
        desc_prclb(prclb_out) = 'total z-vorticity tendency'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'zvtilt'
        desc_prclb(prclb_out) = 'z-vorticity tendency from tilting'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'zvstretch'
        desc_prclb(prclb_out) = 'z-vorticity tendency from stretching'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'zvbcl'
        desc_prclb(prclb_out) = 'z-vorticity from baroclinity'
        unit_prclb(prclb_out) = '1/s/s'

        prclb_out = prclb_out+1
        name_prclb(prclb_out) = 'zvfric'
        desc_prclb(prclb_out) = 'z-vorticity tendency from viscosity'
        unit_prclb(prclb_out) = '1/s/s'
      endif

!-----------------------------------------------------------------------

      end subroutine setup_parcel_budget_vars




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_budget_write(precb,rtime,name_prclb,desc_prclb,unit_prclb,pbdata,pbloc)
      use input
      implicit none

      integer, intent(inout) :: precb
      real, intent(in) :: rtime
      character(len=40), intent(in), dimension(maxvars) :: name_prclb,desc_prclb,unit_prclb
      real, intent(in), dimension(nparcels,npbvals) :: pbdata
      real, intent(inout), dimension(nparcels,3) :: pbloc

      integer :: i,n,np

!----------------------------------------------------------------------
!  write out data

    IF(myid.eq.0)THEN

#ifdef NETCDF
      IF(output_format.eq.2)THEN

        call     writepbdata_nc(precb,rtime,name_prclb,desc_prclb,unit_prclb,pbdata,pbloc(1,1))

#endif
      ENDIF
      if(dowr) write(outfile,*)

    ENDIF   ! endif for myid=0

      return
      end subroutine parcel_budget_write





!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writepbdata_nc(precb,rtime,name_prclb,desc_prclb,unit_prclb,pbdata,pbdata2)
      use input
      use constants
      use netcdf
      use writeout_nc_module, only : disp_err
      implicit none

      integer, intent(inout) :: precb
      real, intent(in) :: rtime
      character(len=40), intent(in), dimension(maxvars) :: name_prclb,desc_prclb,unit_prclb
      real, intent(in), dimension(nparcels,npbvals) :: pbdata
      real, intent(inout), dimension(nparcels) :: pbdata2

      integer :: i,ncid,status,dimid,varid,time_index,n,n2,np,timeid,tfile,xid
      integer :: npid,yhid,zhid
      logical :: allinfo

!-----------------------------------------------------------------------

      if( myid.eq.0 ) print *,'  Entering writepbdata_nc '

      do i=1,maxstring
        string(i:i) = ' '
      enddo

      string(1:22) = 'cm1out_prclvort.nc'

    allinfo = .false.
    IF(precb.eq.1) allinfo=.true.

    ncid = 1

    IF( precb.ne.1 )THEN
      ! cm1r18:  Try to open file.
      !          If error, set precb to 1 and write all info.
      status = nf90_open( path=string , mode=nf90_write , ncid=ncid )
      if( status.eq.nf90_noerr )then
        ! no error, file exists.  Get number of time levels in file:
        call disp_err( nf90_inq_dimid(ncid,'time',timeid) , .true. )
        call disp_err( nf90_inquire_dimension(ncid=ncid,dimid=timeid,len=tfile), .true. )
        if( (tfile+1).lt.precb )then
!!!          if(myid.eq.0) print *,'  tfile,precb = ',tfile,precb
          precb = tfile+1
        endif
      else
        ! if error opening file, then write all info:
        if(myid.eq.0) print *,'  status = ',status
!!!        if(myid.eq.0) print *,nf90_strerror(status)
        allinfo = .true.
        precb = 1
      endif
    ENDIF

    if( myid.eq.0 ) print *,'  pbdata precb = ',precb


  allinfo3:  IF( allinfo )THEN
    ! Definitions/descriptions:

#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create( path=string , cmode=IOR(nf90_netcdf4, nf90_classic_model) , ncid=ncid ) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

    call disp_err( nf90_def_dim(ncid,"xh",nparcels,npid) , .true. )
    call disp_err( nf90_def_dim(ncid,"yh",1,yhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"zh",1,zhid) , .true. )
    call disp_err( nf90_def_dim(ncid,"time",nf90_unlimited,timeid) , .true. )

    call disp_err( nf90_def_var(ncid,"xh",nf90_float,(/npid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","west-east location ... actually, really parcel ID number") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","degree_east") , .true. )

    call disp_err( nf90_def_var(ncid,"yh",nf90_float,(/yhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","south-north location") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","degree_north") , .true. )

    call disp_err( nf90_def_var(ncid,"zh",nf90_float,(/zhid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","height") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","m") , .true. )

    call disp_err( nf90_def_var(ncid,"time",nf90_float,(/timeid/),varid) , .true. )
    call disp_err( nf90_put_att(ncid,varid,"long_name","time") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"units","seconds") , .true. )
    call disp_err( nf90_put_att(ncid,varid,"axis","T") , .true. )

  !---------------------------------------------------------------------------!

    ! new for cm1r19:  use "_prclb" arrays, which are
    !                  defined in parcel_budget_module.F

      prclloop:  &
      do n = 1,prclb_out
        if(myid.eq.0) print *,n,trim(name_prclb(n))
        call disp_err( nf90_def_var(ncid,trim(name_prclb(n)),nf90_float,(/npid,timeid/),varid) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"long_name",trim(desc_prclb(n))) , .true. )
        call disp_err( nf90_put_att(ncid,varid,"units",trim(unit_prclb(n))) , .true. )
!!!#ifdef NCFPLUS
!!!        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/nparcels,1/)) , .true. )
!!!#endif
      enddo  prclloop

  !---------------------------------------------------------------------------!

    call disp_err( nf90_put_att(ncid,NF90_GLOBAL,'CM1 version',trim(cm1version)) , .true. )

    call disp_err( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','CF-1.7') , .true. )
    call disp_err( nf90_put_att(ncid,nf90_global,"missing_value",grads_undef) , .true. )

    call disp_err( nf90_enddef(ncid) , .true. )

  do np=1,nparcels
    call disp_err( nf90_put_var(ncid,npid,float(np),(/np/)) , .true. )
  enddo
    call disp_err( nf90_put_var(ncid,yhid,0.0) , .true. )
    call disp_err( nf90_put_var(ncid,zhid,0.0) , .true. )

!------------------------

  ENDIF  allinfo3

      ! Write data:

      time_index = precb

      call disp_err( nf90_inq_varid(ncid,'time',timeid) , .true. )
      call disp_err( nf90_put_var(ncid,timeid,rtime,(/time_index/)) , .true. )

      call disp_err( nf90_inq_varid(ncid,'x',xid) , .true. )

      DO n=1,prclb_out
        call disp_err( nf90_inq_varid(ncid,trim(name_prclb(n)),varid) , .true. )
        do np=1,nparcels
          pbdata2(np) = pbdata(np,n)
        enddo
        call disp_err( nf90_put_var(ncid,varid,pbdata2,(/1,time_index/),(/nparcels,1/)) , .true. )
      ENDDO

      ! close file

      call disp_err( nf90_close(ncid) , .true. )

      precb = precb + 1

      ! all done

      if( myid.eq.0 ) print *,'  Leaving writepbdata_nc '

      end subroutine writepbdata_nc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tri_interp2(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,pdata)
      use input
      implicit none

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      real :: pdata

      pdata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

      end subroutine tri_interp2


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE parcel_budget_module