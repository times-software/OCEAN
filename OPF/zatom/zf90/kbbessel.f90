! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program kleinman
  implicit none
  !
  integer, parameter :: mr = 4200, nbm = 50, mk = 5, nrvmax = 100
  integer :: nl( mk ), nj2( mk )
  real( kind = kind( 1.0d0 ) ) :: ev( mk ), vps( mr, 0:7 ), cq( mr )
  real( kind = kind( 1.0d0 ) ) :: orb( mr, mk ), pht( mr, 3 )
  real( kind = kind( 1.0d0 ) ) :: ph0( mr ), ph1( mr ), ph2( mr )
  real( kind = kind( 1.0d0 ) ) :: ps0( mr ), ps1( mr ), ps2( mr ), vnl( mr )
  real( kind = kind( 1.0d0 ) ) :: u( 9 ), q( 9 ), z( 9, 9 ), zi( 9, 9 ), zu( 9, 9 )
  !
  integer :: ir, iw, nk, nkm1
  integer :: jrt, nr, i, klocal, nbl, nbh, nx, j, k, l, j2, jdummy
  real( kind = kind( 1.0d0 ) ) :: de, r1, r2, r0, xl, xh
  !
  integer irval, nrval, ileft, iright
  real( kind = kind( 1.0d0 ) ) :: rr, delta, wleft, wright, xlintp
  real( kind = kind( 1.0d0 ) ) :: cs, cp, cd, cf, vadd
  real( kind = kind( 1.0d0 ) ) :: vadds( nrvmax ), vaddp( nrvmax ), vaddd( nrvmax )
  real( kind = kind( 1.0d0 ) ) :: vaddf( nrvmax )
  !
  real( kind = kind( 1.0d0 ) ) :: dum
  real( kind = kind( 1.0d0 ) ), allocatable :: vlocal( : )
  !
  integer :: kk, kklocal
  integer, external :: chunk
  !
  character * 4 :: verb
  !
  ir = 13
  iw = 14
  open ( ir, file = 'sepipt', form = 'formatted', status = 'unknown' )
  open ( iw, file = 'sepopt', form = 'formatted', status = 'unknown' )
  rewind ir
  rewind iw
  read ( ir, * ) nk
  if ( nk.gt.mk ) stop 'increase mk in kleinman'
  read ( ir, * ) ( nl( k ), nj2( k ), ev( k ), k = 1, nk )
  read ( ir, * ) jrt, nr, de
  if ( nr.gt.mr ) stop 'increase mr in kleinman'
  read ( ir, * ) r1, r2
  r0 = r1 * ( ( r2 / r1 ) ** dble( jrt - 1 ) )
  do i = 1, nr
     read ( ir, * ) ( vps( i, k ), k = 1, 7 ), cq( i )
     read ( ir, * ) ( orb( i, j ), j = 1, nk )
  end do
  allocate( vlocal( nr ) )
  read ( 5, * ) klocal, nbl, nbh, xl, xh, nx, verb
  if ( nbh.gt.nbm ) stop 'increase nbm in kleinman'
  if ( klocal.gt.0 ) then
     nkm1 = nk - 1
     kklocal = chunk( nl( klocal ), nj2( klocal ) )
     vlocal( : ) = vps( : , kklocal )
  else
     nkm1 = nk
     open( unit = 99, file = 'ORB', form = 'formatted', status = 'unknown' )
     rewind 99
     do i = 1, nr
        read ( 99, * ) dum, vlocal( i )
     end do
     close( unit = 99 )
     open( 99, file = 'llocal', form = 'formatted', status = 'unknown' )
     rewind 99
     write ( 99, '(1i5)' ) -1
     close( unit = 99 )
  end if
  write ( iw, '(1x,2i8,2(1x,1f20.15))' ) nkm1, nr, r1, r2
  ! 
  ! loop over l, j channels
  ! 
  do k = 1, nk
     l = nl ( k )
     j2 = nj2( k )
     do j = 1, nr
        read ( ir, * ) jdummy, ( pht( jdummy, i ), i = 1, 3 )
     end do
     if ( k .eq. klocal ) then
        open( 99, file = 'llocal', form = 'formatted', status = 'unknown' )
        rewind 99
        write ( 99, '(2x,1i3)' ) l
        close( unit = 99 )
     end if
     if ( k.ne.klocal ) then
        call wvfcs( mr, nr, jrt, de, pht, ph0, ph1, ph2, r1, r2 )
        kk = chunk( nl( k ), nj2( k ) )
        write ( 6, * ) '>>>', k, kk, l, klocal, kklocal
        do j = 1, nr
           ! vnl( j ) = vps( j, kk ) - vps( j, kklocal )
           vnl( j ) = vps( j, kk ) - vlocal( j )
           ps0( j ) = ph0( j ) * vnl( j )
           ps1( j ) = ph1( j ) * vnl( j )
           ps2( j ) = ph2( j ) * vnl( j )
        end do
        ! 
        !
        ! this interlude figures and adds necessary, 
        ! properly j-weighted parts to the s, p, d 
        ! channel potentials
        ! 
        ! get weight coeffs for s-, p-, d-channels
        ! 
        cs = 0
        cp = 0
        cd = 0
        if ( l.eq.0 ) cs = 1
        if ( l.eq.1 ) then
           if ( j2.eq.1 ) cp = 0.33333333d0
           if ( j2.eq.2 ) cp = 1.00000000d0
           if ( j2.eq.3 ) cp = 0.66666667d0
        end if
        if ( l.eq.2 ) then
           if ( j2.eq.3 ) cd = 0.4d0
           if ( j2.eq.4 ) cd = 1
           if ( j2.eq.5 ) cd = 0.6d0
        end if
        if ( l.eq.3 )then
           if( j2.eq.5 ) cf = 3.0d0 / 7.0d0 
           if( j2.eq.6 ) cf = 1.0d0
           if( j2.eq.7 ) cf = 4.0d0 / 7.0d0
        endif

        ! for each r value, 
        ! 
        do irval = 1, nrval
           ! 
           ! set interp parameters for atom prog mesh
           !
           rr = dble( irval - 1 ) * delta
           if ( rr.lt.r1 ) then
              ileft = 4
              iright = 4
              wleft = 0.5d0
              wright = 0.5d0
           else
              xlintp = 1.0d0 + dlog( rr / r1 ) / dlog( r2 / r1 )
              ileft = xlintp
              iright = ileft + 1
              wright = xlintp - dble( ileft )
              wleft = 1.0d0 - wright
           end if
           !
           ! and add the potential according to that...
           ! we convert it here to Rydbergs
           !
           vadd = 2.0d0 * ( wleft * vnl( ileft ) + wright * vnl( iright ) )
           vadds( irval ) = vadds( irval ) + cs * vadd
           vaddp( irval ) = vaddp( irval ) + cp * vadd
           vaddd( irval ) = vaddd( irval ) + cd * vadd
           vaddf( irval ) = vaddf( irval ) + cf * vadd
        end do
        !
        !
        ! resume work
        !
        call getu( r1, r2, ph0, ph1, ph2, ps0, ps1, ps2, jrt, u )
        call getz( u, z )
        call invert( 9, 9, z, zi, zu )
        call getq( q, zi, u )
        call kbopt( l, j2, q, ps0, ps1, ps2, nr, iw )
        if ( nx.ge.0 ) then
!          call kbdiag( r1, r2, nr, jrt, r0, nx, xl, xh, nbl, nbh, k, l, &
!          & ph0, ph1, ph2, ps0, ps1, ps2, q, vnl, vps( 1, kklocal ), orb( 1, k ), verb )
           call kbdiag( r1, r2, nr, jrt, r0, nx, xl, xh, nbl, nbh, k, l, &
           & ph0, ph1, ph2, ps0, ps1, ps2, q, vnl, vlocal, orb( 1, k ), verb )
        end if
     end if
  end do
  !
  write ( 6, * ) 'kleinman terminus achieved'
end program kleinman
!
function chunk( l, j )
  implicit none
  integer :: chunk, l, j
  !
  !
  if ( l.eq.0 ) chunk = 1
  !
  !
  if ( ( l.eq.1 ).and.( j.eq.1 ) ) chunk = 2
  if ( ( l.eq.1 ).and.( j.eq.3 ) ) chunk = 3
  !
  if ( ( l.eq.1 ).and.( j.eq.2 ) ) chunk = 3
  !
  !
  if ( ( l.eq.2 ).and.( j.eq.3 ) ) chunk = 4
  if ( ( l.eq.2 ).and.( j.eq.5 ) ) chunk = 5
  !
  if ( ( l.eq.2 ).and.( j.eq.4 ) ) chunk = 5
  !
  !
  if ( ( l.eq.3 ).and.( j.eq.5 ) ) chunk = 6
  if ( ( l.eq.3 ).and.( j.eq.6 ) ) chunk = 7
  !
  if ( ( l.eq.3 ).and.( j.eq.7 ) ) chunk = 7
  !
  return
end function chunk
!
subroutine kbdiag( r1, r2, nr2, jrt, r0, nx, xlo, xhi, nbl, nbh, k, l, &
  & ph0, ph1, ph2, ps0, ps1, ps2, q, vnl, vlc, orb, verb )
  implicit none
  !
  integer :: nr2, jrt, nx, nbl, nbh, k, l
  real( kind = kind( 1.0d0 ) ) :: r1, r2, r0, xlo, xhi
  real( kind = kind( 1.0d0 ) ) :: ph0( nr2 ), ph1( nr2 ), ph2( nr2 )
  real( kind = kind( 1.0d0 ) ) :: ps0( nr2 ), ps1( nr2 ), ps2( nr2 )
  real( kind = kind( 1.0d0 ) ) :: q( 9 ), vnl( nr2 ), vlc( nr2 ), orb( nr2 )
  character * 4 :: verb
  !
  integer, parameter :: mq = 50, nr = 4200
  integer :: i, j, nbf, ix, mode, iu, ierr
  real( kind = kind( 1.0d0 ) ) :: vloc( nr ), xx
  real( kind = kind( 1.0d0 ) ) :: h( mq, mq ), s( mq, mq ), alfr( mq ), alfi( mq ), beta( mq )
  real( kind = kind( 1.0d0 ) ) :: reval( mq ), qval( mq )
  real( kind = kind( 1.0d0 ) ) :: psii( nr ), psij( nr )
  real( kind = kind( 1.0d0 ) ) :: v1( nr ), v2( nr ), v3( nr ), tmp( nr )
  real( kind = kind( 1.0d0 ) ) :: t, v, vdot1, vdot2, vdot3
  real( kind = kind( 1.0d0 ) ), external :: ekin, olap, pdot
  !
  complex( kind = kind( 1.0d0 ) ) :: ce( mq, mq )
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: ret( :, : ), err( : )
  allocate( ret( nbh, 2 ), err( nbh ) )
  err( : ) = 0
  !
  if ( nr2.gt.nr ) stop 'incr nr in kbdiag'
  if ( nbh.gt.mq ) stop 'incr mq in kbdiag'
  !
  call get( nr2, ph0, ph0, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph0, ph1, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph0, ph2, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph1, ph0, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph1, ph1, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph1, ph2, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph2, ph0, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph2, ph1, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  call get( nr2, ph2, ph2, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  do i = 1, nr2
     vloc( i ) = vlc( i ) + orb( i )
  end do
  do ix = 0, nx
     xx = xlo + ( xhi - xlo ) * dble( ix ) / dble( nx )
     call basis( nbh, r0, xx, l, qval )
     do mode = 1, 2
        do nbf = nbl, nbh
           do i = 1, nbf
              call getps( psii, qval( i ), r1, r2, l, jrt )
              call act( jrt, q, ps0, ps1, ps2, vnl, vloc, r1, r2, psii, v1, v2, v3 )
              do j = 1, nbf
                 call getps( psij, qval( j ), r1, r2, l, jrt )
                 t = ekin( qval( j ), qval( i ), r0, l )
                 vdot1 = pdot( r1, r2, psij, v1, jrt )
                 vdot2 = pdot( r1, r2, psij, v2, jrt )
                 vdot3 = pdot( r1, r2, psij, v3, jrt )
                 if ( mode.eq.1 ) v = vdot2 + vdot1
                 if ( mode.eq.2 ) v = vdot2 + vdot3
                 h( j, i ) = t + v
                 s( j, i ) = olap( qval( j ), qval( i ), r0, l )
              end do
           end do
           call rgdriver( mq, nbf, h, s, ce, ierr, alfr, alfi, beta, reval )
           if ( verb .eq. 'lang' ) then
              write ( 6, '( 20f9.4 )' ) xx, ( reval( i ), i = 1, nbf )
           end if
        end do
        ret( :, mode ) = reval( : )
        iu = 30 + 2 * ( k - 1 ) + mode
        if ( ix.eq.0 ) call gropen( iu )
        write ( iu, '( 1x, 2f18.6 )' ) ( xx, reval( i ), i = 1, nbh )
     end do
     do i = 1, nbh
        err( i ) = max( err( i ), abs( ret( i, 1 ) - ret( i, 2 ) ) )
     end do
  end do
  write ( 6, '(20f9.4)' ) err( : )
  deallocate( ret, err )
  !
  return
end subroutine kbdiag
!
subroutine getps( psi, q, r1, r2, l, jrt )
  implicit none
  integer :: l, jrt, i
  real( kind = kind( 1.0d0 ) ) :: psi( jrt ), q, r1, r2, rfac, r
  real( kind = kind( 1.0d0 ) ), external :: jlof
  rfac = r2 / r1
  r = r1
  do i = 1, jrt
     psi( i ) = jlof( q, r, l )
     r = r * rfac
  end do
  psi( jrt ) = 0.5d0 * psi( jrt )
  return
end subroutine getps
!
subroutine wvfcs( mr, nr, jrt, de, phi, ph0, ph1, ph2, r1, r2 )
  implicit none
  integer :: j, mr, nr, jrt
  real( kind = kind( 1.0d0 ) ) :: phi( mr, 3 ), ph0( mr ), ph1( mr ), ph2( mr )
  real( kind = kind( 1.0d0 ) ) :: de, f1, f2, xn, r1, r2
  real( kind = kind( 1.0d0 ) ), external :: pdot
  f1 = 1.0d0 / ( 2.0d0 * de )
  f2 = 1.0d0 / ( de * de )
  do j = 1, nr
     ph0( j ) = phi( j, 2 )
     ph1( j ) = ( phi( j, 3 ) - phi( j, 1 ) ) * f1
     ph2( j ) = ( phi( j, 3 ) + phi( j, 1 ) - 2.0d0 * phi( j, 2 ) ) * f2
  end do
  xn = 1.0d0 / sqrt( pdot( r1, r2, ph0, ph0, nr ) )
  do j = 1, nr
     ph0( j ) = xn * ph0( j )
  end do
  xn = pdot( r1, r2, ph0, ph1, nr )
  do j = 1, nr
     ph1( j ) = ph1( j ) - xn * ph0( j )
  end do
  xn = 1.0d0 / sqrt( pdot( r1, r2, ph1, ph1, nr ) )
  do j = 1, nr
     ph1( j ) = xn * ph1( j )
  end do
  xn = pdot( r1, r2, ph0, ph2, nr )
  do j = 1, nr
     ph2( j ) = ph2( j ) - xn * ph0( j )
  end do
  xn = pdot( r1, r2, ph1, ph2, nr )
  do j = 1, nr
     ph2( j ) = ph2( j ) - xn * ph1( j )
  end do
  xn = 1.0d0 / sqrt( pdot( r1, r2, ph2, ph2, nr ) )
  do j = 1, nr
     ph2( j ) = xn * ph2( j )
  end do
  return
end subroutine wvfcs
!
subroutine rgdriver( nbm, nb, hh, ss, ce, ierr, ar, ai, be, re )
  implicit none
  integer :: ierr, nbm, nb
  real( kind = kind( 1.0d0 ) ) :: hh( nbm, nbm ), ss( nbm, nbm )
  real( kind = kind( 1.0d0 ) ) :: ar( nbm ), ai( nbm ), be( nbm ), re( nbm )
  complex( kind = kind( 1.0d0 ) ) :: ce( nbm, nbm ) 
  call qzhes( nbm, nb, hh, ss, .true., ce )
  call qzit ( nbm, nb, hh, ss, 0.0d0, .true., ce, ierr )
  call qzval( nbm, nb, hh, ss, ar, ai, be, .true., ce )
  call qzvec( nbm, nb, hh, ss, ar, ai, be, ce )
  call order( nb, ar, ai, be, re )
  return
end subroutine rgdriver
!
subroutine order( nb, alfr, alfi, beta, reval )
  implicit none
  integer :: i, j, nb
  real( kind = kind( 1.0d0 ) ) :: alfr( nb ), alfi( nb ), beta( nb ), reval( nb ), swap
  do i = 1, nb
     if ( dabs( alfi( i ) / beta( i ) ).gt.0.0000000001d0 ) stop 'comp'
     reval( i ) = alfr( i ) / beta( i )
  end do
  do i = 1, nb - 1
     do j = i + 1, nb
        if ( reval( i ).gt.reval( j ) ) then
           swap = reval( i )
           reval( i ) = reval( j )
           reval( j ) = swap
        end if
     end do
  end do
  return
end subroutine order
!
subroutine mapin( w, w00, w01, w02, w11, w12, w22 )
  implicit none
  real( kind = kind( 1.0d0 ) ) :: w( 9 ), w00, w01, w02, w11, w12, w22
  w( 1 ) = w00
  w( 2 ) = w01
  w( 3 ) = w02
  w( 4 ) = w01
  w( 5 ) = w11
  w( 6 ) = w12
  w( 7 ) = w02
  w( 8 ) = w12
  w( 9 ) = w22
  return
end subroutine mapin
!
subroutine basis( nbf, r0, x, l, q )
  implicit none
  real( kind = kind( 1.0d0 ) ) :: q( 50 ), dq, r0, x, qq, xx, qa, xa, qb, xb
  real( kind = kind( 1.0d0 ) ), external :: xget
  integer :: i, l, il, nbf
  dq = 0.001d0
  if ( nbf.eq.0 ) stop 'zero nbf'
  if ( x.gt.( dble( l + 1 ) / r0 ) ) then
     qq = 0.0d0
100  qq = qq - dq
     xx = xget( qq, r0, l )
     if ( xx.lt.x ) go to 100
     qa = qq
     xa = xx
     qb = qq + dq
     xb = xget( qb, r0, l )
     call locate( q( 1 ), x, qa, qb, xa, xb, r0, l )
     il = 2
  else
     il = 1
  endif
  qq = dq / 20.0d0
  qa = qq
  xa = xget( qa, r0, l )
  do i = il, nbf
200  continue
     qq = qq + dq
     xx = xget( qq, r0, l )
     if ( ( xa.gt.x ).and.( xx.gt.xa ) ) then
        qq = ( qq + qa ) / 2.0d0
        go to 200
     end if
     if ( xx.gt.xa ) then
        qa = qq
        xa = xx
        go to 200
     end if
     if ( ( xx - x ) * ( xa - x ).gt.0.0d0 ) then
        qa = qq
        xa = xx
        go to 200
     end if
     qb = qq
     xb = xx
     call locate( q( i ), x, qa, qb, xa, xb, r0, l )
     qa = qq
     xa = xx
  end do
  return
end subroutine basis
!
subroutine locate( q, x, qa, qb, xa, xb, r0, l )
  implicit none
  integer :: l, is
  real( kind = kind( 1.0d0 ) ) :: q, x, qa, xa, qb, xb, r0, qq, xx
  logical :: ready
  real( kind = kind( 1.0d0 ) ), external :: xget
  is = 10
  do
     if ( is.gt.0 ) qq = qa + ( qb - qa ) * ( xa - x ) / ( xa - xb )
     if ( is.le.0 ) qq = 0.5d0 * ( qa + qb )
     is = is - 1
     if ( is.eq. - 10 ) is = 10
     xx = xget( qq, r0, l )
     ready = ( abs( xx - x ) .lt. 0.000001d0 )
     if ( .not. ready ) then
        if ( ( xx - x ) * ( xa - x ) .gt. 0.0d0 ) then
           qa = qq
           xa = xx
        else
           qb = qq
           xb = xx
        endif
     endif
     if ( ready ) exit
  end do
  q = qq
  return
end subroutine locate
!
function olap( q1, q2, r0, l )
  implicit none
  integer :: np, ip, l
  real( kind = kind( 1.0d0 ) ) :: olap, prec, temp, q1, q2, r0, xj, dr, r
  real( kind = kind( 1.0d0 ) ), external :: jlof
  prec = 0.000001d0
  temp = 0.0d0
  if ( dabs( q1 - q2 ).lt.prec ) then
     np = 4000
     dr = r0 / dble( np )
     do ip = 1, np
        r = r0 * ( dble( ip ) - 0.5d0 ) / dble( np )
        xj = jlof( q1, r, l )
        temp = temp + dr * xj * xj
     end do
  end if
  olap = temp
  return
end function olap
!
function ekin( q1, q2, r, l )
  implicit none
  real( kind = kind( 1.0d0 ) ) :: ekin, temp, q1, q2, olap, r
  integer :: l
  temp = 0.5d0 * q1 * q1 * olap( q1, q2, r, l )
  if ( q1.le.0.0d0 ) temp = - temp
  ekin = temp
  return
end function ekin
!
subroutine act( nr, q, ps0, ps1, ps2, vnl, vloc, r1, r2, psi, v1, v2, v3 )
  implicit none
  integer :: nr, j
  real( kind = kind( 1.0d0 ) ) :: q( 9 ), ps0( nr ), ps1( nr ), ps2( nr )
  real( kind = kind( 1.0d0 ) ) :: vnl( nr ), vloc( nr )
  real( kind = kind( 1.0d0 ) ) :: psi( nr ), v1( nr ), v2( nr ), v3( nr )
  real( kind = kind( 1.0d0 ) ) :: d0, d1, d2, r1, r2, pdot
  external pdot
  d0 = pdot( r1, r2, ps0, psi, nr )
  d1 = pdot( r1, r2, ps1, psi, nr )
  d2 = pdot( r1, r2, ps2, psi, nr )
  do j = 1, nr
     v1( j ) = d0 * ( q( 1 ) * ps0( j ) + q( 2 ) * ps1( j ) + q( 3 ) * ps2( j ) )
     v1( j ) = v1( j ) + d1 * ( q( 4 ) * ps0( j ) + q( 5 ) * ps1( j ) + q( 6 ) * ps2( j ) )
     v1( j ) = v1( j ) + d2 * ( q( 7 ) * ps0( j ) + q( 8 ) * ps1( j ) + q( 9 ) * ps2( j ) )
     v2( j ) = vloc( j ) * psi( j )
     v3( j ) = vnl ( j ) * psi( j )
  end do
  return
end subroutine act
!
subroutine get( nr, pl, pr, q, ps0, ps1, ps2, vnl, r1, r2, tmp, verb )
  implicit none
  integer :: nr, i
  real( kind = kind( 1.0d0 ) ) :: pl( nr ), pr( nr ), q( 9 )
  real( kind = kind( 1.0d0 ) ) :: ps0( nr ), ps1( nr ), ps2( nr ), vnl( nr )
  real( kind = kind( 1.0d0 ) ) :: r1, r2, pl0, pl1, pl2, pr0, pr1, pr2
  real( kind = kind( 1.0d0 ) ) :: tmp( nr ), matsl, matnl
  real( kind = kind( 1.0d0 ) ), external :: pdot
  character * 4 :: verb
  do i = 1, nr
     tmp( i ) = pr( i ) * vnl( i )
  end do
  matsl = pdot( r1, r2, pl, tmp, nr )
  pl0 = pdot( r1, r2, pl, ps0, nr )
  pl1 = pdot( r1, r2, pl, ps1, nr )
  pl2 = pdot( r1, r2, pl, ps2, nr )
  pr0 = pdot( r1, r2, pr, ps0, nr )
  pr1 = pdot( r1, r2, pr, ps1, nr )
  pr2 = pdot( r1, r2, pr, ps2, nr )
  matnl = pl0 * ( q( 1 ) * pr0 + q( 2 ) * pr1 + q( 3 ) * pr2 )
  matnl = matnl + pl1 * ( q( 4 ) * pr0 + q( 5 ) * pr1 + q( 6 ) * pr2 )
  matnl = matnl + pl2 * ( q( 7 ) * pr0 + q( 8 ) * pr1 + q( 9 ) * pr2 )
  if ( verb .eq. 'lang' ) write ( 6, '( 1x, 3( 1x, 1f20.10 ) )' ) matsl, matnl, matsl - matnl
  return
end subroutine get
!
function pdot( r1, r2, v1, v2, nr )
  implicit none
  integer :: nr, i
  real( kind = kind( 1.0d0 ) ) :: pdot, temp, r1, r2, rfac, v1( nr ), v2( nr ), rad
  rad = r1
  rfac = r2 / r1
  temp = 0
  do i = 1, nr
     temp = temp + rad * v1( i ) * v2( i )
     rad = rad * rfac
  end do
  pdot = temp * dlog( rfac )
  return
end function pdot
!
subroutine invert( nmax, ndim, smat, sinv, suse )
  implicit none
  integer :: nmax, i, j, k, ii, ndim
  real( kind = kind( 1.0d0 ) ) :: smat( nmax, nmax ), suse( nmax, nmax ), sinv( nmax, nmax )
  real( kind = kind( 1.0d0 ) ) :: rc1, rc2, swap, ratio
  do i = 1, ndim
     do j = 1, ndim
        suse( j, i ) = smat( j, i )
        sinv( j, i ) = 0.0d0
     end do
     sinv( i, i ) = 1.0d0
  end do
  do i = 1, ndim
     ii = i
     do j = i + 1, ndim
        rc1 = dabs( suse( j, i ) )
        rc2 = dabs( suse( ii, i ) )
        if ( rc1.gt.rc2 ) ii = j
     end do
     if ( ii.gt.i ) then
        do j = i, ndim
           swap = suse( i, j )
           suse( i, j ) = suse( ii, j )
           suse( ii, j ) = swap
        end do
        do j = 1, ndim
           swap = sinv( i, j )
           sinv( i, j ) = sinv( ii, j )
           sinv( ii, j ) = swap
        end do
     end if
     if ( suse( i, i ).eq.0.0d0 ) then
        write ( 6, * ) 'ZERO DETERMINANT...'
        stop
     end if
     do j = 1, ndim
        if ( j.ne.i ) then
           ratio = - suse( j, i ) / suse( i, i )
        else
           ratio = 1.0d0 / suse( i, i ) - 1.0d0
        end if
        do k = i, ndim
           suse( j, k ) = suse( j, k ) + ratio * suse( i, k )
        end do
        do k = 1, ndim
           sinv( j, k ) = sinv( j, k ) + ratio * sinv( i, k )
        end do
     end do
  end do
  return
end subroutine invert
!
subroutine getz( u, z )
  implicit none
  real( kind = kind( 1.0d0 ) ) :: u( 9 ), z( 9, 9 )
  integer :: a, b, c, d, ab, cd, ac, bd
  integer, external :: cpd
  do a = 0, 2
     do b = 0, 2
        ab = cpd( a, b )
        do c = 0, 2
           do d = 0, 2
              cd = cpd( c, d )
              ac = cpd( a, c )
              bd = cpd( b, d )
              z( ab, cd ) = u( ac ) * u( bd )
           end do
        end do
     end do
  end do
  return
end subroutine getz
!
function cpd( i, j )
  implicit none
  integer :: i, j, cpd
  cpd = 1 + i + 3 * j
  return
end function cpd
!
subroutine getu( r1, r2, ph0, ph1, ph2, ps0, ps1, ps2, nr, u )
  implicit none
  integer :: nr
  real( kind = kind( 1.0d0 ) ) :: r1, r2, u00, u01, u02, u11, u12, u22, u( 9 ), pdot
  real( kind = kind( 1.0d0 ) ) :: ph0( nr ), ph1( nr ), ph2( nr ), ps0( nr ), ps1( nr ), ps2( nr )
  external pdot
  u00 = pdot( r1, r2, ph0, ps0, nr )
  u01 = pdot( r1, r2, ph0, ps1, nr )
  u02 = pdot( r1, r2, ph0, ps2, nr )
  u11 = pdot( r1, r2, ph1, ps1, nr )
  u12 = pdot( r1, r2, ph1, ps2, nr )
  u22 = pdot( r1, r2, ph2, ps2, nr )
  call mapin( u, u00, u01, u02, u11, u12, u22 )
  return
end subroutine getu
!
subroutine getq( q, zi, u )
  implicit none
  integer :: i, j
  real( kind = kind( 1.0d0 ) ) :: q( 9 ), zi( 9, 9 ), u( 9 )
  do i = 1, 9
     q( i ) = 0
     do j = 1, 9
        q( i ) = q( i ) + zi( i, j ) * u( j )
     end do
  end do
  return
end subroutine getq
!
subroutine kbopt( l, j2, q, ps0, ps1, ps2, nr, iw )
  implicit none
  integer :: l, j2, nr, i, iw
  real( kind = kind( 1.0d0 ) ) :: q( 9 ), ps0( nr ), ps1( nr ), ps2( nr )
  write ( iw, '( 2i5 )' ) l, j2
  write ( iw, '( 1x, 3( 2x, 1f20.10 ) )' ) ( q( i ), i = 1, 9 )
  write ( iw, '( 1x, 3( 2x, 1f20.10 ) )' ) ( ps0( i ), ps1( i ), ps2( i ), i = 1, nr )
  return
end subroutine kbopt
