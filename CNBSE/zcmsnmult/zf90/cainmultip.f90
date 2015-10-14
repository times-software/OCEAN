! Copyright (C) 2013 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program cainmultip
  implicit none
  !
  integer, parameter :: stdin = 5
  !
  include 'cainmusetupnx.h'
  !
  integer :: i, j, ii, jj, zn( 3 )
  integer :: nx, ny, nz, lvl, lvh
  real( kind = kind( 1.0d0 ) ) :: inter, xi, avec( 3, 3 ), bvec( 3, 3 ), bmet( 3, 3 )
  character * 3 :: technique, req, bs, as
! character * 5 :: fnroot, muroot
  character * 5 :: eval
  character * 9 :: ct
  character * 10:: infoname
  ! 
  logical :: conduct, nspn_exist, exst
  integer :: jmin, jmax, ne, nloop, iwrk, i1, i2, need, nspn, iter, inv_loop, ispn
  real( kind = kind( 1.0d0 ) ) :: val, ar, ai, el, eh, ebase, gam0, renorm, kpref, tmp, ener, &
                               f( 2 ), gprc, gres, e_step, e_start, e_stop, core_offset
  real( kind = kind( 1.0d0 ) ) :: relative_error
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: a( : ), b( : ) !, br( : ), bi( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: hv( :, :, : ), newv( :, :, : ), oldv( :, :, : ), rex( :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: v1( : ), v2( : ), cwrk( : ), rhs( : ), bra( : ), pcdiv( : ), x( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: self_energy( :, : )
  !
  call getabb( avec, bvec, bmet )
  open( unit=99, file='xmesh.ipt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nx, ny, nz
  close( unit=99 )
  open( unit=99, file='kmesh.ipt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) zn( : )
  close( unit=99 )
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) atno, ncore, lc
  close( unit=99 )
  write ( add04, '(1a1,1i3.3)' ) 'z', atno
  write ( add10, '(1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'z', atno, 'n', ncore, 'l', lc
  !
  read ( stdin, * ) lc, lvl, lvh
  read ( stdin, * ) xi
! read ( stdin, * ) fnroot, muroot
  !
  include 'cainmusetupx.f90'
  rm1 = -1; rm1 = sqrt( rm1 )
  call soprep( nc, vms, cml, cms, lc, xi, somel )
  somel( :, :, 2 ) = -1.d0 * somel( :, :, 2 )

  !
  allocate( self_energy( n, nspn ) )
  self_energy( :, : ) = 0
  call seprep( n, nspn, e0, self_energy )
  self_energy( :, : ) = conjg( self_energy( :, : ) )

  !
  val = 0
  do ic = 1, nc
     tmp = sum( v( :, ic, 1 ) ** 2 + v( :, ic, 2 ) ** 2 )
     write ( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ic, tmp
     val = val + tmp
  end do
  val = sqrt( val )
  v = v / val
  kpref = 4.0d0 * pi * val ** 2 / ( dble( nq ) * celvol ** 2 ) ! 8 changed to 4 on July 14, 2008
  write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
  open( unit=99, file='mulfile', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1x,1e15.8)' ) kpref
  close( unit=99 )
  open( unit=99, file='nval.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nval
  close( unit=99 )
  !
  read ( stdin, * ) technique
  select case( technique ) 
  case( 'hay' )
     open( unit=99, file='mode', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) inter, jmax
     write(6, *) inter, jmax
     close( unit=99 )
     jmin = 1
!
! we now avoid crashing on May 28, 2008
!    if ( jmax + 5 .gt. n ) stop 'bad jmax'
     if ( jmax + 5 .gt. n ) jmax = n - 5
     allocate( oldv( n, nc, 2 ), newv( n, nc, 2 ) )
     allocate( a( 0 : jmax ), b( 0 : jmax ) )
     read ( stdin, * ) ne, el, eh, gam0, ebase
     oldv = 0
     j = 0
     a = 0
     b = 0
     renorm = 1.0d0
     do while( j .lt. jmax )
        call cainxact( n, nspn, nc, lc, nq, nbd, e0, self_energy, v, hv, pref, &
             iq, qphys, zn, inter, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, &
             vms, cml, cms, xi, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, &
             mhr, mhi, somel, core_offset )
        ar = 0
        ai = 0
        do ic = 1, nc
           ar = ar + dot_product( v( :, ic, 1 ), hv( :, ic, 1 ) ) + dot_product( v( :, ic, 2 ), hv( :, ic, 2 ) )
           ai = ai + dot_product( v( :, ic, 1 ), hv( :, ic, 2 ) ) - dot_product( v( :, ic, 2 ), hv( :, ic, 1 ) )
        end do
        a( j ) = ar
        call redtrid( j, a( 0 ), b( 1 ) )
        newv = hv - a( j ) * v - b( j ) * oldv
        val = 0
        do ic = 1, nc
           val = val + sum( newv( :, ic, 1 ) ** 2 + newv( :, ic, 2 ) ** 2 )
        end do
        b( j + 1 ) = sqrt( val )
        newv = newv / b( j + 1 )
        write ( 6, '(2x,2f10.6,10x,1e11.4)' ) a( j ), b( j + 1 ), ai
        j = j + 1
        if ( j .ge. jmin ) call haydump( ne, el, eh, gam0, nval, eps, j, a( 0 ), b( 1 ), kpref, ebase )
        oldv = v
        v = newv
     end do
  case( 'inv' )
     !
     open( unit=99, file='loopit', form='formatted', status='old' )
     read( 99, * ) e_start, e_stop, e_step
     close( 99 )
     inv_loop = aint( ( e_stop - e_start ) / e_step )
     if (inv_loop .lt. 1 ) inv_loop = 1
     open( unit=99, file='mode', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) inter, jmax
     write(6, *) inter, jmax
     close( unit=99 )
     !  nloop = depth of GMRES, 
     !   gres = resol. gamma, 
     !   gprc = prec. gamma, 
     ! f( 1 ) = tolerance,
     !   ener = E in eV
     open( unit=99, file='vchamp1', form='formatted', status='unknown' )
     rewind 99
     do i = 1, nc
        do j = 1, n
           write ( 99, '(2(1x,1e15.8))' ) v( j, i, : )
        end do 
     end do
     close( unit=99 )
     read ( stdin, * ) nloop, gres, gprc, f( 1 ), ener
     ener = ener / 27.2114d0
     eval = 'zerox'
     iwrk = 1
     allocate( rhs( n * nc ), bra( n * nc ), v1( n * nc ), v2( n * nc ), pcdiv( n * nc ), x( n * nc ), rex( n * nc, 2 ) )
     call vtor( n * nc, v, rhs ) 
     open( unit=99, file='rhschamp1', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(2(1x,1e15.8))' ) rhs
     close( unit=99 )
     bra( : ) = rhs( : )

    ener = e_start / 27.2114d0
    e_step = e_step / 27.2114d0
  do iter = 1, inv_loop
  write( 6, * ) ener * 27.2114d0
    
     rhs( : ) = bra( : )
     v( :, :, 1 ) = 1
     v( :, :, 2 ) = 0
!    call newxact( n, nc, lc, nq, nbd, e0, v, hv, cor, coi, pref, &
!         iq, qphys, nv, zn, 0.0d0, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, &
!         vms, cml, cms, 0.0d0, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, mhr, mhi, somel )
!     call cainxact( n, nc, lc, nq, nbd, e0, v, hv, pref, &
     call cainxact( n, nspn, nc, lc, nq, nbd, e0, self_energy, v, hv, pref, &
             iq, qphys, zn, inter, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, & 
             vms, cml, cms, xi, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, &
             mhr, mhi, somel, core_offset )

     call vtor( n * nc, hv, v1 )
     write ( 6, * ) 'got to here'
     do i = 1, n * nc
        pcdiv( i ) = ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )
     end do
     ct = 'beginning'
     req = '---'
     do while ( req .ne. 'end' )
        call invdrv( x, rhs, n * nc, i1, i2, nloop, need, iwrk, cwrk, v1, v2, bs, as, req, ct, eval, f )
        select case( req )
        case( 'all' ) ! meaning, create an array as shown
           if ( allocated( cwrk ) ) deallocate( cwrk )
           iwrk = need
           allocate( cwrk( need ) )
        case( 'act' ) ! E - S H ... in what follows, v1 must be untouched
           ! v = v1
           call rtov( n * nc, v, v1 )
!          call newxact( n, nc, lc, nq, nbd, e0, v, hv, cor, coi, pref, &
!               iq, qphys, nv, zn, inter, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, &
!               vms, cml, cms, xi, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, mhr, mhi, somel )
!           call cainxact( n, nc, lc, nq, nbd, e0, v, hv, pref, &
           call cainxact( n, nspn, nc, lc, nq, nbd, e0, self_energy, v, hv, pref, &
             iq, qphys, zn, inter, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, & 
             vms, cml, cms, xi, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, &
             mhr, mhi, somel, core_offset )

           ! v2 = hv
           call vtor( n * nc, hv, v2 )
           v2( : ) = ( ener + rm1 * gres ) * v1( : ) - v2( : )
        case( 'prc' ) ! meaning, divide by S(E-H0) ... in what follows, v1 must be untouched
           v2( : ) = v1( : ) * pcdiv( : )
           write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) i1, i2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( bra, x )
           write ( 66, '(1p,2x,3i5,5(1x,1e15.8))' ) i1, i2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( bra, x )
        end select
     end do
     write ( 66, * )
     relative_error = f( 2 ) / ( dimag( dot_product( bra, x ) ) * kpref )
     write ( 76, '(1p,1i5,4(1x,1e15.8))' ) i1, ener, ( 1.0d0 - dot_product( bra, x ) ) * kpref, relative_error
     !
     ener = ener + e_step
   enddo

     call rtov( n * nc, rex, x )
     ! write out electron-core hole amplitudes per channel.
     open( unit=99, file='echamp', form='unformatted', status='unknown' )
     rewind 99
     write ( 99 ) rex
     close( unit=99 )
     open( unit=99, file='xchamp1', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(2(1x,1e15.8))' ) x
     close( unit=99 )
     open( unit=99, file='echamp1', form='formatted', status='unknown' )
     rewind 99
     do i = 1, n * nc
        write ( 99, '(2(1x,1e15.8))' ) rex( i, : )
     end do
     close( unit=99 )
     !
  end select
  !
end program cainmultip
