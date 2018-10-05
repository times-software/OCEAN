! Copyright (C) 2016, 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine spartanfip( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rcocc, lmin, lmax ) 
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: etot, rel, alfa
  real( kind = kind( 1.0d0 ) ) :: dl, zorig, xntot, rcocc( 0 : 3 )
  integer :: nr, njrc( 4 ), nel, isuse, lmin, lmax, skips( 0 : 3 ), iuflag
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  real( kind = kind( 1.0d0 ) ) :: vi( nr, 7 ), cq( nr )
  !
  integer :: l, ne, i, j, irc, idum, nq, nco, ntest, npowr
  real( kind = kind( 1.0d0 ) ) :: emax, dq, qmax, prec, rc, dum
  real( kind = kind( 1.0d0 ) ) :: aexm1( nr, lmin : lmax ), aexm2( nr, lmin : lmax ), aepot( nr, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: psxm1( nr, lmin : lmax ), psxm2( nr, lmin : lmax ), pspot( nr, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: aepl( lmin : lmax ), pspl( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: sig( lmin : lmax ), emin( lmin : lmax ), cush, kappa( lmin : lmax )
  character( len=20 ) :: fnam
  integer, parameter :: prjfile = 93, melfile = 94
  integer, allocatable :: nn( : ), ll( : ), ver( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: eint( : ), frac( : ), phc( : ), coreorb( :, : )
  character( len=20 ), allocatable :: fnl( : )
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  close( unit=99 )
  !
  open( unit=99, file='skip', form='formatted', status='unknown' )
  read ( 99, * )
  do l = 0, 3
     read ( 99, * ) i, skips( l )
  end do
  !
  ntest = 80
  allocate( eint( 0 : ntest ), frac( 0 : ntest ) )
  !
  read ( 5, * ) npowr
  read ( 5, * ) cush, emax, prec
  read ( 5, * ) rc
  read ( 5, * ) dq, qmax
!  nq = 1 + qmax / dq
  nq = 1 + ceiling( qmax / dq )
  ! Require that 'cush' (cushion) be positive
  if( cush .lt. 0 ) then
    write( 6, * ) 'Format error in .fill file. Cushion should be positive'
    write( 6, * ) '   cush     emax     prec'
    write( 6, '(A,F8.4,F9.4,F9.4)' ) '  ', cush, emax, prec
    stop
  endif
  !
  ! reset cutoff radius to radius value on the radial grid
  irc = 0
  do i = 1, nr, 4
     if ( r( i ) .lt. rc ) irc = i
  end do
  write ( fnam, '(1a8,1i3.3)' ) 'radfilez', nint( zorig )
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) r( irc ), nr, irc
  close( unit=99 )
  !
  ! prep for other stuffs
  write(6,*) njrc(:)
  write(6,*) ''
  call escanprep( nr, r, dr, r2, dl, vi, iuflag, cq, njrc, irc, lmin, lmax, sig, emin, cush, kappa, &
       aexm1, aexm2, aepot, psxm1, psxm2, pspot, aepl, pspl )
  !
  write ( fnam, '(1a8,1i3.3)' ) 'prjfilez', nint( zorig )
  open( unit=prjfile, file=fnam, form='formatted', status='unknown' )
  rewind prjfile 
  write ( prjfile, '(3i5,2x,1e15.8)' ) lmin, lmax, nq, dq
  !
  ! loop over values of valence angular momentum ... 
  do l = lmin, lmax
     !
     ! determine suitable projector functions
     call baregrip( l, lmin, lmax, nr, irc, ntest, ne, dl, rel, zorig, emax, prec, &
          r, dr, r2, aexm1, aexm2, aepot, psxm1, psxm2, pspot, aepl, pspl, skips, emin, kappa, nq, dq )
     write ( 6, * )  'projector functions found'
     write ( 6, * ) 'l, ne ... ', l, ne
     !
     write ( prjfile, '(1i5)' ) ne
     !
  end do
  close( unit=prjfile )
  !
  ! loop over all core levels, and output Slater Fk and Gk integrals and radial matrix elements
  allocate( fnl( nco + 1 + lmax - lmin ), nn( nco + 1 + lmax - lmin ), ll( nco + 1 + lmax - lmin ), phc( irc ) )
  write ( fnam, '(1a8,1i3.3)' ) 'radftabz', nint( zorig )
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nco + 1 + lmax - lmin
     read ( 99, * ) fnl( i ), idum, nn( i ), ll( i )
  end do
  close( unit=99 )
  allocate( ver( 1 : maxval( nn ), 0 : maxval( ll ) ) )
  ver( :, : ) = 0
  allocate( coreorb( irc, nco+lmax-lmin+1) )
  do i = 1, nco + 1 + lmax - lmin
     ver( nn( i ), ll( i ) ) = ver( nn( i ), ll( i ) ) + 1
!    write ( fnam, '(1a6,1i3.3,3(1a1,1i2.2))' ) '.corez', nint( zorig ), 'n', nn( i ), 'l', ll( i ), 'v', ver( nn( i ), ll( i ) )
     open( unit=99, file=fnl( i ), form='formatted', status='unknown' )
     rewind 99
     do j = 1, irc
        read ( 99, * ) idum, idum, idum, idum, idum, dum, phc( j )
     end do
     close( unit=99 )
     call getfgnew( nr, nn( i ), ll( i ), dl, zorig, r, phc )
     write ( fnam, '(1a8,1i3.3,2(1a1,1i2.2))' ) 'coreorbz', nint( zorig ), 'n', nn( i ), 'l', ll( i )
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(1a1,1i8)' ) '#', irc
     do j = 1, irc
        if ( abs( phc( j ) ) .lt. 1.d-30 ) phc( j ) = 0.0d0
        write ( 99, '(2(1x,1e22.15))' ) r( j ), phc( j )
        coreorb( j, i ) = phc( j )
     end do
     close( unit=99 )
     call getmeznl( nr, nn( i ), ll( i ), npowr, dl, zorig, r, phc )
     call core2coreme( irc, nco + 1 + lmax - lmin, i, nn, ll, npowr, dl, zorig, r, coreorb )
  end do
  deallocate( coreorb )
  !
  return
end subroutine spartanfip
