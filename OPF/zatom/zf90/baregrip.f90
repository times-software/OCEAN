! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine baregrip( l, lmin, lmax, nr, irc, ntest, north, dl, rel, zorig, emax, prec, &
     r, dr, r2, aexm1, aexm2, aepot, psxm1, psxm2, pspot, aepl, pspl, skips, emin, kappa, nq, dq )
  implicit none
  !
  integer :: l, lmin, lmax, nr, irc, ntest, north
  integer :: skips( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: dl, rel, zorig, emax, prec
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  real( kind = kind( 1.0d0 ) ) :: aexm1( nr, lmin : lmax ), aexm2( nr, lmin : lmax ), aepot( nr, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: psxm1( nr, lmin : lmax ), psxm2( nr, lmin : lmax ), pspot( nr, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: aepl( lmin : lmax ), pspl( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: emin( lmin : lmax ), kappa( lmin : lmax )
  integer, intent( in ) :: nq
  real( kind = kind( 1.0d0 ) ), intent( in ) :: dq

  !
  integer :: i, j, nnode, ief, nnae, nnps
  real( kind = kind( 1.0d0 ) ) :: x
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: et, phi, angae, angps, slpps
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: phiae, phips, pspr, aepr, phirn
  character * 5 :: fnam
  real( kind = kind( 1.0d0 ) ), parameter :: zpseu = 0.0d0, nrl = 0.0d0
  !
  integer :: j1, j2, ii, imin, imax
  real( kind = kind( 1.0d0 ) ) :: e1, e2, a1, a2, ee, abest, slp, slpae, sca, ang1, der0, der1
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: phi1, phibest
  logical :: qual
  !
  allocate( phiae( irc, ntest ), phips( irc, ntest ), pspr( irc, ntest ), aepr( irc, ntest ), phirn( irc, ntest ) )
  allocate( angae( ntest ), angps( ntest ), et( ntest ), phi( nr ), slpps( ntest ) )
  allocate( phi1( irc ), phibest( irc ) )
  write ( fnam, '(1a4,1i1)' ) 'angs', l
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, ntest
     et( i ) = emin( l ) + ( emax - emin( l ) ) * dble( i - 1 ) / dble( ntest - 1 )
     !
     call intego( et( i ), l , kappa( l ), 1000, nnode, irc, ief, x, phi, &
          zorig, aepot( :, l ), aexm1( :, l ), aexm2( :, l ), nr, r, dr, r2, dl, rel, aepl( l ), der0, der1 )
     call normangnodes( irc, phi, phiae( :, i ), nnae, dl, r, 0, angae( i ), slp )
     !
     call intego( et( i ), l , kappa( l ), 1000, nnode, irc, ief, x, phi, &
          zpseu, pspot( :, l ), psxm1( :, l ), psxm2( :, l ), nr, r, dr, r2, dl, nrl, pspl( l ), der0, der1 )
     call normangnodes( irc, phi, phips( :, i ), nnps, dl, r, skips( l ), angps( i ), slpps( i ) )
     write ( 99, '(3(1x,1e15.8),3i5)' ) et( i ), angae( i ), angps( i ), nnae, nnps, skips( l )
  end do
  close( unit=99 )
  !
  do i = 1, ntest
     if ( angps( i ) .gt. angae( ntest ) ) imax = i
  end do
  do i = ntest, 1, -1
     if ( angps( i ) .lt. angae( 1 ) ) imin = i
  end do
  do i = imin, imax
     do j = 1, ntest
        if ( angae( j ) .gt. angps( i ) ) j1 = j
     end do
     e1 = et( j1 ); a1 = angae( j1 )
     do j = ntest, 1, -1
        if ( angae( j ) .lt. angps( i ) ) j2 = j
     end do
     e2 = et( j2 ); a2 = angae( j2 )
     do ii = 1, 5
        ee = ( e1 * ( a2 - angps( i ) ) + e2 * ( angps( i ) - a1 ) ) / ( a2 - a1 )
        call intego( ee, l , kappa( l ), 1000, nnode, irc, ief, x, phi, &
            zorig, aepot( :, l ), aexm1( :, l ), aexm2( :, l ), nr, r, dr, r2, dl, rel, aepl( l ), der0, der1 )
        call normangnodes( irc, phi, phi1, nnae, dl, r, 0, ang1, slp )
        if ( ii .eq. 1 ) then
           qual = .true.
        else
           qual = ( abs( ang1 - angps( i ) ) .lt. abs( abest - angps( i ) ) )
        end if
        if ( qual ) then
           phibest( : ) = phi1
           abest = ang1
           slpae = slp
        end if
        !       write ( 6, '(3i5,10(1x,1e15.8))' ) i, imin, imax, e1, e2, ee, a1, a2, ang1, abest, angps( i )
        if ( ( ang1 - angps( i ) ) * ( angps( i ) - a2 ) .gt. 0.0d0 ) then
           a1 = ang1; e1 = ee
        else
           a2 = ang1; e2 = ee
        end if
     end do
     sca = ( phips( irc, i ) ** 2 + slpps( i ) ** 2 ) / ( phips( irc, i ) * phibest( irc ) + slpps( i ) * slpae )
     phirn( :, i ) = phibest( : ) * sca
  end do
  !  
  call orthred( irc, irc, 1 + imax - imin, north, dl, r, phips( :, imin ), phirn( :, imin ), pspr, aepr, prec )
  !
  call projdumper( l, 'ps', irc, irc, r, dl, north, pspr, zorig, nq, dq )
  call projdumper( l, 'ae', irc, irc, r, dl, north, aepr, zorig, nq, dq )
  !
  write ( fnam, '(1a4, 1i1)' ) 'ldep', l
  open ( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  write (99, *) imin, imax 
  close (unit=99)
  !
  write ( fnam, '(1a4,1i1)' ) 'phrc', l
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1a1,2i8)' ) '#', 1 + imax - imin, irc
  do i = imin, imax
     call reconstruct( 1 + imax - imin, north, phips( :, i ), pspr, aepr, phi, irc, r, dl )
     do j = 1, irc
        write ( 99, '(6(1x,1e15.8))' ) r( j ), et( i ), phiae( j, i ), phirn( j, i ), phips( j, i ), phi( j )
     end do
     write ( 99, * )
  end do
  close( unit=99 )
  !
  return
end subroutine baregrip
