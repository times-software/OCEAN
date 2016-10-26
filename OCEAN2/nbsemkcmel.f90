! Copyright (C) 2015,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine nbsemkcmel( add04, add14 )
  use OCEAN_constants, only : Hartree2eV
  implicit none
  !
  character * 4 :: add04
  character * 14 :: add14
  !
  integer :: lmin, lmax, nr, i, l, ix, mu, nu, ii, j, nrmax
  real( kind = kind( 1.0d0 ) ) :: dl, err, rmax, su, su1, tmp, vrun, x
  real( kind = kind( 1.0d0 ) ) :: v( 100 ), rv( 100 )
  real( kind = kind( 1.0d0 ) ) :: vtrim( 100 ), vdiff( 100 ), bwgt( 0 : 3 )
  character * 5 :: s5
  character * 7 :: s7
  character * 11 :: s11
  character * 20 :: rpot_filename
  integer, allocatable :: nnu( : )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: cmel, nmel, phi, tphi
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rad, dr, val
  !
  write ( s11, '(1a7,1a4)' ) 'prjfile', add04
  open( unit=99, file=s11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lmin, lmax
  allocate( nnu( lmin : lmax ) )
  read ( 99, * ) nnu( : )
  close( unit=99 )
  !
  write ( s11, '(1a7,1a4)' ) 'radfile', add04
  open( unit=99, file=s11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) rmax, nrmax, nr
  if ( rmax .gt. 9.9d0 ) stop 'bad rmax'
  close( unit=99 )
  !
  write(rpot_filename, '(A6,A14)') 'rpot.z', add14
!  open( unit=99, file='rpotfull', form='formatted', status='unknown' )
  open( unit=99, file=rpot_filename, form='formatted', status='unknown' )
  rewind 99
  do i = 1, 100
     read ( 99, * ) v( i ), rv( i )
  end do
  close( unit=99 )
  do i = 100, 1, -1
     if ( rv( i ) .gt. rmax ) then
        vrun = v( i )
        vtrim( i ) = v( i )
        vdiff( i ) = 0
     else
        vtrim( i ) = vrun
        vdiff( i ) = v( i ) - vtrim( i )
     end if     
  end do
  open( unit=99, file='rpottrim', form='formatted', status='unknown' )
  rewind 99
  do i = 1, 100
     write ( 99, '(2f10.5)' ) vtrim( i ), rv( i )
  end do
  close( unit=99 )
  open( unit=99, file='rpotdiff', form='formatted', status='unknown' )
  rewind 99
  do i = 1, 100
     write ( 99, '(2f10.5)' ) vdiff( i ), rv( i )
  end do
  close( unit=99 )
  !
  bwgt( 0 ) = 28.0d0 / 45.0d0 
  bwgt( 1 ) = 64.0d0 / 45.0d0 
  bwgt( 2 ) = 24.0d0 / 45.0d0 
  bwgt( 3 ) = 64.0d0 / 45.0d0 
  do l = lmin, lmax
     allocate( cmel( nnu( l ), nnu( l ) ) )
     allocate( nmel( nnu( l ), nnu( l ) ) )
     write ( s7, '(1a2,1i1,1a4)' ) 'ae', l, add04
     open( unit=99, file=s7, form='formatted', status='unknown' )
     rewind 99
     allocate( rad( nr ), phi( nr, nnu( l ) ), dr( nr ), val( nr ), tphi( nnu( l ), nr ) )
     do i = 1, nr
        read ( 99, * ) rad( i ), tphi( :, i )
     end do
     phi = transpose( tphi )
     deallocate( tphi )
     close( unit=99 )
     dl = log( rad( nr ) / rad( 1 ) ) / dble( nr - 1 )
     x = 1 + nint( log( rmax / rad( 1 ) ) / dl )
     ix = x
     err = abs( dble( ix ) - x ) 
     if ( abs( err ) .gt. 0.0001d0 ) stop 'bad commens.'
     dr( 1 ) = 0.5d0 * dl * rad( 1 )
     do i = 2, ix - 1
        dr( i ) = dl * rad( i ) 
     end do
     dr( ix ) = 0.5d0 * dl * rad( ix )
     do i = 1, ix
        call intval( 100, rv, vdiff, rad( i ), val( i ), 'err', 'err' )
     end do
     do nu = 1, nnu( l )
        do mu = 1, nnu( l )
           su = 0
           su1 = 0
           do i = ix, 5, -4
              do j = 0, 3
                 ii = i - j
                 tmp = phi( ii, mu ) * phi( ii, nu ) * val( ii )
                 su = su + bwgt( j ) * dr( ii ) * rad( ii ) ** 2 * tmp
                 su1 = su1 + bwgt( j ) * dr( ii ) * rad( ii ) ** 2 * phi( ii, mu ) * phi( ii, nu )
              end do 
           end do
           cmel( mu, nu ) = su * Hartree2eV
           nmel( mu, nu ) = su1
        end do
     end do
     write ( s5, '(1a4,1i1)' ) 'cmel', l
     open( unit=99, file=s5, form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(4(1x,1e15.8))' ) cmel( :, : )
     write ( 99, '(4(1x,1e15.8))' ) nmel( :, : )
     close( unit=99 )
     deallocate( cmel, nmel, rad, phi, dr, val )
  end do
  !
  return
end subroutine nbsemkcmel
