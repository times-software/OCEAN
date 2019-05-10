! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine nbsemkcmel( add04, add14 )
  use OCEAN_constants, only : Hartree2eV
  use AI_kinds
  implicit none
  !
  character(len=4) :: add04
  character(len=14) :: add14
  !
  integer :: lmin, lmax, nr, i, l, ix, mu, nu, ii, j, nrmax, k, nv, nvMax, ierr
  real( DP ) :: dl, err, rmax, su, su1, tmp, vrun, x, bwgt( 0 : 3 ), vtmp, rtmp
  real( DP ), allocatable :: v( : ), rv( : ), vtrim( : ), vdiff( : ), tmpArray( : )
!  real( DP ) :: v( 100 ), rv( 100 )
!  real( DP ) :: vtrim( 100 ), vdiff( 100 ), bwgt( 0 : 3 )
  character(len=5) :: s5
  character(len=7) :: s7
  character(len=11) :: s11
  character(len=20) :: rpot_filename
  character(len=24) :: rpotModFilename
  integer, allocatable :: nnu( : )
  real( DP ), allocatable, dimension( :, : ) :: cmel, nmel, phi, tphi
  real( DP ), allocatable, dimension( : ) :: rad, dr, val
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

  nvMax = 1000
  nv = 0
  allocate( v( nvMax ), rv( nvMax ) )
  do
    read( 99, *, iostat=ierr ) vtmp, rtmp
    ! If the read is successful
    if( ierr .eq. 0 ) then
      nv = nv + 1
      ! If we have run out of space, make v and rv bigger
      if( nv .gt. nvMax ) then
        allocate( tmpArray( nvMax ) )
        ! Don't see needing more than 1000, really
        nvMax = nvMax + nvMax
        tmpArray(:) = v(:)
        deallocate( v )
        allocate( v( nvMax ) )
        v( 1 : nv ) = tmpArray(:)
        tmpArray(:) = rv(:)
        deallocate( rv )
        allocate( rv( nvMax ) )
        rv( 1 : nv ) = tmpArray(:)
        deallocate( tmpArray )
      endif
      v( nv ) = vtmp
      rv( nv ) = rtmp
    ! end of file
    elseif( ierr < 0 ) then
      exit
    else
      exit
      ! need to add error handling to this subroutine
    endif
  enddo
  close( 99 )

  allocate( vtrim( nv ), vdiff( nv ) )
  j = 1
  k = nv
  do while ( j .lt. k )
    i = ( j + k ) / 2 
    if( rv( i ) .lt. rmax ) then
      j = i + 1
    elseif( rv( i ) .gt. rmax ) then
      k = i - 1 
    else
      exit
    endif
  enddo
  if( rv( j ) .lt. rmax ) j = j + 1

  write(6,*) 'BINARY'
  write(6,*) i, j, k 
  write(6,*) rv( j - 1 ), rmax, rv( j )
  vrun = v( j )
  do i = 1, j - 1
    vtrim( i ) = vrun
    vdiff( i ) = v( i ) - vrun
  enddo
  do i = j, nv
    vtrim( i ) = v( i )
    vdiff( i ) = 0
  enddo
!  do i = 1, 100
!     read ( 99, * ) v( i ), rv( i )
!  end do
!  close( unit=99 )
!  do i = 100, 1, -1
!     if ( rv( i ) .gt. rmax ) then
!        vrun = v( i )
!        vtrim( i ) = v( i )
!        vdiff( i ) = 0
!     else
!        vtrim( i ) = vrun
!        vdiff( i ) = v( i ) - vtrim( i )
!     end if     
!  end do
  write(rpotModFilename, '(A10,A14)') 'rpottrim.z', add14
  open( unit=99, file=rpotModfilename, form='formatted', status='unknown' )
  rewind 99
  write(99, '(A,I8)' ) '#', nv
  do i = 1, nv !100
     write ( 99, '(E24.16,X,E24.16)' ) vtrim( i ), rv( i )
  end do
  close( unit=99 )
  
  write(rpotModFilename, '(A10,A14)') 'rpotdiff.z', add14
  open( unit=99, file=rpotModfilename, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nv !100
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
        call intval( nv, rv, vdiff, rad( i ), val( i ), 'err', 'err' )
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

  deallocate( v, rv, vdiff, vtrim )
  !
  return
end subroutine nbsemkcmel
