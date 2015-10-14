! Copyright (C) 2013 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cainchact( vr, vi, n, e0, se, hvr, hvi, nq, nbd, zn, inter, amet, nx, ny, nz, ur, ui, & 
                      tau, rcut, rzero, ptab, eps, core_offset )
  implicit none
  !
  integer :: nbd, nq, n, nx, ny, nz, zn( 3 )
  real( kind = kind( 1.0d0 ) ) :: inter, amet( 3, 3 ), tau( 3 )
  real( kind = kind( 1.0d0 ) ) :: vr( n ), vi( n ), e0( n ), hvr( n ), hvi( n )
  complex( kind = kind( 1.0d0 ) ) :: se( n )
  real( kind = kind( 1.0d0 ) ) :: ur( nx * ny * nz, nbd, nq ), ui( nx * ny * nz, nbd, nq )
  real( kind = kind( 1.0d0 ) ) :: rcut, rzero, ptab( 100 )
  real( kind = kind( 1.0d0 ) ), intent( in ) :: core_offset
  !
  integer :: nn1, nfft, jfft
  real( kind = kind( 1.0d0 ) ) :: eps, epsi
  real( kind = kind( 1.0d0 ) ), dimension( n ) :: hvcr, hvci 
  real( kind = kind( 1.0d0 ) ), allocatable :: xwrkr( :, : ), xwrki( :, : ), wrk( : )
  real( kind = kind( 1.0d0 ) ), external :: wav
  logical :: tdlda
  !
  ! diagonal part of channel hamiltonian ...
  hvr( : ) = ( e0( : ) + core_offset ) * vr( : ) + real( se( : ) ) * vr( : ) - aimag( se( : ) ) * vi( : ) 
  hvi( : ) = ( e0( : ) + core_offset ) * vi( : ) + aimag( se( : ) ) * vr( : ) + real( se( : ) ) * vi( : )
!  write(6,*) maxval( real( se ) ), maxval( e0 )
!  write(6,*) minval( aimag( se ) ) , maxval( aimag( se ) )

! if ( 1 .lt. 2 ) return 
  if ( inter .le. 0.01d0 ) return

  inquire( file='tdlda', exist=tdlda )
!  if( tdlda ) then
!    open(unit=99, file='tdlda', form='formatted', status='old' )
!    read( 99, * ) tdlda
!    close( 99 )
!  endif
  if ( tdlda ) return
  !
  ! interaction part of channel hamiltonian ... set up ...
!  open( unit=99, file='epsilon', form='formatted', status='unknown' )
!  rewind 99
!  read ( 99, * ) eps
!  close( unit=99 )
  epsi = 1.d0 / eps
  !
  ! ground work for fft to real space ...
  nn1 = zn( 1 )
  nfft = zn( 1 ) * zn( 2 ) * zn( 3 )
  jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
  allocate( xwrkr( nfft, nx * ny * nz ), xwrki( nfft, nx * ny * nz ), wrk( jfft ) )
  !
  ! fft states to real space ...
  call cainnqtoxr( nx, ny, nz, nfft, nq, nbd, vr, vi, xwrkr, xwrki, tau, nn1, zn, wrk, jfft, ur, ui )
  !
  ! multiply the wave function times the potential in real space ...
  call cainxrvmult( nx, ny, nz, zn, tau, amet, epsi, xwrkr, xwrki, ptab )
  !
  ! go back to ( n, q) space ... then normalize as needed ...
  call cainxrtonq( nx, ny, nz, nfft, nq, nbd, hvcr, hvci, xwrkr, xwrki, tau, nn1, zn, wrk, jfft, ur, ui )
  !
  ! add to result
  hvr( : ) = hvr( : ) - inter * hvcr( : )
  hvi( : ) = hvi( : ) - inter * hvci( : )
  !
  deallocate( xwrkr, xwrki, wrk )
  !
  return
end subroutine cainchact
