! Copyright (C) 2016, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine OCEAN_get_rho( xmesh, celvol, rho, ierr )
  use AI_kinds
  use OCEAN_mpi
  use FFT_wrapper
  use iso_c_binding
  implicit none

!  include 'fftw3.f03'

  integer, intent( in ) :: xmesh( 3 )
  real(dp), intent( in ) :: celvol
  real(dp), intent( out ) :: rho( xmesh(3), xmesh(2), xmesh(1) )
  integer, intent( inout ) :: ierr

  !
  integer :: nfft( 3 ), dumint, i, j, k, ii, jj, kk, ierr_
  complex(dp),  allocatable :: rhoofr(:,:,:), rhoofg(:,:,:)
  real(dp) :: dumr, sumrho !, norm
!  integer*8 :: fftw_plan, fftw_plan2
  type( fft_obj ) :: fo
  integer :: igl(3), igh(3), powlist(3), mul(3), c_nfft(3), iter, offset(3)
  character(len=1) :: fstr

  integer, parameter :: fac(3) = (/ 2, 3, 5 /)
  integer, parameter :: hicap(3) = (/ 20, 8, 4 /)
  integer, parameter :: nfac = 3

  integer, external :: optim




! Lazy start, only MPI master node does work
  if( myid .eq. root ) then
    write(6,*) 'Reading in rho'

    open( unit=99, file='nfft', form='formatted', status='old', IOSTAT=ierr )
    if ( ierr .ne. 0 ) goto 111
    read( 99, * ) nfft( : )
    close( 99 )

    allocate( rhoofr( nfft(1), nfft(2), nfft(3) ) )
!    call fftw_plan_dft_3d( fftw_plan, nfft(1), nfft(2), nfft(3), rhoofr, rhoofr, &
!            FFTW_FORWARD, FFTW_ESTIMATE )
    call FFT_wrapper_init( nfft, fo, rhoofr )
    rhoofr = 0.0_dp


    sumrho = 0.0_dp
    open( unit=99, file='rhoofr', form='formatted', status='old', IOSTAT=ierr )
    if ( ierr .ne. 0 ) goto 111
    read( 99, * ) fstr
    do k = 1, nfft( 3 )
      do j = 1, nfft( 2 )
        do i = 1, nfft( 1 )
          read( 99, * ) dumint, dumint, dumint, dumr
          sumrho = sumrho + dumr
          rhoofr( i, j, k ) = dumr 
        enddo
      enddo
    enddo
    close( 99 )

    write(6,*) 'rhoofr loaded'

    sumrho = sumrho / dble(product( nfft )) * celvol
    write(6,*) '!', sumrho
    write(6,*) '!', minval( real( rhoofr, DP) )
    !
!    call fftw_execute_dft( fftw_plan, rhoofr, rhoofr )
!    call fftw_destroy_plan( fftw_plan )
    call FFT_wrapper_single( rhoofr, OCEAN_FORWARD, fo )
    call FFT_wrapper_delete( fo )
    !
    write( 6, * ) dble(rhoofr( 1, 1, 1 ) ) / dble(product( nfft ))  * celvol , &
          anint( dble(rhoofr( 1, 1, 1 ) ) / dble( product( nfft ) ) * celvol )
!    write( 6, * ) dble(rhoofr( 1, 1, 1 ) ) * celvol , &
!          anint( dble(rhoofr( 1, 1, 1 ) ) * celvol )

    ! need to find compatible grids
    do iter = 1, 3
      call facpowfind( xmesh( iter ), nfac, fac, powlist )
      c_nfft( iter ) = optim( nfft( iter ), nfac, fac, powlist, hicap )
      if( mod( c_nfft( iter ), xmesh( iter ) ) .ne. 0 ) then
        write(6,*) 'you are doing something wrong'
        ierr = 14
        goto 111
      endif
    enddo


    allocate( rhoofg( c_nfft( 1 ), c_nfft( 2 ), c_nfft( 3 ) ), STAT=ierr )
    if( ierr .ne. 0 ) goto 111
!    call  fftw_plan_dft_3d( fftw_plan2, c_nfft(1), c_nfft(2), c_nfft(3), rhoofg, rhoofg, &
!            FFTW_BACKWARD, FFTW_ESTIMATE )
    call FFT_wrapper_init( c_nfft, fo, rhoofg )
    rhoofg = 0.0_dp
    mul( : ) = c_nfft( : ) / nfft( : )
    write(6,*) c_nfft( : )
    write(6,*) nfft( : )
    write(6,*) xmesh( : )
    offset( : ) = c_nfft( : ) - nfft( : )
    write(6,*) offset( : )
    write(6,*) '-----------'
    do k = 1, nfft( 3 ) !/ 2 
      if( k .le. nfft(3)/2 ) then !floor( dble(nfft( 3 )) / 2.0 ) + 1) then
        kk = k
      else
        kk = k + offset( 3 )
      endif
      do j = 1, nfft( 2 ) !/ 2
        if( j .le. nfft(2)/2) then !floor( dble(nfft( 2 ) ) / 2.0) + 1 ) then
          jj = j
        else
          jj  = j + offset( 2 )
        endif
        do i = 1, nfft( 1 ) !/ 2
          if( i .le. nfft(1)/2) then !floor(dble(nfft( 1 ) ) / 2.0) + 1 ) then
            ii = i
          else
            ii = i + offset( 1 )
          endif
          rhoofg( ii, jj, kk ) = rhoofr( i, j, k )
        enddo
      enddo
    enddo
!    call fftw_execute_dft( fftw_plan2, rhoofg, rhoofg )
!    call fftw_destroy_plan( fftw_plan2 )
    call FFT_wrapper_single( rhoofg, OCEAN_BACKWARD, fo )
    call FFT_wrapper_delete( fo )
    !
!    norm = real( product( nfft ), DP )
!    norm = 1.0_dp / norm
!    norm = 1.0_dp
    mul( : ) = c_nfft( : ) / xmesh( : )
    write(6,*) mul( : )
    do i = 0, xmesh( 1 ) - 1
      do j = 0, xmesh( 2 ) - 1
        do k = 0, xmesh( 3 ) - 1
          rho( k+1, j+1, i+1 ) = real(rhoofg( i*mul(1)+1, j*mul(2)+1, k*mul(3)+1 ), DP )
!          rho( k+1, j+1, i+1 ) = norm * real(rhoofg( i*mul(1)+1, j*mul(2)+1, k*mul(3)+1 ), DP )
          if( rho( k+1, j+1, i+1 ) .le. 0.0d0 ) then
            write(6,*) 'low density', rho( k+1, j+1, i+1 ), k+1, j+1, i+1, rhoofg( i*mul(1)+1, j*mul(2)+1, k*mul(3)+1 )
            if( rho( k+1, j+1, i+1 ) .gt. -1.0d-12 ) then
              rho( k+1, j+1, i+1 ) = 0.d0
            else
              write(6,*) 'Negative density found!', rho( k+1, j+1, i+1 ), k+1, j+1, i+1
              ierr = 11
              goto 111
            endif
          endif
        enddo
      enddo
    enddo

    open( unit=99, file='rho2.xpts', form='formatted', status='unknown')
    rewind( 99 )
    sumrho = 0.d0
    do k = 1, xmesh( 3 )
      do j = 1, xmesh( 2 )
        do i = 1, xmesh( 1 )
          write(99,*) i, j, k, rho(k,j,i)
          sumrho = sumrho +  rho(k,j,i)
        enddo
      enddo
    enddo
    close( 99 )
    sumrho = sumrho / dble(product(xmesh)) * celvol
    write(6,*) '!', sumrho
    !
111 continue
    if( allocated( rhoofr ) ) deallocate( rhoofr )
    if( allocated( rhoofg ) ) deallocate( rhoofg )
  endif


#ifdef MPI
  call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
  if( ierr .ne. 0 ) return
  if( ierr_ .ne. 0 ) then
    ierr = ierr_
    return
  endif
  call MPI_BCAST( rho, product(xmesh), MPI_DOUBLE_PRECISION, root, comm, ierr )
  if( ierr .ne. 0 ) return
!  write(6,*) 'Rho shared: ', myid
#endif

end subroutine OCEAN_get_rho




! this finds powers [pow] of up to [nfac] prime factors [fac] 
! to multiply to make [n].
!
subroutine facpowfind( n, nfac, fac, pow )
  implicit none
  !
  integer :: n, nfac
  integer :: fac( nfac ), pow( nfac )
  !
  integer :: ifac, nred, m
  !
  nred = n
  pow( : ) = 0
  ifac = 1
  do while ( nred .gt. 1 )
     m = nred / fac( ifac )
     if ( nred .eq. m * fac( ifac ) ) then
        pow( ifac ) = pow( ifac ) + 1
        nred = m
     else
        ifac = ifac + 1
        if ( ifac .gt. nfac ) stop 'incommensurate grid'
     end if
  end do
  !
  return
end subroutine facpowfind


