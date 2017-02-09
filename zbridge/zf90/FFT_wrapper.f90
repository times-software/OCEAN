! Copyright (C) 2015-2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! The purpose of this module is to wrap all ocean FFT calls.
! Currently we support our legacy mode fft and fftw, but this will provide a 
! centralized location for adding others.
!
module FFT_wrapper

#ifdef __FFTW3
  use iso_c_binding
  include 'fftw3.f03'
#endif

  private
  save

  integer, parameter :: OCEAN_FORWARD = 1
  integer, parameter :: OCEAN_BACKWARD = -1 

  type fft_obj
    integer :: dims(4)
    integer :: jfft

    real(kind=kind(1.0d0)) :: norm
#ifdef __FFTW3
    type(C_PTR) :: fplan, bplan
#endif
  end type fft_obj

  public :: fft_obj
  public :: OCEAN_FORWARD, OCEAN_BACKWARD
  public :: FFT_wrapper_init, FFT_wrapper_delete, FFT_wrapper_split, FFT_wrapper_single

  contains

  subroutine FFT_wrapper_init( zn, fo, io )!, nthread )
    implicit none

    integer, intent(in ) :: zn(3)
    complex(kind=kind(1.0d0)), intent(inout), optional :: io(zn(1),zn(2),zn(3))
!    integer, intent( in ), optional :: nthread
    type( fft_obj ), intent( out ) :: fo
    complex(kind=kind(1.0d0)), allocatable :: cwrk(:)
    !
    fo%dims(1:3) = zn(:)
    fo%dims(4) = product( fo%dims(1:3) )
#ifdef __FFTW3
    fo%norm = 1.0d0 / dble( fo%dims(4) )

!    if( present( nthread ) ) then
!      if( nthread .gt. 0 ) call fftw_plan_with_nthreads( nthreads )
!    endif
    
    if( present( io ) ) then
      fo%fplan = fftw_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), io, io, FFTW_FORWARD, FFTW_PATIENT )
      fo%bplan = fftw_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), io, io, FFTW_BACKWARD, FFTW_PATIENT )
    else
      allocate( cwrk(fo%dims(4)))
      fo%fplan = fftw_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), cwrk, cwrk, FFTW_FORWARD, FFTW_PATIENT )
      fo%bplan = fftw_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), cwrk, cwrk, FFTW_BACKWARD, FFTW_PATIENT )
      deallocate( cwrk )
    endif

!    if( present( nthread ) ) call fftw_plan_with_nthreads( 1 )
#else
    fo%norm = 1.0d0
    fo%jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
#endif
  end subroutine FFT_wrapper_init

  subroutine FFT_wrapper_delete( fo )
    implicit none
    type( fft_obj ), intent( inout ) :: fo

#ifdef __FFTW3
    call fftw_destroy_plan(fo%fplan)
    call fftw_destroy_plan(fo%bplan)
#endif

  end subroutine FFT_wrapper_delete


  subroutine FFT_wrapper_split( r, i, dir, fo )
    implicit none
    type( fft_obj ), intent( inout ) :: fo
    real(kind=kind(1.0d0)), intent( inout ) :: r(fo%dims(4))
    real(kind=kind(1.0d0)), intent( inout ) :: i(fo%dims(4))
    integer, intent( in ) :: dir
    !
#ifdef __FFTW3
    complex(kind=kind(1.0d0)), allocatable :: wrk(:)
#else
    real(kind=kind(1.0d0)), allocatable :: wrk(:)
#endif
    

#ifdef __FFTW3
    allocate( wrk( fo%dims(4) ) )
    wrk(:) = cmplx(r(:),i(:))
    if( dir .eq. OCEAN_FORWARD ) then
      call fftw_execute_dft( fo%fplan, wrk, wrk )
    else
      call fftw_execute_dft( fo%bplan, wrk, wrk )
    endif
    r(:) = real(wrk(:))*fo%norm
    i(:) = aimag(wrk(:))*fo%norm
#else
    allocate( wrk( fo%jfft ) )
    call cfft( r, i, fo%dims(1), fo%dims(1), fo%dims(2), fo%dims(3), dir, wrk, fo%jfft )
#endif

    deallocate( wrk )
  end subroutine FFT_wrapper_split



  subroutine FFT_wrapper_single( io, dir, fo )
    implicit none
    type( fft_obj ), intent( inout ) :: fo
#ifdef __FFTW3
    complex(kind=kind(1.0d0)), intent( inout ) :: io( fo%dims(1), fo%dims(2), fo%dims(3) )
#else
    complex(kind=kind(1.0d0)), intent( inout ) :: io( fo%dims(4))
#endif
    integer, intent( in ) :: dir
    !
#ifndef __FFTW3
    real(kind=kind(1.0d0)), allocatable :: wrk(:), r(:), i(:)
#endif

#ifdef __FFTW3
    if( dir .eq. OCEAN_FORWARD ) then
      call fftw_execute_dft( fo%fplan, io, io )
    else
      call fftw_execute_dft( fo%bplan, io, io )
    endif
    io(:,:,:) = io(:,:,:) * fo%norm
#else
    allocate( r( fo%dims(4) ), i( fo%dims(4) ), wrk( fo%jfft ) )
    r(:) = real(io(:))
    i(:) = aimag(io(:))
    call cfft( r, i, fo%dims(1), fo%dims(1), fo%dims(2), fo%dims(3), dir, wrk, fo%jfft )
    io(:) = cmplx(r(:),i(:))
    deallocate( r, i, wrk )
#endif
  end subroutine FFT_wrapper_single



end module
