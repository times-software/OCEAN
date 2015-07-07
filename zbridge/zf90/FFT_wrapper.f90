module FFT_wrapper

#ifdef __FFTW3
  use iso_c_binding
  include 'fftw3.f03'
#endif

  private
  save
  
  real(kind=kind(1.0d0)), allocatable :: wrk(:)
  complex(kind=kind(1.0d0)), allocatable :: cwrk(:)

  integer :: dims(4)
  integer :: jfft

  real(kind=kind(1.0d0)) :: norm

  integer, parameter :: OCEAN_FORWARD = 1
  integer, parameter :: OCEAN_BACKWARD = -1 

#ifdef __FFTW3
  type(C_PTR) :: fplan, bplan
#endif

  public :: OCEAN_FORWARD, OCEAN_BACKWARD
  public :: FFT_wrapper_init, FFT_wrapper_delete, FFT_wrapper_split, FFT_wrapper_single

  contains

  subroutine FFT_wrapper_init( zn )
    implicit none

    integer, intent(in ) :: zn(3)
    integer :: fftw_flags
    !
    dims(1:3) = zn(:)
    dims(4) = product( dims(1:3) )
#ifdef __FFTW3
    allocate( cwrk(dims(4)))
    norm = 1.0d0 / dble( dims(4) )
    
    fplan = fftw_plan_dft_3d( dims(1), dims(2), dims(3), cwrk, cwrk, FFTW_FORWARD, FFTW_PATIENT )
    bplan = fftw_plan_dft_3d( dims(1), dims(2), dims(3), cwrk, cwrk, FFTW_BACKWARD, FFTW_PATIENT )

#else
    norm = 1.0d0
    jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
    allocate( wrk( jfft ) )
#endif
  end subroutine FFT_wrapper_init

  subroutine FFT_wrapper_delete()
    implicit none

#ifdef __FFTW3
    call dfftw_destroy_plan(fplan)
    call dfftw_destroy_plan(bplan)
    deallocate(cwrk)
#else
    deallocate(wrk)
#endif

  end subroutine FFT_wrapper_delete


  subroutine FFT_wrapper_split( r, i, dir )
    implicit none
    real(kind=kind(1.0d0)), intent( inout ) :: r(dims(4))
    real(kind=kind(1.0d0)), intent( inout ) :: i(dims(4))
    integer, intent( in ) :: dir

#ifdef __FFTW3
    cwrk(:) = cmplx(r(:),i(:))
    if( dir .eq. OCEAN_FORWARD ) then
      call fftw_execute_dft( fplan, cwrk, cwrk )
    else
      call fftw_execute_dft( bplan, cwrk, cwrk )
    endif
    r(:) = real(cwrk(:))*norm
    i(:) = aimag(cwrk(:))*norm
#else
    call cfft( r, i, dims(1), dims(1), dims(2), dims(3), dir, wrk, jfft )
#endif
  end subroutine FFT_wrapper_split



  subroutine FFT_wrapper_single( io, dir )
    implicit none
    complex(kind=kind(1.0d0)), intent( inout ) :: io(dims(4))
    integer, intent( in ) :: dir
    !
#ifndef __FFTW3
    real(kind=kind(1.0d0)) :: r(dims(4)), i(dims(4))
#endif

#ifdef __FFTW3
    if( dir .eq. OCEAN_FORWARD ) then
      call fftw_execute_dft( fplan, io, io )
    else
      call fftw_execute_dft( bplan, io, io )
    endif
    io(:) = io(:) * norm
#else
    r(:) = real(io(:))
    i(:) = aimag(io(:))
    call cfft( r, i, dims(1), dims(1), dims(2), dims(3), dir, wrk, jfft )
    io(:) = cmplx(r(:),i(:))
#endif
  end subroutine FFT_wrapper_single



end module
