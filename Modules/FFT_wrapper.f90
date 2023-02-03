! Copyright (C) 2015-2017 OCEAN collaboration
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
  use ai_kinds

#ifdef __FFTW3
  use iso_c_binding
  include 'fftw3.f03'
#endif

  private
  save

  ! Forward is named to match FFTW convention
  integer, parameter :: OCEAN_FORWARD = -1
  integer, parameter :: OCEAN_BACKWARD = 1 

  type fft_obj
    integer :: dims(4)
    integer :: jfft

    real(DP) :: norm
    real(SP) :: norm_sp
#ifdef __FFTW3
    type(C_PTR) :: fplan, bplan
#endif
    logical :: is_sp
  end type fft_obj

  public :: fft_obj
  public :: OCEAN_FORWARD, OCEAN_BACKWARD
  public :: FFT_wrapper_init, FFT_wrapper_delete, FFT_wrapper_split, FFT_wrapper_single
  public :: FFT_wrapper_init_sp, FFT_wrapper_single_sp

  contains

  subroutine FFT_wrapper_init( zn, fo, io, fh )!, nthread )
    implicit none

    integer, intent(in ) :: zn(3)
    complex(kind=kind(1.0d0)), intent(inout), optional :: io(zn(1),zn(2),zn(3))
    integer, optional :: fh
!    integer, intent( in ), optional :: nthread
    type( fft_obj ), intent( out ) :: fo
    complex(kind=kind(1.0d0)), allocatable :: cwrk(:)
    !
    fo%is_sp = .false.
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
    if( present( fh ) ) then
      write(fh,*) 'Plan using FFTW:', fo%dims(:)
    endif
#else
    fo%norm = 1.0d0
    fo%jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
    if( present( fh ) ) then
      write(fh,*) 'Plan using Legacy:', fo%dims(:)
    endif
#endif
  end subroutine FFT_wrapper_init

  subroutine FFT_wrapper_init_sp( zn, fo, io, fh )!, nthread )
    implicit none

    integer, intent(in ) :: zn(3)
    complex(SP), intent(inout), optional :: io(zn(1),zn(2),zn(3))
    integer, optional :: fh
!    integer, intent( in ), optional :: nthread
    type( fft_obj ), intent( out ) :: fo
    complex(SP), allocatable :: cwrk(:)
#ifndef __FFTW3F
    complex(DP), allocatable :: cwrk_dp(:)
#endif
    !
    fo%dims(1:3) = zn(:)
    fo%dims(4) = product( fo%dims(1:3) )
#ifdef __FFTW3
    fo%is_sp = .true.
    fo%norm = 1.0d0 / dble( fo%dims(4) )
    fo%norm_sp = fo%norm

#ifdef __FFTW3F
    if( present( io ) ) then
      fo%fplan = fftwf_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), io, io, FFTW_FORWARD, FFTW_PATIENT )
      fo%bplan = fftwf_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), io, io, FFTW_BACKWARD, FFTW_PATIENT )
    else
      allocate( cwrk(fo%dims(4)))
      fo%fplan = fftwf_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), cwrk, cwrk, FFTW_FORWARD, FFTW_PATIENT )
      fo%bplan = fftwf_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), cwrk, cwrk, FFTW_BACKWARD, FFTW_PATIENT )
      deallocate( cwrk )
    endif
#else
    fo%is_sp = .false.
    allocate( cwrk_dp(fo%dims(4)))
    fo%fplan = fftw_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), cwrk_dp, cwrk_dp, FFTW_FORWARD, FFTW_PATIENT )
    fo%bplan = fftw_plan_dft_3d( fo%dims(3), fo%dims(2), fo%dims(1), cwrk_dp, cwrk_dp, FFTW_BACKWARD, FFTW_PATIENT )
    deallocate( cwrk_dp )
#endif

#else
    fo%is_sp = .false.
    fo%norm = 1.0_DP
    fo%norm_sp = 1.0_SP
    fo%jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
    if( present( fh ) ) then
      write(fh,*) 'Plan using Legacy:', fo%dims(:)
    endif
#endif
  end subroutine FFT_wrapper_init_sp


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
    type( fft_obj ), intent( in ) :: fo
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
      r(:) = real(wrk(:),kind(1.0d0))*fo%norm
      i(:) = aimag(wrk(:))*fo%norm
    endif
#else
    allocate( wrk( fo%jfft ) )
    call cfft( r, i, fo%dims(1), fo%dims(1), fo%dims(2), fo%dims(3), dir, wrk, fo%jfft )
#endif

    deallocate( wrk )
  end subroutine FFT_wrapper_split



  subroutine FFT_wrapper_single( io, dir, fo, norm )
    implicit none
    type( fft_obj ), intent( in ) :: fo
#ifdef __FFTW3
    complex(kind=kind(1.0d0)), intent( inout ) :: io( fo%dims(1), fo%dims(2), fo%dims(3) )
#else
    complex(kind=kind(1.0d0)), intent( inout ) :: io( fo%dims(4))
#endif
    integer, intent( in ) :: dir
    logical, intent( in ), optional :: norm
    !
#ifndef __FFTW3
    real(kind=kind(1.0d0)), allocatable :: wrk(:), r(:), i(:)
#endif
    logical :: normalize 

    if( present( norm ) ) then
      normalize = norm
    else
      normalize = .true.
    endif

#ifdef __FFTW3
    if( dir .eq. OCEAN_FORWARD ) then
      call fftw_execute_dft( fo%fplan, io, io )
    else
      call fftw_execute_dft( fo%bplan, io, io )
      if( normalize ) then
        io(:,:,:) = io(:,:,:) * fo%norm
      endif
    endif
#else
    allocate( r( fo%dims(4) ), i( fo%dims(4) ), wrk( fo%jfft ) )
    r(:) = real(io(:), kind(1.0d0))
    i(:) = aimag(io(:))
    call cfft( r, i, fo%dims(1), fo%dims(1), fo%dims(2), fo%dims(3), dir, wrk, fo%jfft )
    if( dir .eq. OCEAN_BACKWARD .and. normalize .eqv. .false. ) then
      io(:) = cmplx(r(:),i(:), kind(1.0d0)) * dble( fo%dims(4) )
    else
      io(:) = cmplx(r(:),i(:), kind(1.0d0))
    endif
    deallocate( r, i, wrk )
#endif
  end subroutine FFT_wrapper_single

  subroutine FFT_wrapper_single_sp( io, dir, fo, norm )
    implicit none
    type( fft_obj ), intent( in ) :: fo
#ifdef __FFTW3
    complex(SP), intent( inout ) :: io( fo%dims(1), fo%dims(2), fo%dims(3) )
#else
    complex(SP), intent( inout ) :: io( fo%dims(4))
#endif
    integer, intent( in ) :: dir
    logical, intent( in ), optional :: norm
    !
#ifndef __FFTW3
    real(DP), allocatable :: wrk(:), r(:), i(:)
#endif
#ifndef __FFTW3F
    complex(DP), allocatable :: cwrk(:,:,:)
#endif
    logical :: normalize

    if( present( norm ) ) then
      normalize = norm
    else
      normalize = .true.
    endif

#ifdef __FFTW3
    if( dir .eq. OCEAN_FORWARD ) then
#ifdef __FFTW3F
      call fftwf_execute_dft( fo%fplan, io, io )
#else
      allocate( cwrk( fo%dims(1), fo%dims(2), fo%dims(3) ) )
      cwrk(:,:,:) = io(:,:,:)
      call fftw_execute_dft( fo%fplan, cwrk, cwrk )
      io(:,:,:) = cwrk(:,:,:)
      deallocate( cwrk )
#endif
    else
#ifdef __FFTW3F
      call fftwf_execute_dft( fo%bplan, io, io )
      if( normalize ) then
        io(:,:,:) = io(:,:,:) * fo%norm_sp
      endif
#else
      allocate( cwrk( fo%dims(1), fo%dims(2), fo%dims(3) ) )
      cwrk(:,:,:) = io(:,:,:)
      call fftw_execute_dft( fo%bplan, cwrk, cwrk )
      if( normalize ) then
        io(:,:,:) = cwrk(:,:,:) * fo%norm
      else
        io(:,:,:) = cwrk(:,:,:)
      endif
      deallocate( cwrk )
#endif
    endif
#else
    allocate( r( fo%dims(4) ), i( fo%dims(4) ), wrk( fo%jfft ) )
    r(:) = real(io(:), kind(1.0d0))
    i(:) = aimag(io(:))
    call cfft( r, i, fo%dims(1), fo%dims(1), fo%dims(2), fo%dims(3), dir, wrk, fo%jfft )
    if( dir .eq. OCEAN_BACKWARD .and. normalize .eqv. .false. ) then
      io(:) = cmplx(r(:),i(:), kind(1.0d0)) * dble( fo%dims(4) )
    else
      io(:) = cmplx(r(:),i(:), kind(1.0d0))
    endif
    deallocate( r, i, wrk )
#endif
  end subroutine FFT_wrapper_single_sp

#if 0
  !!! LEGACY
  subroutine cfft(chdr,chdi,nn1,n1,n2,n3,mode,wrk,idwrk)
  !
  !   chdr= ar, chdi= ai, wrk= work
  !   wrk(1)= trigs, ifax= ifax
  !   (inc,jump,n,lot)
  !   mode= isgn
  !
  !      computes a complex 3d fast fourier transform
  !      using the cray2 scilib subroutines
  !      written april 12, 1988. jlm
  !
     implicit real*8 (a-h,o-z)
  !
     dimension chdr(nn1,n2,n3),chdi(nn1,n2,n3),wrk(idwrk)
     dimension ifax(19)
  !
     um = 1.d0
  !
     if (n1.ne.1) then
     call cftfax(n1,ifax,wrk)
     do 10 i=1,n3
       call cfftmlt(chdr(1,1,i),chdi(1,1,i),wrk(2*n1+1),wrk(1),ifax,1,nn1,n1,n2,mode)
  10    continue
     endif
  !
     if (n2.ne.1) then
     call cftfax(n2,ifax,wrk)
     do 20 i=1,n3
       call cfftmlt(chdr(1,1,i),chdi(1,1,i),wrk(2*n2+1),wrk(1),ifax,nn1,1,n2,n1,mode)
  20    continue
     endif
  !
     if (n3.ne.1) then
     call cftfax(n3,ifax,wrk)
     do 30 i=1,n2
       call cfftmlt(chdr(1,i,1),chdi(1,i,1),wrk(2*n3+1),wrk(1),ifax,nn1*n2,1,n3,n1,mode)
  30    continue
     endif
  !
     if (mode .eq. 1) then
       fac = um / real(n1*n2*n3)
       do k=1,n3
         do j=1,n2
           do i=1,n1
             chdr(i,j,k) = fac*chdr(i,j,k)
             chdi(i,j,k) = fac*chdi(i,j,k)
           enddo
         enddo
       enddo
     endif
  !
     return
  end subroutine cfft
  !************************************************************************
  subroutine cfftmlt(ar,ai,work,trigs,ifax,inc,jump,n,lot,isgn)
    implicit real*8 (a-h,o-z)
  !
  ! purpose      performs multiple fast fourier transforms.  this package
  !              will perform a number of simultaneous complex periodic
  !              fourier transforms or corresponding inverse transforms.
  !              that is, given a set of complex gridpoint vectors, the
  !              package returns a set of complex fourier
  !              coefficient vectors, or vice versa.  the length of the
  !              transforms must be a number greater than 1 that has
  !              no prime factors other than 2, 3, and 5.
  !
  !              the package cfft99 contains several user-level routines:
  !
  !            subroutine cftfax
  !                an initialization routine that must be called once
  !                before a sequence of calls to cfft99
  !                (provided that n is not changed).
  !
  !            subroutine cfft99
  !                the actual transform routine routine, cabable of
  !                performing both the transform and its inverse.
  !                however, as the transforms are not normalized,
  !                the application of a transform followed by its
  !                inverse will yield the original values multiplied
  !                by n.
  !
  !
  ! access       *fortran,p=xlib,sn=cfft99
  !
  !
  ! usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 0,
  !              q .ge. 0, and r .ge. 0.  then a typical sequence of
  !              calls to transform a given set of complex vectors of
  !              length n to a set of (unscaled) complex fourier
  !              coefficient vectors of length n is
  !
  !                   dimension ifax(13),trigs(2*n)
  !                   complex a(...), work(...)
  !
  !                   call cftfax (n, ifax, trigs)
  !                   call cfft99 (a,work,trigs,ifax,inc,jump,n,lot,isgn)
  !
  !              the output vectors overwrite the input vectors, and
  !              these are stored in a.  with appropriate choices for
  !              the other arguments, these vectors may be considered
  !              either the rows or the columns of the array a.
  !              see the individual write-ups for cftfax and
  !              cfft99 below, for a detailed description of the
  !              arguments.
  !
  ! history      the package was written by clive temperton at ecmwf in
  !              november, 1978.  it was modified, documented, and tested
  !              for ncar by russ rew in september, 1980.  it was
  !              further modified for the fully complex case by dave
  !              fulker in november, 1980.
  !
  !-----------------------------------------------------------------------
  !
  ! subroutine cftfax (n,ifax,trigs)
  !
  ! purpose      a set-up routine for cfft99.  it need only be
  !              called once before a sequence of calls to cfft99,
  !              provided that n is not changed.
  !
  ! argument     ifax(13),trigs(2*n)
  ! dimensions
  !
  ! arguments
  !
  ! on input     n
  !               an even number greater than 1 that has no prime factor
  !               greater than 5.  n is the length of the transforms (see
  !               the documentation for cfft99 for the definition of
  !               the transforms).
  !
  !              ifax
  !               an integer array.  the number of elements actually used
  !               will depend on the factorization of n.  dimensioning
  !               ifax for 13 suffices for all n less than 1 million.
  !
  !              trigs
  !               a real array of dimension 2*n
  !
  ! on output    ifax
  !               contains the factorization of n.  ifax(1) is the
  !               number of factors, and the factors themselves are stored
  !               in ifax(2),ifax(3),...  if n has any prime factors
  !               greater than 5, ifax(1) is set to -99.
  !
  !              trigs
  !               an array of trigonometric function values subsequently
  !               used by the cft routines.
  !
  !-----------------------------------------------------------------------
  !
  ! subroutine cfft99 (a,work,trigs,ifax,inc,jump,n,lot,isgn)
  !
  ! purpose      perform a number of simultaneous (unnormalized) complex
  !              periodic fourier transforms or corresponding inverse
  !              transforms.  given a set of complex gridpoint
  !              vectors, the package returns a set of
  !              complex fourier coefficient vectors, or vice
  !              versa.  the length of the transforms must be a
  !              number having no prime factors other than
  !              2, 3, and 5.  this routine is
  !              optimized for use on the cray-1.
  !
  ! argument     complex a(n*inc+(lot-1)*jump), work(n*lot)
  ! dimensions   real trigs(2*n), integer ifax(13)
  !
  ! arguments
  !
  ! on input     a
  !               a complex array of length n*inc+(lot-1)*jump containing
  !               the input gridpoint or coefficient vectors.  this array
  !               overwritten by the results.
  !
  !               n.b. although the array a is usually considered to be of
  !               type complex in the calling program, it is treated as
  !               real within the transform package.  this requires that
  !               such type conflicts are permitted in the user"s
  !               environment, and that the storage of complex numbers
  !               matches the assumptions of this routine.  this routine
  !               assumes that the real and imaginary portions of a
  !               complex number occupy adjacent elements of memory.  if
  !               these conditions are not met, the user must treat the
  !               array a as real (and of twice the above length), and
  !               write the calling program to treat the real and
  !               imaginary portions explicitly.
  !
  !              work
  !               a complex work array of length n*lot or a real array
  !               of length 2*n*lot.  see n.b. above.
  !
  !              trigs
  !               an array set up by cftfax, which must be called first.
  !
  !              ifax
  !               an array set up by cftfax, which must be called first.
  !
  !
  !               n.b. in the following arguments, increments are measured
  !               in word pairs, because each complex element is assumed
  !               to occupy an adjacent pair of words in memory.
  !
  !              inc
  !               the increment (in word pairs) between successive element
  !               of each (complex) gridpoint or coefficient vector
  !               (e.g.  inc=1 for consecutively stored data).
  !
  !              jump
  !               the increment (in word pairs) between the first elements
  !               of successive data or coefficient vectors.  on the cray-
  !               try to arrange data so that jump is not a multiple of 8
  !               (to avoid memory bank conflicts).  for clarification of
  !               inc and jump, see the examples below.
  !
  !              n
  !               the length of each transform (see definition of
  !               transforms, below).
  !
  !              lot
  !               the number of transforms to be done simultaneously.
  !
  !              isgn
  !               = -1 for a transform from gridpoint values to fourier
  !                    coefficients.
  !               = +1 for a transform from fourier coefficients to
  !                    gridpoint values.
  !
  ! on output    a
  !               if isgn = -1, and lot gridpoint vectors are supplied,
  !               each containing the complex sequence:
  !
  !               g(0),g(1), ... ,g(n-1)  (n complex values)
  !
  !               then the result consists of lot complex vectors each
  !               containing the corresponding n coefficient values:
  !
  !               c(0),c(1), ... ,c(n-1)  (n complex values)
  !
  !               defined by:
  !                 c(k) = sum(j=0,...,n-1)( g(j)*exp(-2*i*j*k*pi/n) )
  !                 where i = sqrt(-1)
  !
  !
  !               if isgn = +1, and lot coefficient vectors are supplied,
  !               each containing the complex sequence:
  !
  !               c(0),c(1), ... ,c(n-1)  (n complex values)
  !
  !               then the result consists of lot complex vectors each
  !               containing the corresponding n gridpoint values:
  !
  !               g(0),g(1), ... ,g(n-1)  (n complex values)
  !
  !               defined by:
  !                 g(j) = sum(k=0,...,n-1)( g(k)*exp(+2*i*j*k*pi/n) )
  !                 where i = sqrt(-1)
  !
  !
  !               a call with isgn=-1 followed by a call with isgn=+1
  !               (or vice versa) returns the original data, multiplied
  !               by the factor n.
  !
  !
  ! example       given a 64 by 9 grid of complex values, stored in
  !               a 66 by 9 complex array, a, compute the two dimensional
  !               fourier transform of the grid.  from transform theory,
  !               it is known that a two dimensional transform can be
  !               obtained by first transforming the grid along one
  !               direction, then transforming these results along the
  !               orthogonal direction.
  !
  !               complex a(66,9), work(64,9)
  !               real trigs1(128), trigs2(18)
  !               integer ifax1(13), ifax2(13)
  !
  !               set up the ifax and trigs arrays for each direction:
  !
  !               call cftfax(64, ifax1, trigs1)
  !               call cftfax( 9, ifax2, trigs2)
  !
  !               in this case, the complex values of the grid are
  !               stored in memory as follows (using u and v to
  !               denote the real and imaginary components, and
  !               assuming conventional fortran storage):
  !
  !   u(1,1), v(1,1), u(2,1), v(2,1),  ...  u(64,1), v(64,1), 4 nulls,
  !
  !   u(1,2), v(1,2), u(2,2), v(2,2),  ...  u(64,2), v(64,2), 4 nulls,
  !
  !   .       .       .       .         .   .        .        .
  !   .       .       .       .         .   .        .        .
  !   .       .       .       .         .   .        .        .
  !
  !   u(1,9), v(1,9), u(2,9), v(2,9),  ...  u(64,9), v(64,9), 4 nulls.
  !
  !               we choose (arbitrarily) to transorm first along the
  !               direction of the first subscript.  thus we define
  !               the length of the transforms, n, to be 64, the
  !               number of transforms, lot, to be 9, the increment
  !               between elements of each transform, inc, to be 1,
  !               and the increment between the starting points
  !               for each transform, jump, to be 66 (the first
  !               dimension of a).
  !
  !               call cfft99( a, work, trigs1, ifax1, 1, 66, 64, 9, -1)
  !
  !               to transform along the direction of the second subscript
  !               the roles of the increments are reversed.  thus we defin
  !               the length of the transforms, n, to be 9, the
  !               number of transforms, lot, to be 64, the increment
  !               between elements of each transform, inc, to be 66,
  !               and the increment between the starting points
  !               for each transform, jump, to be 1
  !
  !               call cfft99( a, work, trigs2, ifax2, 66, 1, 9, 64, -1)
  !
  !               these two sequential steps results in the two-dimensiona
  !               fourier coefficient array overwriting the input
  !               gridpoint array, a.  the same two steps applied again
  !               with isgn = +1 would result in the reconstruction of
  !               the gridpoint array (multiplied by a factor of 64*9).
  !
  !
  !-----------------------------------------------------------------------
    dimension ar(1),ai(1),work(1),trigs(1),ifax(1)
  !
  !     subroutine "cfft99" - multiple fast complex fourier transform
  !
  !     a is the array containing input and output data
  !     work is an area of size n*lot
  !     trigs is a previously prepared list of trig function values
  !     ifax is a previously prepared list of factors of n
  !     inc is the increment within each data 'vector'
  !         (e.g. inc=1 for consecutively stored data)
  !     jump is the increment between the start of each data vector
  !     n is the length of the data vectors
  !     lot is the number of data vectors
  !     isgn = +1 for transform from spectral to gridpoint
  !           = -1 for transform from gridpoint to spectral
  !
  !
  !     vectorization is achieved on cray by doing the transforms in
  !     parallel.
  !
    nn = n+n
    ink=inc+inc
    jum = jump+jump
    nfax=ifax(1)
    jnk = 2
    jst = 2
    if (isgn.ge.0) go to 30
  !
  !     the innermost temperton routines have no facility for the
  !     forward (isgn = -1) transform.  therefore, the input must be
  !     rearranged as follows:
  !
  !     the order of each input vector,
  !
  !     g(0), g(1), g(2), ... , g(n-2), g(n-1)
  !
  !     is reversed (excluding g(0)) to yield
  !
  !     g(0), g(n-1), g(n-2), ... , g(2), g(1).
  !
  !     within the transform, the corresponding exponential multiplier
  !     is then precisely the conjugate of that for the normal
  !     ordering.  thus the forward (isgn = -1) transform is
  !     accomplished
  !
  !     for nfax odd, the input must be transferred to the work array,
  !     and the rearrangement can be done during the move.
  !
    jnk = -2
    jst = nn-2
    if (mod(nfax,2).eq.1) goto 40
  !
  !     for nfax even, the rearrangement must be applied directly to
  !     the input array.  this can be done by swapping elements.
  !
    ibase = 1
    ilast = (n-1)*inc
    nh = n/2
    do 20 l=1,lot
    i1 = ibase+inc
    i2 = ibase+ilast
    do 10 m=1,nh
  !     swap real and imaginary portions
    hreal = ar(i1)
    himag = ai(i1)
    ar(i1) = ar(i2)
    ai(i1) = ai(i2)
    ar(i2) = hreal
    ai(i2) = himag
    i1 = i1+inc
    i2 = i2-inc
  10 continue
    ibase = ibase+jump
  20 continue
    goto 100
  !
  30 continue
    if (mod(nfax,2).eq.0) goto 100
  !
  40 continue
  !
  !     during the transform process, nfax steps are taken, and the
  !     results are stored alternately in work and in a.  if nfax is
  !     odd, the input data are first moved to work so that the final
  !     result (after nfax steps) is stored in array a.
  !
  !      write(*,*)'Cheng'      

    ibase=1
    jbase=1
    do 60 l=1,lot
  !     move real and imaginary portions of element zero
    work(jbase) = ar(ibase)
    work(jbase+1) = ai(ibase)
    i=ibase+inc
    j=jbase+jst
    do 50 m=2,n
  !     move real and imaginary portions of other elements (possibly in
  !     reverse order, depending on jst and jnk)
    work(j) = ar(i)
    work(j+1) = ai(i)
    i=i+inc
    j=j+jnk
  50 continue
    ibase=ibase+jump
    jbase=jbase+nn
  60 continue
  !
  !      write(*,*)'Yinghua'
  100 continue
  !
  !     perform the transform passes, one pass for each factor.  during
  !     each pass the data are moved from a to work or from work to a.
  !
  !     for nfax even, the first pass moves from a to work
    igo = 110
  !     for nfax odd, the first pass moves from work to a
    if (mod(nfax,2).eq.1) igo = 120
    la=1
    do 140 k=1,nfax
    if (igo.eq.120) go to 120
  110 continue
    call vpassm(ar,ai,work(1),work(2),trigs,inc,2,jump,nn,lot,n,ifax(k+1),la)
    igo=120
    go to 130
  120 continue
    call vpassm(work(1),work(2),ar,ai,trigs,2,inc,nn,jump,lot,n,ifax(k+1),la)
    igo=110
  130 continue
    la=la*ifax(k+1)
  140 continue
  !
  !     at this point the final transform result is stored in a.
  !
    return
  end
  !****************************
  subroutine cftfax(n,ifax,trigs)
    implicit real*8 (a-h,o-z)
    dimension ifax(13),trigs(1)
  !
  !     this routine was modified from temperton"s original
  !     by dave fulker.  it no longer produces factors in ascending
  !     order, and there are none of the original 'mode' options.
  !
  ! on input     n
  !               the length of each complex transform to be performed
  !
  !               n must be greater than 1 and contain no prime
  !               factors greater than 5.
  !
  ! on output    ifax
  !               ifax(1)
  !                 the number of factors chosen or -99 in case of error
  !               ifax(2) thru ifax( ifax(1)+1 )
  !                 the factors of n in the followin order:  appearing
  !                 first are as many factors of 4 as can be obtained.
  !                 subsequent factors are primes, and appear in
  !                 ascending order, except for multiple factors.
  !
  !              trigs
  !               2n sin and cos values for use by the transform routine
  !
    call fact(n,ifax)
    k = ifax(1)
    if (k .lt. 1 .or. ifax(k+1) .gt. 5) ifax(1) = -99
    if (ifax(1) .le. 0 )then
      write(*,1900)n
  1900    format(' fftfax - invalid n=',i20)
      stop 'bad numbers...game over, man!'
      endif
    call cftrig (n, trigs)
    return
  end
  !*******************************
  subroutine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
    implicit real*8 (a-h,o-z)
    dimension a(n),b(n),c(n),d(n),trigs(n)
  !
  !     "vpassm" - multiple version of "vpassa"
  !     performs one pass through data
  !     as part of multiple complex fft routine
  !     a is first real input vector
  !     b is first imaginary input vector
  !     c is first real output vector
  !     d is first imaginary output vector
  !     trigs is precalculated table of sines " cosines
  !     inc1 is addressing increment for a and b
  !     inc2 is addressing increment for c and d
  !     inc3 is addressing increment between a"s & b"s
  !     inc4 is addressing increment between c"s & d"s
  !     lot is the number of vectors
  !     n is length of vectors
  !     ifac is current factor of n
  !     la is product of previous factors
  !
    data sin36/0.587785252292473d0/,cos36/0.809016994374947d0/, &
         sin72/0.951056516295154d0/,cos72/0.309016994374947d0/, &
         sin60/0.866025403784437d0/,hlf/0.5d0/
  !
    m=n/ifac
    iink=m*inc1
    jink=la*inc2
    jump=(ifac-1)*jink
    ibase=0
    jbase=0
    igo=ifac-1
    if (igo.gt.4) return
    go to (10,50,90,130),igo
  !
  !     coding for factor 2
  !
  10 ia=1
    ja=1
    ib=ia+iink
    jb=ja+jink
    do 20 l=1,la
    i=ibase
    j=jbase
    do 15 ijk=1,lot
    c(ja+j)=a(ia+i)+a(ib+i)
    d(ja+j)=b(ia+i)+b(ib+i)
    c(jb+j)=a(ia+i)-a(ib+i)
    d(jb+j)=b(ia+i)-b(ib+i)
    i=i+inc3
    j=j+inc4
  15 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  20 continue
    if (la.eq.m) return
    la1=la+1
    jbase=jbase+jump
    do 40 k=la1,m,la
    kb=k+k-2
    c1=trigs(kb+1)
    s1=trigs(kb+2)
    do 30 l=1,la
    i=ibase
    j=jbase
    do 25 ijk=1,lot
    c(ja+j)=a(ia+i)+a(ib+i)
    d(ja+j)=b(ia+i)+b(ib+i)
    c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
    d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
    i=i+inc3
    j=j+inc4
  25 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  30 continue
    jbase=jbase+jump
  40 continue
    return
  !
  !     coding for factor 3
  !
  50 ia=1
    ja=1
    ib=ia+iink
    jb=ja+jink
    ic=ib+iink
    jc=jb+jink
    do 60 l=1,la
    i=ibase
    j=jbase
    do 55 ijk=1,lot
    c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
    d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
    c(jb+j)=(a(ia+i)-hlf*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
    c(jc+j)=(a(ia+i)-hlf*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
    d(jb+j)=(b(ia+i)-hlf*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
    d(jc+j)=(b(ia+i)-hlf*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
    i=i+inc3
    j=j+inc4
  55 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  60 continue
    if (la.eq.m) return
    la1=la+1
    jbase=jbase+jump
    do 80 k=la1,m,la
    kb=k+k-2
    kc=kb+kb
    c1=trigs(kb+1)
    s1=trigs(kb+2)
    c2=trigs(kc+1)
    s2=trigs(kc+2)
    do 70 l=1,la
    i=ibase
    j=jbase
    do 65 ijk=1,lot
    c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
    d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
    c(jb+j)= & 
        c1*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
       -s1*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
    d(jb+j)= &
        s1*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
       +c1*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
    c(jc+j)= & 
        c2*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
       -s2*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
    d(jc+j)= &
        s2*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
       +c2*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
    i=i+inc3
    j=j+inc4
  65 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  70 continue
    jbase=jbase+jump
  80 continue
    return
  !
  !     coding for factor 4
  !
  90 ia=1
    ja=1
    ib=ia+iink
    jb=ja+jink
    ic=ib+iink
    jc=jb+jink
    id=ic+iink
    jd=jc+jink
    do 100 l=1,la
    i=ibase
    j=jbase
    do 95 ijk=1,lot
    c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
    c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
    d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
    d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
    c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
    c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
    d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
    d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
    i=i+inc3
    j=j+inc4
  95 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  100 continue
    if (la.eq.m) return
    la1=la+1
    jbase=jbase+jump
    do 120 k=la1,m,la
    kb=k+k-2
    kc=kb+kb
    kd=kc+kb
    c1=trigs(kb+1)
    s1=trigs(kb+2)
    c2=trigs(kc+1)
    s2=trigs(kc+2)
    c3=trigs(kd+1)
    s3=trigs(kd+2)
    do 110 l=1,la
    i=ibase
    j=jbase
    do 105 ijk=1,lot
    c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
    d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
    c(jc+j)= & 
        c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) & 
       -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
    d(jc+j)= & 
        s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
       +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
    c(jb+j)= &
        c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
       -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
    d(jb+j)= &
        s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
       +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
    c(jd+j)= &
        c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
       -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
    d(jd+j)= &
        s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
       +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
    i=i+inc3
    j=j+inc4
  105 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  110 continue
    jbase=jbase+jump
  120 continue
    return
  !
  !     coding for factor 5
  !
  130 ia=1
    ja=1
    ib=ia+iink
    jb=ja+jink
    ic=ib+iink
    jc=jb+jink
    id=ic+iink
    jd=jc+jink
    ie=id+iink
    je=jd+jink
    do 140 l=1,la
    i=ibase
    j=jbase
    do 135 ijk=1,lot
    c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
    d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
    c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
    c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
    d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
    d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
    c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
    c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
    d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
    d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) & 
      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
    i=i+inc3
    j=j+inc4
  135 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  140 continue
    if (la.eq.m) return
    la1=la+1
    jbase=jbase+jump
    do 160 k=la1,m,la
    kb=k+k-2
    kc=kb+kb
    kd=kc+kb
    ke=kd+kb
    c1=trigs(kb+1)
    s1=trigs(kb+2)
    c2=trigs(kc+1)
    s2=trigs(kc+2)
    c3=trigs(kd+1)
    s3=trigs(kd+2)
    c4=trigs(ke+1)
    s4=trigs(ke+2)
    do 150 l=1,la
    i=ibase
    j=jbase
    do 145 ijk=1,lot
    c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
    d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
    c(jb+j)= &
        c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
          -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) & 
       -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
          +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
    d(jb+j)= &
        s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
          -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
       +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
          +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
    c(je+j)= &
        c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
          +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
       -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
          -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
    d(je+j)= &
        s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
          +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
       +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
          -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
    c(jc+j)= &
        c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
          -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) & 
       -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
          +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
    d(jc+j)= &
        s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
          -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
       +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
          +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
    c(jd+j)= &
        c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
          +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
       -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
          -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
    d(jd+j)= & 
        s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
          +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
       +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) & 
          -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
    i=i+inc3
    j=j+inc4
  145 continue
    ibase=ibase+inc1
    jbase=jbase+inc2
  150 continue
    jbase=jbase+jump
  160 continue
    return
  end
  !********************************************
  subroutine fact(n,ifax)
    implicit real*8 (a-h,o-z)
  !     factorization routine that first extracts all factors of 4
    dimension ifax(13)
    if (n.gt.1) go to 10
    ifax(1) = 0
    if (n.lt.1) ifax(1) = -99
    return
  10 nn=n
    k=1
  !     test for factors of 4
  20 if (mod(nn,4).ne.0) go to 30
    k=k+1
    ifax(k)=4
    nn=nn/4
    if (nn.eq.1) go to 80
    go to 20
  !     test for extra factor of 2
  30 if (mod(nn,2).ne.0) go to 40
    k=k+1
    ifax(k)=2
    nn=nn/2
    if (nn.eq.1) go to 80
  !     test for factors of 3
  40 if (mod(nn,3).ne.0) go to 50
    k=k+1
    ifax(k)=3
    nn=nn/3
    if (nn.eq.1) go to 80
    go to 40
  !     now find remaining factors
  50 l=5
    max = dsqrt(dfloat(nn))
    inc=2
  !     inc alternately takes on values 2 and 4
  60 if (mod(nn,l).ne.0) go to 70
    k=k+1
    ifax(k)=l
    nn=nn/l
    if (nn.eq.1) go to 80
    go to 60
  70 if (l.gt.max) go to 75
    l=l+inc
    inc=6-inc
    go to 60
  75 k = k+1
    ifax(k) = nn
  80 ifax(1)=k-1
  !     ifax(1) now contains number of factors
    return
  end
  !*************************************
  subroutine cftrig(n,trigs)
    implicit real*8 (a-h,o-z)
    dimension trigs(1)
    pi=4.0d0*datan(1.0d0)
    del=(pi+pi)/dfloat(n)
    l=n+n
    do 10 i=1,l,2
    angle=0.5d0*dfloat(i-1)*del
    trigs(i  )=dcos(angle)
    trigs(i+1)=dsin(angle)
  10 continue
    return
  end
#endif
end module
