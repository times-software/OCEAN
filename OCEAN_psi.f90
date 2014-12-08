module OCEAN_psi
  use AI_kinds
  use iso_c_binding

  implicit none
  save
  private

!  COMPLEX(DP),  POINTER :: psi(:,:,:)
!!DEC$ ATTRIBUTES ALIGN: 32 :: psi

!  COMPLEX(DP), ALLOCATABLE, TARGET :: old_psi(:,:,:)
!!DEC$ ATTRIBUTES ALIGN: 32 :: old_psi
!  COMPLEX(DP), ALLOCATABLE, TARGET :: new_psi(:,:,:)
!!DEC$ ATTRIBUTES ALIGN: 32 :: new_psi

!  COMPLEX(DP), ALLOCATABLE, TARGET :: lr_psi(:,:,:)
!!DEC$ ATTRIBUTES ALIGN: 32 :: lr_psi
!  COMPLEX(DP), ALLOCATABLE, TARGET :: mult_psi(:,:,:)
!!DEC$ ATTRIBUTES ALIGN: 32 :: mult_psi

!  COMPLEX(DP), ALLOCATABLE, TARGET :: hpsi(:,:,:)
!!DEC$ ATTRIBUTES ALIGN: 32 :: hpsi


  REAL(DP), public :: kpref
  INTEGER :: psi_bands_pad
  INTEGER :: psi_kpts_pad
  INTEGER :: psi_length
  INTEGER :: psi_val_bands
  INTEGER :: psi_con_bands
  INTEGER :: psi_val_length

  LOGICAL :: core = .true.
  LOGICAL :: valence = .false.


  type OCEAN_vector
    REAL(DP), POINTER :: r(:,:,:) => null()
    REAL(DP), POINTER :: i(:,:,:) => null()

    COMPLEX(DP), POINTER :: val(:,:,:) => null()

#ifdef CONTIGUOUS
    CONTIGUOUS :: r, i, val
#endif

#ifdef __INTEL
!dir$ attributes align:64 :: r, i, val
#endif

    TYPE(C_PTR) :: rcptr
    TYPE(C_PTR) :: icptr
    
    INTEGER :: bands_pad
    INTEGER :: kpts_pad
  end type


  public :: OCEAN_psi_init, OCEAN_psi_kill, OCEAN_psi_load, OCEAN_psi_sum_lr,  &
            OCEAN_psi_sum, OCEAN_psi_ab, OCEAN_psi_set_prec, OCEAN_psi_write, &
            OCEAN_psi_swap, OCEAN_psi_dot, OCEAN_psi_nrm, OCEAN_psi_scal, &
            OCEAN_psi_axpy

  public :: OCEAN_vector
  contains

  subroutine OCEAN_psi_swap( a, b, c )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a, b, c
    type( OCEAN_vector ) :: d

    call OCEAN_psi_assign( a, d )
    call OCEAN_psi_assign( b, a )
    call OCEAN_psi_assign( c, b )
    call OCEAN_psi_assign( d, c )

  end subroutine

  subroutine OCEAN_psi_assign( a, b )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a, b

    b%r => a%r
    b%i => a%i
    
    b%bands_pad = a%bands_pad
    b%kpts_pad = b%kpts_pad

  end subroutine OCEAN_psi_assign

  real(dp) function OCEAN_psi_nrm( a )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a
    real(dp) :: val
    real(dp), external :: DDOT, DZNRM2

    if( core ) then
      OCEAN_psi_nrm = DDOT( psi_length, a%r, 1, a%r, 1 )
      OCEAN_psi_nrm = DDOT( psi_length, a%i, 1, a%i, 1 ) + OCEAN_psi_nrm
    else
      OCEAN_psi_nrm = 0.0_dp
    endif

    if( valence ) then
      val = DZNRM2( psi_val_length, a%val, 1 )
      OCEAN_psi_nrm = OCEAN_psi_nrm + val * val
    endif

    OCEAN_psi_nrm = sqrt( OCEAN_psi_nrm )

  end function OCEAN_psi_nrm

  complex(dp) function OCEAN_psi_dot( a, b )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a, b
    real(dp) :: r, i
    real(dp), external :: DDOT
    compleX(dp), external :: ZDOTC

    if( core ) then
      r = DDOT( psi_length, a%r, 1, b%r, 1 ) &
        + DDOT( psi_length, a%i, 1, b%i, 1 )
      i = DDOT( psi_length, a%r, 1, b%i, 1 ) &
        - DDOT( psi_length, a%i, 1, b%r, 1 )
    else
      r = 0.0_dp
      i = 0.0_dp
    endif

    if( valence ) then
      OCEAN_psi_dot = ZDOTC( psi_val_length, a%val, 1, b%val, 1 ) + cmplx( r, i )
    else
      OCEAN_psi_dot = cmplx( r, i )
    endif

  end function OCEAN_psi_dot

  subroutine OCEAN_psi_scal( a, x )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( inout ) :: x


    if( core ) then
      call DSCAL( psi_length, a, x%r, 1 )
      call DSCAL( psi_length, a, x%i, 1 )
    endif

    if( valence ) then
      call ZDSCAL( psi_val_length, a, x%val, 1 )
    endif

  end subroutine OCEAN_psi_scal

  subroutine OCEAN_psi_axpy( a, x, y )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( in ) :: x
    type( OCEAN_vector ), intent( inout ) :: y
    !
    complex(dp) :: za

    if( core ) then
      call DAXPY( psi_length, a, x%r, 1, y%r, 1 )
      call DAXPY( psi_length, a, x%i, 1, y%i, 1 )
    endif

    if( valence ) then
      za = a
      call ZAXPY( psi_val_length, za, x%val, 1, y%val )
    endif


  end subroutine OCEAN_psi_axpy
  

  subroutine OCEAN_psi_set_prec( sys, energy, gprc, psi_in, psi_out )
    use OCEAN_system
    implicit none
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: energy, gprc
    type(OCEAN_vector), intent(in) :: psi_in
    type(OCEAN_vector), intent(inout) :: psi_out
    !
    real( DP ) :: gprc_sqd, denom
    integer :: ialpha, ikpt, ibnd

    gprc_sqd = gprc * gprc
!        pcdiv( i ) = ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )
    do ialpha = 1, sys%nalpha
      do ikpt = 1, sys%nkpts
        do ibnd = 1, sys%num_bands
          denom = ( energy - psi_in%r( ibnd, ikpt, ialpha ) ) ** 2 & 
                + gprc_sqd + psi_in%i( ibnd, ikpt, ialpha ) ** 2
          psi_out%r( ibnd, ikpt, ialpha ) = ( energy - psi_in%r( ibnd, ikpt, ialpha ) ) / denom
          psi_out%i( ibnd, ikpt, ialpha ) = -( psi_in%i( ibnd, ikpt, ialpha ) + gprc ) / denom
        enddo
      enddo
    enddo

  end subroutine OCEAN_psi_set_prec

  subroutine OCEAN_psi_ab( sys, a, b1, b2, imag_a, psi, hpsi, old_psi, ierr )
    use OCEAN_system
    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    real(DP), intent( out ) :: a 
    real(DP), intent( out ) :: imag_a 
    real( DP ), intent( in ) :: b1
    real(DP), intent( out ) :: b2
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    
    type(OCEAN_vector) :: temp_psi

    complex(DP) :: temp_a
    real(DP), external :: DZNRM2
    complex(DP), external :: ZDOTC

    integer :: ialpha, ikpt

    temp_a = cmplx( 0.0_DP )
    a = 0.0_DP
    imag_a = 0.0_DP
    do ialpha = 1, sys%nalpha
      do ikpt = 1, sys%nkpts
!        temp_a = temp_a + dot_product( hpsi(:,ikpt,ialpha), psi(:,ikpt,ialpha) )
        a = a + sum( hpsi%r(:,ikpt,ialpha) * psi%r(:,ikpt,ialpha) ) &
                   + sum(hpsi%i(:,ikpt,ialpha) * psi%i(:,ikpt,ialpha) )
        imag_a = imag_a - sum( hpsi%r(:,ikpt,ialpha) * psi%i(:,ikpt,ialpha) &
                          - hpsi%i(:,ikpt,ialpha) * psi%r(:,ikpt,ialpha ) )
      enddo
    enddo
!    temp_a = ZDOTC( psi_bands_pad * psi_kpts_pad * sys%nalpha, psi, 1, hpsi, 1 )
!    a = real(temp_a)
!    imag_a = aimag(temp_a)

    hpsi%r( :, :, : ) = hpsi%r( :, :, : ) - a * psi%r( :, :, : ) - b1 * old_psi%r( :, :, : )
    hpsi%i( :, :, : ) = hpsi%i( :, :, : ) - a * psi%i( :, :, : ) - b1 * old_psi%i( :, :, : )
!#ifdef BLAS2
!    b2 = DZNRM2( psi_bands_pad * psi_kpts_pad * sys%nalpha, new_psi, 1 )
!#else
    b2 = 0.0_DP
    do ialpha = 1, sys%nalpha
      do ikpt = 1, sys%nkpts
        b2 = b2 + sum( hpsi%r( :, ikpt, ialpha )**2 + hpsi%i( :, ikpt, ialpha )**2 )
      enddo
    enddo
    b2 = sqrt( b2 )
!#endif
    hpsi%r( :, :, : ) = hpsi%r( :, :, : ) / b2
    hpsi%i( :, :, : ) = hpsi%i( :, :, : ) / b2


    !JTV Need to figure out fortran pointers instead of mem copies, mem copies are dumb
!    temp_psi%r => old_psi%r
!    temp_psi%i => old_psi%i
!    old_psi%r => psi%r
!    old_psi%i => psi%i
!    psi%r => hpsi%r
!    psi%i => hpsi%i
!    hpsi%r => temp_psi%r
!    hpsi%i => temp_psi%i
    old_psi%r( :, :, : ) = psi%r( :, :, : )
    old_psi%i( :, :, : ) = psi%i( :, :, : )
    psi%r( :, :, : ) = hpsi%r( :, :, : )
    psi%i( :, :, : ) = hpsi%i( :, :, : )

  end subroutine OCEAN_psi_ab

  subroutine OCEAN_psi_sum_lr( sys, p, ierr ) 
    use mpi
    use OCEAN_system
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: p

#ifdef MPI
! Doing all_reduce here so that the multiplet part can hide the latency maybe
!    call MPI_ALLREDUCE( MPI_IN_PLACE, lr_psi, psi_bands_pad * psi_kpts_pad * sys%nalpha, &
!                        MPI_DOUBLE_COMLPEX, MPI_SUM, comm, ierr )
      call MPI_ALLREDUCE( MPI_IN_PLACE, p%r, p%bands_pad*p%kpts_pad*sys%nalpha, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
      call MPI_ALLREDUCE( MPI_IN_PLACE, p%i, p%bands_pad*p%kpts_pad*sys%nalpha, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
#endif

  end subroutine OCEAN_psi_sum_lr


  subroutine OCEAN_psi_sum( sys, hpsi, p, q, ierr )
    use OCEAN_system
    use ocean_mpi
    use mpi
    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys

    type(OCEAN_vector), intent(inout) :: p, q, hpsi

!    if( sys%long_range .and. sys%mult ) then
!      hpsi(:,:,:) = hpsi(:,:,:) + lr_psi(:,:,:) + mult_psi(:,:,:)
!    elseif( sys%long_range ) then
!      hpsi(:,:,:) = hpsi(:,:,:) + lr_psi(:,:,:)
!    elseif( sys%mult ) then
!      hpsi(:,:,:) = hpsi(:,:,:) + mult_psi(:,:,:)
!    endif

!    hpsi%r(:,:,:) = hpsi%r(:,:,:) + p%r(:,:,:) - q%r(:,:,:)
!    hpsi%i(:,:,:) = hpsi%i(:,:,:) + p%i(:,:,:) - q%i(:,:,:)

! OMP?
    p%r(:,:,:) = p%r(:,:,:) - q%r(:,:,:)
    p%i(:,:,:) = p%i(:,:,:) - q%i(:,:,:)

#ifdef MPI
    call MPI_ALLREDUCE( MPI_IN_PLACE, p%r, p%bands_pad*p%kpts_pad*sys%nalpha, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
    call MPI_ALLREDUCE( MPI_IN_PLACE, p%i, p%bands_pad*p%kpts_pad*sys%nalpha, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
#endif

    hpsi%r(:,:,:) = hpsi%r(:,:,:) + p%r(:,:,:)
    hpsi%i(:,:,:) = hpsi%i(:,:,:) + p%i(:,:,:)


  end subroutine OCEAN_psi_sum


  subroutine OCEAN_psi_init( sys, p, ierr )
    use OCEAN_system
!    use OCEAN_mpi, only : myid, root
    implicit none
    include 'fftw3.f03'
    
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( out ) :: p
  
    integer, parameter :: cacheline_by_Z = 1
!    type(C_PTR) :: cptr

    if( associated( p%r ) ) then
      call fftw_free( p%rcptr )
      p%r => null()
    endif
    if( associated( p%i ) ) then
      call fftw_free( p%icptr )
      p%i => null()
    endif
    
    if( mod( sys%num_bands, cacheline_by_Z ) == 0 ) then
      p%bands_pad = sys%num_bands
    else
      p%bands_pad =  cacheline_by_Z * ( sys%num_bands / cacheline_by_Z + 1 )
    endif

    if( mod( sys%nkpts, cacheline_by_Z ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      p%kpts_pad = sys%nkpts
    else 
      p%kpts_pad =  cacheline_by_Z * ( sys%nkpts / cacheline_by_Z + 1 )
    endif

    


    p%rcptr = fftw_alloc_real( int(p%bands_pad * p%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( p%rcptr, p%r, [p%bands_pad, p%kpts_pad, sys%nalpha ] )
    p%icptr = fftw_alloc_real( int(p%bands_pad * p%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( p%icptr, p%i, [p%bands_pad, p%kpts_pad, sys%nalpha ] )
!    allocate( p%r(p%bands_pad, p%kpts_pad, sys%nalpha), &
!              p%i(p%bands_pad, p%kpts_pad, sys%nalpha ) )

    p%r(:,:,:) = 0.0_DP
    p%i(:,:,:) = 0.0_DP

    
!    allocate(  & !psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
!              old_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
!              new_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
!              hpsi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
!              STAT=ierr )
!    if(sys%long_range) allocate( lr_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ) )
!    if(sys%mult) allocate( mult_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ) )

!    psi = 0_DP
!    if( myid .eq. root ) write(6,*) associated( p%r ), associated( p%i )



    psi_bands_pad = p%bands_pad
    psi_kpts_pad = p%kpts_pad
    psi_length = psi_bands_pad * psi_kpts_pad * sys%nalpha
    psi_val_bands = sys%val_bands
    psi_con_bands = sys%num_bands
    psi_val_length = sys%val_bands * sys%num_bands * sys%nkpts
    core = .true.
    valence = .false.


  end subroutine OCEAN_psi_init


  subroutine OCEAN_psi_kill( p, ierr )
    implicit none 
    include 'fftw3.f03'
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent(inout ) :: p

!    deallocate( psi )
    if( associated( p%r ) ) then
      call fftw_free( p%rcptr )
      p%r => null()
    endif 
    if( associated( p%i ) ) then
      call fftw_free( p%icptr )
      p%i => null()
    endif

    p%bands_pad = 0
    p%kpts_pad = 0
    
  end subroutine OCEAN_psi_kill


  subroutine OCEAN_psi_load_old( sys, p, ierr )
    use OCEAN_mpi 
    use OCEAN_system
    use mpi
    use AI_kinds

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: ialpha, icms, icml, ivms, ikpt, iband, lc
    real(DP) :: val, nrm, tmpr, tmpi, pi
    real(DP) :: tau( 3 )
    character(LEN=8) :: str 
    real(DP), external :: DZNRM2
    lc = sys%ZNL(3)
    pi = 4.0_DP * ATAN( 1.0_DP )

    if( myid .eq. root ) then
      write(6,*) 'Reading in projector coefficients'
      write(6,*) lc

      ialpha = 0
      if( sys%nspn == 1 ) then
        do icms = -1, 1, 2
          do icml = -lc, lc
            write ( str, '(1a4,1i2.2)' ) 'beff', 1 + icml + lc
            open( unit=99, file=str, form='unformatted', status='old' )
            write(6,*) str
            rewind 99
            read ( 99 ) tau( : )
            do ivms = -1, 1, 2 
              ialpha = ialpha + 1
              write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
              if ( icms .eq. ivms ) then
                do ikpt = 1, sys%nkpts
                  do iband = 1, sys%num_bands
                    read( 99 ) p%r(iband,ikpt,ialpha), p%i(iband,ikpt,ialpha)
!                    read ( 99 ) tmpr, tmpi
!                    psi( iband, ikpt, ialpha ) = cmplx( tmpr, tmpi )
                  enddo
                end do  
              end if
            end do
            close( unit=99 )
          end do
        end do
      else
        do icms = 1, 2
          do icml = -lc, lc
            do ivms = 1, 2
              write ( str, '(1a4,1i2.2,1a1,1i1.1)' ) 'beff', 1 + icml + lc, '.', ivms
              open( unit=99, file=str, form='unformatted', status='old' )
              rewind 99
              read ( 99 ) tau( : )
              ialpha = ialpha + 1
              write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
              if ( icms .eq. ivms ) then
                do ikpt = 1, sys%nkpts
                  do iband = 1, sys%num_bands
                    read( 99 ) p%r(iband,ikpt,ialpha), p%i(iband,ikpt,ialpha)
!                    read ( 99 ) tmpr, tmpi
!                    psi( iband, ikpt, ialpha ) = cmplx( tmpr, tmpi )
                  enddo 
                end do  
              end if
              close( unit=99 )
            end do
          end do
        end do
      endif

      write (6,*) 'band states have been read in'


      val = 0.0_DP
      do ialpha = 1, sys%nalpha
! !#ifdef BLAS2
!        nrm = DZNRM2( psi(1,1,ialpha), psi_bands_pad*psi_kpts_pad, 1 )
!        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm**2
!        val = val + nrm
! !#else
        nrm = 0.0_DP
        do ikpt = 1, sys%nkpts
!          nrm = nrm + dot_product( psi( 1 : sys%num_bands, ikpt, ialpha ), psi( 1 : sys%num_bands, ikpt, ialpha ) )
          nrm = nrm + sum(p%r(:,ikpt,ialpha)**2 + p%i(:,ikpt,ialpha)**2)
        enddo
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
        val = val +  nrm 
! !#endif
      enddo
      val = sqrt(val)
      kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 ) 
      val = 1.0_DP / val
      !psi = psi / val
      p%r = p%r * val
      p%i = p%i * val
      write(6,*) pi, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) kpref
      close( unit=99 )
    endif
        
#ifdef MPI
    if( myid .eq. root ) write(6,*) p%bands_pad, p%kpts_pad, sys%nalpha
    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, p%bands_pad*p%kpts_pad*sys%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, p%bands_pad*p%kpts_pad*sys%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif



  111 continue

  end subroutine OCEAN_psi_load_old


  subroutine OCEAN_psi_load( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use mpi
    use AI_kinds

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: ialpha, icms, icml, ivms, ikpt, iband, lc
    real(DP) :: val, nrm, tmpr, tmpi, pi
    real(DP) :: tau( 3 )
    character(LEN=8) :: str
    real(DP), external :: DZNRM2
    lc = sys%ZNL(3)
    pi = 4.0_DP * ATAN( 1.0_DP )

    if( myid .eq. root ) then

      write(6,*) 'Reading in projector coefficients'
      call OCEAN_psi_dotter( sys, p, ierr )
      if( ierr .ne. 0 ) goto 111
      write (6,*) 'band states have been read in',  sys%nalpha

    
      val = 0.0_DP
      do ialpha = 1, sys%nalpha
        nrm = 0.0_DP
        do ikpt = 1, sys%nkpts 
          nrm = nrm + sum(p%r(:,ikpt,ialpha)**2 + p%i(:,ikpt,ialpha)**2)
        enddo
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
        val = val +  nrm
      enddo
      val = sqrt(val)
      kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
      val = 1.0_DP / val
      p%r = p%r * val
      p%i = p%i * val 
      write(6,*) pi, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) kpref
      close( unit=99 )
    endif
        
#ifdef MPI
    if( myid .eq. root ) write(6,*) p%bands_pad, p%kpts_pad, sys%nalpha
    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, p%bands_pad*p%kpts_pad*sys%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, p%bands_pad*p%kpts_pad*sys%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

  111 continue

  end subroutine OCEAN_psi_load

  subroutine OCEAN_psi_write( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use mpi
    use AI_kinds
    
    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    integer, intent( inout ) :: ierr

    character( LEN = 17 ) :: rhs_filename

    complex(DP), allocatable :: out_vec(:,:,:)


    if( myid .ne. root ) return

    allocate( out_vec(sys%num_bands, sys%nkpts, sys%nalpha ) )
    out_vec(:,:,:) = cmplx( p%r(1:sys%num_bands,1:sys%nkpts,:), p%i(1:sys%num_bands,1:sys%nkpts,:) )

    write(rhs_filename,'(A4,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'rhs_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', '1s', '_', sys%cur_run%photon
    open(unit=99,file=rhs_filename,form='unformatted',status='unknown')
    rewind( 99 )
    write( 99 ) out_vec
    close( 99 )

    deallocate(out_vec) 

  end subroutine OCEAN_psi_write



  subroutine OCEAN_psi_dotter( sys, p, ierr )
    use OCEAN_system
 
    implicit none
 
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    real(DP) :: tau( 3 ), rr, ri, ir, ii
    real(DP), allocatable, dimension(:,:) :: pcr, pci
    real(DP), allocatable, dimension(:,:) :: mer, mei
    integer :: nptot, ntot, ialpha, icms, ivms, icml, ikpt, iband, iter

    character (LEN=127) :: cks_filename
    character (LEN=5) :: cks_prefix
    character (LEN=18) :: mel_filename

    !select case (runtype)
    select case ( sys%cur_run%calc_type)
    case( 'XES' )
      cks_prefix = 'cksv.'
    case( 'XAS' )
      cks_prefix = 'cksc.'
    case default
      cks_prefix = 'cksc.'
    end select
    write(cks_filename,'(A5,A2,I4.4)' ) cks_prefix, sys%cur_run%elname, sys%cur_run%indx

    open(unit=99,file=cks_filename,form='unformatted',status='old')
    rewind( 99 )
    read ( 99 ) nptot, ntot
    read ( 99 ) tau( : )
    allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
    read ( 99 ) pcr
    read ( 99 ) pci
    close( unit=99 )



    allocate( mer( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3) ),  &
              mei( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3) ) )

    write(mel_filename,'(A5,A1,I3.3,A1,I2.2,A1,I2.2,A1,I2.2)' ) 'mels.', 'z', sys%cur_run%ZNL(1), &
            'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3), 'p', sys%cur_run%photon
    open( unit=99, file=mel_filename, form='formatted', status='old' )
    rewind( 99 )
    do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
      do iter = 1, nptot
        read( 99, * ) mer( iter, icml ), mei( iter, icml )
      enddo
    enddo 
    close( 99 )

    ialpha = 0
    if( sys%nspn == 1 ) then
      do icms = -1, 1, 2
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = -1, 1, 2
            ialpha = ialpha + 1
            if( icms .eq. ivms ) then
              iter = 0
              do ikpt = 1, sys%nkpts
                do iband = 1, sys%num_bands
                  iter = iter + 1
                  rr = dot_product( mer( :, icml ), pcr( :, iter ) )
                  ri = dot_product( mer( :, icml ), pci( :, iter ) )
                  ir = dot_product( mei( :, icml ), pcr( :, iter ) )
                  ii = dot_product( mei( :, icml ), pci( :, iter ) )
                  p%r(iband,ikpt,ialpha) = rr - ii
                  p%i(iband,ikpt,ialpha) = -ri - ir
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    else
      ierr = -1
      return
    endif

    deallocate( pcr, pci, mer, mei )
    

  end subroutine OCEAN_psi_dotter
end module OCEAN_psi
