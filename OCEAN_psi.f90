module OCEAN_psi
  use AI_kinds
  use iso_c_binding

  implicit none
  save
  private

  REAL(DP), public :: kpref

  INTEGER :: psi_bands_pad
  INTEGER :: psi_kpts_pad
  INTEGER :: psi_core_alpha
  INTEGER :: psi_core_beta
  INTEGER :: psi_val_alpha
  INTEGER :: psi_val_beta
  INTEGER :: psi_val_bands

  INTEGER :: CACHE_DOUBLE = 8

  LOGICAL :: is_init = .false.

  INTEGER, PARAMETER, PUBLIC :: core_vector = 1
  INTEGER, PARAMETER, PUBLIC :: val_vector = 2

  type OCEAN_vector
    REAL(DP), ALLOCATABLE :: r(:,:,:,:,:) 
    REAL(DP), ALLOCATABLE :: i(:,:,:,:,:) 

#ifdef CONTIGUOUS
    CONTIGUOUS :: r, i
#endif

#ifdef __INTEL
!dir$ attributes align:64 :: r, i, val
#endif

    INTEGER :: nband
    INTEGER :: vband ! valence or 1
    INTEGER :: nkpts
    INTEGER :: alpha
    INTEGER :: beta  ! always 1 for right now
    INTEGER :: full_size
    INTEGER :: async_size
    INTEGER :: my_type

    ! MPI requests for sharing OCEAN_vector
    INTEGER, ALLOCATABLE :: r_requests(:)
    INTEGER, ALLOCATABLE :: i_requests(:)

  end type

  type OCEAN_pvector
    TYPE(OCEAN_vector), POINTER :: p => null
  end type OCEAN_pvector


  public :: OCEAN_psi_init, OCEAN_psi_kill, OCEAN_psi_load, OCEAN_psi_sum_lr,  &
            OCEAN_psi_sum, OCEAN_psi_set_prec, OCEAN_psi_write, &
            OCEAN_psi_swap, OCEAN_psi_dot, OCEAN_psi_nrm, OCEAN_psi_scal, &
            OCEAN_psi_axpy, OCEAN_psi_new

  public :: OCEAN_vector, OCEAN_pvector
  contains

  subroutine OCEAN_psi_swap( a, b, c )
    implicit none
    type( OCEAN_pvector ), intent( inout ) :: a, b, c
    type( OCEAN_pvector ) :: d

    d%p => a%p
    a%p => b%p
    b%p => c%p
    c%p => d%p

  end subroutine

  real(dp) function OCEAN_psi_nrm( a )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a
    real(dp), external :: DDOT

    OCEAN_psi_nrm = DDOT( psi_length, a%r, 1, a%r, 1 )
    OCEAN_psi_nrm = DDOT( psi_length, a%i, 1, a%i, 1 ) + OCEAN_psi_nrm

    OCEAN_psi_nrm = sqrt( OCEAN_psi_nrm )

  end function OCEAN_psi_nrm

  complex(dp) function OCEAN_psi_dot( a, b )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a, b
    real(dp) :: r, i
    real(dp), external :: DDOT

    if( a%my_type .ne. b%my_type ) then
      OCEAN_psi_dot = cmplx( 0.0_dp )
      return
    endif

    r = DDOT( a%full_size, a%, a%r, 1, b%r, 1 ) &
      + DDOT( a%full_size, a%i, 1, b%i, 1 )
    i = DDOT( a%full_size, a%r, 1, b%i, 1 ) &
      - DDOT( a%full_size, a%i, 1, b%r, 1 )

    OCEAN_psi_dot = cmplx( r, i )

  end function OCEAN_psi_dot

  subroutine OCEAN_psi_scal( a, x )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( inout ) :: x


    call DSCAL( x%full_size, a, x%r, 1 )
    call DSCAL( x%full_size, a, x%i, 1 )

  end subroutine OCEAN_psi_scal

  subroutine OCEAN_psi_axpy( a, x, y )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( in ) :: x
    type( OCEAN_vector ), intent( inout ) :: y

    if( x%my_type .ne. y%my_type ) then
!      ierr = -1
      return
    endif

    call DAXPY( x%full_size, a, x%r, 1, y%r, 1 )
    call DAXPY( x%full_size, a, x%i, 1, y%i, 1 )

  end subroutine OCEAN_psi_axpy
  

  subroutine OCEAN_psi_set_prec( energy, gprc, psi_in, psi_out )
    implicit none
    real( DP ), intent( in ) :: energy, gprc
    type(OCEAN_vector), intent(in) :: psi_in
    type(OCEAN_vector), intent(inout) :: psi_out
    !
    real( DP ) :: gprc_sqd, denom
    integer :: ibeta, ialpha, ikpt, ibnd1, ibnd2

!        pcdiv( i ) = ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )
    if( psi_in%my_type .ne. psi_out%my_type ) then
      ierr = -1
      return
    endif

    gprc_sqd = gprc * gprc

    do ibeta = 1, psi_in%nbeta
      do ialpha = 1, psi_in%nalpha
        do ikpt = 1, psi_in%nkpts
          do ibnd2 = 1, psi_in%vband
            do ibnd1 = 1, psi_in%nband
              denom = ( energy - psi_in%r( ibnd1, ibnd2, ikpt, ialpha, ibeta ) ) ** 2 & 
                    + gprc_sqd + psi_in%i( ibnd1, ibnd2, ikpt, ialpha, ibeta ) ** 2
              psi_out%r( ibnd1, ibnd2, ikpt, ialpha, ibeta ) = ( energy - psi_in%r(ibnd1,ibnd2,ikpt,ialpha,ibeta) ) / denom
              psi_out%i( ibnd1, ibnd2, ikpt, ialpha, ibeta ) = -( psi_in%i(ibnd1,ibnd2,ikpt,ialpha,ibeta) + gprc ) / denom
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine OCEAN_psi_set_prec

#ifdef FALSE
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
#endif

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


  subroutine OCEAN_psi_sum( hpsi, p, q, ierr )
    use ocean_mpi
    use mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent(in) :: q
    type(OCEAN_vector), intent(inout) :: p, hpsi

    integer :: ialpha, ibeta, ireq

    if( hpsi%my_type .ne. p%my_type .or. hpsi%my_type .ne. q%my_type ) then
      ierr = -1
      return
    endif

    ireq = 0
    do ibeta = 1, q%nbeta
      do ialpha = 1, q%nalpha
        ireq = ireq + 1
        p%r(:,:,:,ialpha,ibeta) = p%r(:,:,:,ialpha,ibeta) - q%r(:,:,:,ialpha,ibeta)
        p%i(:,:,:,ialpha,ibeta) = p%i(:,:,:,ialpha,ibeta) - q%i(:,:,:,ialpha,ibeta)

#ifdef MPI
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%r, p%bands_pad*p%kpts_pad*sys%nalpha, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%r_request(ireq), ierr )
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%i, p%bands_pad*p%kpts_pad*sys%nalpha, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%i_request(ireq), ierr )
#endif
      enddo
    enddo

    ireq = 0
    do ibeta = 1, q%nbeta
      do ialpha = 1, q%nalpha
        ireq = ireq + 1
#ifdef MPI
        call MPI_WAIT( p%r_request( ireq ), MPI_STATUS_IGNORE, ierr )
#endif
        hpsi%r(:,:,:,ialpha,ibeta) = hpsi%r(:,:,:,ialpha,ibeta) + p%r(:,:,:,ialpha,ibeta)
#ifdef MPI
        call MPI_WAIT( p%i_request( ireq ), MPI_STATUS_IGNORE, ierr )
#endif
        hpsi%i(:,:,:,ialpha,ibeta) = hpsi%i(:,:,:,ialpha,ibeta) + p%i(:,:,:,ialpha,ibeta)
      enddo
    enddo

  end subroutine OCEAN_psi_sum


  subroutine OCEAN_psi_init( sys, ierr )
    use OCEAN_system
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    if( mod( sys%num_bands, CACHE_DOUBLE ) == 0 ) then
      psi_bands_pad = sys%num_bands
    else
      psi_bands_pad = CACHE_DOUBLE * ( sys_num_bands / CACHE_DOUBLE + 1 )
    endif

    if( sys%val_bands .le. 1 ) then
      psi_val_bands = 1
    elseif( mod( sys%val_bands, CACHE_DOUBLE ) == 0 ) then
      psi_val_bands = sys%val_bands
    else
      psi_val_bands = CACHE_DOUBLE * ( sys%val_bands / CACHE_DOUBLE + 1 )
    endif

    if( mod( sys%nkpts, CACHE_DOUBLE ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      psi_kpts_pad = sys%nkpts
    else
      psi_kpts_pad =  CACHE_DOUBLE * ( sys%nkpts / CACHE_DOUBLE + 1 )
    endif

    psi_core_alpha = sys%nalpha
    psi_core_beta = 1
 
    psi_val_alpha = sys%nspn ** 2
    psi_val_beta = 1

    is_init = .true.

  end subroutine OCEAN_psi_init

  subroutine OCEAN_psi_new( p, mytype, ierr, q )
    use OCEAN_system
    implicit none
    
    integer, intent(inout) :: ierr
    integer, intent( in ) :: my_type
    type(OCEAN_vector), intent( out ) :: p
    type(OCEAN_vector), intent(in), optional :: q

    if( .not. is_init ) then
      ierr = -1
      return
    endif
  
    p%my_type = mytype
    
    select case( p%my_type )

    case( core_vector )
      p%nband = psi_bands_pad
      p%vband = 1
      p%nkpts = psi_kpts_pad
      p%alpha = psi_core_alpha
      p%beta  = psi_core_beta
    case( val_vector )
      p%nband = psi_bands_pad
      p%vband = psi_val_bands
      p%nkpts = psi_kpts_pad
      p%alpha = psi_val_alpha
      p%beta  = psi_val_beta
    case default
      ierr = -1
      return
    end select


    p%full_size = p%nband*p%vband*p%nkpts*p%alpha*p%beta
    p%async_size = p%nband*p%vband*p%nkpts

    allocate( p%r(p%nband,p%vband,p%nkpts,p%alpha,p%beta), &
              p%i(p%nband,p%vband,p%nkpts,p%alpha,p%beta), &
              p%r_request(p%alpha*p%beta), p%i_request(p%alpha*p%beta), STAT=ierr )

    ! initialize
    if( present( q ) ) then
      if( q%my_type .ne. p%my_type ) then
        ierr = -2
        return
      endif
      p%r = q%r
      p%i = q%i
    else
      p%r = 0.0_DP
      p%i = 0.0_DP
    endif

  end subroutine OCEAN_psi_new



  subroutine OCEAN_psi_kill( p, ierr )
    implicit none 
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

!    deallocate( psi )
    if( allocated( p%r ) deallocate( p%r, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%i ) deallocate( p%i, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%r_requests ) deallocate( p%r_requests, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%i_requests ) deallocate( p%i_requests, STAT=ierr )
    
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
                    read( 99 ) p%r(iband,1,ikpt,ialpha,1), p%i(iband,1,ikpt,ialpha,1)
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
                    read( 99 ) p%r(iband,1,ikpt,ialpha,1), p%i(iband,1,ikpt,ialpha,1)
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
          nrm = nrm + sum(p%r(:,1,ikpt,ialpha,1)**2 + p%i(:,1,ikpt,ialpha,1)**2)
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

    call MPI_BCAST( p%r, p%full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, p%full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
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
          nrm = nrm + sum(p%r(:,1,ikpt,ialpha,1)**2 + p%i(:,1,ikpt,ialpha,1)**2)
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

    call MPI_BCAST( p%r, p%full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, p%full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
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
    out_vec(:,:,:) = cmplx( p%r(1:sys%num_bands,1,1:sys%nkpts,:,1), p%i(1:sys%num_bands,1,1:sys%nkpts,:,1) )

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
    real(DP), allocatable, dimension(:,:,:) :: mer, mei
    integer :: nptot, ntot, ialpha, icms, ivms, icml, ikpt, iband, iter, is

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



    allocate( mer( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3), 2 ),  &
              mei( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3), 2 ) )

    write(mel_filename,'(A5,A1,I3.3,A1,I2.2,A1,I2.2,A1,I2.2)' ) 'mels.', 'z', sys%cur_run%ZNL(1), &
            'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3), 'p', sys%cur_run%photon
    open( unit=99, file=mel_filename, form='formatted', status='old' )
    rewind( 99 )
!    do is = 1, 2
      do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
        do iter = 1, nptot
          read( 99, * ) mer( iter, icml, 1 ), mei( iter, icml, 1 )
        enddo
      enddo 
!    enddo
    close( 99 )

    ialpha = 0
    is = 0
    is = 1
    if( sys%nspn == 1 ) then
      do icms = -1, 1, 2
!        is = is + 1
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = -1, 1, 2
            ialpha = ialpha + 1
            if( icms .eq. ivms ) then
              iter = 0
              do ikpt = 1, sys%nkpts
                do iband = 1, sys%num_bands
                  iter = iter + 1
                  rr = dot_product( mer( :, icml, is ), pcr( :, iter ) )
                  ri = dot_product( mer( :, icml, is ), pci( :, iter ) )
                  ir = dot_product( mei( :, icml, is ), pcr( :, iter ) )
                  ii = dot_product( mei( :, icml, is ), pci( :, iter ) )
                  p%r(iband,1,ikpt,ialpha,1) = rr - ii
                  p%i(iband,1,ikpt,ialpha,1) = -ri - ir
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
