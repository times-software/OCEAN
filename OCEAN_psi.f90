module OCEAN_psi
  use AI_kinds

  implicit none
  save
  private

  COMPLEX(DP), ALLOCATABLE, TARGET :: psi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: psi

  COMPLEX(DP), ALLOCATABLE, TARGET :: old_psi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: old_psi
  COMPLEX(DP), ALLOCATABLE, TARGET :: new_psi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: new_psi

  COMPLEX(DP), ALLOCATABLE, TARGET :: lr_psi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: lr_psi
  COMPLEX(DP), ALLOCATABLE, TARGET :: mult_psi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: mult_psi

  COMPLEX(DP), ALLOCATABLE, TARGET :: hpsi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: hpsi


  REAL(DP) :: kpref
  INTEGER :: psi_bands_pad
  INTEGER :: psi_kpts_pad


  public :: OCEAN_psi_init, OCEAN_psi_kill, OCEAN_psi_load, OCEAN_psi_sum_lr,  &
            OCEAN_psi_sum, OCEAN_psi_ab
  contains

  subroutine OCEAN_psi_ab( sys, a, b1, b2, ierr )
    use OCEAN_system
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_system), intent( in ) :: sys
    complex(DP), intent( out ) :: a 
    real( DP ), intent( in ) :: b1
    real(DP), intent( out ) :: b2

    a = dot_product( hpsi(:,:,:), psi(:,:,:) )

    new_psi( :, :, : ) = hpsi( :, :, : ) - real(a) * psi( :, :, : ) - b1 * old_psi( :, :, : )
    b2 = DZNRM2( psi_bands_pad * psi_kpts_pad * sys%nalpha, new_psi, 1 )
    new_psi( :, :, : ) = new_psi( :, :, : ) / b2


    !JTV Need to figure out fortran pointers instead of mem copies, mem copies are dumb
    old_psi( :, :, : ) = psi( :, :, : )
    psi( :, :, : ) = new_psi( :, :, : )

  end subroutine OCEAN_psi_ab

  subroutine OCEAN_psi_sum_lr( sys, ierr )
    use mpi
    use OCEAN_system
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_system), intent( in ) :: sys

#ifdef MPI
! Doing all_reduce here so that the multiplet part can hide the latency maybe
    call MPI_ALLREDUCE( MPI_IN_PLACE, lr_psi, psi_bands_pad * psi_kpts_pad * sys%nalpha, &
                        MPI_DOUBLE_COMLPEX, MPI_SUM, comm, ierr )
#endif

  end subroutine OCEAN_psi_sum_lr


  subroutine OCEAN_psi_sum( sys, ierr )
    use OCEAN_system
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_system), intent( in ) :: sys

    if( sys%long_range .and. sys%mult ) then
      hpsi(:,:,:) = hpsi(:,:,:) + lr_psi(:,:,:) + mult_psi(:,:,:)
    elseif( sys%long_range ) 
      hpsi(:,:,:) = hpsi(:,:,:) + lr_psi(:,:,:)
    elseif( sys%mult )
      hpsi(:,:,:) = hpsi(:,:,:) + mult_psi(:,:,:)
    endif

  end subroutine OCEAN_psi_sum


  subroutine OCEAN_psi_init( sys, ierr )
    use OCEAN_system

    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_system), intent( in ) :: sys
  
    integer, parameter :: cacheline_by_Z = 4
    
    if( mod( sys%num_bands, cacheline_by_Z ) == 0 ) then
      psi_bands_pad = sys%num_bands
    else
      psi_bands_pad =  cacheline_by_Z * ( sys%num_bands / cache_line_by_Z + 1 )
    endif

    if( mod( sys%nkpts, cacheline_by_Z ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      psi_kpts_pad = sys%nkpts
    else 
      psi_kpts_pad =  cacheline_by_Z * ( sys%nkpts / cache_line_by_Z + 1 )
    endif


    allocate( psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
              old_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
              new_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
              hpsi( psi_bands_pad, psi_kpts_pad, sys%nalpha ),  &
              STAT=ierr )
    if(sys%lr) allocate( lr_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ) )
    if(sys%mult) allocate( mult_psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ) )

    psi = 0_DP

  end subroutine OCEAN_psi_init


  subroutine OCEAN_psi_kill( ierr )
    implicit none 
    integer, intent(inout) :: ierr

    deallocate( psi )
  end subroutine OCEAN_psi_kill( ierr )


  subroutine OCEAN_psi_load( sys, ierr )
    use OCEAN_mpi 
    use mpi

    implicit none

    type(OCEAN_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: ialpha, icms, icml, ivms, ikpt, iband
    real(DP) :: val, nrm, tmpr, tmpi
    real(DP) :: tau( 3 )

    if( myid .eq. root ) then
      write(6,*) 'Reading in projector coefficients'

      ialpha = 0
      if( sys%nspn == 1 ) then
        do icms = -1, 1, 2
          do icml = -lc, lc
            write ( str, '(1a4,1i2.2)' ) 'beff', 1 + icml + lc
            open( unit=99, file=str, form='unformatted', status='old' )
            rewind 99
            read ( 99 ) tau( : )
            do ivms = -1, 1, 2 
              ialpha = ialpha + 1
              write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
              if ( icms .eq. ivms ) then
                do ikpt = 1, sys%nkpts
                  do iband = 1, sys%num_bands
                    read ( 99 ) tmpr, tmpi
                    psi( iband, ikpt, ialpha ) = cmplx( tmpr, tmpi )
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
            ivms = 1, 2
              write ( str, '(1a4,1i2.2,1a1,1i1.1)' ) 'beff', 1 + icml + lc, '.', ivms
              open( unit=99, file=str, form='unformatted', status='old' )
              rewind 99
              read ( 99 ) tau( : )
              ialpha = ialpha + 1
              write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
              if ( icms .eq. ivms ) then
                do ikpt = 1, sys%nkpts
                  do iband = 1, sys%num_bands
                    read ( 99 ) tmpr, tmpi
                    psi( iband, ikpt, ialpha ) = cmplx( tmpr, tmpi )
                  enddo 
                end do  
              end if
              close( unit=99 )
            end do
          end do
        end do
      endif

      write (6,*) 'band states have been read in'


      val = 0_DP
      do ialpha = 1, sys%nalpha
#ifdef BLAS
        nrm = DZNRM2( psi(1,1,ialpha), psi_bands_pad*psi_kpts_pad, 1 )
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm**2
        val = val + nrm
#else
        nrm = 0_DP
        do ikpt = 1, sys%nkpts
          nrm = nrm + sum( psi( :, ikpt, ialpha ) * conjg( psi( :, ikpt, ialpha ) ) )
        enddo
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
        val = val + sqrt( nrm )
#endif
      enddo
      psi = psi / val
      kpref = 4_DP * pi * val**2 / (dble(sys%nkpts) * sys%celvol ** 2 ) 
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) kpref
      close( unit=99 )
    endif
        
#ifdef MPI
    call MPI_BCAST( kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( psi, psi_bands_pad*psi_kpts_pad*sys%nalpha, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif



  111 continue

  end subroutine OCEAN_psi_load
