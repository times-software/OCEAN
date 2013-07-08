module OCEAN_psi
  use AI_kinds

  implicit none
  save

  INTEGER :: psi_bands_pad
  INTEGER :: psi_kpts_pad
  REAL(DP) :: kpref

  COMPLEX(DP), ALLOCATABLE, TARGET :: psi(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: psi


  contains


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

    if( mod( sys%nkpts, cacheline_by_Z ) == 0 ) then
      psi_kpts_pad = sys%nkpts
    else 
      psi_kpts_pad =  cacheline_by_Z * ( sys%nkpts / cache_line_by_Z + 1 )
    endif


    allocate( psi( psi_bands_pad, psi_kpts_pad, sys%nalpha ), STAT=ierr )
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

    integer :: ialpha
    real(DP) :: val, nrm

    if( myid .eq. root ) then
      write(6,*) 'Reading in projector coefficients'

      


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

    

#endif






  111 continue

  end subroutine OCEAN_psi_load
