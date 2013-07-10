module OCEAN_energies
  use AI_kinds

  implicit none
  save

  REAL(DP), ALLOCATABLE, TARGET :: energies(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: psi

  real(DP), ALLOCATABLE, TARGET :: imag_selfenergy(:,:,:)
!DEC$ ATTRIBUTES ALIGN: 32 :: imag_selfenergy

  INTEGER :: energy_bands_pad
  INTEGER :: energy_kpts_pad

  LOGICAL :: have_selfenergy

  
  contains

  subroutine OCEAN_energies_init(  sys, ierr )
    use OCEAN_system

    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_system), intent( in ) :: sys

    integer, parameter :: cacheline_by_Z = 4

    if( mod( sys%num_bands, cacheline_by_Z ) == 0 ) then
      energy_bands_pad = sys%num_bands
    else
      energy_bands_pad =  cacheline_by_Z * ( sys%num_bands / cache_line_by_Z + 1 )
    endif

    if( mod( sys%nkpts, cacheline_by_Z ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      energy_kpts_pad = sys%nkpts
    else
      energy_kpts_pad =  cacheline_by_Z * ( sys%nkpts / cache_line_by_Z + 1 )
    endif


    allocate( energies( psi_bands_pad, psi_kpts_pad, sys%nspn ), STAT=ierr )
    energies = 0_DP


    allocate( imag_selfenergy( psi_bands_pad, psi_kpts_pad, sys%nspn ), STAT=ierr )
    imag_selfenergy = 0_DP
    

  end subroutine OCEAN_energies_init


  subroutine OCEAN_energies_kill( ierr )
    implicit none


    deallocate( energies )
    deallocate( imag_selfenergy )

  end subroutine OCEAN_energies_kill

  subroutine OCEAN_energies_load( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi

    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_system), intent( in ) :: sys


    integer :: nbd, nq, nspn
    character, len=9 :: infoname

    if( sys%conduct ) then
      infoname = 'wvfcninfo'
    else
      infoname = 'wvfvainfo'
    endif

    if( myid .eq. root ) then
      write(6,*) 'Loading energies ...'
      open( unit=99, file=infoname, form='unformatted', status='old' )
      rewind(99)
      if( sys%nspn .eq. 2 ) then
         read( 99 ) nbd, nq, nspn
      else
        read ( 99 ) nbd, nq
        nspn = 1
      endif
      if( nbd .ne. sys%num_bands ) then
        write(6,*) 'Band mismatch!', sys%num_bands, nbd
        ierr = 1
        goto 111
      elseif( nq .ne. sys%nkpts ) then
        write(6,*) 'K-point mismatch', sys%nkpts, nq
        ierr = 1
        goto 111
      elseif( nspn .ne. sys%nspn ) then
        write(6,*) 'Spin mismatch', sys%nspn, nspn
        ierr = 1
        goto 111
      endif

!      allocate( tmp_e0( sys%num_bands * sys%nkpts, sys%nspn ) )
      read(99) energies( 1 : sys%num_bands, 1 : sys%nkpts, : )

      close( 99 )

    endif
! Need to sychronize to test for ierr from above?
#ifdef MPI
    call MPI_BARRIER( comm )

    call MPI_BCAST( energies, energy_bands_pad*energy_kpts_pad*sys%nspn, MPI_DOUBLE, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

  111 continue

end module OCEAN_energies
