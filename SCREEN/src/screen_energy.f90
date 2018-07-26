! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 04-2017
!
!
module screen_energy
  use ai_kinds, only : DP

  implicit none
  private
  save

  ! Energies are the same for all the sites
  !   Should be able to duplicate across all MPI procs
  !   1MB = 131k bands @ 1 k-point, 16k bands @ 2x2x2
  real( DP ), protected, allocatable :: energies( :, :, : )
  real( DP ), protected :: mu_ryd    ! midpoint between HOMO/LUMO of states used in screening
  real( DP ), protected :: efermi    ! input eFermi from DFT calc
  real( DP ), protected :: mindiff   ! min distance between \mu and HOMO or LUMO
  real( DP ), protected :: maxdiff   ! max distance between lowest occupied/highest unocc and \mu
  real( DP ), protected :: geodiff   ! sqrt( mindiff * maxdiff )

  integer :: nspin
  integer :: nbands
  integer :: nkpts

  public :: energies, mu_ryd, geodiff

  public :: screen_energy_init
  public :: screen_energy_kill
  public :: screen_energy_find_fermi
  public :: screen_energy_load

  contains


  ! To replicate current routines find highest energy below read-in Fermi level and
  !   lowest energy above read-in Fermi level. Set `efermi' to be between these two
  !   to make it as far as possible from any poles. 
  subroutine screen_energy_find_fermi( ierr )
    use ocean_mpi, only : myid, root
    use ocean_constants, only : Rydberg2eV
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: vlryd, vhryd, clryd, chryd
    integer :: padder, ispn, ik, ib
    !
    call screen_energy_read_fermi( ierr )
    if( ierr .ne. 0 ) return
    !
    ! Can't be bothered to be clever here.
    ! Only crawling through a few 1000s, probably fastest to repeat for everyone

    ! Initialize valence to lowest band, conduction to highest
    vlryd = energies( 1, 1, 1 )
    vhryd = energies( 1, 1, 1 )
    clryd = energies( nbands, 1, 1 )
    chryd = energies( nbands, 1, 1 )
    !
!    padder = min( 16, nbv )
    do ispn = 1, nspin
      do ik = 1, nkpts
        do ib = 1,  min( 16, nbands )
          if( energies( ib, ik, ispn ) .lt. vlryd ) vlryd = energies( ib, ik, ispn )
        enddo
        do ib = 1, nbands !max( 1, nbv - padder ), min( nbands, nbv + padder )
          if( energies( ib, ik, ispn ) .gt. efermi ) then
            if( energies( ib, ik, ispn ) .lt. clryd ) clryd = energies( ib, ik, ispn )
          else
            if( energies( ib, ik, ispn ) .gt. vhryd ) vhryd = energies( ib, ik, ispn )
          endif
        enddo
        do ib = max( 1, nbands - 16 ), nbands
          if( energies( ib, ik, ispn ) .gt. chryd ) chryd = energies( ib, ik, ispn )
        enddo
      enddo
    enddo

    ! If we have a degeneracy at the Fermi level then there will be problems
    if( abs( clryd - vhryd ) .lt. 0.000001_DP ) then
      ierr = -1
      if( myid .eq. root ) write( 6, * ) 'ERROR! Degeneracy at Fermi level'
      return
    endif

    mu_ryd = ( clryd + vhryd ) / 2.0_DP
    mindiff = min( mu_ryd - vhryd, clryd - mu_ryd )
    maxdiff = max( mu_ryd - vlryd, chryd - mu_ryd )
    geodiff = sqrt( mindiff * maxdiff )

    if( myid .eq. root ) then
      write( 6, '(A,1F12.6)' ) 'Input Fermi energy is: ', efermi*Rydberg2eV
      write( 6, '(A,1F12.6)' ) 'Mid-gap/Chemical pot.: ', mu_ryd*Rydberg2eV
      write( 6, '(A)' ) "    #### Energy summary (eV) ####"
      write( 6, '(A)' ) " Valence minimum   Valence maximum   Conduction min.   Conduction max."
      !                 " -X.XXXXXXXXE+YY   -X.XXXXXXXXE+YY   -X.XXXXXXXXE+YY   -X.XXXXXXXXE+YY
      write( 6, '(4(1x,1e15.8,2x))' ) vlryd*Rydberg2eV, vhryd*Rydberg2eV, clryd*Rydberg2eV, chryd*Rydberg2eV
      write( 6, '(A)' ) "   Fermi/midgap    Min. difference   Max. difference  Geo. mean of diffs"
      !                 " -X.XXXXXXXXE+YY   -X.XXXXXXXXE+YY   -X.XXXXXXXXE+YY   -X.XXXXXXXXE+YY
      write ( 6, '(4(1x,1e15.8,2x))' ) mu_ryd*Rydberg2eV, mindiff*Rydberg2eV, maxdiff*Rydberg2eV, geodiff*Rydberg2eV
      write( 6, '(A)' ) "    #############################"
    endif
    
  end subroutine screen_energy_find_fermi


  subroutine screen_energy_read_fermi( ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER, MPI_SUCCESS, MPI_DOUBLE_PRECISION
    integer, intent( inout ) :: ierr
#ifdef MPI
    integer :: ierr_
#endif
    !
    if( myid .eq. root ) then
      open( unit=99, file='efermiinrydberg.ipt', form='formatted', status='old', iostat=ierr, err=10 )
      read( 99, *, iostat=ierr ) efermi
      close( 99 )
    endif
10  continue
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 .or. ierr_ .ne. MPI_SUCCESS ) then
      if( ierr .eq. 0 ) ierr = ierr_
      return
    endif
    call MPI_BCAST( efermi, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif
  end subroutine screen_energy_read_fermi


  subroutine screen_energy_init( ierr )
    use screen_system, only : system_parameters, params
    !
    integer, intent( inout ) :: ierr 
    !
    nspin  = params%nspin
    nkpts  = product( params%kmesh( : ) )
    nbands = params%nbands
    !
    allocate( energies( nbands, nkpts, nspin ), STAT=ierr )
    !
    ! If we need some clever OMP NUMA then first touch goes here
    energies( :, :, : ) = 0.0_DP
    !
  end subroutine screen_energy_init

  subroutine screen_energy_kill()
    deallocate( energies )
  end subroutine


  subroutine screen_energy_load( ierr )
    use OCEAN_dft_files, only : odf_read_energies_single
    use OCEAN_mpi, only : myid, root, comm
    integer, intent( inout ) :: ierr
    
    call odf_read_energies_single( myid, root, comm, energies, ierr )

  end subroutine screen_energy_load

end module screen_energy
