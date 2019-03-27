! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2019
! Large sections taken from SCREEN/src/screen_system.f90
!
module prep_system
  use ai_kinds, only : DP
  implicit none
  private


  type atoms
    real( DP ) :: reduced_coord( 3 )
    character( len=2 ) :: el_name
  end type atoms

  type physical_system
    real( DP ) :: avecs( 3, 3 )
    real( DP ) :: celvol
    real( DP ) :: bvecs( 3, 3 )
    integer :: natoms
    type( atoms ), allocatable :: atom_list( : )
  end type physical_system

  type system_parameters
    integer :: brange( 4 )
    integer :: kmesh( 3 )
    integer :: xmesh( 3 )
    integer :: nkpts
    integer :: nspin
    integer :: nbands
    real( DP ) :: kshift( 3 )
    logical :: isSplit !is the DFT run separated into k and k+q?
    logical :: haveShift
  end type system_parameters

  type calculation_parameters
    logical :: tmels
    logical :: makeU2
    logical :: makeCKS

    ! This should almost certainly be thrown over to the dft reader
    character(len=128) :: outdir
    character(len=128) :: shiftOutdir

    character(len=2), allocatable :: nameCKS( : )
    integer, allocatable :: indxCKS( : )
  end type calculation_parameters

  type( physical_system ), save :: psys
  type( system_parameters ), save :: params
  type( calculation_parameters ), save :: calcParams

  public :: physical_system, atoms, system_parameters
  public :: psys, params


  public :: prep_system_load, prep_system_summarize

  contains


  subroutine prep_system_load( ierr )
    use OCEAN_mpi, only : myid, root, comm, nproc, MPI_INTEGER, MPI_SUCCESS
    integer, intent( inout ) :: ierr
    !
#ifdef MPI
    integer :: ierr_
#endif
!    logical :: loud


    if( myid .eq. root ) then
      call load_xyz( ierr )
      if( ierr .ne. 0 ) goto 111

      call load_abvecs( ierr )
      if( ierr .ne. 0 ) goto 111

      call load_params( ierr )
      if( ierr .ne. 0 ) goto 111

      call load_calcParams( psys, ierr )
      if( ierr .ne. 0 ) goto 111
    endif

111 continue

#ifdef MPI
    if( nproc .gt. 1 ) then

      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. MPI_SUCCESS ) return

      call share_xyz( ierr )
      if( ierr .ne. 0 ) return

      call share_abvecs( ierr )
      if( ierr .ne. 0 ) return

      call share_params( ierr )
      if( ierr .ne. 0 ) return

      call share_calcParams( ierr )
      if( ierr .ne. 0 ) return

    endif
#endif
    params%nkpts = product( params%kmesh(:) )
    params%nbands = params%brange(4) - params%brange(3 ) & 
                  + params%brange(2) - params%brange(1) + 2

!    loud = ( root == myid )
!    call reconcile_inputs( loud, ierr )

  end subroutine prep_system_load


  subroutine share_calcParams( ierr )
    use OCEAN_mpi
    integer, intent( inout ) :: ierr
    !
    integer :: nCKS
    if( nproc .eq. 1 ) return
#ifdef MPI
    call MPI_BCAST( calcParams%tmels, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( calcParams%makeU2, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( calcParams%makeCKS, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( calcParams%outdir, 128, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( calcParams%shiftOutdir, 128, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    if( calcParams%makeCKS ) then
      if( myid .eq. 0 ) nCKS = size( calcParams%indxCKS )

      call MPI_BCAST( nCKS, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return

      if( myid .ne. root ) allocate( calcParams%indxCKS( nCKS ), calcParams%nameCKS( nCKS ) )

      call MPI_BCAST( calcParams%indxCKS, nCKS, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return

      call MPI_BCAST( calcParams%nameCKS, 2*nCKS, MPI_CHARACTER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return

    else
      if( myid .ne. root ) allocate( calcParams%indxCKS( 0 ), calcParams%nameCKS( 0 ) )
    endif
#endif
  end subroutine share_calcParams

  subroutine share_params( ierr )
    use OCEAN_mpi
    integer, intent( inout ) :: ierr
    !
    if( nproc .eq. 1 ) return
    !
#ifdef MPI
    call MPI_BCAST( params%kshift, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%brange, 4, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%kmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%xmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%nspin, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%isSplit, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%haveShift, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif


  end subroutine share_params


  subroutine share_abvecs( ierr )
    use OCEAN_mpi
    integer, intent( inout ) :: ierr
    !
    if( nproc .eq. 1 ) return
    !
    call MPI_BCAST( psys%avecs, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( psys%bvecs, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( psys%celvol, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

  end subroutine


  subroutine share_xyz( ierr )
    use OCEAN_mpi
    integer, intent( inout ) :: ierr
    !
    real( DP ), allocatable :: temp_coords( :, : )
    character( len=2 ), allocatable :: temp_names( : )
    integer :: ii
    !
    if( nproc .eq. 1 ) return

    call MPI_BCAST( psys%natoms, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    !
    if( myid .ne. root ) allocate( psys%atom_list( psys%natoms ) )
    allocate( temp_names( psys%natoms ), temp_coords( 3, psys%natoms ) )

    if( myid .eq. root ) then
      do ii = 1, psys%natoms
        temp_names( ii ) = psys%atom_list( ii )%el_name
        temp_coords( :, ii ) = psys%atom_list( ii )%reduced_coord( : )
      enddo
    endif

    call MPI_BCAST( temp_names, 2 * psys%natoms, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( temp_coords, 3 * psys%natoms, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    if( myid .ne. root ) then
      do ii = 1, psys%natoms
        psys%atom_list( ii )%el_name = temp_names( ii )
        psys%atom_list( ii )%reduced_coord( : ) = temp_coords( :, ii )
      enddo
    endif

    deallocate( temp_names, temp_coords )

    return

  end subroutine share_xyz



  subroutine prep_system_summarize( ierr )
    use OCEAN_mpi, only : myid, root
    integer, intent( inout ) :: ierr

    if( myid .eq. root ) then
      write( 6, * ) "-------------------------------"
    endif
  end subroutine prep_system_summarize


  ! NOT MPI SAFE ( in so much as it will let every process hit the filesystem )
  subroutine load_calcParams( psys_, ierr )
    type( physical_system ), intent( in ) :: psys_
    integer, intent( inout ) :: ierr

    integer :: nCKS, i, elnum
    character(len=2) :: elname
    character(len=64) :: prefix, workDir
    logical :: ex

    inquire( file='prep.tmels', exist=ex )
    if( ex ) then
      open( unit=99, file='prep.tmels', form='formatted', status='old' )
      read( 99, * ) calcParams%tmels
      close( 99 )
    else
      calcParams%tmels = .true.
    endif

    inquire( file='prep.u2', exist=ex )
    if( ex ) then
      open( unit=99, file='prep.u2', form='formatted', status='old' )
      read( 99, * ) calcParams%makeU2
      close( 99 )
    else
      calcParams%makeU2 = .true.
    endif

    inquire( file='prep.cks', exist=ex )
    if( ex ) then
      open( unit=99, file='prep.cks', form='formatted', status='old' )
      read( 99, * ) nCKS
      if( nCKS .gt. 0 ) then
        calcParams%makeCKS = .true.
        allocate( calcParams%nameCKS( nCKS ), calcParams%indxCKS( nCKS ) )
        do i = 1, nCKS
          read( 99, * ) calcParams%nameCKS( i ), calcParams%indxCKS( i )
        enddo
      else
        calcParams%makeCKS = .false.
        allocate( calcParams%nameCKS( 0 ), calcParams%indxCKS( 0 ) )
      endif
      close( 99 )
    else
      calcParams%makeCKS = .false.
      allocate( calcParams%nameCKS( 0 ), calcParams%indxCKS( 0 ) )
    endif

    inquire( file='work_dir', exist=ex )
    if( ex ) then
      open( unit=99, file='work_dir', form='formatted', status='old' )
      read( 99, * ) workDir
      close( 99 )
    else
      workDir = './Out'
    endif

    inquire( file='prefix', exist=ex )
    if( ex ) then
      open( unit=99, file='prefix', form='formatted', status='old' )
      read( 99, * ) prefix
      close( 99 )
    else
      prefix = 'system'
    endif

    write( calcParams%outdir, '(A,A,A,A,A)' ) trim(workDir), '/', trim(prefix), '.save', '/'
    write( calcParams%shiftOutdir, '(A,A,A,A,A)' ) trim(workDir), '/', trim(prefix), '_shift.save', '/'


  end subroutine load_calcParams

  ! NOT MPI SAFE ( in so much as it will let every process hit the filesystem )
  subroutine load_params( ierr )
    integer, intent( inout ) :: ierr
    !
    real(dp), parameter :: tol = 0.000000001
    integer :: brange( 4 )
    logical :: ex
    !
    open( unit=99, file='kmesh.ipt', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open kmesh.ipt', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) params%kmesh( : )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read kmesh.ipt', ierr
      return
    endif
    close( 99, IOSTAT=ierr)
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close kmesh.ipt', ierr
      return
    endif
    params%nkpts = product( params%kmesh(:) )

    open( unit=99, file='nspin', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open nspin', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) params%nspin
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read nspin', ierr
      return
    endif
    close( 99, IOSTAT=ierr)
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close nspin', ierr
      return
    endif

    inquire( file='k0.ipt', exist=ex ) 
    if( ex ) then
      open( unit=99, file='k0.ipt', form='formatted', status='old', IOSTAT=ierr )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to open k0.ipt', ierr
        return
      endif
      read( 99, *, IOSTAT=ierr ) params%kshift( : )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to read k0.ipt', ierr
        return
      endif
      close( 99, IOSTAT=ierr)
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to close k0.ipt', ierr
        return
      endif
    else
      params%kshift( 1 ) = 0.125_DP
      params%kshift( 2 ) = 0.25_DP
      params%kshift( 3 ) = 0.375_DP
    endif


    if( abs( params%kshift( 1 ) ) + abs( params%kshift( 2 ) ) + abs( params%kshift( 3 ) ) &
        < tol ) then
      params%haveShift = .false.
    else
      params%haveShift = .true.
    endif

    open( unit=99, file='brange.ipt', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open kmesh.ipt', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) params%brange( : )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read bands.ipt', ierr
      return
    endif
    close( 99, IOSTAT=ierr)
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close bands.ipt', ierr
      return
    endif

    inquire( file='dft.split', exist=ex )
    if( ex ) then
      open( unit=99, file='dft.split', form='formatted', status='old', IOSTAT=ierr )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to open dft.split', ierr
        return
      endif
      read( 99, *, IOSTAT=ierr ) params%isSplit
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to read dft.split', ierr
        return
      endif
      close( 99, IOSTAT=ierr)
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to close dft.split', ierr
        return
      endif
    else
      params%isSplit = .false.
    endif

    open(unit=99, file='xmesh.ipt', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open xmesh.ipt', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) params%xmesh(:)
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read xmesh.ipt', ierr
      return
    endif
    close( 99, IOSTAT=ierr)
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close xmesh.ipt', ierr
      return
    endif


  end subroutine load_params

  ! NOT MPI SAFE ( in so much as it will let every process hit the filesystem )
  subroutine load_abvecs( ierr )
    integer, intent( inout ) :: ierr
    !
    open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open avecsinbohr.ipt', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) psys%avecs( :, : )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read avecsinbohr.ipt', ierr
      return
    endif
    close( 99, IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close avecsinbohr.ipt', ierr
      return
    endif
    call getomega( psys%avecs, psys%celvol )

    open( unit=99, file='bvecs', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open bvecs', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) psys%bvecs( :, : )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read bvecs', ierr
      return
    endif
    close( 99, IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close bvecs', ierr
      return
    endif

  end subroutine load_abvecs

  ! NOT MPI SAFE ( in so much as it will let every process hit the filesystem )
  subroutine load_xyz( ierr )
    integer, intent( inout ) :: ierr
    !
    integer :: ii
    logical :: ex

    inquire( file='xyz.wyck', exist=ex )
    if( ex .eqv. .false. ) then
      psys%natoms = 0
      allocate( psys%atom_list( 0 ) )
      return
    endif

    open( unit=99, file='xyz.wyck', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open xyz.wyck', ierr
      goto 111
    endif
    ii = 0
    read( 99, *, IOSTAT=ierr, ERR=10 ) psys%natoms
    if( psys%natoms .gt. 0 ) then
      allocate( psys%atom_list( psys%natoms ), STAT=ierr )
    else
      write( 6, * ) 'FATAL ERROR: Negative number of atoms in xyz.wyck', psys%natoms
      ierr = -1
    endif
    if( ierr .ne. 0 ) goto 111

    do ii = 1, psys%natoms
      read( 99, *, IOSTAT=ierr, ERR=10 ) psys%atom_list( ii )%el_name, &
                                         psys%atom_list( ii )%reduced_coord( : )
    enddo
    close( 99, IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Trouble closing xyz.wyck'
      goto 111
    endif

10    continue
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Trouble reading xyz.wyck', ii
      goto 111
    endif

111 continue
    return
  end subroutine load_xyz

  subroutine getomega( avec, omega )
    !
    real(DP), intent( in ) :: avec( 3, 3 )
    real(DP), intent( out ) :: omega
    integer :: i, j, k
    !
    omega = 0.0_DP
    do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     omega = omega + avec( i, 1 ) * avec( j, 2 ) * avec( k, 3 )
     omega = omega - avec( i, 1 ) * avec( k, 2 ) * avec( j, 3 )
    end do
    omega = abs( omega )
    !
    return
  end subroutine getomega

end module prep_system
