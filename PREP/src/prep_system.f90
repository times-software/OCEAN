! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2019
! Large sections taken from SCREEN/src/screen_system.'formatted'0
!
module prep_system
  use ai_kinds, only : DP
  implicit none
  private


  type atoms
    real( DP ) :: reduced_coord( 3 )
    real( DP ) :: xcoord( 3 )
    character( len=2 ) :: el_name
    integer :: indx
    integer :: z
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
    real( DP ) :: k0( 3 )
    real( DP ) :: qshift( 3 )
    logical :: isSplit !is the DFT run separated into k and k+q?
    logical :: haveShift
  end type system_parameters

  type calculation_parameters
    logical :: legacyFiles

    ! Need to standardize these and then acccess them from the other routines
    logical :: makeTmels
    logical :: makeU2
    logical :: makeCKS

    ! This should almost certainly be thrown over to the dft reader
    character(len=128) :: outdir
    character(len=128) :: shiftOutdir

    ! I'm not sure about these 
    character(len=2), allocatable :: nameCKS( : )
    integer, allocatable :: indxCKS( : )
  end type calculation_parameters

  type( physical_system ), save :: psys
  type( system_parameters ), save :: params
  type( calculation_parameters ), save :: calcParams

  public :: physical_system, atoms, system_parameters, calculation_parameters
  public :: psys, params, calcParams


  public :: prep_system_load, prep_system_summarize, prep_system_snatch

  public ::  prep_system_ikpt2kvec, prep_system_umklapp

  contains

  subroutine prep_system_ikpt2kvec( ikpt, addShift, kvec, kvecCart )
    use ai_kinds, only : QP
    integer, intent( in ) :: ikpt
    logical, intent( in ) :: addShift
    real(DP), intent( out ) :: kvec(3), kvecCart(3)
    !
    integer :: ikx, iky, ikz, ik, i

    ikx = ( ikpt - 1 ) / ( params%kmesh( 3 ) * params%kmesh( 2 ) )
    ik = ikpt - ikx * params%kmesh( 3 ) * params%kmesh( 2 ) 

    iky = ( ik - 1 ) / params%kmesh( 3 )
    ikz = ik - iky * params%kmesh( 3 ) - 1

    
    kvec( 1 ) = ( params%k0( 1 ) + real( ikx, QP ) ) / real( params%kmesh( 1 ), QP )
    kvec( 2 ) = ( params%k0( 2 ) + real( iky, QP ) ) / real( params%kmesh( 2 ), QP )
    kvec( 3 ) = ( params%k0( 3 ) + real( ikz, QP ) ) / real( params%kmesh( 3 ), QP )

!    if( addShift ) kvec(:) = kvec(:) + params%qshift(:)
    ! Moving the occupied to the shifted grid by shifting them back with respect to the unocc
    if( .not. addShift ) kvec(:) = kvec(:) - params%qshift(:)

    kvecCart(:) =  0
    do i = 1, 3
      kvecCart( : ) = kvecCart( : ) + psys%bvecs( :, i ) * kvec( i )
    enddo

  end subroutine prep_system_ikpt2kvec

  subroutine prep_system_umklapp( ikpt, addshift, umklapp )
    integer, intent( in ) :: ikpt
    logical, intent( in ) :: addShift
    integer, intent( out ) :: umklapp( 3 )

    real(DP) :: kvec(3), kvecCart(3)
    integer :: i

    call prep_system_ikpt2kvec( ikpt, addShift, kvec, kvecCart )

    umklapp(:) = 0
    do i = 1, 3

      do while( kvec(i) .lt. -1.0_DP )
        kvec(i) = kvec(i) + 1.0_DP
        umklapp(i) = umklapp(i) - 1
      enddo

      do while( kvec(i) .gt. 1.0_DP )
        kvec(i) = kvec(i) - 1.0_DP
        umklapp(i) = umklapp(i) + 1
      enddo
    enddo

  end subroutine prep_system_umklapp

  subroutine prep_system_load( ierr )
    use OCEAN_mpi, only : myid, root, comm, nproc, MPI_INTEGER, MPI_SUCCESS
    integer, intent( inout ) :: ierr
    !
#ifdef MPI
    integer :: ierr_
#endif
!    logical :: loud


    if( myid .eq. root ) then
      call load_abvecs( ierr )
      if( ierr .ne. 0 ) goto 111

      call load_atomList( ierr )
!      call load_xyz( ierr )
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
    call MPI_BCAST( calcParams%makeTmels, 1, MPI_LOGICAL, root, comm, ierr )
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
    call MPI_BCAST( params%k0, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
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

    call MPI_BCAST( params%qshift, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
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
    real( DP ), allocatable :: temp_reduced( :, : ), temp_xcoord( :, : )
    character( len=2 ), allocatable :: temp_names( : )
    integer, allocatable :: zeeAndIndx( :, : )
    integer :: ii
    !
    if( nproc .eq. 1 ) return

    call MPI_BCAST( psys%natoms, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    !
    if( myid .ne. root ) allocate( psys%atom_list( psys%natoms ) )
    allocate( temp_names( psys%natoms ), temp_reduced( 3, psys%natoms ), &
              zeeAndIndx( 2, psys%natoms ), temp_xcoord( 3, psys%natoms ) )

    if( myid .eq. root ) then
      do ii = 1, psys%natoms
        temp_names( ii ) = psys%atom_list( ii )%el_name
        temp_reduced( :, ii ) = psys%atom_list( ii )%reduced_coord( : )
        temp_xcoord( :, ii ) = psys%atom_list( ii )%xcoord( : )
        zeeAndIndx( 1, ii ) = psys%atom_list( ii )%z
        zeeAndIndx( 2, ii ) = psys%atom_list( ii )%indx
      enddo
    endif

    call MPI_BCAST( temp_names, 2 * psys%natoms, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( temp_reduced, 3 * psys%natoms, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( temp_xcoord, 3 * psys%natoms, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( zeeAndIndx, 2 * psys%natoms, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    if( myid .ne. root ) then
      do ii = 1, psys%natoms
        psys%atom_list( ii )%el_name = temp_names( ii )
        psys%atom_list( ii )%reduced_coord( : ) = temp_reduced( :, ii )
        psys%atom_list( ii )%xcoord( : ) = temp_xcoord( :, ii )
        psys%atom_list( ii )%z = zeeAndIndx( 1, ii )
        psys%atom_list( ii )%indx = zeeAndIndx( 2, ii )
        write(1000+myid,*) psys%atom_list( ii )%el_name, psys%atom_list( ii )%z, psys%atom_list( ii )%indx
        write(1000+myid,*) psys%atom_list( ii )%reduced_coord( : )
        write(1000+myid,*) psys%atom_list( ii )%xcoord( : )
      enddo
    endif

    deallocate( temp_names, temp_reduced, zeeAndIndx )

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
      read( 99, * ) calcParams%makeTmels
      close( 99 )
    else
      calcParams%makeTmels = .true.
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
      read( 99, *, IOSTAT=ierr ) params%k0( : )
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
      ierr = 124712
      return
      params%k0( 1 ) = 0.125_DP
      params%k0( 2 ) = 0.25_DP
      params%k0( 3 ) = 0.375_DP
    endif


    inquire( file='qinunitsofbvectors.ipt', exist=ex )
    if( ex ) then
      open( unit=99, file='qinunitsofbvectors.ipt', form='formatted', status='old', IOSTAT=ierr )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to open qinunitsofbvectors.ipt', ierr
        return
      endif
      read( 99, *, IOSTAT=ierr ) params%qshift( : )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to read qinunitsofbvectors.ipt', ierr
        return
      endif
      close( 99, IOSTAT=ierr)
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Failed to close qinunitsofbvectors.ipt', ierr
        return
      endif
    else
      params%qshift( : ) = 0.0_DP
    endif


    if( abs( params%qshift( 1 ) ) + abs( params%qshift( 2 ) ) + abs( params%qshift( 3 ) ) &
        < tol ) then
      params%haveShift = .false.
    else
      params%haveShift = .true.
    endif

    write(6,*) params%qshift( : ), params%haveShift

    open( unit=99, file='brange.ipt', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open brange.ipt', ierr
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
    use ocean_phys, only : ophys_getBvecs
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

    call ophys_getBvecs( psys%avecs, psys%bvecs, psys%celvol )

#if 0
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
#endif

  end subroutine load_abvecs

  ! NOT MPI SAFE ( in so much as it will let every process hit the filesystem )
  ! must be called after abvec are loaded
  subroutine load_atomList( ierr )
    use periodic, only : getsymbol_underscore 
    use ocean_phys, only : ophys_fixCoords
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: xinput(:,:), xred(:,:), xcoord(:,:)
    integer :: ntypat, natom, i
    integer, allocatable :: znucl(:), typat(:), tempIndx(:)
    character(len=8) :: coord

    open(unit=99,file='ntype',form='formatted',status='old')
    read(99,*) ntypat
    close(99)
    !
    open(unit=99,file='natoms',form='formatted',status='old')
    read(99,*) natom
    close(99)
    !
    allocate(znucl(ntypat),typat(natom),xinput(3,natom))

    open(unit=99,file='znucl',form='formatted',status='old')
    read(99,*) znucl(:)
    close(99)

    open(unit=99,file='typat',form='formatted',status='old')
    read(99,*) typat(:)
    close(99)

    open(unit=99,file='taulist',form='formatted',status='old')
    read(99,*) xinput(:,:)
    close(99)

    open(unit=99,file='coord',form='formatted',status='old')
    read(99,*) coord
    close(99)

    allocate( xred( 3,natom), xcoord(3,natom) )
    call ophys_fixCoords( psys%avecs, coord, xinput, xred, xcoord, ierr )
    if( ierr .ne. 0 ) return

    ! must be less than 160 elements
    allocate( tempIndx( 160 ) )
    tempIndx(:) = 0

    psys%natoms = natom
    allocate( psys%atom_list( natom ) )

    do i = 1, psys%natoms
      psys%atom_list( i )%reduced_coord( : ) = xred( :, i )
      psys%atom_list( i )%xcoord( : ) = xcoord( :, i )
!      psys%atom_list( i )%el_name = elements(znucl(typat( i ) ) )
      call getsymbol_underscore( znucl(typat( i ) ), psys%atom_list( i )%el_name )
      psys%atom_list( i )%z = znucl( typat( i ) )
      tempIndx( psys%atom_list( i )%z ) = tempIndx( psys%atom_list( i )%z ) + 1
      psys%atom_list( i )%indx = tempIndx( psys%atom_list( i )%z )
    enddo

    deallocate( tempIndx, znucl, typat, xred, xcoord, xinput )

  end subroutine load_atomList
  
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

    subroutine prep_system_snatch( element, indx, tau, xcoord, ierr )
    character( len=2 ), intent( in ) :: element
    integer, intent( in ) :: indx
    real( DP ), intent( out ) :: tau( 3 )
    real( DP ), intent( out ) :: xcoord( 3 )
    integer, intent( inout ) :: ierr
    !
    integer :: ii, nmatch

    nmatch = 0
    do ii = 1, psys%natoms
      if( element .eq. psys%atom_list( ii )%el_name ) then
        nmatch = nmatch + 1
        if( nmatch .eq. indx ) then
          tau( : ) = psys%atom_list( ii )%reduced_coord(: )
          xcoord( : ) = psys%atom_list( ii )%xcoord(: )
          goto 111
        endif
      endif
    enddo

    write( 6, * ) 'Atom coord not found!'
    ierr = -1
    return

111 continue

!    xcoord = matmul( psys%avecs, tau )

  end subroutine prep_system_snatch

end module prep_system
