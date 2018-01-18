! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2017
!
!
module screen_system
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
    integer :: nkpts
    integer :: nspin
    integer :: nbands
    real( DP ) :: kshift( 3 )
  end type system_parameters

  type calculation_parameters
    character(len=4) :: convertStyle = 'real'
    character(len=4) :: chi0Integrand = 'half'
    character(len=6) :: inversionStyle = 'sinqr'
    integer :: QuadOrder = 16
  end type calculation_parameters

  type( physical_system ), save :: psys
  type( system_parameters ), save :: params
  type( calculation_parameters ), save :: calcParams

  public :: physical_system, atoms, system_parameters
  public :: psys, params
  public :: screen_system_load, screen_system_snatch, screen_system_summarize
  public :: screen_system_returnKvec
  public :: screen_system_convertStyle, screen_system_chi0Integrand
  public :: screen_system_invStyle, screen_system_QuadOrder

  contains 

  pure function screen_system_chi0Integrand() result( integrand )
    character(len=4) :: integrand
    integrand = calcParams%chi0Integrand
  end function screen_system_chi0Integrand

  pure function screen_system_QuadOrder() result( i )
    integer :: i
    i = calcParams%QuadOrder
  end function screen_system_QuadOrder

  pure function screen_system_invStyle() result (is)
    character(len=6) :: is
    is = calcParams%inversionStyle
  end function screen_system_invStyle

  pure function screen_system_convertStyle() result (cs)
    character(len=4) :: cs
    cs = calcParams%convertStyle
  end function screen_system_convertStyle

  pure function screen_system_returnKvec( sp, ikpt ) result( Kvec )
    type( system_parameters ), intent( in ) :: sp
    integer, intent( in ) :: ikpt
    real(DP) :: Kvec(3)
    !
    integer :: i, j

    j = ( ikpt - 1 ) / ( sp%kmesh(2) * sp%kmesh(3) ) 
    i = ikpt - 1 - j * sp%kmesh(2) * sp%kmesh(3)

    Kvec(1) = ( sp%kshift(1) + real( j, DP ) ) / real( sp%kmesh(1), DP )

    Kvec(2) = ( sp%kshift(2) + real( i/sp%kmesh(3), DP) ) / real( sp%kmesh(2), DP )

    Kvec(3) = ( sp%kshift(3) + real( mod( i, sp%kmesh(3) ), DP ) ) / real( sp%kmesh(3), DP )
    
  end function screen_system_returnKvec

  subroutine screen_system_summarize( ierr )
    use OCEAN_mpi, only : myid, root
    integer, intent( inout ) :: ierr
    !
    integer :: ii
  
    if( myid .eq. root ) then
      write( 6, * ) "-------------------------------"
      write( 6, * ) "N atoms: ", psys%natoms
      do ii = 1, psys%natoms
        write( 6, '(A2,A1,1X,3F16.8)' ) psys%atom_list( ii )%el_name, ':', psys%atom_list( ii )%reduced_coord( : )
      enddo
      write( 6, * ) "-------------------------------"
      write( 6, * ) "Expected K-point mesh: "
      do ii = 1, params%nkpts
        write(6,*) screen_system_returnKvec( params, ii )
      enddo
      write( 6, * ) "-------------------------------"
  
    endif

  end subroutine screen_system_summarize

  subroutine screen_system_snatch( element, indx, tau, xcoord, ierr )
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
          goto 111
        endif
      endif
    enddo

    write( 6, * ) 'Atom coord not found!'
    ierr = -1
    return

111 continue

    xcoord = matmul( psys%avecs, tau )

  end subroutine screen_system_snatch

  subroutine screen_system_load( ierr )
    use OCEAN_mpi, only : myid, root, comm, nproc, MPI_BCAST, MPI_INTEGER, MPI_SUCCESS
    integer, intent( inout ) :: ierr
    !
#ifdef MPI
    integer :: ierr_
#endif
    

    if( myid .eq. root ) then
      call load_xyz( ierr )
      if( ierr .ne. 0 ) goto 111

      call load_abvecs( ierr )
      if( ierr .ne. 0 ) goto 111

      call load_params( ierr ) 
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

    endif
#endif
    params%nkpts = product( params%kmesh(:) )
    params%nbands = params%brange(4) - params%brange(3) + params%brange(2) - params%brange(1) + 2

  end subroutine screen_system_load


  subroutine share_params( ierr )
    use OCEAN_mpi
    integer, intent( inout ) :: ierr
    !
    if( nproc .eq. 1 ) return
    !
    call MPI_BCAST( params%kshift, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    
    call MPI_BCAST( params%brange, 4, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%kmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return

    call MPI_BCAST( params%nspin, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return


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

  ! NOT MPI SAFE ( in so much as it will let every process hit the filesystem )
  subroutine load_params( ierr )
    integer, intent( inout ) :: ierr
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

    open( unit=99, file='brange.ipt', form='formatted', status='old', IOSTAT=ierr )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to open brange.ipt', ierr
      return
    endif
    read( 99, *, IOSTAT=ierr ) params%brange( : )
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to read brange.ipt', ierr
      return
    endif
    close( 99, IOSTAT=ierr)
    if( ierr .ne. 0 ) then
      write( 6, * ) 'FATAL ERROR: Failed to close brange.ipt', ierr
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


end module screen_system
