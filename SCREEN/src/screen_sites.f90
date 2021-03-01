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
module screen_sites
  use ai_kinds, only : DP
  use screen_grid, only : sgrid
  use screen_paral, only : site_parallel_info
  use screen_wavefunction, only : screen_wvfn

  implicit none
  private


  type site_info
    character( len=2 ) :: elname
    integer :: indx
    integer :: z
  end type site_info

!  type site_parallel_info
!    integer :: total_procs  ! this is a copy of ocean_mpi : nproc
!    integer :: nprocs   ! how many procs assigned to this site
!    integer :: num_groups   ! how many groups exist
!    integer :: mygroup  ! group index
!    integer :: myid   ! mpi index w/i site group
!    integer :: comm
!    integer :: root
!  end type site_parallel_info


  type site
    type( sgrid ) :: grid
    type( site_info ) :: info
!    type( site_parallel_info ) :: pinfo
    type( screen_wvfn ) :: wvfn
    real(DP), allocatable :: shells(:)
  end type site


  type( site_parallel_info ), save, public :: pinfo

!  type( site ), allocatable, save :: all_sites( : )
!  integer, save :: n_sites

!  integer, allocatable :: tmp_indices( : )
!  character( len=2 ), allocatable :: tmp_elnames( : )

  public :: site, site_info
  public :: screen_sites_load, screen_sites_prep, screen_sites_load_atom
  public :: screen_sites_returnWavefunctionDims, screen_sites_returnWavefunctionBK

  contains

  pure function screen_sites_returnWavefunctionDims( a ) result( dims )
    use screen_wavefunction, only : screen_wvfn_returnWavefunctionDims
    type( site ), intent( in ) :: a
    integer :: dims(2)
    
    dims(:) = screen_wvfn_returnWavefunctionDims( a%wvfn )
  end function screen_sites_returnWavefunctionDims

  pure function screen_sites_returnWavefunctionBK( a ) result( dims )
    use screen_wavefunction, only : screen_wvfn_returnWavefunctionBK
    type( site ), intent( in ) :: a
    integer :: dims(2)

    dims(:) = screen_wvfn_returnWavefunctionBK( a%wvfn )
  end function screen_sites_returnWavefunctionBK

  subroutine screen_sites_prep( nsites, ierr )
    integer, intent( out ) :: nsites
    integer, intent( inout ) :: ierr 
    !
    call count_sitelist( nsites, ierr )

  end subroutine

!> @brief load the info for all the sites. Redirects to either the core or grid versions
  subroutine screen_sites_load( n_sites, all_sites, ierr )
    use screen_system, only : screen_system_mode
    !
    integer, intent( in ) :: n_sites
    type( site ), intent( inout ) :: all_sites( : )
    integer, intent( inout ) :: ierr
    !
    select case( screen_system_mode() )

      case( 'core' )
        call screen_sites_load_core( n_sites, all_sites, ierr )
      case( 'grid' )
        call screen_sites_load_grid( n_sites, all_sites, ierr )
      case default
        ierr = 1
    end select 

  end subroutine screen_sites_load

!> @brief Sets up sites on the FFT grid xmesh which is used by the BSE
!! @todo Move the xmesh file read up somewhere else so it is done only once (currently twice)
  subroutine screen_sites_load_grid( n_sites, all_sites, ierr )
    use screen_system, only : tau2xcoord
    use screen_grid, only : screen_grid_init
    use screen_paral, only : screen_paral_init
    use screen_wavefunction, only : screen_wvfn_init
    !
    integer, intent( in ) :: n_sites
    type( site ), intent( inout ) :: all_sites( : )
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: tau( 3 ), xcoord( 3 )
    integer :: i, ix, iy, iz, nshells, xmesh(3)
    real(DP), allocatable :: shells(:)

    call count_shells( nshells, ierr )
    if( ierr .ne. 0 ) return
    allocate( shells( nshells ) )
    call read_shells( nshells, shells, ierr )
    if( ierr .ne. 0 ) return

    call get_xmesh( xmesh, ierr )
    if( ierr .ne. 0 ) return

    i = 0
    do ix = 1, xmesh( 1 )
      tau( 1 ) = real( ix - 1, DP ) / real( xmesh( 1 ), DP )
      do iy = 1, xmesh( 2 )
        tau( 2 ) = real( iy - 1, DP ) / real( xmesh( 2 ), DP )
        do iz = 1, xmesh( 3 )
          tau( 3 ) = real( iz - 1, DP ) / real( xmesh( 3 ), DP )
          i = i + 1
!         do i = 1, n_sites
          all_sites( i )%info%indx = i
          all_sites( i )%info%elname = '__'
          all_sites( i )%info%z = 0
      
          allocate( all_sites( i )%shells( nshells ) )
          all_sites( i )%shells( : ) = shells( 1 : nshells )

          ! convert tau to xcoord
!          xcoord = matmul( psys%avecs, tau )
          xcoord = tau2xcoord( tau )
          ! After the first grid just copy all the initialization data instead of reading from file
          if( i .eq. 1 ) then
            call screen_grid_init( all_sites( i )%grid, xcoord, ierr )
          else
            call screen_grid_init( all_sites( i )%grid, xcoord, ierr, all_sites( 1 )%grid )
          endif
          if( ierr .ne. 0 ) return
        enddo
      enddo
    enddo

    call screen_paral_init( n_sites, pinfo, ierr )

    do i = 1, n_sites
      call screen_wvfn_init( pinfo, all_sites( i )%grid, all_sites( i )%wvfn, i, ierr )
      if( ierr .ne. 0 ) return
    enddo

    deallocate( shells )

  end subroutine screen_sites_load_grid

!> @brief loads xmesh, which might be better moved to system at some point
  subroutine get_xmesh( xmesh, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    !
    integer, intent( out ) :: xmesh( 3 )
    integer, intent( inout ) :: ierr
    !

    if( myid .eq. root ) then
      open( unit=99, file='xmesh.ipt', form='formatted', status='old' )
      read( 99, * ) xmesh( : )
      close( 99 )
    endif

#ifdef MPI
    call MPI_BCAST( xmesh, 3, MPI_INTEGER, root, comm, ierr )
#endif

  end subroutine get_xmesh

  subroutine screen_sites_load_core( n_sites, all_sites, ierr )
    use OCEAN_mpi
    use screen_system, only : screen_system_snatch
    use screen_grid, only : screen_grid_init
    use screen_paral, only : screen_paral_init
    use screen_wavefunction, only : screen_wvfn_init
    use periodic, only : get_atom_number
    !
    integer, intent( in ) :: n_sites
    type( site ), intent( inout ) :: all_sites( : )
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: tau( 3 ), xcoord( 3 )
    integer :: i, atomNumber, nshells
    real(DP), allocatable :: shells(:)
    integer, allocatable :: tmp_indices( : )
    character( len=2 ), allocatable :: tmp_elnames( : )

    allocate( tmp_indices( n_sites ), tmp_elnames( n_sites ) )
    call load_sitelist( n_sites, tmp_indices, tmp_elnames, ierr )
    if( ierr .ne. 0 ) return

    ! currently every site gets the same shell list
    call count_shells( nshells, ierr )
    if( ierr .ne. 0 ) return
    allocate( shells( nshells ) )
    call read_shells( nshells, shells, ierr )
    if( ierr .ne. 0 ) return

    ! Store into sites
    do i = 1, n_sites
      all_sites( i )%info%indx = tmp_indices( i )
      all_sites( i )%info%elname = tmp_elnames( i )
      call get_atom_number( atomNumber, tmp_elnames( i ) )
      if( atomNumber .gt. 0 ) then
        all_sites( i )%info%z = atomNumber
      else
        ierr = 10
        return
      endif

      allocate( all_sites( i )%shells( nshells ) )
      all_sites( i )%shells( : ) = shells( 1 : nshells )
    enddo
    deallocate( tmp_indices, tmp_elnames, shells )

    do i = 1, n_sites
      call screen_system_snatch( all_sites( i )%info%elname, all_sites( i )%info%indx, tau, xcoord, ierr )
      if( ierr .ne. 0 ) then
        if( myid .eq. root ) write( 6, * ) 'Problems fetching location of site:', i, & 
                                           all_sites( i )%info%elname, all_sites( i )%info%indx
        return
      endif

      ! After the first grid just copy all the initialization data instead of reading from file
      if( i .eq. 1 ) then
        call screen_grid_init( all_sites( i )%grid, xcoord, ierr )
      else
        call screen_grid_init( all_sites( i )%grid, xcoord, ierr, all_sites( 1 )%grid )
      endif
      if( ierr .ne. 0 ) return

    enddo

    call screen_paral_init( n_sites, pinfo, ierr )

    do i = 1, n_sites 
      call screen_wvfn_init( pinfo, all_sites( i )%grid, all_sites( i )%wvfn, i, ierr )
      if( ierr .ne. 0 ) return
    enddo

  end subroutine screen_sites_load_core

  subroutine screen_sites_load_atom( n_sites, all_sites, ierr )
    use OCEAN_mpi
    use screen_system, only : screen_system_grab
    use screen_grid, only : screen_grid_init
    use periodic, only : get_atom_number
    !
    integer, intent( in ) :: n_sites
    type( site ), intent( inout ) :: all_sites( : )
    integer, intent( inout ) :: ierr

    real( DP ) :: tau( 3 ), xcoord( 3 )
    integer :: i

    do i = 1, n_sites
      call screen_system_grab( i, all_sites( i )%info%elname, all_sites( i )%info%indx, tau, xcoord, ierr )
      if( ierr .ne. 0 ) return

      call get_atom_number( all_sites( i )%info%z, all_sites( i )%info%elname )
      
      if( i .eq. 1 ) then
        call screen_grid_init( all_sites( i )%grid, xcoord, ierr, Z=all_sites( i )%info%z )
      elseif( all_sites( i )%info%z .eq. all_sites( i-1 )%info%z ) then
        call screen_grid_init( all_sites( i )%grid, xcoord, ierr, all_sites( i-1 )%grid )
      else
        call screen_grid_init( all_sites( i )%grid, xcoord, ierr, Z=all_sites( i )%info%z )
      endif

    enddo

  end subroutine screen_sites_load_atom

      

!> @brief Returns a count of unique sites for calculating the screening. Either the number of
!> atomic sites that are being screened, or the count of xmesh for valence code
  subroutine count_sitelist( nsites, ierr )
    use OCEAN_mpi, only : myid, root, comm, MPI_INTEGER
    use screen_system, only : screen_system_mode
    !
    integer, intent( out ) :: nsites
    integer, intent( inout ) :: ierr
    !
    integer :: ierr_, xmesh(3)
    !
    if( myid .eq. root ) then

      select case( screen_system_mode() )

        case( 'core' ) 

          open( unit=99, file='sitelist', form='formatted', status='old', iostat=ierr, action='read' )
          if( ierr .ne. 0 ) then
            write(6,*) 'Failed to open sitelist'
            goto 10
          endif

          rewind(99)
          read( 99, *, iostat=ierr ) nsites
          if( ierr .ne. 0 ) write(6,*) 'Failed to read sitelist'
          close( 99 )

        case( 'grid' )

          open( unit=99, file='xmesh.ipt', form='formatted', status='old', iostat=ierr, action='read' )
          if( ierr .ne. 0 ) then
            write(6,*) 'Failed to open xmesh.ipt'
            goto 10
          endif

          rewind(99)
          read( 99, *, iostat=ierr ) xmesh(:)
          if( ierr .ne. 0 ) then 
            write(6,*) 'Failed to read xmesh'
          else
            nsites = product( xmesh )
          endif
          close( 99 )

        case default
          ierr = 1
      end select
    endif

10  continue

#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to sync ierr in screen_sites_count_sitelist'
      ierr = ierr_
      return
    endif

    call MPI_BCAST( nsites, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine count_sitelist

  subroutine load_sitelist( nsites, tmp_indices, tmp_elnames, ierr )
    use OCEAN_mpi, only : myid, comm, root, MPI_INTEGER, MPI_CHARACTER
    !
    integer, intent( in ) :: nsites
    integer, intent( inout ) :: ierr
    integer, intent( out ) :: tmp_indices( : )
    character( len=2 ), intent( out ) :: tmp_elnames( : )

    integer :: ierr_, nsites_, i, idum
    !
    if( myid .eq. root ) then
      open( unit=99, file='sitelist', form='formatted', status='old', iostat=ierr, action='read' )
      if( ierr .ne. 0 ) then
        write(6,*) 'Failed to open sitelist'
        goto 10
      endif

      rewind(99)
      read( 99, *, iostat=ierr ) nsites_
      if( ierr .ne. 0 ) then
        write(6,*) 'Failed to read sitelist'
        goto 10
      endif
      
      if( nsites .ne. nsites_ ) then
        write(6,*) 'Mismatch in nsites. File sitelist must have changed during execution!'
        ierr = -6
        goto 10
      endif

      do i = 1, nsites
        read( 99, *, iostat=ierr ) tmp_elnames( i ), idum, tmp_indices( i )
        if( ierr .ne. 0 ) then
          write( 6, * ) 'Failed reading sitelist. Site #', i
          goto 10
        endif
      enddo

      close( 99 )

    endif

10  continue

#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to sync ierr in screen_sites_count_sitelist'
      ierr = ierr_
      return
    endif

    call MPI_BCAST( tmp_elnames, 2*nsites, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( tmp_indices, nsites, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine load_sitelist

  subroutine count_shells( nshells, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    integer, intent( out ) :: nshells
    integer, intent( inout ) :: ierr

    if( myid .eq. root ) then
      open( unit=99, file='shells', form='formatted', status='old' )
      rewind( 99 )
      read( 99, * ) nshells
      close( 99 )
    endif
#ifdef MPI
    call MPI_BCAST( nshells, 1, MPI_INTEGER, root, comm, ierr )
#endif
  end subroutine count_shells

  subroutine read_shells( nshells, shells, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER, MPI_DOUBLE_PRECISION
    integer, intent( inout ) :: nshells
    real(DP), intent( out ) :: shells( nshells )
    integer, intent( inout ) :: ierr

    integer :: nshells_, ierr_

    ierr_ = 0
    nshells_ = nshells
    if( myid .eq. root ) then
      open( unit=99, file='shells', form='formatted', status='old' )
      rewind( 99 )
      read( 99, * ) nshells

      if( nshells .gt. nshells_ ) then
        ierr = 1
        goto 11
      endif
      read(99,*) shells(:)
      close( 99 )
    endif

11  continue
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
#endif
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif
#ifdef MPI
    call MPI_BCAST( nshells, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( shells, nshells, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif
  end subroutine read_shells

end module screen_sites
