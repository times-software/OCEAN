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

  private


  type site_info
    character( len=2 ) :: elname
    integer :: indx
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
  end type site


  type( site_parallel_info ), save, public :: pinfo

!  type( site ), allocatable, save :: all_sites( : )
!  integer, save :: n_sites

!  integer, allocatable :: tmp_indices( : )
!  character( len=2 ), allocatable :: tmp_elnames( : )

  public :: site
  public :: screen_sites_load, screen_sites_prep
  public :: screen_sites_returnWavefunctionDims, screen_sites_returnWavefunctionBK

  contains

  pure function screen_sites_returnWavefunctionDims( a ) result( dims )
    use screen_wavefunction, only : screen_wvfn
    type( site ), intent( in ) :: a
    integer :: dims(2)
    
    dims(1) = a%wvfn%mypts
    dims(2) = a%wvfn%npts
  end function screen_sites_returnWavefunctionDims

  pure function screen_sites_returnWavefunctionBK( a ) result( dims )
    use screen_wavefunction, only : screen_wvfn
    type( site ), intent( in ) :: a
    integer :: dims(2)

    dims(1) = a%wvfn%mybands
    dims(2) = a%wvfn%mykpts
  end function screen_sites_returnWavefunctionBK

  subroutine screen_sites_prep( nsites, ierr )
    integer, intent( out ) :: nsites
    integer, intent( inout ) :: ierr 
    !
    call count_sitelist( nsites, ierr )

  end subroutine

  subroutine screen_sites_load( n_sites, all_sites, ierr )
    use OCEAN_mpi
    use screen_system, only : screen_system_snatch
    use screen_grid, only : screen_grid_init
    use screen_paral, only : screen_paral_init
    use screen_wavefunction, only : screen_wvfn_init
    !
    integer, intent( in ) :: n_sites
    type( site ), intent( inout ) :: all_sites( : )
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: tau( 3 ), xcoord( 3 )
    integer :: i
    integer, allocatable :: tmp_indices( : )
    character( len=2 ), allocatable :: tmp_elnames( : )

    allocate( tmp_indices( n_sites ), tmp_elnames( n_sites ) )
    call load_sitelist( n_sites, tmp_indices, tmp_elnames, ierr )
    if( ierr .ne. 0 ) return
    ! Store into sites
    do i = 1, n_sites
      all_sites( i )%info%indx = tmp_indices( i )
      all_sites( i )%info%elname = tmp_elnames( i )
    enddo
    deallocate( tmp_indices, tmp_elnames )

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

  end subroutine screen_sites_load


  subroutine count_sitelist( nsites, ierr )
    use OCEAN_mpi, only : myid, root, comm, MPI_INTEGER
    !
    integer, intent( out ) :: nsites
    integer, intent( inout ) :: ierr
    !
    integer :: ierr_
    !
    if( myid .eq. root ) then
      open( unit=99, file='sitelist', form='formatted', status='old', iostat=ierr, action='read' )
      if( ierr .ne. 0 ) then
        write(6,*) 'Failed to open sitelist'
        goto 10
      endif

      rewind(99)
      read( 99, *, iostat=ierr ) nsites
      if( ierr .ne. 0 ) write(6,*) 'Failed to read sitelist'
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

    integer :: ierr_, nsites_, i
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
        read( 99, *, iostat=ierr ) tmp_elnames( i ), tmp_indices( i )
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

end module screen_sites
