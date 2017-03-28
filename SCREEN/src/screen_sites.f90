module screen_sites
  use ai_kinds, only : DP
  use screen_grid, only : sgrid
  use screen_paral, only : site_parallel_info

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
  end type site


  type( site_parallel_info ), save :: pinfo

  type( site ), allocatable, save :: all_sites( : )
  integer, save :: n_sites


  public :: screen_sites_load

  contains

  subroutine screen_sites_load( ierr )
    use OCEAN_mpi
    use screen_system, only : screen_system_snatch
    !
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: tau( 3 )
    integer :: i

    call screen_sites_load_hfinlist( ierr )
    if( ierr .ne. 0 ) return

    do i = 1, n_sites
      call screen_system_snatch( all_sites( i )%info%elname, all_sites( i )%info%indx, tau, ierr )
      if( ierr .ne. 0 ) then
        if( myid .eq. root ) write( 6, * ) 'Problems fetching location of site:', i, & 
                                           all_sites( i )%info%elname, all_sites( i )%info%indx
        return
      endif

      ! After the first grid just copy all the initialization data instead of reading from file
      if( i .eq. 1 ) then
        call screen_grid_init( all_sites( i )%grid, tau, ierr )
      else
        call screen_grid_init( all_sites( i )%grid, tau, ierr, all_sites( i - 1 )%grid )
      endif

    enddo

    call screen_paral_init( n_sites, pinfo, ierr )

  end subroutine screen_sites_load

  subroutine screen_sites_load_hfinlist( ierr )
    use OCEAN_mpi
    integer, intent( inout ) :: ierr
    !
    integer, allocatable :: tmp_indices( : ), tmp_indices2( : )
    character( len=2 ), allocatable :: tmp_elnames( : ), tmp_elnames2( : )
    character( len=50) :: pspname
    integer :: tmp_nsites, i, j
    integer :: ZNL(3)
    logical :: dup

    if( myid .eq. root ) then
      allocate( tmp_indices2( 10000 ), tmp_elnames2( 10000 ), STAT=ierr )
      if( ierr .ne. 0 ) goto 10

      tmp_nsites = 0
      open( unit=99, file='hfinlist', form='formatted', status='old' )
      do while( ierr .eq. 0 )
        tmp_nsites = tmp_nsites + 1
        read( 99, *, IOSTAT=ierr ) pspname, ZNL(:), tmp_elnames2( tmp_nsites ), & 
                                                    tmp_indices2( tmp_nsites )
      enddo
      close( 99 )

      tmp_nsites = tmp_nsites - 1
      if( tmp_nsites .lt. 1 ) then
        ierr  = -1
        goto 10
      endif

      allocate( tmp_indices( tmp_nsites ), tmp_elnames( tmp_nsites ) )

      ! hfinlist will have duplicate sites iff the edges are different, e.g. 1s vs 2p vd 3d
      !   The calculation of \chi is the same regardless of edge details
      n_sites = 1
      tmp_elnames( 1 ) = tmp_elnames2( 1 )
      tmp_indices( 1 ) = tmp_indices2( 1 )

      do i = 2, tmp_nsites
        dup = .false.
        do j = 1, n_sites
          if( ( tmp_indices( j ) .eq. tmp_indices2( i ) ) .or. &
              ( tmp_elnames( j ) .eq. tmp_elnames2( i ) )  ) then
            dup = .true.
            exit
          endif
        enddo
        if( .not. dup ) then
          n_sites = n_sites + 1
          tmp_indices( n_sites ) = tmp_indices2( i ) 
          tmp_elnames( n_sites ) = tmp_elnames2( i )
        endif
      enddo
      deallocate( tmp_elnames2, tmp_indices2 )

      write( 6, * ) 'N sites: ', n_sites
      do i = 1, n_sites
        write( 6, * ) '  ', tmp_elnames( i ), tmp_indices( i )
      enddo

    endif
10  continue


    ! Share with all the threads
#ifdef MPI
    if( nproc .gt. 1 ) then
      call MPI_BCAST( n_sites, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return
      if( myid .ne. root ) allocate( tmp_elnames( n_sites ), tmp_indices( n_sites ) )
      call MPI_BCAST( tmp_elnames, 2*n_sites, MPI_CHARACTER, root, comm, ierr )
      if( ierr .ne. 0 ) return
      call MPI_BCAST( tmp_indices, n_sites, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return
    endif
#endif 

    ! Store into sites
    allocate( all_sites( n_sites ) )
    do i = 1, n_sites
      all_sites( i )%info%indx = tmp_indices( i )
      all_sites( i )%info%elname = tmp_elnames( i )
    enddo
    deallocate( tmp_indices, tmp_elnames )


  end subroutine screen_sites_load_hfinlist

end module screen_sites
