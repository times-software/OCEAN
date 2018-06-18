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
! The idea here is to supplant the PREP directory (old DENDIP)
!   This module will contain the wrappers to interface with various DFT style files
!
!
! 
module ocean_prep
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM
#endif


  implicit none
  save
  private

  integer, parameter :: LEGACY = 0
  integer, parameter :: ABINIT = 1
  integer, parameter :: QESPRESSO = 2


  type op_parallel_info
    integer :: total_procs  ! this is (often) copy of ocean_mpi : nproc
    integer :: nprocs   ! how many procs assigned to this site
    integer :: num_groups   ! how many groups exist
    integer :: mygroup  ! group index
    integer :: myid   ! mpi index w/i site group
    integer :: root = 0
#ifdef MPI_F08
    type( MPI_COMM ) :: comm
#else
    integer :: comm
#endif
  end type op_parallel_info

  integer :: op_nkpts
  integer :: op_nbands
!  integer :: kpts( 3 )
  integer :: op_brange( 4 )

!  real( DP ) :: k0( 3 )

  type( op_parallel_info ) :: pinfo

  character :: op_selector  ! Can be 'C'onduction, 'V'alence, or 'A'll
  logical :: op_skips  ! Determines if the number of kpts is actually doubled 
                       ! with unshifted/shifted grids interleaved 
  logical :: is_init = .false.
  logical :: is_read_init = .false.

  public :: ocean_prep_init
  
  contains


  subroutine ocean_prep_read_init( flavor, ierr )
    use OCEAN_mpi, only : myid, root, comm
    integer, intent( in ) :: flavor
    integer, intent( inout ) :: ierr
    !

    select case( flavor )

    case( LEGACY )
      call olf_read_init( myid, root, comm, ierr )

    case default
      ierr = 11
    end select
    if( ierr .ne. 0 ) return

    is_read_init = .true.

  end subroutine ocean_prep_read_init


  ! depending on file format this may use group comms
  subroutine ocean_prep_test_kpt( flavor, ikpt, is_my_kpt, my_ngvecs, my_bands, ierr )
    use ocean_legacy_files, only : olf_get_ngvecs_at_kpt
    integer, intent( in ) :: flavor
    integer, intent( in ) :: ikpt
    logical, intent( out ) :: is_my_kpt
    integer, intent( out ) :: my_ngvecs, my_bands
    integer, intent( inout ) :: ierr
    !
    if( .not. is_init ) then
      ierr = 11
      return
    endif

    if( .not. is_read_init ) then
      call ocean_prep_read_init( flavor, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( mod( ikpt - 1, pinfo%num_groups ) .eq. pinfo%mygroup ) then
      is_my_kpt = .true.
    else
      is_my_kpt = .false.
      my_ngvecs = 0
      my_bands = 0
      return
    endif

    ! For those whose k-point it is determine bands
    ! everyone gets at least int division bands
    my_bands = op_nbands / pinfo%nprocs
    ! if there are uneven bands, the lowest remainder get an additional
    if( pinfo%myid .lt. mod( op_nbands, pinfo%nprocs ) ) my_bands = my_bands + 1

    ! Figure out gvecs
    select case( flavor )

    case( LEGACY )
      call olf_get_ngvecs_at_kpt( ikpt, pinfo%myid, pinfo%root, pinfo%comm, my_ngvecs, ierr )

    case default
      ierr = 11
    end select
    if( ierr .ne. 0 ) return


  end subroutine ocean_prep_test_kpt

  subroutine ocean_prep_read_at_kpt( flavor, ikpt, my_bands, my_ngvecs, gvecs, wavefunctions, ierr )
    use ocean_legacy_files, only : olf_read_at_kpt
    integer, intent( in ) :: flavor, ikpt, my_ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, my_ngvecs )
    complex( DP ), intent( out ) :: wavefunctions( my_ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
    select case ( flavor )
      
    case( LEGACY )
      call olf_read_at_kpt( ikpt, pinfo%myid, pinfo%root, pinfo%nprocs, pinfo%comm, op_nbands, &
                            my_bands, my_ngvecs, gvecs, wavefunctions, ierr )
    case default
      ierr =12
    end select    

  end subroutine ocean_prep_read_at_kpt

  subroutine ocean_prep_paral_init( ierr )
    use ocean_mpi, only : myid, root, comm, nproc, & 
      MPI_SUCCESS, MPI_COMM_SPLIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_INTEGER
    !
    integer, intent( inout ) :: ierr
    !
    integer :: test_size

    if( .not. is_init ) then
      ierr = 10
      return
    endif

    ! Need to determine the size of the processor groups 
    !  Each group will tackle a single k-point and divide the bands

    call how_many_groups( op_nkpts, nproc, pinfo%num_groups )
    !
    pinfo%nprocs = nproc / pinfo%num_groups
    pinfo%total_procs = pinfo%num_groups * pinfo%nprocs
    !
    ! integer division will give us pinfo%nprocs for group 0, and etc
    pinfo%mygroup = myid / pinfo%nprocs
    call MPI_COMM_SPLIT( comm, pinfo%mygroup, myid, pinfo%comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    !
    call MPI_COMM_RANK( pinfo%comm, pinfo%myid, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    !
    call MPI_COMM_SIZE( pinfo%comm, test_size, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    if( test_size .ne. pinfo%nprocs ) then
      write( 6, * ) 'Problems creating site comms', myid, test_size, pinfo%nprocs
      ierr = -101
      return
    endif

  end subroutine ocean_prep_paral_init

  subroutine ocean_prep_init( ierr, nkpts, brange, selector, skips )
    use ocean_mpi, only : myid, root, comm, MPI_BCAST, MPI_INTEGER
    !
    integer, intent( inout ) :: ierr
    integer, intent( in ), optional :: nkpts
    integer, intent( in ), optional :: brange( 4 )
    character, intent( in ), optional :: selector
    logical, intent( in ), optional :: skips
    !
    integer :: kmesh(3), ierr_

    if( present( skips ) ) then
      op_skips = skips
    else
      op_skips = .false.
    endif

    if( present( nkpts ) ) then
      op_nkpts = nkpts
    else
      if( myid .eq. root ) then
        open( unit=99, file='kmesh.ipt', form='formatted', status='old', IOSTAT=ierr )
        if( ierr .ne. 0 ) then
          write(6,*) 'FATAL ERROR: Failed to open kmesh.ipt', ierr
          goto 10
        endif
        read( 99, *, IOSTAT=ierr ) kmesh( : )
        if( ierr .ne. 0 ) then
          write( 6, * ) 'FATAL ERROR: Failed to read kmesh.ipt', ierr
          goto 10
        endif
        close( 99 )

        op_nkpts = product( kmesh( : ) )
      endif
10    continue 
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) return
      call MPI_BCAST( op_nkpts, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return
    endif
      
    if( present( brange ) ) then
      op_brange( : ) = brange( : )
    else
      if( myid .eq. root ) then
        open( unit=99, file='brange.ipt', form='formatted', status='old', IOSTAT=ierr )
        if( ierr .ne. 0 ) then
          write( 6, * ) 'FATAL ERROR: Failed to open brange.ipt', ierr
          goto 11
        endif
        read( 99, *, IOSTAT=ierr ) op_brange( : )
        if( ierr .ne. 0 ) then
          write( 6, * ) 'FATAL ERROR: Failed to read brange.ipt', ierr
          goto 11
        endif
        close( 99, IOSTAT=ierr)
        if( ierr .ne. 0 ) then
          write( 6, * ) 'FATAL ERROR: Failed to close brange.ipt', ierr
          goto 11
        endif
      endif
11    continue
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) return
      call MPI_BCAST( op_brange, 4, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( present( selector ) ) then
      op_selector = selector
    else
      op_selector = 'A'
    endif
    if( op_selector .eq. 'c' ) op_selector = 'C'
    if( op_selector .eq. 'v' ) op_selector = 'V'
    if( op_selector .ne. 'C' .and. op_selector .ne. 'V' ) op_selector = 'A'
  

    select case( op_selector )

      case( 'C' )
        op_nbands = brange( 4 ) - brange( 3 ) + 1
      
      case( 'V' )
        op_nbands = brange( 2 ) - brange( 1 ) + 1
  
      case default
        op_selector = 'A'
        op_nbands = brange( 4 ) - brange( 1 ) + 1

    end select

    is_init = .true.

    call ocean_prep_paral_init( ierr )

  end subroutine ocean_prep_init


  ! At the moment just a copy from screen_paral
  !   In the future we'll probably want to tune (such as nbands awareness )
  subroutine how_many_groups( n_sites, nproc, ngroups )
    use OCEAN_mpi, only : myid, root
    integer, intent( in ) :: n_sites, nproc
    integer, intent( out ) :: ngroups
    !
    real( DP ) :: score, best_score
    integer :: i, j, best_ngroup

    if( n_sites .eq. 1 .or. nproc .eq. 1 ) then
      ngroups = 1
      return
    endif

    if( mod( nproc, n_sites ) .eq. 0 ) then
      ngroups = n_sites
      return
    endif

    if( mod( n_sites, nproc ) .eq. 0 ) then
      ngroups = nproc
      return
    endif

    best_ngroup = 1
    best_score = 0.0_DP

    do i = 1, n_sites
      ! i groups, j per group, i * j being used
      j = nproc / i
      score = real( j * i, DP ) / real( nproc, DP )
      ! and then if mod( n_sites, ngroups ) .ne. 0 then more idle
      if( mod( n_sites, i ) .ne. 0 ) then
        score = score * real( mod( n_sites, i ), DP ) / real( i, DP )
      endif

      if( score .gt. best_score ) then
        best_score = score
        best_ngroup = i
      endif
    enddo


    ngroups = best_ngroup
    if( myid .eq. root ) write( 6, * ) 'Using groups: ', ngroups, best_score

    return

  end subroutine how_many_groups


end module ocean_prep
