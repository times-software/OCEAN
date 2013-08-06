module OCEAN_system
  use AI_kinds
  use mpi
  implicit none


  type, public :: o_system
    real(DP)         :: celvol
    real(DP)         :: avec(3,3)
    real(DP)         :: bvec(3,3)
    real(DP)         :: bmet(3,3)
    integer( S_INT ) :: nkpts
    integer( S_INT ) :: nxpts
    integer( S_INT ) :: nalpha
    integer( S_INT ) :: num_bands
    integer( S_INT ) :: xmesh( 3 )
    integer( S_INT ) :: kmesh( 3 )
    integer( S_INT ) :: ZNL(3)
    integer( S_INT ) :: nspn = 1
    integer( S_INT ) :: nobf = 0
    integer( S_INT ) :: nruns

    logical          :: e0
    logical          :: mult
    logical          :: long_range
    logical          :: obf
    logical          :: conduct
    character(len=5) :: calc_type

    type(o_run), pointer :: cur_run => null()

  end type o_system


  type :: o_run
    real(DP) :: tau(3)
    integer( S_INT ) :: nalpha
    integer( S_INT ) :: ZNL(3)
    
    integer :: indx
    integer :: photon
    character(len=2) :: elname
    character(len=2) :: corelevel
    character(len=255) :: basename
    character(len=255) :: filename

    character(len=5) :: calc_type

    logical          :: e0
    logical          :: mult
    logical          :: long_range
    logical          :: obf
    logical          :: conduct
    
    type(o_run), pointer :: prev_run => null()
    type(o_run), pointer :: next_run => null()

  end type o_run


  

  contains 

  subroutine OCEAN_sys_init( sys, ierr )
    use OCEAN_mpi, ONLY : myid, comm, root
    implicit none
     

    type( o_system ), intent( out ) :: sys
    integer, intent( inout ) :: ierr

    real( DP ) :: inter
    real( DP ), parameter :: inter_min = 0.000001
    integer :: nruns 

    if( myid .eq. root ) then

      open(unit=99,file='xmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%xmesh(:)
      close(99)
      sys%nxpts = product( sys%xmesh(:) )

      open(unit=99,file='kmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%kmesh(:)
      close(99)
      sys%nkpts = product( sys%kmesh(:) )

      open(unit=99,file='ZNL',form='formatted',status='old')
      rewind(99) 
      read(99,*) sys%ZNL(:)
      close(99) 
      ! nalpha is ( nspin valence ) * ( nspin core ) * ( 2 * l_core + 1 )
      sys%nalpha = 4 * ( 2 * sys%ZNL(3) + 1 )

      open(unit=99,file='nbuse.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%num_bands
      close(99)

      call getabb( sys%avec, sys%bvec, sys%bmet )
      call getomega( sys%avec, sys%celvol )     


      sys%mult = .true.
      sys%long_range = .true.

      open(unit=99,file='cks.normal',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%conduct
      close(99)
      if( .not. sys%conduct ) then
        sys%mult = .false.
        sys%long_range = .false.
      endif

      open(unit=99,file='mode',form='formatted',status='old')
      rewind(99)
      read(99,*) inter
      close(99)
      if( inter .lt. inter_min ) then
        sys%mult = .false.
        sys%long_range = .false.
      endif
      

      sys%e0 = .true.
      sys%obf = .false.
      sys%calc_type = 'NaN'
!      sys%conduct = .true.

      
      
    endif
#ifdef MPI
! Could create an mpi_datatype, but probably not worth it
    call MPI_BCAST( sys%celvol, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%avec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bvec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bmet, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nkpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nxpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nalpha, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%num_bands, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%xmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%kmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%ZNL, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111


    call MPI_BCAST( sys%e0, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%mult, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%long_range, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%obf, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%conduct, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%calc_type, 5, MPI_CHAR, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111


111 continue

#endif
!    allocate( sys%cur_run, STAT=ierr )
!    if( ierr .ne. 0 ) return
!
!    sys%cur_run%nalpha = sys%nalpha
!    sys%cur_run%ZNL(:) = sys%ZNL(:)
!    sys%cur_run%indx = 1
!    sys%cur_run%photon = 1
!    sys%cur_run%elname = 'F_'
!    sys%cur_run%corelevel = '1s'
!    sys%cur_run%basename = 'xas'
!!    write(sys%cur_run%filename,'(A,A,A,I2.2,A,A,A,I2.2)') sys%cur_run%basename, &
!!          '_', sys%cur_run%elname, '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, &
!!          '_', sys%cur_run%photon
!    sys%cur_run%tau(:) = 0.0_DP
    call OCEAN_runlist_init( sys, nruns, ierr )
    sys%nruns = nruns

  end subroutine OCEAN_sys_init

  subroutine OCEAN_runlist_init( sys, running_total, ierr )
    use OCEAN_mpi, ONLY : myid, comm, root
    implicit none

    type( o_system ), intent( out ) :: sys
    integer, intent( out ) :: running_total
    integer, intent( inout ) :: ierr

    real(DP) :: tau(3)
    integer :: ZNL(3), indx, photon
    character(len=2) :: elname, corelevel, ein
    character(len=5) :: calc_type
    type(o_run), pointer :: temp_prev_run, temp_cur_run

    integer :: ntot, nmatch, iter
    logical :: found
    real(DP) :: tmp(3)


    running_total = 0 
    temp_prev_run => sys%cur_run

    if( myid .eq. root ) then
      open(unit=99,file='runlist',form='formatted',status='old')
      rewind(99)
    endif

    do
      if( myid .eq. root ) read(99,*,END=111)  ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type

      if( myid .eq. root ) then
        open(unit=98,file='xyz.wyck',form='formatted',status='old')
        rewind(98)
        read(98,*) ntot
        nmatch = 0
        found = .false.
        do iter = 1, ntot
          read ( 98, * ) ein, tmp( : )
          if ( ein .eq. elname ) then
            nmatch = nmatch + 1
            if ( nmatch .eq. indx ) then
              tau( : ) = tmp( : )
              found = .true.
            end if
          end if
          if ( found ) goto 112
        end do
        if( .not. found ) then
          ierr = -1
          return
        endif
112     continue 
        write ( 6, '(1a15,3f10.5)' ) 'snatched alpha=', tau( : )
        close(98)
      endif


#ifdef MPI
      call MPI_BCAST( tau, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( ZNL, 3, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( elname, 2, MPI_CHAR, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( corelevel, 2, MPI_CHAR, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( indx, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( photon, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( calc_type, 5, MPI_CHAR, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

      


      allocate( temp_cur_run, STAT=ierr )
      if( ierr .ne. 0 ) return
      if( running_total .eq. 0 ) then
        sys%cur_run => temp_cur_run
      else
        temp_cur_run%prev_run => temp_prev_run
        temp_prev_run%next_run => temp_cur_run
      endif
      temp_cur_run%tau(:) = tau(:)
      temp_cur_run%ZNL(:) = ZNL(:)
      temp_cur_run%elname = elname
      temp_cur_run%indx = indx
      temp_cur_run%corelevel = corelevel
      temp_cur_run%calc_type = calc_type
      temp_cur_run%photon = photon
      
      temp_cur_run%basename = 'abs'
      write(temp_cur_run%filename,'(A3,A1,A2,A1,I2.2,A1,A2,A1,I2.2)' ) temp_cur_run%basename, '_', temp_cur_run%elname, &
            '.', temp_cur_run%indx, '_', temp_cur_run%corelevel, '_', temp_cur_run%photon

      temp_prev_run => temp_cur_run
      running_total = running_total + 1
    enddo

111 continue

    if( running_total .lt. 1 ) then
      ierr = -1
      if(myid .eq. root ) write(6,*) 'Failed to read in any runs'
    endif
    if( myid .eq. root ) then
      close( 99 )
      write(6,*) 'Number of calcs to complete: ', running_total
    endif

    
  end subroutine OCEAN_runlist_init


    
end module OCEAN_system
