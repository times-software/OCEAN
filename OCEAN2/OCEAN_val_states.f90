module OCEAN_val_states
  use AI_kinds

  implicit none
  save

  real(dp), allocatable :: re_val( :, :, :, : )
  real(dp), allocatable :: im_val( :, :, :, : )
  real(dp), allocatable :: re_con( :, :, :, : )
  real(dp), allocatable :: im_con( :, :, :, : )

  integer :: nkpts
  integer :: nxpts
  integer :: nbc
  integer :: nbv
  integer :: nspn
  integer :: nxpts_pad
  integer :: val_pad
  integer :: con_pad

  integer :: max_nxpts
  integer :: startx

  integer, allocatable :: nxpts_by_mpiID( : )
  integer, allocatable :: startx_by_mpiID( : )

  integer, parameter :: cache_double = 8
  logical, private :: is_init = .false.
  logical, private :: is_loaded = .false.

#ifdef CONTIGUOUS
  CONTIGUOUS :: re_val, im_val, re_con, im_con
#endif

#ifdef __INTEL
!dir$ attributes align:64 :: re_val, im_val, re_con, im_con
#endif

  
  public :: OCEAN_val_states_load, OCEAN_val_states_init

  contains

  subroutine OCEAN_val_states_load( sys, ierr )
    use OCEAN_system
    implicit none
    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    if( .not. is_init ) then
      ierr = -1
      return
    endif
 
    if( is_loaded ) return


    allocate( re_val( nxpts_pad, val_pad, nkpts, nspn ), &
              im_val( nxpts_pad, val_pad, nkpts, nspn ), &
              re_con( nxpts_pad, con_pad, nkpts, nspn ), &
              im_con( nxpts_pad, con_pad, nkpts, nspn ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_val_states_read( sys, ierr )
  end subroutine OCEAN_val_states_load

  subroutine OCEAN_val_states_init( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, nproc
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    ! 
    integer :: nx_remain, i

    if( is_init ) then
      return
    endif

    nkpts = sys%nkpts
    nspn  = sys%nspn
    nbc   = sys%cur_run%num_bands
    nbv   = sys%cur_run%val_bands

    nxpts = 0
    startx = 1
    nx_remain = sys%nxpts

    allocate( nxpts_by_mpiID( 0:nproc-1 ), startx_by_mpiID( 0:nproc-1 ) )

    do i = 0, myid
      startx = startx + nxpts
      nxpts = nx_remain / ( nproc - i )
      ! round up to cache line
      if( nxpts .ge. cache_double ) then
        nxpts = nxpts + mod( nxpts, cache_double )
      endif
      nx_remain = nx_remain - nxpts
      nxpts_by_mpiID( i ) = nxpts
    enddo

    do i = myid + 1, nproc - 1
      nxpts_by_mpiID( i ) = nx_remain / ( nproc - i )
      nx_remain = nx_remain - nxpts_by_mpiID( i )
    enddo

    startx_by_mpiID( 0 ) = 1
    do i = 1, nproc - 1
      startx_by_mpiID( i ) = startx_by_mpiID( i - 1 ) + nxpts_by_mpiID( i - 1 )
    enddo


    max_nxpts = maxval( nxpts_by_mpiID )
    if( nxpts .lt. cache_double ) then
      nxpts_pad = cache_double
      max_nxpts = max( max_nxpts, cache_double ) 
    else
      nxpts_pad = nxpts
    endif

    if( nxpts .lt. 1 ) then
      ierr = -1
      return
    endif

    con_pad = nbc + mod( nbc, cache_double )
    val_pad = nbv + mod( nbv, cache_double )

    is_init = .true.

  end subroutine
  

  subroutine OCEAN_val_states_read( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, nproc, root, comm
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
!    logical :: want_val = .true.
!    logical :: want_con = .true.
!    logical :: want_core = .false.
    logical :: invert_xmesh, ex
    character(len=3) :: bloch_type

    integer :: brange(4), err

    if( myid .eq. root ) then
      open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) brange(:)
      close( unit=99 )

      inquire(file="bloch_type",exist=ex)
      if( ex ) then
        open(unit=99,file="bloch_type", form='formatted', status='old' )
        read(99,*) bloch_type
        invert_xmesh = .true.
        close(99)
      else
        bloch_type = 'old'
        invert_xmesh = .true.
      endif


      if( nbv .gt. (1+brange(2)-brange(1)) ) then
        write(6,*) 'Not enough valence bands!'
        ierr = -1
      endif
      if( nbc .gt. (1+brange(4)-brange(3)) ) then
        write(6,*) 'Not enough conduction bands!'
        ierr = ierr - 2
      endif
    endif

#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, err )
    if( err .ne. MPI_SUCCESS ) then
      ierr = -4
      return
    endif
#endif
    if( ierr .ne. 0 ) return

#ifdef MPI
    call MPI_BCAST( bloch_type, 3, MPI_CHARACTER, root, comm, ierr )
    call MPI_BCAST( invert_xmesh, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    select case( bloch_type )
      case( 'new' ) 
        call load_new_u2( sys, brange, ierr )
      case default
        return
    end select 

  end subroutine OCEAN_val_states_read

  subroutine load_new_u2( sys, brange, ierr )
    use OCEAN_system
    use OCEAN_mpi
    implicit none
      
    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: brange( 4 )
    integer, intent( inout ) :: ierr

    complex(dp), allocatable :: u2_buf(:,:,:,:)
    real(dp), allocatable :: re_share_buffer(:,:,:), im_share_buffer(:,:,:)

    logical :: io_group = .false.
    integer :: width(3), io_comm, fmode, fhu2
    integer :: iq, nelement, ibd, ncount, iproc, xiter, iz, iy, ix
#ifdef MPI
    integer(MPI_OFFSET_KIND) :: offset
    integer :: u2_status(MPI_STATUS_SIZE)
#endif

    if( myid .eq. root ) then
      write(6,*) 'New U2 format'

      open(unit=99,file='obf_control',form='formatted',status='old')
      rewind(99)
      read(99,*) width(1)
      read(99,*) width(2)
      read(99,*) width(3)
      close(99)
    endif

#ifdef MPI
    call MPI_BCAST( width, 3, MPI_INTEGER, root, comm, ierr )
#endif
    if( ierr .ne. 0 ) return


    if( myid .eq. root ) then
      write(6,*) width(:)
    endif

    if( myid .eq. root ) then
      io_group = .true.
    endif

#ifdef MPI
    if( io_group ) then
      call MPI_COMM_SPLIT( comm, 1, myid, io_comm, ierr )
      if( ierr .ne. 0 ) return
    else
      call MPI_COMM_SPLIT( comm, MPI_UNDEFINED, myid, io_comm, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( myid .eq. root ) then
      fmode = MPI_MODE_RDONLY
      call MPI_FILE_OPEN( io_comm, 'u2par.dat', fmode, MPI_INFO_NULL, fhu2, ierr )
      if( ierr/=0 ) then
        return
      endif
      offset=0
      call MPI_FILE_SET_VIEW( fhu2, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr )
      if( ierr/=0 ) then
        return
      endif 

    endif
#else
    if( myid .eq. root ) then
      open(unit=fhu2,file='u2par.dat',form='binary',status='old',ERR=ierr)
      if( ierr .ne. 0 ) return
    endif
#endif

    if( myid .eq. root ) then
      allocate( u2_buf( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), max(brange(2), sys%cur_run%num_bands ) ) )
      allocate( re_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), nproc ), &
                im_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), nproc )  ) 
    else
      allocate( u2_buf( 1, 1, 1, 1 ) )
      allocate( re_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 1 ), &
                im_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 1 )  )
    endif

    do iq = 1, nkpts

      if( myid .eq. root ) then
        ! read valence
#ifdef MPI
        !JTV need to do some 2GB chunking at some point. Blerg.
        nelement = nbv * sys%nxpts
        call MPI_FILE_READ_AT( fhu2, offset, u2_buf, nelement, MPI_DOUBLE_COMPLEX, &
                                 u2_status, ierr )
        if( ierr .ne. 0 ) then
          write(6,*) 'u2 failed', ierr, iq
          return
        endif
        call MPI_GET_COUNT( u2_status, MPI_DOUBLE_COMPLEX, ncount, ierr )
        if( ncount .ne. nelement) then
          write(6,*) 'u2 read failed.', ncount, nelement, iq
        endif
        offset = offset + INT( (brange(2)-brange(1) + 1 ), MPI_OFFSET_KIND ) * INT( sys%nxpts, MPI_OFFSET_KIND )
#else
        read(fhu2) u2_buf(:,:,:,1:brange(2)-brange(1)+1)
#endif        
        do ibd = 1, nbv
          iproc = 0
          xiter = 0
          do iz = 1, sys%xmesh(3)
            do iy = 1, sys%xmesh(2)
              do ix = 1, sys%xmesh(1)
                xiter = xiter + 1
                if( xiter .gt. nxpts_by_mpiID( iproc ) ) then
                  iproc = iproc + 1
                  xiter = 1
                endif
                re_share_buffer( xiter, ibd, iproc ) = real(u2_buf( iz, iy, ix, ibd ))
                im_share_buffer( xiter, ibd, iproc ) = aimag(u2_buf( iz, iy, ix, ibd ))
              end do
            end do
          end do
        end do
      endif
!
#ifdef MPI
      do iproc = 0, nproc-1
        if( myid .eq. root .and. iproc .ne. myid) then
          call MPI_SEND( re_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, iproc, iproc, comm, ierr )
          call MPI_SEND( im_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, iproc, iproc+nproc, comm, ierr )
        elseif( myid .eq. iproc ) then
          call MPI_RECV( re_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, root, iproc, comm, MPI_STATUS_IGNORE, ierr ) 
          call MPI_RECV( im_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
        endif
      enddo
#endif

      iproc = myid
      if( myid .ne. root ) iproc = 1

      re_val( :, :, iq, 1 ) = re_share_buffer( 1:nxpts, 1:nbv, iproc )
      im_val( :, :, iq, 1 ) = im_share_buffer( 1:nxpts, 1:nbv, iproc )



      if( myid .eq. root ) then
      ! read conduction
      
#ifdef MPI
        !JTV need to do some 2GB chunking at some point. Blerg.
        nelement = nbc * sys%nxpts
        call MPI_FILE_READ_AT( fhu2, offset, u2_buf, nelement, MPI_DOUBLE_COMPLEX, &
                                 u2_status, ierr )
        if( ierr .ne. 0 ) then
          write(6,*) 'u2 failed', ierr, iq
          return
        endif
        call MPI_GET_COUNT( u2_status, MPI_DOUBLE_COMPLEX, ncount, ierr )
        if( ncount .ne. nelement) then
          write(6,*) 'u2 read failed.', ncount, nelement, iq
        endif
        offset = offset + INT( (brange(4)-brange(3) + 1 ), MPI_OFFSET_KIND ) * INT( sys%nxpts, MPI_OFFSET_KIND )
#else
        read(fhu2) u2_buf(:,:,:,1:brange(4)-brange(3)+1)
#endif        
        do ibd = 1, nbc
          iproc = 0
          xiter = 0
          do iz = 1, sys%xmesh(3)
            do iy = 1, sys%xmesh(2)
              do ix = 1, sys%xmesh(1)
                xiter = xiter + 1
                if( xiter .gt. nxpts_by_mpiID( iproc ) ) then
                  iproc = iproc + 1
                  xiter = 1
                endif
                re_share_buffer( xiter, ibd, iproc ) = real(u2_buf( iz, iy, ix, ibd ))
                im_share_buffer( xiter, ibd, iproc ) = aimag(u2_buf( iz, iy, ix, ibd ))
              end do
            end do
          end do
        end do
      endif
!
#ifdef MPI
      do iproc = 0, nproc-1
        if( myid .eq. root .and. iproc .ne. myid) then
          call MPI_SEND( re_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, iproc, iproc, comm, ierr )
          call MPI_SEND( im_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, iproc, iproc+nproc, comm, ierr )
        elseif( myid .eq. iproc ) then
          call MPI_RECV( re_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, root, iproc, comm, MPI_STATUS_IGNORE, ierr )
          call MPI_RECV( im_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
        endif
      enddo
#endif

      iproc = myid
      if( myid .ne. root ) iproc = 1

      re_con( :, :, iq, 1 ) = re_share_buffer( 1:nxpts, 1:nbc, iproc )
      im_con( :, :, iq, 1 ) = im_share_buffer( 1:nxpts, 1:nbc, iproc )

    enddo ! iq

    if( myid .eq. root ) then
#ifdef MPI
      call MPI_FILE_CLOSE( fhu2, ierr )
#else
      close( fhu2 )
#endif
    endif

    deallocate( re_share_buffer, im_share_buffer, u2_buf )

  end subroutine load_new_u2

end module OCEAN_val_states
