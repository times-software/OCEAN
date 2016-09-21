module OCEAN_val_states
  use AI_kinds

  implicit none
  save
  private


  real(dp), public, protected, allocatable :: re_val( :, :, :, : )
  real(dp), public, protected, allocatable :: im_val( :, :, :, : )
  real(dp), public, protected, allocatable :: re_con( :, :, :, : )
  real(dp), public, protected, allocatable :: im_con( :, :, :, : )

  integer, public, protected :: nkpts
  integer, public, protected :: nxpts
  integer, public, protected :: nbc
  integer, public, protected :: nbv
  integer, public, protected :: nspn
  integer, public, protected :: nxpts_pad
  integer, public, protected :: val_pad
  integer, public, protected :: con_pad

  integer, public, protected :: max_nxpts
  integer, public, protected :: startx

  integer, public, protected, allocatable :: nxpts_by_mpiID( : )
  integer, public, protected, allocatable :: startx_by_mpiID( : )

  integer, public, parameter :: cache_double = 1
  logical, private :: is_init = .false.
  logical, private :: is_loaded = .false.

!#ifdef CONTIGUOUS
!  CONTIGUOUS :: re_val, im_val, re_con, im_con
!#endif

#ifdef __INTEL
!dir$ attributes align:64 :: re_val, im_val, re_con, im_con
#endif

  
  public :: OCEAN_val_states_load, OCEAN_val_states_init, OCEAN_val_states_returnPadXpts

  contains

  subroutine val_states_add_phase( re_phase, im_phase )
    implicit none
    
    real(dp), intent( in ), dimension( nxpts_pad, nkpts ) :: re_phase, im_phase
!    integer, intent( inout ) :: ierr
    !
    integer :: ik, ix, ispn

    write(6,*) 'PHASES!'
    write(6,*) cmplx( re_val(1,1,1,1), im_val(1,1,1,1), DP )
    write(6,*) cmplx( re_val(1,4,1,1), im_val(1,4,1,1), DP )
    write(6,*) cmplx( re_val(2,1,2,1), im_val(2,1,2,1), DP )

    ! DROT is "backwards" from how we want the phases to go hence minus sign in definition of im_phase
    do ispn = 1, nspn
      do ik = 1, nkpts
        do ix = 1, nxpts
          call DROT( nbv, re_val( ix, 1, ik, ispn ), nxpts_pad, & 
                          im_val( ix, 1, ik, ispn ), nxpts_pad, &
                     re_phase( ix, ik ), im_phase( ix, ik ) )
          call DROT( nbc, re_con( ix, 1, ik, ispn ), nxpts_pad, &
                          im_con( ix, 1, ik, ispn ), nxpts_pad, &
                     re_phase( ix, ik ), im_phase( ix, ik ) )
        enddo
      enddo
    enddo
    write(6,*) cmplx( re_val(1,1,1,1), im_val(1,1,1,1), DP )
    write(6,*) cmplx( re_val(1,4,1,1), im_val(1,4,1,1), DP )
    write(6,*) cmplx( re_val(2,1,2,1), im_val(2,1,2,1), DP )

  end subroutine val_states_add_phase

  subroutine val_states_generate_phases( sys, re_phase, im_phase, ierr )
    use OCEAN_system
    use OCEAN_constants, only : PI_dp
    implicit none

    type(O_system), intent( in ) :: sys
    real(dp), intent( out ), dimension( nxpts_pad, nkpts ) :: re_phase, im_phase
    integer, intent( inout ) :: ierr
    ! 
    real(dp) :: c( 3 ), dotProduct !, xk( 3 )
    real(dp), allocatable :: kfac( :, : ), xfac( :, : )

    integer :: ix, iy, iz, iix, iik, ikx, iky, ikz

    allocate( kfac( 3, nkpts ), xfac( 3, nxpts ), STAT=ierr )
    if( ierr .ne. 0 ) return
  
    c( : ) = 2.0_dp * PI_dp / sys%kmesh( : )
    
    if( sys%nkpts .ne. nkpts ) then
      ierr = 2000
      return
    endif

    !
    
    !
    iik = 0
    do ikx = 0, sys%kmesh( 1 ) - 1
      do iky = 0, sys%kmesh( 2 ) - 1
        do ikz = 0, sys%kmesh( 3 ) - 1
          iik = iik + 1
!          xk( : ) = c( : ) * real( (/ ikx, iky, ikz /), dp )
          kfac( :, iik ) = c( : ) * real( (/ ikx, iky, ikz /), dp )
        enddo
      enddo
    enddo


    if( .false. ) then
      ix = 1 + ( startx - 1 ) / (sys%xmesh(2)*sys%xmesh(3))
      iy = 1 + mod( ( startx - 1 )/sys%xmesh(3), sys%xmesh(2 ) )
      iz = 1 + mod( startx - 1, sys%xmesh(2)*sys%xmesh(3) ) 
      do iix = startx, nxpts+startx-1
!        xfac( 1, iix - startx + 1 ) = real( 1 + ( iix - 1) / (sys%xmesh(2)*sys%xmesh(3)), dp )/ real( sys%xmesh( 1 ), dp )
!        xfac( 2, iix - startx + 1 ) = real( 1 + mod( ( iix - 1 )/sys%xmesh(3), sys%xmesh(2 )), dp )/ real( sys%xmesh( 2 ), dp ) 
!        xfac( 3, iix - startx + 1 ) = real( 1 + mod( iix - 1, sys%xmesh(2)*sys%xmesh(3) ) / real( sys%xmesh( 3 ), dp )
      enddo
    else

      iix = 0
      do ix = 1, sys%xmesh( 1 )
        do iy = 1, sys%xmesh( 2 )
          do iz = 1, sys%xmesh( 3 )
            iix = iix + 1
            xfac( 1, iix ) = dble( ix - 1 ) / dble( sys%xmesh(1) )
            xfac( 2, iix ) = dble( iy - 1 ) / dble( sys%xmesh(2) )
            xfac( 3, iix ) = dble( iz - 1 ) / dble( sys%xmesh(3) )
          enddo
        enddo
      enddo
    endif


    ! DROT is "backwards" from how we want the phases to go hence minus sign in definition of im_phase
    do iik = 1, nkpts
      do iix = 1, nxpts
        dotProduct = dot_product( xfac( :, iix ), kfac( :, iik ) )
        re_phase( iix, iik ) = dcos( dotProduct )
        im_phase( iix, iik ) = -dsin( dotProduct )
      enddo
    enddo

    deallocate( kfac, xfac )


  end subroutine val_states_generate_phases

  subroutine OCEAN_val_states_returnPadXpts( x_pad, ierr )
    implicit none
    integer, intent( out ) :: x_pad
    integer, intent( inout ) :: ierr
    !
    if( .not. is_init ) then
      ierr = 201
      return
    endif
    !
    x_pad = nxpts_pad
  end subroutine OCEAN_val_states_returnPadXpts


  ! The loading routines are redundantly done (currently)
  ! While the core and valence routines have different 
  !  storage orderings, there should be a more generic
  !  set of routines that then branches out for DFT-specific versions
  subroutine OCEAN_val_states_load( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    implicit none
    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    real(dp), allocatable :: re_phase( :, : ), im_phase( :, : )

    if( .not. is_init ) then
      ierr = -1
      return
    endif
 
    if( is_loaded ) return

    if( myid .eq. root ) write( 6, * )'Loading valence states'


    allocate( re_val( nxpts_pad, val_pad, nkpts, nspn ), &
              im_val( nxpts_pad, val_pad, nkpts, nspn ), &
              re_con( nxpts_pad, con_pad, nkpts, nspn ), &
              im_con( nxpts_pad, con_pad, nkpts, nspn ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_val_states_read( sys, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) write( 6, * ) 'Adding phase'
    allocate( re_phase( nxpts_pad, nkpts ), im_phase( nxpts_pad, nkpts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call val_states_generate_phases( sys, re_phase, im_phase, ierr )
    if( ierr .ne. 0 ) return 

    call val_states_add_phase( re_phase, im_phase )

    deallocate( re_phase, im_phase )

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

    !JTV Ladder needs modification to deal with padding
    nxpts_pad = nxpts

    if( nxpts .lt. 1 ) then
      ierr = -1
      return
    endif

    con_pad = nbc + mod( nbc, cache_double )
    val_pad = nbv + mod( nbv, cache_double )

    nxpts_pad = nxpts
    con_pad = nbc
    val_pad = nbv

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
!    logical :: invert_xmesh, ex
!    character(len=3) :: bloch_type

    integer :: brange(4), err

    if( myid .eq. root ) then
      open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) brange(:)
      close( unit=99 )

!      inquire(file="bloch_type",exist=ex)
!      if( ex ) then
!        open(unit=99,file="bloch_type", form='formatted', status='old' )
!        read(99,*) bloch_type
!        invert_xmesh = .true.
!        close(99)
!      else
!        bloch_type = 'old'
!        invert_xmesh = .true.
!      endif


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

!#ifdef MPI
!    call MPI_BCAST( bloch_type, 3, MPI_CHARACTER, root, comm, ierr )
!    call MPI_BCAST( invert_xmesh, 1, MPI_LOGICAL, root, comm, ierr )
!    if( ierr .ne. 0 ) return
!#endif

    select case( sys%bloch_selector )
      case( 1 ) 
        call load_new_u2( sys, brange, ierr )
      case( 0 )
        call load_old_u2( sys, brange, ierr )
      case default
        ierr = 500
        if( myid .eq. root ) write(6,*) 'Unsupported bloch_selector:', sys%bloch_selector
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
      allocate( re_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 0:nproc-1 ), &
                im_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 0:nproc-1 )  ) 
    else
      allocate( u2_buf( 1, 1, 1, 1 ) )
      allocate( re_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 1 ), &
                im_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 1 )  )
    endif

    do iq = 1, nkpts

      if( myid .eq. root ) then
!        write(6,*) iq, nkpts
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

!JTV -- might need to check on xyz inversion
        do ibd = 1, nbv
          iproc = 0
          xiter = 0
          do ix = 1, sys%xmesh(1)
            do iy = 1, sys%xmesh(2)
              do iz = 1, sys%xmesh(3)
                xiter = xiter + 1
                if( xiter .gt. nxpts_by_mpiID( iproc ) ) then
                  iproc = iproc + 1
                  xiter = 1
                endif
                re_share_buffer( xiter, ibd, iproc ) = real(u2_buf( iz, iy, ix, ibd ), DP)
                im_share_buffer( xiter, ibd, iproc ) = real(aimag(u2_buf( iz, iy, ix, ibd )),DP)
              end do
            end do
          end do
        end do
      endif
!
!      if( myid .eq. root ) write(6,*) 'Val read', myid, root, nproc-1
#ifdef MPI
      do iproc = 0, nproc-1
        if( myid .eq. root ) then
          if( iproc .ne. myid) then
!            write(6,*) 'Val send'
            call MPI_SEND( re_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, iproc, iproc, comm, ierr )
            call MPI_SEND( im_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, iproc, iproc+nproc, comm, ierr )
          endif
        elseif( myid .eq. iproc ) then
          call MPI_RECV( re_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, root, iproc, comm, MPI_STATUS_IGNORE, ierr ) 
          call MPI_RECV( im_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
        endif
      enddo
#endif
!      if( myid .eq. root ) write(6,*) 'Val shared'


      iproc = myid
      if( myid .ne. root ) iproc = 1

      re_val( 1:nxpts, 1:nbv, iq, 1 ) = re_share_buffer( 1:nxpts, 1:nbv, iproc )
      im_val( 1:nxpts, 1:nbv, iq, 1 ) = im_share_buffer( 1:nxpts, 1:nbv, iproc )



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
!JTV -- ditto comment above re: xyz ordering
        do ibd = 1, nbc
          iproc = 0
          xiter = 0
          do ix = 1, sys%xmesh(1)
            do iy = 1, sys%xmesh(2)
              do iz = 1, sys%xmesh(3)
                xiter = xiter + 1
                if( xiter .gt. nxpts_by_mpiID( iproc ) ) then
                  iproc = iproc + 1
                  xiter = 1
                endif
                re_share_buffer( xiter, ibd, iproc ) = real(u2_buf( iz, iy, ix, ibd ), DP)
                im_share_buffer( xiter, ibd, iproc ) = real(aimag(u2_buf( iz, iy, ix, ibd )), DP)
              end do
            end do
          end do
        end do
      endif
!
#ifdef MPI
      do iproc = 0, nproc-1
        if( myid .eq. root ) then
          if( iproc .ne. myid) then
            call MPI_SEND( re_share_buffer(:,:,iproc), max_nxpts*nbc, MPI_DOUBLE_PRECISION, iproc, iproc, comm, ierr )
            call MPI_SEND( im_share_buffer(:,:,iproc), max_nxpts*nbc, MPI_DOUBLE_PRECISION, iproc, iproc+nproc, comm, ierr )
          endif
        elseif( myid .eq. iproc ) then
          call MPI_RECV( re_share_buffer, max_nxpts*nbc, MPI_DOUBLE_PRECISION, root, iproc, comm, MPI_STATUS_IGNORE, ierr )
          call MPI_RECV( im_share_buffer, max_nxpts*nbc, MPI_DOUBLE_PRECISION, root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
        endif
      enddo
#endif

      iproc = myid
      if( myid .ne. root ) iproc = 1

      re_con( 1:nxpts, 1:nbc, iq, 1 ) = re_share_buffer( 1:nxpts, 1:nbc, iproc )
      im_con( 1:nxpts, 1:nbc, iq, 1 ) = im_share_buffer( 1:nxpts, 1:nbc, iproc )

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


  ! This is quite redudant 
  subroutine load_old_u2( sys, brange, ierr )
    use OCEAN_system
    use OCEAN_mpi
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: brange( 4 )
    integer, intent( inout ) :: ierr

    real(dp), allocatable :: re_share_buffer(:,:,:), im_share_buffer(:,:,:), ur(:,:,:,:), ui(:,:,:,:)

    logical :: io_group = .false.
    integer :: fhu2, idum(3)
    integer :: iq, nelement, ibd, ncount, iproc, xiter, iz, iy, ix, ispn


    if( myid .eq. root ) then
      write(6,*) 'Old U2 format'
    endif

    if( myid .eq. root ) then
      open(unit=fhu2,file='u2.dat',form='unformatted',status='old')
      if( ierr .ne. 0 ) return
    endif

    if( myid .eq. root ) then
      allocate( re_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 0:nproc-1 ), &
                im_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 0:nproc-1 ), &
                ur( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), max( nbv, nbc ) ), &
                ui( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), max( nbv, nbc ) ), STAT=ierr )
    else
      allocate( re_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 1 ), &
                im_share_buffer( max_nxpts, max(brange(2), sys%cur_run%num_bands ), 1 ), &
                ur( 1, 1, 1, 1 ), ui( 1, 1, 1, 1 ), STAT=ierr )
    endif
    if( ierr .ne. 0 ) return


    do ispn = 1, sys%nspn
      do iq = 1, sys%nkpts

      
        if( myid .eq. root ) then
          do ibd = 1, brange(2)-brange(1) + 1
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  read ( fhu2 ) idum( 1 : 3 ), ur( iz, iy, ix, ibd ), ui( iz, iy, ix, ibd )
                end do
              end do
            end do
          enddo
          do ibd = 1, nbv
            iproc = 0
            xiter = 0
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  xiter = xiter + 1
                  if( xiter .gt. nxpts_by_mpiID( iproc ) ) then
                    iproc = iproc + 1
                    xiter = 1
                  endif
                  re_share_buffer( xiter, ibd, iproc ) = ur( iz, iy, ix, ibd )
                  im_share_buffer( xiter, ibd, iproc ) = ui( iz, iy, ix, ibd )
                end do
              end do
            end do
          end do
        endif

#ifdef MPI
        do iproc = 0, nproc-1
          if( myid .eq. root ) then
            if( iproc .ne. myid) then
  !            write(6,*) 'Val send'
              call MPI_SEND( re_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, &
                             iproc, iproc, comm, ierr )
              call MPI_SEND( im_share_buffer(:,:,iproc), max_nxpts*nbv, MPI_DOUBLE_PRECISION, &
                             iproc, iproc+nproc, comm, ierr )
            endif
          elseif( myid .eq. iproc ) then
            call MPI_RECV( re_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, &
                           root, iproc, comm, MPI_STATUS_IGNORE, ierr )
            call MPI_RECV( im_share_buffer, max_nxpts*nbv, MPI_DOUBLE_PRECISION, &
                           root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
          endif
        enddo
#endif
  !      if( myid .eq. root ) write(6,*) 'Val shared'


        iproc = myid
        if( myid .ne. root ) iproc = 1

        re_val( 1:nxpts, 1:nbv, iq, ispn ) = re_share_buffer( 1:nxpts, 1:nbv, iproc )
        im_val( 1:nxpts, 1:nbv, iq, ispn ) = im_share_buffer( 1:nxpts, 1:nbv, iproc )


        ! Conduction bands are stacked on top of valence
        if( myid .eq. root ) then
          do ibd = 1, brange(4)-brange(3) + 1
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  read ( fhu2 ) idum( 1 : 3 ), ur( iz, iy, ix, ibd ), ui( iz, iy, ix, ibd )
                end do
              end do
            end do
          enddo
          do ibd = 1, nbc
            iproc = 0
            xiter = 0
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  xiter = xiter + 1
                  if( xiter .gt. nxpts_by_mpiID( iproc ) ) then
                    iproc = iproc + 1
                    xiter = 1
                  endif
                  re_share_buffer( xiter, ibd, iproc ) = ur( iz, iy, ix, ibd )
                  im_share_buffer( xiter, ibd, iproc ) = ui( iz, iy, ix, ibd )
                end do
              end do
            end do
          end do
        endif

!
#ifdef MPI
        do iproc = 0, nproc-1
          if( myid .eq. root ) then
            if( iproc .ne. myid) then
              call MPI_SEND( re_share_buffer(:,:,iproc), max_nxpts*nbc, MPI_DOUBLE_PRECISION, &
                             iproc, iproc, comm, ierr )
              call MPI_SEND( im_share_buffer(:,:,iproc), max_nxpts*nbc, MPI_DOUBLE_PRECISION, &
                             iproc, iproc+nproc, comm, ierr )
            endif
          elseif( myid .eq. iproc ) then
            call MPI_RECV( re_share_buffer, max_nxpts*nbc, MPI_DOUBLE_PRECISION, &
                           root, iproc, comm, MPI_STATUS_IGNORE, ierr )
            call MPI_RECV( im_share_buffer, max_nxpts*nbc, MPI_DOUBLE_PRECISION, &
                           root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
          endif
        enddo
#endif

        iproc = myid
        if( myid .ne. root ) iproc = 1

        re_con( 1:nxpts, 1:nbc, iq, 1 ) = re_share_buffer( 1:nxpts, 1:nbc, iproc )
        im_con( 1:nxpts, 1:nbc, iq, 1 ) = im_share_buffer( 1:nxpts, 1:nbc, iproc )

      enddo ! iq
    enddo  ! ispn

    deallocate( re_share_buffer, im_share_buffer, ur, ui )

    if( myid .eq. root ) then
      close( fhu2 )
    endif

  end subroutine load_old_u2

end module OCEAN_val_states
