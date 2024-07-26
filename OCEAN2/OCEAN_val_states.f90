module OCEAN_val_states
  use AI_kinds

  implicit none
  save
  private


  real(dp), public, protected, allocatable :: re_val( :, :, :, :, : )
  real(dp), public, protected, allocatable :: im_val( :, :, :, :, : )
  real(dp), public, protected, allocatable :: re_con( :, :, :, :, : )
  real(dp), public, protected, allocatable :: im_con( :, :, :, :, : )

  real(sp), public, protected, allocatable :: re_val_sp( :, :, :, : )
  real(sp), public, protected, allocatable :: im_val_sp( :, :, :, : )
  real(sp), public, protected, allocatable :: re_con_sp( :, :, :, : )
  real(sp), public, protected, allocatable :: im_con_sp( :, :, :, : )


  integer, public, protected :: nkpts
  integer, public, protected :: nxpts
  integer, public, protected :: nbc
  integer, public, protected :: nbv
  integer, public, protected :: nspn
  integer, public, protected :: nxpts_pad
  integer, public, protected :: val_pad
  integer, public, protected :: con_pad
  integer, public, protected :: nbw

  integer, public, protected :: max_nxpts
  integer, public, protected :: startx

  integer, public, protected, allocatable :: nxpts_by_mpiID( : )
  integer, public, protected, allocatable :: startx_by_mpiID( : )

  integer, public, parameter :: cache_double = 1
  logical, private :: is_init = .false.
  logical, private :: is_loaded = .false.

!  logical, public, parameter :: use_sp = .false.


#ifdef __INTEL
!dir$ attributes align:64 :: re_val, im_val, re_con, im_con
#endif

  
  public :: OCEAN_val_states_load, OCEAN_val_states_init, OCEAN_val_states_returnPadXpts

  contains

  subroutine val_states_add_phase( re_phase, im_phase )
    use OCEAN_mpi, only : myid
    implicit none
    
    real(dp), intent( in ), dimension( nxpts_pad, nkpts ) :: re_phase, im_phase
!    integer, intent( inout ) :: ierr
    !
    integer :: ik, ix, ispn, ibw

    if( nxpts .lt. 1 ) return

    if( .false. ) then
      write(6,*) 'PHASES!'
  !    write(6,*) myid, cmplx( re_val(1,1,1,1), im_val(1,1,1,1), DP )
      write(6,*) myid, startx+1
      write(6,*) myid, cmplx( re_val(2,1,2,1,1), im_val(2,1,2,1,1), DP )
      write(6,*) myid, startx+min(110,nxpts)-1, startx, nxpts
      write(6,*) myid, cmplx( re_val(min(110,nxpts),1,2,1,1), im_val(min(110,nxpts),1,2,1,1), DP )
    endif

    ! DROT is "backwards" from how we want the phases to go hence minus sign in definition of im_phase
    do ibw = 1, nbw
      do ispn = 1, nspn
        do ik = 1, nkpts
          do ix = 1, nxpts
            call DROT( nbv, re_val( ix, 1, ik, ispn, ibw ), nxpts_pad, & 
                            im_val( ix, 1, ik, ispn, ibw ), nxpts_pad, &
                       re_phase( ix, ik ), im_phase( ix, ik ) )
            call DROT( nbc, re_con( ix, 1, ik, ispn, ibw ), nxpts_pad, &
                            im_con( ix, 1, ik, ispn, ibw ), nxpts_pad, &
                       re_phase( ix, ik ), im_phase( ix, ik ) )
          enddo
        enddo
      enddo
    enddo

    if( .false. ) then
  !    write(6,*) myid, cmplx( re_val(1,1,1,1), im_val(1,1,1,1), DP )
      write(6,*) myid, cmplx( re_val(2,1,2,1,1), im_val(2,1,2,1,1), DP )
      write(6,*) myid, cmplx( re_val(min(110,nxpts),1,2,1,1), im_val(min(110,nxpts),1,2,1,1), DP )
    endif

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
        xfac( 1, iix - startx + 1 ) = real( 1 + ( iix - 1) / (sys%xmesh(2)*sys%xmesh(3)), dp )/ real( sys%xmesh( 1 ), dp )
        xfac( 2, iix - startx + 1 ) = real( 1 + mod( ( iix - 1 )/sys%xmesh(3), sys%xmesh(2 )), dp )/ real( sys%xmesh( 2 ), dp ) 
        xfac( 3, iix - startx + 1 ) = real( 1 + mod( iix - 1, sys%xmesh(2)*sys%xmesh(3) ), dp ) / real( sys%xmesh( 3 ), dp )
      enddo
    else

      ! This is the "dumb" way, but should work for now. Only called once, so not a terrible thing, 
      !  but it will become more burdensome as the total number of x-points grows
      iix = 0
      do ix = 1, sys%xmesh( 1 )
        do iy = 1, sys%xmesh( 2 )
          do iz = 1, sys%xmesh( 3 )
            iix = iix + 1
            if( iix .ge. startx .and. ( iix - startx .lt. nxpts ) ) then
              xfac( 1, iix - startx + 1 ) = dble( ix - 1 ) / dble( sys%xmesh(1) )
              xfac( 2, iix - startx + 1 ) = dble( iy - 1 ) / dble( sys%xmesh(2) )
              xfac( 3, iix - startx + 1 ) = dble( iz - 1 ) / dble( sys%xmesh(3) )
            endif
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


    allocate( re_val( max( 1, nxpts_pad), val_pad, nkpts, nspn, nbw ), &
              im_val( max( 1, nxpts_pad), val_pad, nkpts, nspn, nbw ), &
              re_con( max( 1, nxpts_pad), con_pad, nkpts, nspn, nbw ), &
              im_con( max( 1, nxpts_pad), con_pad, nkpts, nspn, nbw ), STAT=ierr )
    if( ierr .ne. 0 ) return
!    re_val(:,:,:,:) = 0.0_DP
!    im_val(:,:,:,:) = 0.0_DP
!    re_con(:,:,:,:) = 0.0_DP
!    im_con(:,:,:,:) = 0.0_DP

    call OCEAN_val_states_read( sys, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) write( 6, * ) 'Adding phase'
    allocate( re_phase( nxpts_pad, nkpts ), im_phase( nxpts_pad, nkpts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call val_states_generate_phases( sys, re_phase, im_phase, ierr )
    if( ierr .ne. 0 ) return 

    call val_states_add_phase( re_phase, im_phase )

    deallocate( re_phase, im_phase )

    !TODO make either or with dp?
    if( sys%use_sp ) then
      if( nbw .ne. 1 ) then
        write(6,*) 'SP and BWFLAG not enabled'
        ierr = 101112
        return
      endif
      allocate( re_val_sp( max( 1, nxpts_pad), val_pad, nkpts, nspn ), &
                im_val_sp( max( 1, nxpts_pad), val_pad, nkpts, nspn ), &
                re_con_sp( max( 1, nxpts_pad), con_pad, nkpts, nspn ), &
                im_con_sp( max( 1, nxpts_pad), con_pad, nkpts, nspn ), STAT=ierr )
      re_val_sp(:,:,:,:) = re_val(:,:,:,:,1)
      im_val_sp(:,:,:,:) = im_val(:,:,:,:,1)
      re_con_sp(:,:,:,:) = re_con(:,:,:,:,1)
      im_con_sp(:,:,:,:) = im_con(:,:,:,:,1)
    else
      allocate( re_val_sp(0,0,0,0), im_val_sp(0,0,0,0), re_con_sp(0,0,0,0), im_con_sp(0,0,0,0) )
    endif
    if( ierr .ne. 0 ) return

    is_loaded = .true.

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

    if( sys%bwflg ) then
      nbw = 2
    else
      nbw = 1
    endif

    nxpts = 0
    startx = 1
    nx_remain = sys%nxpts

    allocate( nxpts_by_mpiID( 0:nproc-1 ), startx_by_mpiID( 0:nproc-1 ) )

    if( .false. ) then
      max_nxpts = ceiling( dble( sys%nxpts ) / dble( nproc ) )
      
      nxpts = max_nxpts
      nxpts_by_mpiID( 0 ) = nxpts
      startx_by_mpiID( 0 ) = 1
      do i = 1, nproc - 1
        nx_remain = nx_remain - nxpts
        nxpts = min( max_nxpts, nx_remain )
        nxpts_by_mpiID( i ) = nxpts
        startx_by_mpiID( i ) = startx_by_mpiID( i - 1 ) + nxpts_by_mpiID( i - 1 )
        if( nxpts .eq. 0 ) startx_by_mpiID( i ) = 1
      enddo

      nxpts = nxpts_by_mpiID( myid )
      startx = startx_by_mpiID( myid )

    else
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
    endif

    !JTV Ladder needs modification to deal with padding
    nxpts_pad = nxpts

    if( nxpts .lt. 0 ) then
      ierr = -1
      return
    endif

    con_pad = nbc + mod( nbc, cache_double )
    val_pad = nbv + mod( nbv, cache_double )

    nxpts_pad = nxpts
    con_pad = nbc
    val_pad = nbv

    is_init = .true.

    write(1000+myid,*) 'con_pad  ', con_pad, nbc
    write(1000+myid,*) 'val_pad  ', val_pad, nbv
    write(1000+myid,*) 'nxpts_pad', nxpts_pad, nxpts
    write(1000+myid,*) 'max_nxpts', max_nxpts
    write(1000+myid,*) 'start_x', startx
    write(1000+myid,*) '------------------------------'
    do i = 0, nproc-1
      write(1000+myid,*) startx_by_mpiID( i ), nxpts_by_mpiID( i )
    enddo
    flush(1000+myid)
    

  end subroutine
  

  subroutine OCEAN_val_states_read( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi!, only : myid, nproc, root, comm
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer :: err

    if( myid .eq. root ) then

!JTV is this all redundant?
      if( nbv .gt. (1+sys%brange(2)-sys%brange(1)) ) then
        write(6,*) 'Not enough valence bands!'
        ierr = -1
      endif
      if( nbc .gt. (1+sys%brange(4)-sys%brange(3)) ) then
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


    select case( sys%bloch_selector )
      case( 1 ) 
        call load_new_u2( sys, sys%brange, ierr )
      case( 0 )
        call load_old_u2( sys, sys%brange, ierr )
      case( 2 )
        call load_raw( sys, ierr )
      case( 3 )
        call load_single_prefixu2dat( sys, ierr )
      case default
        ierr = 500
        if( myid .eq. root ) write(6,*) 'Unsupported bloch_selector:', sys%bloch_selector
        return
    end select 

  end subroutine OCEAN_val_states_read

!> @brief Read in con.u2.dat and val.u2.dat from a single master processor and share to all
  subroutine load_single_prefixu2dat( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: readU2(:), transposeU2(:,:,:), share_buffer(:,:,:)
    integer :: ispn, iq, ibd, ii, ix, iy, iz, iproc, file_brange(4), max_band, min_band, allbands
    logical :: ex
#ifdef MPI__F08
    type( MPI_REQUEST ), allocatable :: request(:)
#else
    integer, allocatable :: request(:)
#endif



!TODO: need to figure out nbv+nbc, allowing for overlapping bands and a minimum that isn't 1
! Then allocate share_buffer to be of size allbands
! read in all the bands, store to the correct con and val sections
! con backwards stores the unocc from val.u2.dat
!

    if( myid .eq. root ) then
!      inquire(file='val.control.txt', exist=ex )
!      if( ex ) then
!        open(unit=99,file='val.control.txt', form='formatted', status='old')
!        read(99,*) file_brange(1:2)
!        read(99,*) file_brange(3:4)
!        close(99)
!      else
        file_brange(1:4) = sys%file_brange(1:4)
!      endif
      open(unit=99, file='val.u2.dat', form='unformatted', status='old', access='stream' )
    endif

    min_band = sys%brange(1)
    if( sys%bwflg ) then
      max_band = sys%brange(4)
      allbands = sys%brange(4) - sys%brange(1) + 1
    else
      allbands = nbv
      max_band = sys%brange(2)
    endif

    if( myid .eq. root ) then
      allocate( readU2( sys%nxpts ), transposeU2( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1) ), &
                share_buffer( max_nxpts, max(allbands,nbc), 0:nproc-1 ), request(0:nproc), STAT=ierr )
      request(:) = MPI_REQUEST_NULL
    else
      allocate( share_buffer( max_nxpts, max(allbands,nbc), 1 ) )
    endif
    if( ierr .ne. 0 ) return


    do ispn = 1, sys%nspn
      do iq = 1, sys%nkpts

        if( myid .eq. root ) then
          ! currently the valence code expects the real-space to be stored (z,y,x)!
          ! con.u2.dat and val.u2.dat store it (x,y,z)
!          do ibd = file_brange(1), sys%brange(1) - 1
!            read( 99 ) readU2
!          enddo

!          do ibd = 1, nbv
          do ibd = file_brange(1), file_brange(2)
            read( 99 ) readU2
            if( ibd .gt. max_band .or. ibd .lt. min_band ) cycle

            ii = 0
            do iz = 1, sys%xmesh(3)
              do iy = 1, sys%xmesh(2)
                do ix = 1, sys%xmesh(1)
                  ii = ii + 1
                  transposeU2( iz, iy, ix ) = readU2( ii )
                enddo
              enddo
            enddo

            iproc = 0
            ii = 0
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  ii = ii + 1
                  if( ii .gt. nxpts_by_mpiID( iproc ) ) then
                    iproc = iproc + 1
                    ii = 1
                  endif 
                  share_buffer( ii, ibd - min_band + 1, iproc ) = transposeU2( iz, iy, ix )
                enddo
              enddo
            enddo

          enddo ! ibd

          
        endif

        if( myid .eq. root ) then
          do iproc = 0, nproc-1
            if( iproc .ne. myid ) then
              call MPI_ISEND( share_buffer(:,:,iproc), max_nxpts*allbands, MPI_DOUBLE_COMPLEX, &
                              iproc, 1, comm, request( iproc ), ierr )
            endif
          enddo
          do ibd = 1, nbv
            re_val( 1:nxpts, ibd, iq, ispn, 1 ) = real( share_buffer( 1:nxpts, ibd, myid ), DP )
            im_val( 1:nxpts, ibd, iq, ispn, 1 ) = aimag( share_buffer( 1:nxpts, ibd, myid ) )
          enddo
          if( sys%bwflg ) then
            do ii = 1, nbc
              ibd = ii + sys%brange(3) - sys%brange(1)
              re_con( 1:nxpts, ii, iq, ispn, 2 ) = real( share_buffer( 1:nxpts, ibd, myid ), DP )
              im_con( 1:nxpts, ii, iq, ispn, 2 ) = aimag( share_buffer( 1:nxpts, ibd, myid ) )
            enddo
          endif
          call MPI_WAITALL( nproc, request, MPI_STATUSES_IGNORE, ierr )
        else
          call MPI_RECV( share_buffer, max_nxpts*allbands, MPI_DOUBLE_COMPLEX, &
                         root, 1, comm, MPI_STATUS_IGNORE, ierr )
          do ibd = 1, nbv
            re_val( 1:nxpts, ibd, iq, ispn, 1 ) = real( share_buffer( 1:nxpts, ibd, 1 ), DP )
            im_val( 1:nxpts, ibd, iq, ispn, 1 ) = aimag( share_buffer( 1:nxpts, ibd, 1 ) )
          enddo
          if( sys%bwflg ) then
            do ii = 1, nbc
              ibd = ii + sys%brange(3) - sys%brange(1)
              re_con( 1:nxpts, ii, iq, ispn, 2 ) = real( share_buffer( 1:nxpts, ibd, 1 ), DP )
              im_con( 1:nxpts, ii, iq, ispn, 2 ) = aimag( share_buffer( 1:nxpts, ibd, 1 ) )
            enddo
          endif
        endif

      enddo
    enddo

    if( myid .eq. root ) then
      close(99)

      open(unit=99, file='con.u2.dat', form='unformatted', status='old', access='stream' )
    endif

    if( sys%bwflg ) then
      min_band = sys%brange(1)
      max_band = sys%brange(4)
      allbands = sys%brange(4) - sys%brange(1) + 1
    else    
      min_band = sys%brange(3)
      allbands = nbc
      max_band = sys%brange(4)
    endif  

    do ispn = 1, sys%nspn
      do iq = 1, sys%nkpts

        if( myid .eq. root ) then
          ! currently the valence code expects the real-space to be stored (z,y,x)!
          ! con.u2.dat and val.u2.dat store it (x,y,z)
          do ibd = file_brange(3), file_brange(4)
            read( 99 ) readU2
            if( ibd .gt. max_band .or. ibd .lt. min_band ) cycle

            ii = 0
            do iz = 1, sys%xmesh(3)
              do iy = 1, sys%xmesh(2)
                do ix = 1, sys%xmesh(1)
                  ii = ii + 1
                  transposeU2( iz, iy, ix ) = readU2( ii )
                enddo
              enddo
            enddo

            iproc = 0
            ii = 0
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  ii = ii + 1
                  if( ii .gt. nxpts_by_mpiID( iproc ) ) then
                    iproc = iproc + 1
                    ii = 1
                  endif
                  share_buffer( ii, ibd-min_band+1, iproc ) = transposeU2( iz, iy, ix )
                enddo
              enddo
            enddo

          enddo ! ibd
        endif

        if( myid .eq. root ) then
          do iproc = 0, nproc-1
            if( iproc .ne. myid ) then
              call MPI_ISEND( share_buffer(:,:,iproc), max_nxpts*allbands, MPI_DOUBLE_COMPLEX, &
                              iproc, 1, comm, request( iproc ), ierr )
            endif
          enddo
          if( sys%bwflg ) then
            do ibd = 1, nbv
              re_val( 1:nxpts, ibd, iq, ispn, 2 ) = real( share_buffer( 1:nxpts, ibd, myid ), DP )
              im_val( 1:nxpts, ibd, iq, ispn, 2 ) = aimag( share_buffer( 1:nxpts, ibd, myid ) )
            enddo
          endif
          do ibd = 1, nbc
            re_con( 1:nxpts, ibd, iq, ispn, 1 ) = real( share_buffer( 1:nxpts, ibd+sys%brange(3)-min_band, myid ), DP )
            im_con( 1:nxpts, ibd, iq, ispn, 1 ) = aimag( share_buffer( 1:nxpts, ibd+sys%brange(3)-min_band, myid ) )
          enddo
          call MPI_WAITALL( nproc, request, MPI_STATUSES_IGNORE, ierr )
        else
          call MPI_RECV( share_buffer, max_nxpts*allbands, MPI_DOUBLE_COMPLEX, &
                         root, 1, comm, MPI_STATUS_IGNORE, ierr )
          if( sys%bwflg ) then
            do ibd = 1, nbv
              re_val( 1:nxpts, ibd, iq, ispn, 2 ) = real( share_buffer( 1:nxpts, ibd, 1 ), DP )
              im_val( 1:nxpts, ibd, iq, ispn, 2 ) = aimag( share_buffer( 1:nxpts, ibd, 1 ) )
            enddo
          endif
          do ibd = 1, nbc
            re_con( 1:nxpts, ibd, iq, ispn, 1 ) = real( share_buffer( 1:nxpts, ibd+sys%brange(3)-min_band, 1 ), DP )
            im_con( 1:nxpts, ibd, iq, ispn, 1 ) = aimag( share_buffer( 1:nxpts, ibd+sys%brange(3)-min_band, 1 ) )
          enddo
        endif
      enddo
    enddo

    deallocate( share_buffer )
    if( myid .eq. 0 ) then
      close(99)
      deallocate( readU2, transposeU2, request )
    endif
  

  end subroutine load_single_prefixu2dat

!> @brief Read in DFT states directly from u(G) form instead of u(x)
  subroutine load_raw( sys, ierr )
    use OCEAN_system, only : O_system
    
    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

  end subroutine load_raw

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

      re_val( 1:nxpts, 1:nbv, iq, 1, 1 ) = re_share_buffer( 1:nxpts, 1:nbv, iproc )
      im_val( 1:nxpts, 1:nbv, iq, 1, 1 ) = im_share_buffer( 1:nxpts, 1:nbv, iproc )



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

      re_con( 1:nxpts, 1:nbc, iq, 1, 1 ) = re_share_buffer( 1:nxpts, 1:nbc, iproc )
      im_con( 1:nxpts, 1:nbc, iq, 1, 1 ) = im_share_buffer( 1:nxpts, 1:nbc, iproc )

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

    fhu2 = 99

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
          do ibd = 1, brange(1) - 1
            do ix = 1, sys%xmesh(1)
              do iy = 1, sys%xmesh(2)
                do iz = 1, sys%xmesh(3)
                  read ( fhu2 ) 
                end do
              end do
            end do
          enddo
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

        re_val( 1:nxpts, 1:nbv, iq, ispn, 1 ) = re_share_buffer( 1:nxpts, 1:nbv, iproc )
        im_val( 1:nxpts, 1:nbv, iq, ispn, 1 ) = im_share_buffer( 1:nxpts, 1:nbv, iproc )


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

        re_con( 1:nxpts, 1:nbc, iq, ispn, 1 ) = re_share_buffer( 1:nxpts, 1:nbc, iproc )
        im_con( 1:nxpts, 1:nbc, iq, ispn, 1 ) = im_share_buffer( 1:nxpts, 1:nbc, iproc )

      enddo ! iq
    enddo  ! ispn

    deallocate( re_share_buffer, im_share_buffer, ur, ui )

    if( myid .eq. root ) then
      close( fhu2 )
    endif

  end subroutine load_old_u2

end module OCEAN_val_states
