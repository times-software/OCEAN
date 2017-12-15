! Copyright (C) 2015, 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_bloch
  use AI_kinds

  save
  private


  ! these guys have no position phasing
  real(DP), pointer :: re_bloch_state( :, :, :, : )
  real(DP), pointer :: im_bloch_state( :, :, :, : )


  ! these guys do. There is no reason for them not to
  real(DP), pointer :: re_v_bloch_state( :, :, : )
  real(DP), pointer :: im_v_bloch_state( :, :, : )

#ifdef __INTEL
!dir$ attributes align:64 :: re_bloch_state, im_bloch_state
!dir$ attributes align:64 :: re_v_bloch_state, im_v_bloch_state
#endif


  
  INTEGER :: my_xpts
  INTEGER :: my_kpts
  INTEGER :: my_num_bands
  INTEGER :: my_val_bands
  INTEGER :: my_start_nx

  INTEGER :: xshift_override( 3 )
  
  LOGICAL :: is_init = .false.
  LOGICAL :: is_val_init = .false.
  LOGICAL :: is_loaded = .false.
  LOGICAL :: is_val_loaded = .false.


  public :: OCEAN_bloch_init, OCEAN_bloch_load, OCEAN_bloch_lrINIT, OCEAN_bloch_lrLOAD, OCEAN_bloch_is_loaded

  contains

  logical function OCEAN_bloch_is_loaded()
    OCEAN_bloch_is_loaded = is_loaded
    return
  end function OCEAN_bloch_is_loaded

  subroutine OCEAN_bloch_lrLOAD( sys, tau, xshift, rbs_out, ibs_out, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    real(DP), intent(out) :: rbs_out(my_num_bands,my_kpts,my_xpts,sys%nspn) 
    real(DP), intent(out) :: ibs_out(my_num_bands,my_kpts,my_xpts,sys%nspn) 
    real(DP), intent(inout) :: tau(3)
    integer, intent( out ) :: xshift( 3 )

    real(DP) :: cphs, sphs, pi, phsx, phsy, phsz
    integer :: iq, iq1, iq2, iq3
    integer :: ix, iy, iz, xtarg, ytarg, ztarg, xph, yph, zph
    integer :: xiter, ibd, ispn


    ! Change phase as if things were flipped around. Don't actually re-distribute
    rbs_out = 0.0_dp
    ibs_out = 0.0_dp
    pi = 4.0d0 * atan( 1.0d0 )
    !
    ! Reasoning for the below;
    !  If the kmesh is odd, then the FT of the kmesh will give an equal number of
    !  cells on either side of the central one -- the one with the core hole. 
    !  This means that the core hole should be placed as near as (0.5,0.5,0.5)
    !
    !  On the other hand, if the kmesh is even, the bonus cell will be tacked on 
    !  over toward the negative side. This means that we should center the core 
    !  hole as near as possible to (0,0,0)
    !
    !  We test each dimension individually of course. 
    !
    xshift(:) = 0
    if( mod( sys%kmesh(1), 2 ) .eq. 0 ) then
      xshift( 1 ) = floor( real(sys%xmesh(1), DP ) * tau(1) )
    else
      xshift( 1 ) = floor( real(sys%xmesh(2), DP ) * (tau(1)-0.5d0 ) )
    endif
    if( mod( sys%kmesh(2), 2 ) .eq. 0 ) then
      xshift( 2 ) = floor( real(sys%xmesh(2), DP ) * tau(2) )
    else
      xshift( 2 ) = floor( real(sys%xmesh(2), DP ) * (tau(2)-0.5d0 ) )
    endif
    if( mod( sys%kmesh(3), 2 ) .eq. 0 ) then
      xshift( 3 ) = floor( real(sys%xmesh(3), DP ) * tau(3) )
    else
      xshift( 3 ) = floor( real(sys%xmesh(3), DP ) * (tau(3)-0.5d0 ) )
    endif
    ! 
    xshift(:) = xshift(:) * xshift_override(:)
    !
    if( myid .eq. root ) write(6,*) 'Shifting X-grid by ', xshift(:)
    if( myid .eq. root ) write(6,*) 'Original tau ', tau(:)
    tau( 1 ) = tau(1) - real(xshift(1), DP )/real(sys%xmesh(1), kind( 1.0d0 ))
    tau( 2 ) = tau(2) - real(xshift(2), DP )/real(sys%xmesh(2), kind( 1.0d0 ))
    tau( 3 ) = tau(3) - real(xshift(3), DP )/real(sys%xmesh(3), kind( 1.0d0 ))
    if( myid .eq. root ) write(6,*) 'New tau      ', tau(:)


!    nx_left = sys%nxpts
!    nx_start = 1
!    do i = 0, nproc - 1
!      nx_tmp = nx_left / ( nproc - i )
!      nx_left = nx_left - nx_tmp
!      if( myid .eq. i ) my_start = nx_start
!      nx_start = nx_start + nx_tmp
!    enddo

    do ispn = 1, sys%nspn

      xiter = 0
      do iz = 1, sys%xmesh(3)
        ztarg = iz - xshift( 3 )
        if( ztarg .gt. sys%xmesh(3) ) then
          ztarg = ztarg - sys%xmesh(3)
          zph = -sys%xmesh(3)
        elseif( ztarg .lt. 1 ) then
          ztarg = ztarg + sys%xmesh(3)
          zph = sys%xmesh(3)
        else
          zph = 0
        endif
        do iy = 1, sys%xmesh(2)
          ytarg = iy - xshift( 2 )
          if( ytarg .gt. sys%xmesh(2) ) then
            ytarg = ytarg - sys%xmesh(2)
            yph = -sys%xmesh(2)
          elseif( ytarg .lt. 1 ) then
            ytarg = ytarg + sys%xmesh(2)
            yph = sys%xmesh(2)
          else
            yph = 0
          endif 
          do ix = 1, sys%xmesh(1)
            xiter = xiter + 1
            if( xiter .lt. my_start_nx ) cycle
  !JTV ???
            if( xiter .gt. my_start_nx + my_xpts - 1 ) goto 111
            xtarg = ix - xshift( 1 )
            if( xtarg .gt. sys%xmesh(1) ) then
              xtarg = xtarg - sys%xmesh(1)
              xph = -sys%xmesh(1)
            elseif( xtarg .lt. 1 ) then
              xtarg = xtarg + sys%xmesh(1) 
              xph = sys%xmesh(1)
            else
            xph = 0
            endif 
            
            iq = 0
            do iq1 = 1, sys%kmesh( 1 )
              do iq2 = 1, sys%kmesh( 2 )
                do iq3 = 1, sys%kmesh( 3 )
                  iq = iq + 1

                  phsx = 2.0d0 * pi * dble( ( xph + ix - 1 ) * ( iq1 - 1 ) ) / dble( sys%xmesh(1) * sys%kmesh( 1 ) )
                  phsy = 2.0d0 * pi * dble( ( yph + iy - 1 ) * ( iq2 - 1 ) ) / dble( sys%xmesh(2) * sys%kmesh( 2 ) )
                  phsz = 2.0d0 * pi * dble( ( zph + iz - 1 ) * ( iq3 - 1 ) ) / dble( sys%xmesh(3) * sys%kmesh( 3 ) )
                  cphs = dcos( phsx + phsy + phsz ) 
                  sphs = dsin( phsx + phsy + phsz )
                  !JTV replace this stupid thing with a rotation thingy
                  ! rbs_out = re_bloch_state
                  ! ibs_out = im_bloch_state
                  ! call DROT( nbd, rbs_out, 1, ibs_out, 1, cphs, -sphs )
                  do ibd = 1, my_num_bands
                    rbs_out( ibd, iq, xiter - my_start_nx + 1, ispn )  & 
                          = cphs * re_bloch_state( ibd, iq, xiter - my_start_nx + 1, ispn )  &
                          - sphs * im_bloch_state( ibd, iq, xiter - my_start_nx + 1, ispn ) 
                    ibs_out( ibd, iq, xiter - my_start_nx + 1, ispn )  &
                          = cphs * im_bloch_state( ibd, iq, xiter - my_start_nx + 1, ispn )  &
                          + sphs * re_bloch_state( ibd, iq, xiter - my_start_nx + 1, ispn )
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
      enddo

111 continue
  enddo ! spin

          


  end subroutine OCEAN_bloch_lrLOAD


  subroutine OCEAN_bloch_lrINIT( xpts, kpts, num_bands, start_nx, ierr )
    implicit none
    integer, intent( out ) :: xpts, kpts, num_bands, start_nx
    integer, intent( inout ) :: ierr

    if( .not. is_loaded ) then 
      ierr = -1
      return
    endif
    xpts = my_xpts
    kpts = my_kpts
    num_bands = my_num_bands
    start_nx = my_start_nx


  end subroutine OCEAN_bloch_lrINIT

  subroutine OCEAN_bloch_init( sys, ierr )
    use OCEAN_system
    implicit none

    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

!    if( sys%conduct ) then
      call OCEAN_bloch_con_init( sys, ierr )
      if( ierr .ne. 0 ) return
!    endif

    if( sys%valence ) then
      call OCEAN_bloch_val_init( sys, ierr )
      if( ierr .ne. 0 ) return
    endif

  end subroutine OCEAN_bloch_init
     

  subroutine OCEAN_bloch_con_init( sys, ierr )
    use OCEAN_mpi!, only : myid, nproc, root, comm
    use OCEAN_system
    use iso_c_binding
    implicit none
!    include 'fftw3.f03'

    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr


    integer( S_INT ) :: nx_left, nx_start, nx_tmp, i
    type(C_PTR) :: cptr

    logical :: havefile


    if( is_init ) return

    ! check on xshift override
    xshift_override( : ) = 1
    if( myid .eq. root ) then
      inquire(file='xshift_override.ipt',exist=havefile )
      if( havefile ) then
        open( unit=99,file='xshift_override.ipt',form='formatted',status='old')
        read(99,*) xshift_override(:)
        close(99)
        if( maxval( xshift_override ) .gt. 1 .or. minval( xshift_override ) .lt. 0 ) then
          write(6,*) 'Invalid value in xshift_override.ipt. Reset to 1'
          xshift_override(:) = 1
        endif
      endif
    endif
    call MPI_BCAST( xshift_override, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    ! SOP
    my_num_bands = sys%num_bands
    my_kpts = sys%nkpts
    nx_left = sys%nxpts
    nx_start = 1
    do i = 0, nproc - 1
      nx_tmp = nx_left / ( nproc - i )
      nx_left = nx_left - nx_tmp
      if( myid .eq. root ) write(6,*) i, nx_tmp, nx_left, nx_start
      if( i .eq. myid ) then
        my_xpts = nx_tmp
        my_start_nx = nx_start
      endif
      nx_start = nx_start + nx_tmp
    enddo
    if( myid .eq. root ) then
      write(6,*) 'LR: NB, NK, NX_START, NX'
      write(6,*) my_num_bands, my_kpts, myid
      write(6,*) my_start_nx, my_xpts
    endif
    

!#ifdef __INTEL
    allocate( re_bloch_state( my_num_bands, my_kpts, my_xpts, sys%nspn ), &
              im_bloch_state( my_num_bands, my_kpts, my_xpts, sys%nspn ), &
              STAT=ierr )
    if( ierr .ne. 0 ) return
!#else
!    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts * sys%nspn, C_SIZE_T) )
!    call c_f_pointer( cptr, re_bloch_state, [my_num_bands, my_kpts, my_xpts, sys%nspn] )
!    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts * sys%nspn, C_SIZE_T) )
!    call c_f_pointer( cptr, im_bloch_state, [my_num_bands, my_kpts, my_xpts, sys%nspn] )
!#endif

    is_init = .true.

  end subroutine OCEAN_bloch_con_init

  subroutine OCEAN_bloch_val_init( sys, ierr )
    use OCEAN_mpi!, only : myid, nproc, root, comm
    use OCEAN_system
    use iso_c_binding
    implicit none
!    include 'fftw3.f03'
    
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    

    integer( S_INT ) :: nx_left, nx_start, nx_tmp, i
    type(C_PTR) :: cptr 
    

    if( is_val_init ) return

    ! SOP
    my_val_bands = sys%val_bands

    if( .not. is_init ) then
      my_num_bands = sys%num_bands
      my_kpts = sys%nkpts
      nx_left = sys%nxpts
      nx_start = 1
      do i = 0, nproc - 1
        nx_tmp = nx_left / ( nproc - i )
        nx_left = nx_left - nx_tmp
        if( myid .eq. root ) write(6,*) i, nx_tmp, nx_left, nx_start
        if( i .eq. myid ) then
          my_xpts = nx_tmp
          my_start_nx = nx_start
        endif
        nx_start = nx_start + nx_tmp
      enddo
      if( myid .eq. root ) then
        write(6,*) 'LR: NB, NK, NX_START, NX'
        write(6,*) my_num_bands, my_kpts, myid
        write(6,*) my_start_nx, my_xpts
      endif
    endif


!#ifdef __INTEL
    allocate( re_v_bloch_state( my_num_bands, my_kpts, my_xpts ), &
              im_v_bloch_state( my_num_bands, my_kpts, my_xpts ), &
              STAT=ierr )
    if( ierr .ne. 0 ) return
!#else
!    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts, C_SIZE_T) )
!    call c_f_pointer( cptr, re_v_bloch_state, [my_num_bands, my_kpts, my_xpts] )
!    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts, C_SIZE_T) )
!    call c_f_pointer( cptr, im_v_bloch_state, [my_num_bands, my_kpts, my_xpts] )
!#endif

    is_val_init = .true.

  end subroutine OCEAN_bloch_val_init




  subroutine OCEAN_bloch_load( sys, ierr )
    use OCEAN_mpi!, only : myid, nproc, root, comm
    use OCEAN_system
!    use mpi
    implicit none

    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr


    integer, parameter :: u2dat = 35
    !
    integer :: nx, ny, nz, nbd, nq
!    real( kind = DP  ), dimension( nx, ny, nz, nbd ) :: ur, ui
    !
    integer :: iq, ibd, ig, idum( 3 ), ix, iy, iz, ivl, ivh, icl, ich, ispn, i
    integer :: iq1, iq2, iq3, dumint, icl2, ich2, ivh2, xshift(3)
    integer :: xtarg, ytarg, ztarg, xph, yph, zph
    real( kind = DP  ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
    real( kind = DP  ) :: su, sul, suh
    real( kind = DP  ), allocatable, dimension( :, :, :, : ) :: tmp_ur, tmp_ui, ur, ui
    real( DP ), allocatable :: re_transpose( :, : ), im_transpose( :, : )
    complex( DP ), allocatable :: u2_buf( :, :, :, : )
!    complex( kind = DP  ), allocatable, dimension( :, : ) :: tmp_bloch
    logical :: metal, normal, ex
    integer :: nx_left, nx_start, nx_tmp, xiter, ii, width(3)
    character(len=3) :: bloch_type
    logical :: invert_xmesh

    !
    integer :: fmode, fhu2
    integer(kind=MPI_OFFSET_KIND) :: offset, offset_extra, offset_nx, offset_start
    integer(8) :: time1, time2, tics_per

    integer :: ncount, u2_status(MPI_STATUS_SIZE), u2_type

    if( is_loaded ) return

    if( myid .eq. root ) then
      write(6,*) 'OCEAN_bloch_load'
    endif

    call OCEAN_bloch_init( sys, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) "Bloch init failed"
      return
    endif

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. 0 ) write(6,*) "Finished Bloch init"

    if( ( .not. associated( re_bloch_state ) ) .or. ( .not. associated( im_bloch_state ) ) ) then
      ierr = -1
      if( myid .eq. root ) write( 6, * ) 'Trying to read in band states before allocation'
      return
    endif

    

    if( myid .eq. 0 ) then
      sul = 1.0d0
      suh = 1.0d0
      open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
      rewind 99 
      read ( 99, * ) ivl, ivh, icl, ich
      close( unit=99 )
      ivh2 = ivh 
    !  icl2 = icl
      open( unit=99, file='metal', form='formatted', status='old')
      read( 99, * ) metal
      close( 99 )
      if( sys%legacy_ibeg .and. ( metal .or. ( sys%nspn .eq. 2 ) ) ) then
        open( unit=36, file='ibeg.h', form='formatted', status='old' )
        rewind( 36 )
      endif


      inquire(file="bloch_type",exist=ex)
      if( ex ) then
        open(unit=99,file="bloch_type", form='formatted', status='old' )
        read(99,*) bloch_type
!        read(99,*) invert_xmesh
        invert_xmesh = .true.
        close(99)
      else
        bloch_type = 'old'
        invert_xmesh = .true.
      endif

      nbd = sys%num_bands 
      if ( nbd .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'

      write(6,*) 'nspn: ', sys%nspn
!      if( sys%nspn .ne. 1 ) then
!        write(6,*) 'No spin yet'
!        ierr = 1
!        goto 111 
!      endif

    endif

#ifdef MPI
    call MPI_BCAST( bloch_type, 3, MPI_CHARACTER, root, comm, ierr )
    call MPI_BCAST( invert_xmesh, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    nx = sys%xmesh(1)
    ny = sys%xmesh(2)
    nz = sys%xmesh(3)

    select case( bloch_type )
      case('old')

    allocate( ur(sys%xmesh(1),sys%xmesh(2),sys%xmesh(3),sys%num_bands), &
              ui(sys%xmesh(1),sys%xmesh(2),sys%xmesh(3),sys%num_bands) )
    ! As per usual, do this the dumbest way first, 
    ! 1) generic copy of original serial method 
    ! 2) MPI

      
    if( myid .eq. 0 ) then 
      write(6,*) size(re_bloch_state,1), size(re_bloch_state,2),size(re_bloch_state,3), size(re_bloch_state,4)
      allocate( tmp_ur( nx, ny, nz, nbd ), tmp_ui( nz, ny, nz, nbd ), &
                re_transpose( sys%num_bands, sys%nxpts ), &
                im_transpose( sys%num_bands, sys%nxpts ) )
      open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
      rewind u2dat 
    else
      nx_left = sys%nxpts
      do i = 0, myid
        nx_tmp = nx_left / ( nproc - i )
        nx_left = nx_left - nx_tmp
      enddo
      allocate( re_transpose( my_num_bands, nx_tmp ), im_transpose( my_num_bands, nx_tmp ) )
    endif

    do ispn = 1, sys%nspn
      iq = 0
      do iq1 = 1, sys%kmesh( 1 )
       do iq2 = 1, sys%kmesh( 2 )
        do iq3 = 1, sys%kmesh( 3 )
        iq = iq + 1
  !    do iq = 1, sys%nkpts

          if( myid .eq. 0 ) then
            open( unit=99, file='gumatprog', form='formatted', status='unknown' )
            rewind 99
            write ( 99, '(2i8)' ) iq, sys%nkpts
            close( unit=99 )
            if( sys%legacy_ibeg .and. (metal .or. ( sys%nspn .eq. 2 ) ) ) then
              read( 36, * ) dumint, ivh2
              ivh2 = ivh2 - 1
            endif

      !  Skip all of the occupied bands (and for metals)
            do ibd = ivl, ivh2
              do ig = 1, nx * ny * nz
                 read ( u2dat )
              end do
            end do
            do ibd = 1, nbd
              do ix = 1, nx
                 do iy = 1, ny
                    do iz = 1, nz
                       read ( u2dat ) idum( 1 : 3 ), ur( ix, iy, iz, ibd ), ui( ix, iy, iz, ibd )
                    end do
                 end do
              end do
              su = sum( ur( :, :, :, ibd ) ** 2 + ui( :, :, :, ibd ) ** 2 )
              sul = min( su, sul )
              suh = max( su, suh )
            end do
      ! Adding 22 nov 2010 get rid of the un-used wfns at the top
            do ibd = ivh2 + nbd + 1, ivh - ivl + ich - icl + 2
              do ig = 1, nx * ny * nz
                 read ( u2dat )
              end do
            enddo

          xiter = 0
          do iz = 1, nz
              do iy = 1, ny
                 do ix = 1, nx
                    xiter = xiter + 1
                    do ibd = 1, nbd
                      re_transpose( ibd, xiter ) = ur( ix, iy, iz, ibd )
                      im_transpose( ibd, xiter ) = ui( ix, iy, iz, ibd )
                    end do
                 end do
              end do
           end do
         endif
         if( myid .eq. root .and. mod(iq,10) .eq. 0 ) write(6,*) iq

         nx_left = sys%nxpts
         nx_start = 1
         do i = 0, nproc - 1
            nx_tmp = nx_left / ( nproc - i )
            nx_left = nx_left - nx_tmp
            if( myid .eq. root .and. iq .eq. 1 ) write(6,*) i, nx_start, nx_tmp
            if( i .eq. root .and. myid .eq. root ) then
              re_bloch_state( :, iq, :, ispn ) = re_transpose( :, nx_start : nx_start + nx_tmp - 1 )
              im_bloch_state( :, iq, :, ispn ) = im_transpose( :, nx_start : nx_start + nx_tmp - 1 )
#ifdef MPI
            elseif( myid .eq. root ) then
              call MPI_SEND( re_transpose(1,nx_start), my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i, comm, ierr )
              call MPI_SEND( im_transpose(1,nx_start), my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i+nproc, comm, ierr )
            elseif( myid .eq. i ) then
!              if( iq .eq. 1 .and. ispn .eq. 1 ) then
!                allocate( re_transpose( my_num_bands, nx_tmp ), im_transpose( my_num_bands, nx_tmp ) )
!              endif

              call MPI_RECV( re_transpose, my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                            i, comm, MPI_STATUS_IGNORE, ierr )
              call MPI_RECV( im_transpose, my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                            i+nproc, comm, MPI_STATUS_IGNORE, ierr )
              re_bloch_state( :, iq, :, ispn ) = re_transpose( :, : )
              im_bloch_state( :, iq, :, ispn ) = im_transpose( :, : )
#endif
            endif
           nx_start = nx_start + nx_tmp
         enddo



        enddo
       enddo
      enddo
    enddo

    if( myid .eq. 0 ) close(u2dat)

    deallocate( re_transpose, im_transpose, ur, ui )
    if( myid .eq. 0 ) deallocate( tmp_ur, tmp_ui )

!!!!!!!!!!!!!!!!!!!!
    case( 'new' )

#ifndef MPI
      if( myid .eq. 0 ) write(6,*) "New requires MPI. Attempting old"
      bloch_type = 'old'
      ierr = -1
!      goto 112
#else
    call load_new_u2(  sys, ierr, 36, metal, ivh2 )
    if( ierr .ne. 0 ) then
      if( myid .eq. 0 ) write(6,*) 'new_u2 failed'
      return
    endif

#endif
    case( 'obf' )

      ierr = -2

#ifndef MPI
      if( myid .eq. root ) write(6,*) 'OBF requires MPI'
      ierr = -1
      return
#else

! Called by obf_load
!      call OCEAN_obf_init( sys, ierr )
!      if( ierr .ne. 0 ) return


! Want to change to have dynamic FFT here I think, or at least the option
!      call OCEAN_obf_load( sys, ierr )
!      if( ierr .ne. 0 ) return
!
!      call OCEAN_obf_obf2bloch( sys, re_bloch_state, im_bloch_state, &
!                                my_num_bands, my_kpts, my_xpts, ierr )
      

#endif

!!!!!!!!!!!!!!!
    case default
      if( myid .eq. 0 ) write(6,*) 'Unrecognized BLOCH type. Attempting type = old'
      bloch_type = 'old'

    end select

!111 continue


    if( myid .eq. 0 ) write ( 6, '(1a16,2f20.15)' ) 'norm bounds ... ', sul, suh

    if( myid .eq. 0 .and. ( metal .or. ( sys%nspn .eq. 2 ) ) ) close( 36 )

    is_loaded = .true.

  end subroutine OCEAN_bloch_load


  subroutine load_new_u2( sys, ierr, fh_ibeg, metal, ivh2 )
    use OCEAN_mpi!, only : myid, nproc, root, comm
    use OCEAN_system
!    use mpi
    implicit none
    
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: fh_ibeg !ibeg( sys%nkpts, sys%nspn )
    logical, intent( in ) :: metal
    integer, intent( inout ) :: ivh2
    integer, intent( inout ) :: ierr

    real( DP ), allocatable :: re_transpose( :, : ), im_transpose( :, : )
    complex( DP ), allocatable :: u2_buf( :, :, :, : )
    
    integer :: ix, iy, iz, nx_left, nx_start, nx_tmp, xiter, iq, ibd, i, ispn

    integer :: width(3), iq_ten
    integer :: fmode, fhu2, u2_type, u2_status(MPI_STATUS_SIZE), ncount, io_comm, nelement
    integer(kind=MPI_OFFSET_KIND) :: offset, offset_extra, offset_nx, offset_start
    integer(8) :: time1, time2, tics_per
    
    logical :: io_group = .false.

    logical :: io_style = .false.

    
    call MPI_BARRIER( comm, ierr )

    iq_ten = max(1,sys%nkpts/10)

    
    if( myid .eq. root ) then
      write(6,*) 'New U2 format'

      open(unit=99,file='obf_control',form='formatted',status='old')
      rewind(99)
      read(99,*) width(1)
      read(99,*) width(2)
      read(99,*) width(3)
      close(99)


      write(6,*) my_start_nx, my_xpts

      if( sys%kshift ) then
        write(6,*) 'Two k-grids in u2'
      endif
    endif

    call MPI_BCAST( width, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return


    if( myid .eq. root ) then
      write(6,*) width(:)
    endif


    ! Create i/o comm
    if( myid .eq. root ) then
      io_group = .true.
    endif
  
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
      !JTV At this point it would be good to create a custom MPI_DATATYPE
      !  so that we can get optimized file writing
      call MPI_TYPE_CONTIGUOUS( sys%nxpts, MPI_DOUBLE_COMPLEX, u2_type, ierr )
      if( ierr/=0 ) then
        return
      endif
      call MPI_TYPE_COMMIT( u2_type, ierr )
      if( ierr/=0 ) then
        return
      endif

      if( io_style) then
        call MPI_FILE_SET_VIEW( fhu2, offset, u2_type, u2_type, 'native', MPI_INFO_NULL, ierr )
      else
        call MPI_FILE_SET_VIEW( fhu2, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr )
      endif
      if( ierr/=0 ) then
        return
      endif
    endif



!JTV change up for valence/conduction choice

    call system_clock( time1, tics_per, time2 )

    if( myid .eq. root ) then
      allocate( re_transpose( sys%num_bands, sys%nxpts ), &
                im_transpose( sys%num_bands, sys%nxpts ) )
      allocate( u2_buf( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), sys%num_bands ) )
    endif

!  u2 is stored by k-point. If there is a shift then
!    we store bands 1 - max_val :/ unshifted
!             bands 1 - max :/ shifted
!    
    if( sys%kshift ) then
      if( sys%conduct ) then
        offset_start = width(2)
        offset_extra = 0
      else 
        offset_start = 0
        offset_extra = width(2)
      endif
!            offset_extra = width(2)
!          else
!            offset_extra = width3
!          endif
    else
      offset_extra = 0
      offset_start = 0
    endif

    if( myid .eq. root ) then
      write(6,*) offset_start, offset_extra
    endif


    do ispn = 1, sys%nspn
     do iq = 1, sys%nkpts

      if( myid .eq. root ) then

        if( sys%legacy_ibeg .and. (metal .or. ( sys%nspn .eq. 2 ) ) ) then
          read(fh_ibeg, * ) i, ivh2
!          ivh2 = ibeg( iq, 1 ) - 1
          ivh2 = ivh2 - 1
        endif
!        if( mod(iq,10) .eq. 0 ) write(6,*) iq, ivh2
        if( mod(iq,iq_ten) .eq. 0 ) write(6,*) iq, ivh2


        offset = int( ( iq - 1 ), MPI_OFFSET_KIND ) * int( width( 3 ), MPI_OFFSET_KIND ) &
               + int( (ispn - 1)*sys%nkpts, MPI_OFFSET_KIND ) * int( width( 3 ), MPI_OFFSET_KIND )
        if( .not. io_style ) then
          offset = offset * int( sys%nxpts, MPI_OFFSET_KIND )
        endif

        if( sys%valence ) then
          exit ! I need to fix

          do ibd = 1, ivh2
            call MPI_FILE_READ_AT( fhu2, offset, u2_buf(1,1,1,ibd), 1, u2_type, u2_status, ierr )
            if( ierr .ne. 0 ) then
              write(6,*) 'u2 failed', ierr, iq, ispn
              return
            endif
            call MPI_GET_COUNT( u2_status, u2_type, ncount, ierr )
            if( ncount .ne. 1 ) then
              write(6,*) 'u2 read failed.', ncount, sys%num_bands, iq, ispn
              ierr = -100
              return
            endif
            offset = offset + 1
          enddo

          ! zero out bonus bands
          do ibd = ivh2 + 1, sys%val_bands
            u2_buf(:,:,:,ibd) = 0.0_dp
          enddo

        else
          if( io_style ) then
            offset = offset + int(ivh2,MPI_OFFSET_KIND)
          else
            offset = offset + int(ivh2,MPI_OFFSET_KIND) * int(sys%nxpts, MPI_OFFSET_KIND )
          endif
        endif

      if( io_style .and. .true. ) then
        call mpi_file_read_at( fhu2, offset, u2_buf, sys%num_bands, u2_type, u2_status, ierr )
        if( ierr .ne. 0 ) then
          write(6,*) 'u2 failed', ierr, iq, ispn
          return
        endif
        call MPI_GET_COUNT( u2_status, u2_type, ncount, ierr )
        if( ncount .ne. sys%num_bands ) then
          write(6,*) 'u2 read failed.', ncount, sys%num_bands, iq, ispn
        endif
      elseif( io_style ) then
        do ibd = 1, sys%num_bands
          call mpi_file_read_at( fhu2, offset, u2_buf(1,1,1,ibd), 1, u2_type, u2_status, ierr )
          if( ierr .ne. 0 ) then
            write(6,*) 'u2 failed', ierr, iq, ispn
            return
          endif
          call MPI_GET_COUNT( u2_status, u2_type, ncount, ierr )
          if( ncount .ne. 1 ) then
            write(6,*) 'u2 read failed.', ncount, sys%num_bands, iq, ispn
            ierr = -100
            return
          endif
          offset = offset + 1
        enddo
      else
        do ibd = 1, sys%num_bands, 512
          nelement = min( sys%num_bands - ibd + 1, 512 )
          nelement = nelement * sys%nxpts

          call MPI_FILE_READ_AT( fhu2, offset, u2_buf(1,1,1,ibd), nelement, MPI_DOUBLE_COMPLEX, &
                                 u2_status, ierr )
          if( ierr .ne. 0 ) then
            write(6,*) 'u2 failed', ierr, iq, ispn
            return
          endif
          call MPI_GET_COUNT( u2_status, MPI_DOUBLE_COMPLEX, ncount, ierr )
          if( ncount .ne. nelement) then
            write(6,*) 'u2 read failed.', ncount, nelement, iq, ispn
  !          ierr = -100
  !          return
          endif

          offset = offset + nelement

        enddo

      endif

        offset = offset + ( width(3) - ivh2 - sys%num_bands) 
        xiter = 0
        do iz = 1, sys%xmesh(3)
          do iy = 1, sys%xmesh(2)
            do ix = 1, sys%xmesh(1)
              xiter = xiter + 1
              do ibd = 1, sys%num_bands
                re_transpose( ibd, xiter ) = real( u2_buf( iz, iy, ix, ibd ), DP )
                im_transpose( ibd, xiter ) = aimag( u2_buf( iz, iy, ix, ibd ) )
              end do
            end do
          end do
        end do

      endif

     nx_left = sys%nxpts
     nx_start = 1 
     do i = 0, nproc - 1
       nx_tmp = nx_left / ( nproc - i )
       nx_left = nx_left - nx_tmp
       if( myid .eq. root .and. iq .eq. 1 ) write(6,*) i, nx_start, nx_tmp
       if( i .eq. root .and. myid .eq. root ) then
         re_bloch_state( :, iq, :, ispn ) = re_transpose( :, nx_start : nx_start + nx_tmp - 1 )
         im_bloch_state( :, iq, :, ispn ) = im_transpose( :, nx_start : nx_start + nx_tmp - 1 )
#ifdef MPI
       elseif( myid .eq. root ) then
         call MPI_SEND( re_transpose(1,nx_start), my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i, comm, ierr )
         call MPI_SEND( im_transpose(1,nx_start), my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i+nproc, comm, ierr )
       elseif( myid .eq. i ) then
         if( iq .eq. 1 .and. ispn .eq. 1 ) then
           allocate( re_transpose( my_num_bands, nx_tmp ), im_transpose( my_num_bands, nx_tmp ) )
         endif

         call MPI_RECV( re_transpose, my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                        i, comm, MPI_STATUS_IGNORE, ierr )
         call MPI_RECV( im_transpose, my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                        i+nproc, comm, MPI_STATUS_IGNORE, ierr )
         re_bloch_state( :, iq, :, ispn ) = re_transpose( :, : )
         im_bloch_state( :, iq, :, ispn ) = im_transpose( :, : )
#endif
       endif
       nx_start = nx_start + nx_tmp
      enddo
     enddo
    enddo 

  deallocate( re_transpose, im_transpose )

  call system_clock( time2 )

  if( io_group ) then
    call MPI_FILE_CLOSE( fhu2, ierr )
  endif


  if( myid .eq. root ) deallocate( u2_buf )

  if( myid .eq. root ) write(6,*) 'Read-in took ', dble( time2-time1 ) / dble( tics_per )


  end subroutine load_new_u2

end module OCEAN_bloch
