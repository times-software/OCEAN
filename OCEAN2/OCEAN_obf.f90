! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_obf
  use AI_kinds

  save
  private

  real(DP), pointer, public :: re_obf( :, : )
  real(DP), pointer, public :: im_obf( :, : )

  real(DP), pointer, public :: re_obf2u( :, :, : )
  real(DP), pointer, public :: im_obf2u( :, :, : )

  real(DP), pointer, public :: re_obf_phs( :, : )
  real(DP), pointer, public :: im_obf_phs( :, : )

#ifdef HAVE_CONTIGUOUS
  CONTIGUOUS :: re_obf, im_obf, re_obf2u, im_obf2u, re_obf_phs, im_obf_phs
#endif


  INTEGER :: my_xpts
  INTEGER :: my_kpts
  INTEGER :: my_num_bands
  INTEGER :: my_start_nx
  INTEGER :: my_obf
  INTEGER, public :: num_obf
  INTEGER :: val_num_bands
  INTEGER :: con_num_bands

  LOGICAL :: is_init = .false.
  LOGICAL :: is_loaded = .false.


  public :: OCEAN_obf_init, OCEAN_obf_load, OCEAN_obf_is_loaded, OCEAN_obf_lrINIT, &
            OCEAN_obf_lrLOAD, OCEAN_obf_make_phase, OCEAN_obf_obf2bloch

  contains

  subroutine OCEAN_obf_lrLOAD( sys, tau, xshift, rbs_out, ibs_out, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    real(DP), intent(out) :: rbs_out(my_num_bands,my_kpts,my_xpts)
    real(DP), intent(out) :: ibs_out(my_num_bands,my_kpts,my_xpts)
    real(DP), intent(inout) :: tau(3)
    integer, intent( out ) :: xshift( 3 )

    real(DP) :: cphs, sphs, pi, phsx, phsy, phsz, su2, su3
!    complex(DP) :: su2
    integer :: iq, iq1, iq2, iq3, iobf
    integer :: ix, iy, iz, xtarg, ytarg, ztarg, xph, yph, zph
    integer :: xiter, ibd, ibd2
    real(DP), allocatable :: re_bloch_state( :, :, : ), im_bloch_state( :,:,:), su(:,:)


    allocate(  re_bloch_state( my_num_bands, my_kpts, my_xpts ), &
               im_bloch_state( my_num_bands, my_kpts, my_xpts ) )
    re_bloch_state(:,:,:) = 0.0_DP
    im_bloch_state(:,:,:) = 0.0_DP


    allocate( su( num_obf, 1 ) ) 
    su(:,:) = 0.0_DP
!    re_obf(:,:) = re_obf(:,:) / sqrt( real( sys%nxpts, DP ) )
!    im_obf(:,:) = im_obf(:,:) / sqrt( real( sys%nxpts, DP ) )

    do xiter = 1, my_xpts
      do iobf = 1, num_obf
        su(iobf,1) = su(iobf,1) + re_obf(iobf,xiter)**2 + im_obf(iobf,xiter)**2
      enddo
    enddo
    if( myid .eq. 0 ) write(6,*) 'Orthog test obf: ', MINVAL(su), MAXVAL(su)
    deallocate(su)

    allocate( su( sys%num_bands, sys%nkpts ) )
    su = 0.0_DP


    if( myid .eq. root ) write(6,*) my_xpts, sys%nkpts, num_obf
    do xiter = 1, my_xpts
      do iq = 1, sys%nkpts
        do iobf = 1, num_obf
          re_bloch_state( :, iq, xiter ) = re_bloch_state(:,iq,xiter) &
                                         + re_obf2u(:,iobf,iq) * re_obf( iobf,xiter ) &
                                         - im_obf2u(:,iobf,iq) * im_obf( iobf,xiter )
          im_bloch_state( :, iq, xiter ) = im_bloch_state(:,iq,xiter) &
                                         + im_obf2u(:,iobf,iq) * re_obf( iobf,xiter ) &
                                         + re_obf2u(:,iobf,iq) * im_obf( iobf,xiter )
          do ibd = 1, sys%num_bands
            su(ibd,iq) = su(ibd,iq) + re_bloch_state(ibd,iq,xiter)**2 + im_bloch_state(ibd,iq,xiter)**2
          enddo
        enddo
      enddo
    enddo
#ifdef MPI
    call MPI_ALLREDUCE( MPI_IN_PLACE, su, sys%nkpts*sys%num_bands, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, comm, ierr )
#endif
    if( myid .eq. root ) then 
      write(6,*) 'bloch created'
      write(6,*) 'Orthog test: ', MINVAL(su), MAXVAL(su)
    endif

    su(:,:) = 1.0_DP / sqrt(su(:,:))
    do xiter = 1, my_xpts
      do iq = 1, sys%nkpts
        do ibd = 1, sys%num_bands
          re_bloch_state(ibd,iq,xiter) = re_bloch_state(ibd,iq,xiter) * su(ibd,iq)
          im_bloch_state(ibd,iq,xiter) = im_bloch_state(ibd,iq,xiter) * su(ibd,iq)
        enddo
      enddo
    enddo


!    do iq = 1, sys%nkpts
!      do ibd = 1, sys%num_bands
!        do ibd2 = 1, ibd - 1
!          su2 = sum( re_bloch_state(ibd,iq,:) * re_bloch_state(ibd2,iq,:) &
!                   + im_bloch_state(ibd,iq,:) * im_bloch_state(ibd2,iq,:) )
!          su3 = sum( - re_bloch_state(ibd,iq,:) * im_bloch_state(ibd2,iq,:) &
!                     + im_bloch_state(ibd,iq,:) * re_bloch_state(ibd2,iq,:) )
!          if( iq .eq. 1 ) then
!            write(6,*) ibd, ibd2, su2, su3
!
!          endif
!          re_bloch_state(ibd,iq,:) = re_bloch_state(ibd,iq,:) - su2 * re_bloch_state(ibd2,iq,:) + su3 * im_bloch_state(ibd2,iq,:)
!          im_bloch_state(ibd,iq,:) = im_bloch_state(ibd,iq,:) - su2 * im_bloch_state(ibd2,iq,:) + su3 * re_bloch_state(ibd2,iq,:)
!        enddo
!        su2 = sum( re_bloch_state(ibd,iq,:)**2+im_bloch_state(ibd,iq,:)**2)
!        re_bloch_state(ibd,iq,:) = re_bloch_state(ibd,iq,:) / sqrt(su2)
!        im_bloch_state(ibd,iq,:) = im_bloch_state(ibd,iq,:) / sqrt(su2)
!      enddo
!    enddo

!    re_bloch_state = 0.0_DP
!    im_bloch_state = 0.0_DP
    




    
    deallocate(su)


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
      xshift( 1 ) = floor( real(sys%xmesh(1), kind( 1.0d0 )) * tau(1) )
    else
      xshift( 1 ) = floor( real(sys%xmesh(2), kind( 1.0d0 )) * (tau(1)-0.5d0 ) )
    endif
    if( mod( sys%kmesh(2), 2 ) .eq. 0 ) then
      xshift( 2 ) = floor( real(sys%xmesh(2), kind( 1.0d0 )) * tau(2) )
    else
      xshift( 2 ) = floor( real(sys%xmesh(2), kind( 1.0d0 )) * (tau(2)-0.5d0 ) )
    endif
    if( mod( sys%kmesh(3), 2 ) .eq. 0 ) then
      xshift( 3 ) = floor( real(sys%xmesh(3), kind( 1.0d0 )) * tau(3) )
    else
      xshift( 3 ) = floor( real(sys%xmesh(3), kind( 1.0d0 )) * (tau(3)-0.5d0 ) )
    endif
    ! 
    if( myid .eq. root ) write(6,*) 'Shifting X-grid by ', xshift(:)
    if( myid .eq. root ) write(6,*) 'Original tau ', tau(:)
    tau( 1 ) = tau(1) - real(xshift(1), kind( 1.0d0 ))/real(sys%xmesh(1), kind( 1.0d0 ))
    tau( 2 ) = tau(2) - real(xshift(2), kind( 1.0d0 ))/real(sys%xmesh(2), kind( 1.0d0 ))
    tau( 3 ) = tau(3) - real(xshift(3), kind( 1.0d0 ))/real(sys%xmesh(3), kind( 1.0d0 ))
    if( myid .eq. root ) write(6,*) 'New tau      ', tau(:)
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
!          if( myid .eq. root ) write(6,*) xiter
          if( xiter .lt. my_start_nx ) cycle
!JTV ???
          if( xiter .ge. my_start_nx + my_xpts ) goto 111
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
                  rbs_out( ibd, iq, xiter - my_start_nx + 1 )  &
                        = cphs * re_bloch_state( ibd, iq, xiter - my_start_nx + 1 )  &
                        - sphs * im_bloch_state( ibd, iq, xiter - my_start_nx + 1 )
                  ibs_out( ibd, iq, xiter - my_start_nx + 1 )  &
                        = cphs * im_bloch_state( ibd, iq, xiter - my_start_nx + 1 )  &
                        + sphs * re_bloch_state( ibd, iq, xiter - my_start_nx + 1 )
                enddo
              enddo
            enddo
          enddo

        enddo
      enddo
    enddo

111 continue


  deallocate(  re_bloch_state, im_bloch_state )


  end subroutine OCEAN_obf_lrLOAD

  subroutine OCEAN_obf_make_phase( sys, tau, xshift, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    integer, intent( out ) :: xshift( 3 )
    real(DP), intent( inout ) :: tau( 3 )

    real(DP) :: cphs, sphs, pi, phsx, phsy, phsz
    integer :: iq, iq1, iq2, iq3
    integer :: ix, iy, iz, xtarg, ytarg, ztarg, xph, yph, zph
    integer :: xiter


    re_obf_phs(:,:) = -100.0_DP
    im_obf_phs(:,:) = -100.0_DP
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
      xshift( 1 ) = floor( real(sys%xmesh(1), kind( 1.0d0 )) * tau(1) )
    else
      xshift( 1 ) = floor( real(sys%xmesh(2), kind( 1.0d0 )) * (tau(1)-0.5d0 ) )
    endif
    if( mod( sys%kmesh(2), 2 ) .eq. 0 ) then
      xshift( 2 ) = floor( real(sys%xmesh(2), kind( 1.0d0 )) * tau(2) )
    else
      xshift( 2 ) = floor( real(sys%xmesh(2), kind( 1.0d0 )) * (tau(2)-0.5d0 ) )
    endif
    if( mod( sys%kmesh(3), 2 ) .eq. 0 ) then
      xshift( 3 ) = floor( real(sys%xmesh(3), kind( 1.0d0 )) * tau(3) )
    else
      xshift( 3 ) = floor( real(sys%xmesh(3), kind( 1.0d0 )) * (tau(3)-0.5d0 ) )
    endif
    ! 
    if( myid .eq. root ) write(6,*) 'Shifting X-grid by ', xshift(:)
    if( myid .eq. root ) write(6,*) 'Original tau ', tau(:)
    tau( 1 ) = tau(1) - real(xshift(1), kind( 1.0d0 ))/real(sys%xmesh(1), kind( 1.0d0 ))
    tau( 2 ) = tau(2) - real(xshift(2), kind( 1.0d0 ))/real(sys%xmesh(2), kind( 1.0d0 ))
    tau( 3 ) = tau(3) - real(xshift(3), kind( 1.0d0 ))/real(sys%xmesh(3), kind( 1.0d0 ))
    if( myid .eq. root ) write(6,*) 'New tau      ', tau(:)
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
!          if( myid .eq. root ) write(6,*) xiter
          if( xiter .lt. my_start_nx ) cycle
!JTV ???
          if( xiter .ge. my_start_nx + my_xpts ) goto 111
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
                re_obf_phs( iq, xiter - my_start_nx + 1 ) = cphs
                im_obf_phs( iq, xiter - my_start_nx + 1 ) = sphs
!                if( iq1 .eq. 2 .and. xiter .eq. 2 ) then
!                  write(6,*) ix, iy, iz
!                  write(6,*) iq1, iq2, iq3
!                  write(6,*) phsx, phsy, phsz
!                  write(6,*) dcos( phsx + phsy + phsz )
!                endif
              enddo
            enddo
          enddo

        enddo
      enddo
    enddo

111 continue

!  write(6,*) 'Phased', sys%xmesh(1),  sys%kmesh( 1 )
!  write(6,*) re_obf_phs( 1, 1 ), im_obf_phs( 1, 1 )
!  write(6,*) re_obf_phs( 2, 2 ), im_obf_phs( 2, 2 )
!  write(6,*) re_obf_phs( 3, 3 ), im_obf_phs( 3, 3 )
  

  end subroutine





  subroutine OCEAN_obf_lrINIT( xpts, kpts, num_bands, start_nx, ierr )
    implicit none
    integer, intent( out ) :: xpts, kpts, num_bands, start_nx
    integer, intent( inout ) :: ierr

    if( .not. is_loaded ) then
      write(6,*) 'ERROR'
      ierr = -1
      return
    endif
    xpts = my_xpts
    kpts = my_kpts
    num_bands = my_num_bands
    start_nx = my_start_nx


  end subroutine OCEAN_obf_lrINIT

  logical function OCEAN_obf_is_loaded()
   OCEAN_obf_is_loaded = is_loaded
    return
  end function OCEAN_obf_is_loaded 

  subroutine OCEAN_obf_init( sys, ierr )
    use OCEAN_mpi, only : myid, nproc, root, comm
    use OCEAN_system
    use iso_c_binding
    use mpi
    implicit none
    include 'fftw3.f03'

    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr


    integer( S_INT ) :: nx_left, nx_start, nx_tmp, i
    type(C_PTR) :: cptr


    if( is_init ) return

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


    if( myid .eq. root ) then
      open(unit=99,file='obf_control',form='formatted',status='old')
      read(99,*) num_obf
      read(99,*) val_num_bands
      read(99,*) con_num_bands
      close(99)
      write(6,*) num_obf, val_num_bands, con_num_bands
    endif
#ifdef MPI
    call MPI_BCAST( num_obf, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) then 
      write(6,*) 'MPI_BCAST FAILED', ierr
      return
    endif
    call MPI_BCAST( val_num_bands, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( con_num_bands, 1, MPI_INTEGER, root, comm, ierr )
#endif
    my_obf = num_obf


    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_obf, C_SIZE_T) )
    call c_f_pointer( cptr, re_obf2u, [my_num_bands, my_obf, my_kpts] )
    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_obf, C_SIZE_T) )
    call c_f_pointer( cptr, im_obf2u, [my_num_bands, my_obf, my_kpts] )
  
    cptr = fftw_alloc_real( int(my_obf * my_xpts, C_SIZE_T) )
    call c_f_pointer( cptr, re_obf, [my_obf, my_xpts] )
    cptr = fftw_alloc_real( int(my_obf * my_xpts, C_SIZE_T) )
    call c_f_pointer( cptr, im_obf, [my_obf, my_xpts] )

    cptr = fftw_alloc_real( int(my_xpts*sys%nkpts, C_SIZE_T) )
    call c_f_pointer( cptr, re_obf_phs, [sys%nkpts,my_xpts] )
    cptr = fftw_alloc_real( int(my_xpts*sys%nkpts, C_SIZE_T) )
    call c_f_pointer( cptr, im_obf_phs, [sys%nkpts,my_xpts] )

    is_init = .true.

  end subroutine OCEAN_obf_init

  subroutine OCEAN_obf_load( sys, ierr )
    use OCEAN_mpi, only : myid, nproc, root, comm
    use OCEAN_system
    use mpi
    implicit none

    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr


!    integer, parameter :: eig = 35

      !
!    integer :: nx, ny, nz
!    real( kind = kind( 1.0d0 ) ), dimension( nx, ny, nz, nbd ) :: ur, ui
    !
    integer :: iq, ibd, ix, iy, iz, ivl, ivh, icl, ich, i
    integer :: dumint, ivh2, iobf, fheig, finfo, obf_width
    integer(kind=MPI_OFFSET_KIND) :: offset
!    integer :: xtarg, ytarg, ztarg, xph, yph, zph
!    real( kind = kind( 1.0d0 ) ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
!    real( kind = kind( 1.0d0 ) ) :: su, sul, suh
!    real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: tmp_ur, tmp_ui, ur, ui
!    real( DP ), allocatable :: re_transpose( :, : ), im_transpose( :, : )
    complex( DP ), allocatable, dimension( :, :, : ) :: temp_obf
    real( DP ), allocatable, dimension( : ) ::  re_temp_obf, im_temp_obf
    complex( DP ), allocatable, dimension( :, : ) :: tmp_eigvec
    logical :: metal
    integer :: nx_left, nx_start, nx_tmp, xiter

    if( is_loaded ) return

    if( is_init .eqv. .false. ) then
      call OCEAN_obf_init( sys, ierr )
      if( ierr .ne. 0 ) return
    endif

    !!! (1)
    ! Load up the obfs
    ! Currently writen as a complex
    if( myid .eq. root ) then
      open( unit=99,file='u1.dat',form='unformatted',status='old')
      rewind(99)
      allocate( temp_obf( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1) ), &
                re_temp_obf( sys%nxpts ), im_temp_obf( sys%nxpts ) )
    endif

    do iobf = 1, num_obf  
      if( myid .eq. root ) then 
        read( 99 ) temp_obf(:,:,:)
!        re_trans_temp_obf = transpose( real( temp_obf ) )
!        im_trans_temp_obf = transpose( aimag( temp_obf ) ) 

        xiter = 0
        do iz = 1, sys%xmesh(3)
          do iy = 1, sys%xmesh(2)
            do ix = 1, sys%xmesh(1)
              xiter = xiter + 1
              re_temp_obf( xiter ) = real( temp_obf( iz, iy, ix ) )
              im_temp_obf( xiter ) = aimag( temp_obf( iz, iy, ix ) )
            enddo
          enddo
        enddo
      endif
      
      nx_left = sys%nxpts
      nx_start = 1
      do i = 0, nproc - 1 
        nx_tmp = nx_left / ( nproc - i )
        nx_left = nx_left - nx_tmp
        if( myid .eq. root .and. iobf .eq. 1 ) write(6,*) i, nx_start, nx_tmp
        if( i .eq. root .and. myid .eq. root ) then
          re_obf( iobf, : ) = re_temp_obf( nx_start : nx_start + nx_tmp - 1 )
          im_obf( iobf, : ) = im_temp_obf( nx_start : nx_start + nx_tmp - 1 )
#ifdef MPI
        elseif( myid .eq. root ) then
          call MPI_SEND( re_temp_obf(nx_start), nx_tmp, MPI_DOUBLE_PRECISION, i, i, comm, ierr )
          call MPI_SEND( im_temp_obf(nx_start), nx_tmp, MPI_DOUBLE_PRECISION, i, i+nproc, comm, ierr )
        elseif( myid .eq. i ) then
          if( iobf .eq. 1 ) then
            allocate( re_temp_obf( nx_tmp + 1), im_temp_obf( nx_tmp + 1) )
          endif

!JTV probably want to restructure to call all the recieves earlier on
          call MPI_RECV( re_temp_obf, nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                        i, comm, MPI_STATUS_IGNORE, ierr )
          call MPI_RECV( im_temp_obf, nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                        i+nproc, comm, MPI_STATUS_IGNORE, ierr )
          re_obf( iobf, : ) = re_temp_obf( 1:nx_tmp )
          im_obf( iobf, : ) = im_temp_obf( 1:nx_tmp )
#endif
        endif
        nx_start = nx_start + nx_tmp
      enddo
#ifdef MPI
      call MPI_BARRIER( comm, ierr )
#endif
    enddo

!    re_obf(:,:) = re_obf(:,:) / sqrt( real( sys%nxpts, DP ) )
!    im_obf(:,:) = im_obf(:,:) / sqrt( real( sys%nxpts, DP ) )

    deallocate( re_temp_obf, im_temp_obf )

    if( myid .eq. 0 ) deallocate( temp_obf )




!JTV here we should decided valence or conduction; hardwired for absorption right now
    obf_width = con_num_bands
#ifdef MPI    
    call MPI_INFO_CREATE( finfo, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to create info'
      return
    endif
    call MPI_FILE_OPEN( comm, 'con_eigvecs.dat', MPI_MODE_RDONLY, finfo, fheig, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to open file con_eigvecs.dat'
      return
    endif
    offset = 0
    call MPI_FILE_SET_VIEW( fheig, offset, MPI_DOUBLE_COMPLEX, &
                          MPI_DOUBLE_COMPLEX, 'native', finfo, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to set file view for con_eigvecs.dat'
      return
    endif
#else
    if( myid .eq. 0 ) write(6,*) 'OBF only works when compiled with MPI support!!'
    ierr = -1
    return
#endif

    
      

    if( myid .eq. 0 ) then

      !if absorption need con, emission need val
      

      open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) ivl, ivh, icl, ich
      close( unit=99 )
      ivh2 = ivh
      open( unit=99, file='metal', form='formatted', status='old')
      read( 99, * ) metal
      close( 99 )
      if( metal ) then
        open( unit=36, file='ibeg.h', form='formatted', status='old' )
      endif

      if ( sys%num_bands .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'
!      open( unit=eig, file='con_eigvecs.dat', form='unformatted', status='unknown' )
!      rewind eig
      write(6,*) 'nspn: ', sys%nspn
      if( sys%nspn .ne. 1 ) then
        write(6,*) 'No spin yet'
        ierr = 1  
        return
      endif

    endif

#ifdef MPI
    call MPI_BCAST( metal, 1, MPI_LOGICAL, root, comm, ierr )
    call MPI_BCAST( ivl, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( ivh, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( icl, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( ich, 1, MPI_INTEGER, root, comm, ierr )
#endif
    ivh2 = ivh
  
    allocate( tmp_eigvec( num_obf, sys%num_bands ) )

    do iq = 1, sys%nkpts
      if( myid .eq. 0 ) then 
        open( unit=99, file='gumatprog', form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(2i8)' ) iq, sys%nkpts
        close( unit=99 )
        if( metal ) then
          read( 36, * ) dumint, ivh2
          ivh2 = ivh2 - 1 
        endif
      endif

      if( metal ) call MPI_BCAST( ivh2, 1, MPI_INTEGER, root, comm, ierr )


#ifdef MPI
    !  Skip all of the occupied bands (and for metals)
    ! obf_width is the number of bands in total stored per k
    ! need to offset by an additional num_obf * (valence bands )
      offset = (iq-1)*num_obf*obf_width + num_obf*ivh2
      call MPI_FILE_READ_AT_ALL( fheig, offset, tmp_eigvec, num_obf*sys%num_bands, &
                                 MPI_DOUBLE_COMPLEX, &
                                 MPI_STATUS_IGNORE, ierr )
#endif
      do ibd = 1, sys%num_bands
        re_obf2u( ibd, :, iq ) = real( tmp_eigvec( :, ibd ) )
        im_obf2u( ibd, :, iq ) = aimag( tmp_eigvec( :, ibd ) )
      end do
      
    enddo

    if( myid .eq. root ) then
      close(36)
    endif

    deallocate( tmp_eigvec )
    call MPI_FILE_CLOSE( fheig, ierr )
    

    is_loaded = .true.

    if( myid .eq. root ) write(6,*) 'OBFs loaded'


  end subroutine OCEAN_obf_load


  subroutine OCEAN_obf_obf2bloch( sys, rbs, ibs, nbands, nkpts, nxpts, ierr )
    use OCEAN_mpi, only : myid, nproc, root, comm
    use OCEAN_system
    use mpi
    implicit none
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    integer, intent( in ) :: nbands, nkpts, nxpts
    real(DP), intent( out ), dimension(nbands, nkpts, nxpts) :: rbs, ibs
    !
    ! 
    

    


  end subroutine OCEAN_obf_obf2bloch

end module OCEAN_obf
