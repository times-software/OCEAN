module OCEAN_bloch
  use AI_kinds

  save
  private


  ! these guys have no position phasing
  real(DP), pointer, contiguous :: re_bloch_state( :, :, : )
  real(DP), pointer, contiguous :: im_bloch_state( :, :, : )


  
  INTEGER :: my_xpts
  INTEGER :: my_kpts
  INTEGER :: my_num_bands
  INTEGER :: my_start_nx
  
  LOGICAL :: is_init = .false.
  LOGICAL :: is_loaded = .false.


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
    real(DP), intent(out) :: rbs_out(my_num_bands,my_kpts,my_xpts) 
    real(DP), intent(out) :: ibs_out(my_num_bands,my_kpts,my_xpts) 
    real(DP), intent(inout) :: tau(3)
    integer, intent( out ) :: xshift( 3 )

    real(DP) :: cphs, sphs, pi, phsx, phsy, phsz
    integer :: iq, iq1, iq2, iq3
    integer :: ix, iy, iz, xtarg, ytarg, ztarg, xph, yph, zph
    integer :: xiter, ibd


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


!    nx_left = sys%nxpts
!    nx_start = 1
!    do i = 0, nproc - 1
!      nx_tmp = nx_left / ( nproc - i )
!      nx_left = nx_left - nx_tmp
!      if( myid .eq. i ) my_start = nx_start
!      nx_start = nx_start + nx_tmp
!    enddo

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
    use OCEAN_mpi, only : myid, nproc, root, comm
    use OCEAN_system
    use iso_c_binding
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
    

    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts, C_SIZE_T) )
    call c_f_pointer( cptr, re_bloch_state, [my_num_bands, my_kpts, my_xpts] )
    cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts, C_SIZE_T) )
    call c_f_pointer( cptr, im_bloch_state, [my_num_bands, my_kpts, my_xpts] )

    is_init = .true.


  end subroutine OCEAN_bloch_init

  subroutine OCEAN_bloch_load( sys, ierr )
    use OCEAN_mpi, only : myid, nproc, root, comm
    use OCEAN_system
    use mpi
    implicit none

    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr


    integer, parameter :: u2dat = 35
    !
    integer :: nx, ny, nz, nbd, nq
!    real( kind = kind( 1.0d0 ) ), dimension( nx, ny, nz, nbd ) :: ur, ui
    !
    integer :: iq, ibd, ig, idum( 3 ), ix, iy, iz, ivl, ivh, icl, ich, ispn, i
    integer :: iq1, iq2, iq3, dumint, icl2, ich2, ivh2, xshift(3)
    integer :: xtarg, ytarg, ztarg, xph, yph, zph
    real( kind = kind( 1.0d0 ) ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
    real( kind = kind( 1.0d0 ) ) :: su, sul, suh
    real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: tmp_ur, tmp_ui, ur, ui
    real( DP ), allocatable :: re_transpose( :, : ), im_transpose( :, : )
!    complex( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: tmp_bloch
    logical :: metal, normal
    integer :: nx_left, nx_start, nx_tmp, xiter, ii

    if( is_loaded ) return

    if( is_init .eqv. .false. ) then 
      call OCEAN_bloch_init( sys, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( ( .not. associated( re_bloch_state ) ) .or. ( .not. associated( im_bloch_state ) ) ) then
      ierr = -1
      if( myid .eq. root ) write( 6, * ) 'Trying to read in band states before allocation'
      return
    endif

    nx = sys%xmesh(1)
    ny = sys%xmesh(2)
    nz = sys%xmesh(3)
    nbd = sys%num_bands 
    !  tmp_bloch( nbd, nx*ny*nz )
    allocate( ur(sys%xmesh(1),sys%xmesh(2),sys%xmesh(3),sys%num_bands), &
              ui(sys%xmesh(1),sys%xmesh(2),sys%xmesh(3),sys%num_bands) )
    ! As per usual, do this the dumbest way first, 
    ! 1) generic copy of original serial method 
    ! 2) MPI

      

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
      if( metal ) then
        open( unit=36, file='ibeg.h', form='formatted', status='old' )
      endif

      if ( nbd .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'
      open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
      rewind u2dat 
      write(6,*) 'nspn: ', sys%nspn
      if( sys%nspn .ne. 1 ) then
        write(6,*) 'No spin yet'
        ierr = 1
        goto 111 
      endif

      allocate( tmp_ur( nx, ny, nz, nbd ), tmp_ui( nz, ny, nz, nbd ), &
                re_transpose( sys%num_bands, sys%nxpts ), im_transpose( sys%num_bands, sys%nxpts ) )
    

    endif

      
    if( myid .eq. 0 ) write(6,*) size(re_bloch_state,1), size(re_bloch_state,2),size(re_bloch_state,3)

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
          if( metal ) then
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
!                  phsx = 2.0d0 * pi * dble( ( ix - 1 ) * ( iq1 - 1 ) ) / dble( nx * sys%kmesh( 1 ) )
!                  phsy = 2.0d0 * pi * dble( ( iy - 1 ) * ( iq2 - 1 ) ) / dble( ny * sys%kmesh( 2 ) )
!                  phsz = 2.0d0 * pi * dble( ( iz - 1 ) * ( iq3 - 1 ) ) / dble( nz * sys%kmesh( 3 ) )
!                  cphs = dcos( phsx + phsy + phsz )
!                  sphs = dsin( phsx + phsy + phsz )
                  do ibd = 1, nbd
!                     psir = cphs * ur( ix, iy, iz, ibd ) - sphs * ui( ix, iy, iz, ibd )
!                     psii = cphs * ui( ix, iy, iz, ibd ) + sphs * ur( ix, iy, iz, ibd )
!!                     tmp_ur( xtarg, ytarg, ztarg, ibd ) = psir
!!                     tmp_ui( xtarg, ytarg, ztarg, ibd ) = psii
!                    re_transpose( ibd, xiter ) = psir
!                    im_transpose( ibd, xiter ) = psii
                    re_transpose( ibd, xiter ) = ur( ix, iy, iz, ibd )
                    im_transpose( ibd, xiter ) = ui( ix, iy, iz, ibd )
                  end do
               end do
            end do
         end do
!!         ur( :, :, :, : ) = tmp_ur( :, :, :, : )
!!         ui( :, :, :, : ) = tmp_ui( :, :, :, : )
!         xiter = 0
!         do iz = 1, nz
!           do iy = 1, ny
!             do ix = 1, nx
!               xiter = xiter + 1
!!               tmp_bloch( :, xiter ) = cmplx( ur( ix, iy, iz, : ), ui( ix, iy, iz, : ) )
!!                write(6,*) ix, iy, iz, iq, xiter
!!                re_bloch_state( :, iq, xiter ) = ur( ix, iy, iz, : )
!!                im_bloch_state( :, iq, xiter ) = ui( ix, iy, iz, : )
!                re_transpose( :, xiter ) = tmp_ur( ix, iy, iz, : )
!                im_transpose( :, xiter ) = tmp_ui( ix, iy, iz, : )
!             enddo
!           enddo
!         enddo
       endif
       if( myid .eq. root .and. mod(iq,10) .eq. 0 ) write(6,*) iq

       nx_left = sys%nxpts
       nx_start = 1
       do i = 0, nproc - 1
          nx_tmp = nx_left / ( nproc - i )
          nx_left = nx_left - nx_tmp
          if( myid .eq. root .and. iq .eq. 1 ) write(6,*) i, nx_start, nx_tmp
          if( i .eq. root .and. myid .eq. root ) then
            re_bloch_state( :, iq, : ) = re_transpose( :, nx_start : nx_start + nx_tmp - 1 )
            im_bloch_state( :, iq, : ) = im_transpose( :, nx_start : nx_start + nx_tmp - 1 )
#ifdef MPI
          elseif( myid .eq. root ) then
            call MPI_SEND( re_transpose(1,nx_start), my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i, comm, ierr )
            call MPI_SEND( im_transpose(1,nx_start), my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i+nproc, comm, ierr )
          elseif( myid .eq. i ) then
            if( iq .eq. 1 ) then
              allocate( re_transpose( my_num_bands, nx_tmp ), im_transpose( my_num_bands, nx_tmp ) )
            endif

            call MPI_RECV( re_transpose, my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                          i, comm, MPI_STATUS_IGNORE, ierr )
            call MPI_RECV( im_transpose, my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                          i+nproc, comm, MPI_STATUS_IGNORE, ierr )
            re_bloch_state( :, iq, : ) = re_transpose( :, : )
            im_bloch_state( :, iq, : ) = im_transpose( :, : )
#endif
          endif
         nx_start = nx_start + nx_tmp
       enddo



      enddo
     enddo
    enddo

111 continue


    if( myid .eq. 0 ) write ( 6, '(1a16,2f20.15)' ) 'norm bounds ... ', sul, suh

    if( myid .eq. 0 ) close(u2dat)

    deallocate( re_transpose, im_transpose, ur, ui )
    if( myid .eq. 0 ) deallocate( tmp_ur, tmp_ui )

    is_loaded = .true.

  end subroutine OCEAN_bloch_load


end module OCEAN_bloch
