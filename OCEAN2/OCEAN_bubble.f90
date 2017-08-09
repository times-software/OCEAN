! Copyright (C) 2016 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! This is the RPA bubble diagrams for the valence level

module OCEAN_bubble

  use AI_kinds
  use iso_c_binding
  use FFT_wrapper, only : fft_obj
  implicit none
  private
  save

  logical :: MATCH_AI2NBSE = .true.

  complex( C_DOUBLE_COMPLEX ), allocatable :: scratch( : )
  real( DP ), allocatable  :: bubble( : )
  real(DP), allocatable :: re_scratch(:,:), im_scratch(:,:)

!  integer, allocatable :: xstart_by_mpiID(:)

!  type(C_PTR)        :: fplan
!  type(C_PTR)        :: bplan

  type( fft_obj ) :: fo

  logical :: is_init = .false.

  public :: AI_bubble_prep, AI_bubble_act, AI_bubble_clean

  contains

  real(dp) function gvec_length( gvec, bvec) result( length )
    implicit none
    integer :: gvec( 3 )
    real(dp) :: bvec(3,3)
    !
    !write(71,*) gvec(:) 
    length = ( bvec(1,1) * dble(gvec(1)) + bvec(2,1) * dble(gvec(2)) + bvec(3,1) * dble(gvec(3)) ) ** 2.0_dp &
           + ( bvec(1,2) * dble(gvec(1)) + bvec(2,2) * dble(gvec(2)) + bvec(3,2) * dble(gvec(3)) ) ** 2.0_dp &
           + ( bvec(1,3) * dble(gvec(1)) + bvec(2,3) * dble(gvec(2)) + bvec(3,3) * dble(gvec(3)) ) ** 2.0_dp
  end function gvec_length

  subroutine AI_bubble_clean( )
    use OCEAN_mpi, only : myid, root
!    use iso_c_binding
    use FFT_wrapper, only : FFT_wrapper_delete
    !
    implicit none
!    include 'fftw3.f03'
    !
!    integer, intent( inout ) :: ierr
    !
    if( myid .eq. root ) then
      call FFT_wrapper_delete( fo )
!      call dfftw_destroy_plan( bplan )
!      call dfftw_destroy_plan( fplan )
    endif

  end subroutine AI_bubble_clean

  subroutine AI_bubble_prep( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi
    use OCEAN_val_states, only : nxpts_by_mpiID
    use FFT_wrapper, only : FFT_wrapper_init
!    use iso_c_binding
    !
    implicit none
!    include 'fftw3.f03'
    !
    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, nthreads, xmesh( 3 )
!    type(fftw_iodim) :: guru_dims( 3 )
!$  integer, external :: omp_get_max_threads
    real( DP ) :: length
    !
    nthreads = 1
!$  nthreads = omp_get_max_threads()

    if( is_init ) then
      if( myid .eq. root ) write(6,*) '  bubble cached'
      return
    endif

    ! The data for bubble gets saved down to a single array for FFT
    !  In the future we will want two! different sorts for spin up and spin down
    if( myid .eq. root ) then
      allocate( bubble(  sys%nxpts ), &
                scratch( sys%nxpts ), re_scratch( sys%nxpts, sys%valence_ham_spin ), &
                im_scratch( sys%nxpts, sys%valence_ham_spin ), &
                STAT=ierr )
      if( ierr .ne. 0 ) then
        write( 1000+myid,* ) 'Failed to allocate bubble and scratch'
        goto 111
      endif
!      call dfftw_plan_with_nthreads( nthreads )
!!      call dfftw_plan_dft_3d( fplan, sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), &
!!                              scratch, scratch, FFTW_FORWARD, FFTW_PATIENT )
!!      call dfftw_plan_dft_3d( bplan, sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), &
!!                              scratch, scratch, FFTW_BACKWARD, FFTW_PATIENT )

      

      ! For valence states the ordering of bands the array is ( z, y, x )
      ! Therefore we need to invert before passing it in
      xmesh( 1 ) = sys%xmesh( 3 )
      xmesh( 2 ) = sys%xmesh( 2 )
      xmesh( 3 ) = sys%xmesh( 1 )
      call FFT_wrapper_init( xmesh, fo, scratch )
! Also works (the same as above)
!      fplan = fftw_plan_dft_3d( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), & 
!                                scratch, scratch, FFTW_FORWARD, FFTW_PATIENT )
!      bplan = fftw_plan_dft_3d( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), & 
!                                scratch, scratch, FFTW_BACKWARD, FFTW_PATIENT )


!      fplan = fftw_plan_guru_split_dft( 3,  


!      call dfftw_plan_with_nthreads( 1 )
!      call AI_bubble_wall_search( sys, length, ierr )
      call wall_search_old_way( sys, length, ierr )
      if( ierr .ne. 0 ) goto 111
      call AI_bubble_populate( sys, length, ierr )
      if( ierr .ne. 0 ) goto 111

!      allocate( xstart_by_mpiID( 0:nproc-1 ), STAT=ierr )
!      if( ierr .ne. 0 ) goto 111
!      j = 1
!      do i = 0, nproc - 1
!        xstart_by_mpiID( i ) = j
!        j = j + nxpts_by_mpiID( i )
!      enddo

    else ! To allow aggressive error checking
      allocate( bubble( 1 ), &
                scratch( 1 ), &
                STAT=ierr )
      if( ierr .ne. 0 ) then
        write( 1000+myid ) 'Failed to allocate bubble and scratch'
        goto 111
      endif
    endif
    !

    is_init = .true.
111 continue

  end subroutine AI_bubble_prep
!
!
!
  subroutine AI_bubble_wall_search( sys, length, ierr )
    use OCEAN_system
    !
    implicit none
    !
    type( o_system ), intent( in ) :: sys
    real( DP ), intent( out ) :: length
    integer( S_INT ), intent( inout ) :: ierr
    !
    integer :: ix, iy, iz, gvec( 3 )
    real( DP ) :: min_edge_length, prev_min_edge
    !
!    real( DP ), external :: gvec_length
    !
    ! exclude anything greater than or equal to the min distance we missed.
    !  ie. if x goes from -5 to 5 than the |x| = 6 surface will tell us what we missed
    !  this may be slightly more accepting than taking everything equal or less than the |x|=5 surface
    !  since it may include some corners that would otherwise be thrown out.
    gvec = (/ floor( dble(sys%xmesh( 1 )) / 2.0_dp ) + 1, floor(dble(sys%xmesh( 2 )) / 2.0_dp )+ 1, &
              floor(dble(sys%xmesh( 3 )) / 2.0_dp ) + 1/)
    min_edge_length = gvec_length( gvec, sys%bvec )
    gvec = (/ floor( dble(sys%xmesh( 1 )) / 2.0_dp ) + 2, floor(dble(sys%xmesh( 2 )) / 2.0_dp )+ 1, &
              floor(dble(sys%xmesh( 3 )) / 2.0_dp ) + 1/)

    prev_min_edge = gvec_length( gvec, sys%bvec )
    do ix = ceiling(dble(-sys%xmesh(1))/2.0_dp) , floor(dble(sys%xmesh(1))/2.0_dp  )
      do iy = ceiling(dble(-sys%xmesh(2)) / 2.0_dp) , floor(dble(sys%xmesh( 2 )) / 2.0_dp )
        gvec = (/ ix, iy, ceiling(dble(-sys%xmesh(3))/2.0_dp)  /)
        length = gvec_length( gvec, sys%bvec )
        if ( length .lt. min_edge_length ) then 
          prev_min_edge = min_edge_length
          min_edge_length = length
        endif
        gvec = (/ ix, iy, floor(dble(sys%xmesh(3))/2.0_dp) /)
        length = gvec_length( gvec, sys%bvec )
        if ( length .lt. min_edge_length ) then 
          prev_min_edge = min_edge_length
          min_edge_length = length
        endif
      enddo
      do iz = ceiling(dble(-sys%xmesh( 3 )) / 2.0_dp) , floor(dble(sys%xmesh( 3 )) / 2.0_dp  )
        gvec = (/ ix, ceiling(dble(-sys%xmesh(2))/2.0_dp) , iz /)
        length = gvec_length( gvec, sys%bvec )
        if ( length .lt. min_edge_length ) then 
          prev_min_edge = min_edge_length
          min_edge_length = length
        endif
        gvec = (/ ix, floor(dble(sys%xmesh(2))/2.0_dp) , iz /)
        length = gvec_length( gvec, sys%bvec )
        if ( length .lt. min_edge_length ) then 
          prev_min_edge = min_edge_length
          min_edge_length = length
        endif
      enddo
    enddo
    do iy = ceiling(dble(-sys%xmesh(2)) / 2.0_dp) , floor(dble(sys%xmesh( 2 )) / 2.0_dp )
      do iz = ceiling(dble(-sys%xmesh( 3 )) / 2.0_dp) , floor(dble(sys%xmesh( 3 )) / 2.0_dp )
        gvec = (/ ceiling(dble(-sys%xmesh(1))/2.0_dp) , iy, iz /)
        length = gvec_length( gvec, sys%bvec )
        if ( length .lt. min_edge_length ) then 
          prev_min_edge = min_edge_length
          min_edge_length = length
        endif
        gvec = (/ floor(dble(sys%xmesh(2))/2.0_dp) , iy, iz /)
        length = gvec_length( gvec, sys%bvec )
        if ( length .lt. min_edge_length ) then
          prev_min_edge = min_edge_length
          min_edge_length = length
        endif
      enddo
    enddo
    length = min_edge_length
    length = length + 1.0d-8
    write( 6, * ) 'Max gvec length: ', length, prev_min_edge
  end subroutine AI_bubble_wall_search
!
!
  subroutine wall_search_old_way( sys, length, ierr )
    use OCEAN_system
    !
    implicit none
    !
    type( o_system ), intent( in ) :: sys
    real( DP ), intent( out ) :: length
    integer( S_INT ), intent( inout ) :: ierr
    !
    integer :: ix, iy, iz, gvec( 3 )
    real( DP ) :: temp_length, next_length
    !

    length = huge( 0.0_dp )

    do ix = floor(dble(-sys%xmesh(1))/2.0_dp) , ceiling(dble(sys%xmesh(1))/2.0_dp  )
      do iy = floor(dble(-sys%xmesh(2)) / 2.0_dp) , ceiling(dble(sys%xmesh( 2 )) / 2.0_dp )
        do iz = floor(dble(-sys%xmesh(3)) / 2.0_dp) , ceiling(dble(sys%xmesh( 3 )) / 2.0_dp )

          if( 2 * abs( iz ) .ge. sys%xmesh(3) .or. &
              2 * abs( iy ) .ge. sys%xmesh(2) .or. &
              2 * abs( ix ) .ge. sys%xmesh(1) ) then

            gvec = (/ ix, iy, iz /)
            temp_length = gvec_length( gvec, sys%bvec )

            if( temp_length .lt. length ) then
!              write(6,*)  ix, iy, iz , temp_length, length
              length = temp_length
            endif
          endif
        enddo
      enddo
    enddo
!    write(6,*)  ''

    next_length = 0.0_dp
    do ix = floor(dble(-sys%xmesh(1))/2.0_dp) , ceiling(dble(sys%xmesh(1))/2.0_dp  )
      do iy = floor(dble(-sys%xmesh(2)) / 2.0_dp) , ceiling(dble(sys%xmesh( 2 )) / 2.0_dp )
        do iz = floor(dble(-sys%xmesh(3)) / 2.0_dp) , ceiling(dble(sys%xmesh( 3 )) / 2.0_dp )

            gvec = (/ ix, iy, iz /)
            temp_length = gvec_length( gvec, sys%bvec )

            if( temp_length .lt. length .and. temp_length .gt. next_length) then
!              write(6,*)  ix, iy, iz , temp_length, next_length
              next_length = temp_length
            endif
        enddo
      enddo
    enddo



    write( 6, * ) 'Max gvec length: ', length, next_length, 0.5_dp * ( length + next_length )
    length = 0.5_dp * ( length + next_length )

    

  end subroutine wall_search_old_way

! 
  subroutine AI_bubble_populate( sys, length, ierr )
    use ai_kinds
    use OCEAN_system
    use OCEAN_constants, only : pi_dp, Hartree2eV
    !
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    real( DP ), intent( in ) :: length
    integer( S_INT ), intent( inout ) :: ierr
    !
    integer :: ix, iy, iz, ij, ik, igvec, izz, iyy, ixx, temp_gvec( 3 ), iter
    real( DP ) :: mul, gsqd, qq( 3 ), bmet( 3, 3 ) !, pi
  
    real( DP ), allocatable :: TdBubble( :, :, : )
!    real( DP ), external :: gvec_length
    ! 
!    pi = 4.0d0 * atan( 1.0d0 )
    allocate( TdBubble( sys%xmesh( 3 ), sys%xmesh( 2 ), sys%xmesh( 1 ) ) )
    TdBubble(:,:,:) = 0.0_DP

    igvec = 0
!    pi = 4.0d0 * atan( 1.0d0 )
    do ij = 1, 3
      do ik = 1, 3
        bmet(ij,ik) = dot_product( sys%bvec(:,ij), sys%bvec( :, ik ) )
!        write(6,*) bmet(ij,ik), sys%bvec(ij, ik)
      enddo
!      write(6,*) sys%bvec( :, ij ) 
    enddo
    !
    iter = 0
    bubble( : ) = 0.0_dp
    write( 6, * ) 'qinb: ', sys%qinunitsofbvectors( : )
!    write(6,* ) sys%nkpts , sys%celvol
!    write(6,* )  floor( dble(-sys%xmesh(1))/2.0_dp ), ceiling( dble(sys%xmesh(1))/2.0_dp )

    !
    do ix = ceiling( dble(-sys%xmesh(1))/2.0_dp ), floor( dble(sys%xmesh(1))/2.0_dp )
      temp_gvec( 1 ) = ix 
      do iy =  ceiling( dble(-sys%xmesh(2))/2.0_dp ), floor( dble(sys%xmesh(2))/2.0_dp )
        temp_gvec( 2 ) = iy 
        do iz = ceiling( dble(-sys%xmesh(3))/2.0_dp ), floor( dble(sys%xmesh(3))/2.0_dp )
          temp_gvec( 3 ) = iz 
          qq( : ) = sys%qinunitsofbvectors( : ) + real( temp_gvec( : ), DP )
          gsqd = 0.0_dp
          do ij = 1, 3
            do ik = 1, 3
              gsqd = gsqd + qq( ij ) * qq( ik ) * bmet( ij, ik )
            enddo
          enddo
          ixx = 1 + ix
          iyy = 1 + iy
          izz = 1 + iz
          if ( ixx .le. 0 ) ixx = ixx + sys%xmesh( 1 )
          if ( iyy .le. 0 ) iyy = iyy + sys%xmesh( 2 )
          if ( izz .le. 0 ) izz = izz + sys%xmesh( 3 )
          iter = izz + sys%xmesh(3)*(iyy-1) + sys%xmesh(3)*sys%xmesh(2)*(ixx-1) 
          if( gvec_length( temp_gvec, sys%bvec ) .ge. length ) then
            mul = 0.0_dp
!            write(6,*) ix, iy, iz, gvec_length( temp_gvec, sys%bvec ), .false.
          else
            mul = 4.0_dp * pi_dp / ( sys%celvol * gsqd * dble(sys%nkpts) )
            ! fake to better match old
!            mul = mul * 27.2114d0 / Hartree2eV
            igvec = igvec + 1
!            write(6,*) ix, iy, iz, gvec_length( temp_gvec, sys%bvec ), .true.
          endif
          Tdbubble( izz, iyy, ixx ) = mul
!          write(103,'(6I5,X,E22.7,E22.7)') ix, iy,iz, ixx, iyy, izz, mul * 27.2114d0, gsqd
          !  bubble( iter ) = mul
        enddo
      enddo 
    enddo 
    igvec = igvec - 1
!    bubble( 1 ) = 0.d0 ! skip G = 0
    Tdbubble( 1, 1, 1 ) = 0.0_dp
    bubble = reshape( Tdbubble, (/ sys%nxpts /) )
    write( 6, * )'Num gvecs retained: ', igvec, sys%nxpts
    write( 6, * ) maxval( bubble ) * Hartree2eV !27.2114d0
    deallocate( Tdbubble )

!    write(6,*) 'BUBBLE:'
!    do ix = 1, sys%nxpts
!      write(6,*) ix, bubble( ix ) * 27.2114d0
!    enddo
!    write(6,*) 'END BUBBLE:'

  end subroutine AI_bubble_populate
!
!
!
!> @author John Vinson, NIST
!
!> @brief Calculates the action of the exchange (bubble) operator on the valence vector
!
!> @details The exchange operator requires that both the electron and hole have 
!! the same spin, but can mix up-up with down-down. We therefore, when spin is 
!! being used, loop over the first section of the calculation and accumulate 
!! the result. Then we loop over the second half to apply it to both up-up and 
!! down-down.
  subroutine AI_bubble_act( sys, psi, psiout, ierr )
    use OCEAN_val_states, only : nkpts, nxpts, nbc, nbv, nxpts_pad, &
                                 re_val, im_val, re_con, im_con, &
                                 cache_double, startx_by_mpiID, nxpts_by_mpiID
    use OCEAN_psi
    use OCEAN_mpi!, only : nproc, myid, root, comm
    use OCEAN_system
    use iso_c_binding
    use FFT_wrapper, only : OCEAN_BACKWARD, OCEAN_FORWARD, FFT_wrapper_single
    implicit none
!    include 'fftw3.f03'
    !
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( in ) :: psi
    type( OCEAN_vector ), intent( inout ) :: psiout
    integer, intent( inout ) :: ierr
    !
    integer :: bciter, bviter, ik, ispn, dft_spin, psi_spin, jspn, dft_spin2, psi_spin2
    integer :: ix, iix, xstop, xwidth, nthreads, iproc
    integer :: psi_con_pad!, nxpts_pad
!    integer :: nbv, nbc, nkpts
    integer, external :: omp_get_max_threads
    real( DP ), allocatable :: re_l_bubble( : ), im_l_bubble( : )
    real( DP ), allocatable :: re_amat( :, :, : ), im_amat( :, :, : )
    !
    real( DP ), parameter :: one = 1.0_dp
    real( DP ), parameter :: zero = 0.0_dp
    real( DP ), parameter :: minusone = -1.0_dp
    real( DP ) :: spin_prefac, minus_spin_prefac
    !
    !
    ! Populate sizes from psi and val_states
    call OCEAN_psi_returnBandPad( psi_con_pad, ierr )
    if( ierr .ne. 0 ) return
!    call OCEAN_val_states_returnPadXpts( nxpts_pad, ierr )
!    if( ierr .ne. 0 ) return

!    do bviter = 1, nbv
!      do bciter = 1, nbc
!        write(105,*) cmplx( psi%valr( bciter, bviter, 1, 1 ), psi%vali( bciter, bviter, 1, 1 ), DP )
!      enddo
!    enddo
!    write(105,*)
!
    !
    if( sys%valence_ham_spin .eq. 1 ) then
     spin_prefac = 2.0_dp
    else
      spin_prefac = 1.0_dp
    endif
    minus_spin_prefac = -spin_prefac
      
    allocate( re_l_bubble( max(1,nxpts) ), im_l_bubble( max(1,nxpts) ) )
    re_l_bubble( : ) = 0.0_dp
    im_l_bubble( : ) = 0.0_dp


    nthreads = 1
!$  nthreads = omp_get_max_threads()
!    xwidth = max( (nxpts_pad/cache_double) * nkpts / nthreads, 1 )
!    xwidth = xwidth * cache_double

    ! If more threads than kpts
    xwidth = (nxpts_pad/cache_double) * nkpts / nthreads
    xwidth = max( xwidth / nkpts, 1 )

!    xwidth = max( nxpts_pad / cache_double, 1 )
    xwidth = min( xwidth * cache_double, nxpts_pad )
    
!    write(6,*) nxpts_pad, xwidth, startx_by_mpiID( myid )

    allocate( re_amat( max(1,nxpts_pad), nbv, nkpts ), im_amat( max(nxpts_pad,1), nbv, nkpts ) )

!$OMP PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( sys, psi, re_amat, im_amat, nxpts, nkpts, nbv, nbc, re_con, im_con ) &
!$OMP& SHARED( nxpts_pad, re_val, im_val, psi_con_pad, im_l_bubble, re_l_bubble, nproc, myid ) &
!$OMP& SHARED( re_scratch, im_scratch, nxpts_by_mpiID, comm, MPI_STATUS_IGNORE, ierr, nthreads ) &
!$OMP& SHARED( scratch, fo, bubble, spin_prefac, psiout, minus_spin_prefac, startx_by_mpiID ) &
!$OMP& PRIVATE( bciter, bviter, xstop, xwidth, ix, iix, ik, ispn, dft_spin, psi_spin, dft_spin2, psi_spin2 )




  do ispn = 1, sys%valence_ham_spin
    psi_spin = ispn * ispn  ! Either 1 or 4 
    dft_spin = min( ispn, sys%nspn )   ! The DFT states need not be different for spin up/down

    if( nxpts .gt. 0 ) then
    xwidth = nxpts
    ix = 1
!$OMP DO COLLAPSE(1)
      do ik = 1, nkpts
!      do ix = 1, nxpts, xwidth
          call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik, dft_spin), nxpts_pad, &
                      psi%valr( 1, 1, ik, psi_spin ), psi_con_pad, zero, re_amat( ix, 1, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', xwidth, nbv, nbc, minusone, im_con(ix, 1, ik,dft_spin), nxpts_pad, &
                      psi%vali( 1, 1, ik, psi_spin ), psi_con_pad, one, re_amat( ix, 1, ik ), nxpts_pad )

          call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, im_con(ix, 1, ik,dft_spin), nxpts_pad, &
                      psi%valr( 1, 1, ik, psi_spin ), psi_con_pad, zero, im_amat( ix, 1, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik,dft_spin), nxpts_pad, &
                      psi%vali( 1, 1, ik, psi_spin ), psi_con_pad, one, im_amat( ix, 1, ik ), nxpts_pad )
!      enddo
      enddo
!$OMP END DO
    else 
      re_amat( :,:,: ) = 0.0_dp
      im_amat( :,:,: ) = 0.0_dp
    endif

!  write(6,*) 'A:', maxval( re_amat ), maxval( im_amat )

! Should move zero-ing of re/im_l_bubble here to get omp lined up correctly

    xwidth = nxpts / ( nthreads * 8 )
    xwidth = max( 1, xwidth )
    xwidth = xwidth * 8

!$OMP DO 
  do ix = 1, nxpts, xwidth
    xstop = min( nxpts, ix + xwidth - 1 )
!      ix = 1
!      xstop = nxpts
          do ik = 1, nkpts
            do bviter = 1, nbv
              do iix = ix, xstop
                re_l_bubble( iix ) = re_l_bubble( iix ) &
                                   + re_val( iix, bviter, ik, dft_spin ) * re_amat( iix, bviter, ik ) &
                                   + im_val( iix, bviter, ik, dft_spin ) * im_amat( iix, bviter, ik ) 
                im_l_bubble( iix ) = im_l_bubble( iix ) &
                                   + re_val( iix, bviter, ik, dft_spin ) * im_amat( iix, bviter, ik ) &
                                   - im_val( iix, bviter, ik, dft_spin ) * re_amat( iix, bviter, ik ) 
              enddo
            enddo
        enddo
  enddo
!$OMP END DO

!    write(6,*) maxval( re_l_bubble ), maxval( im_l_bubble )

!    write(6,*) startx_by_mpiID( 0 ), nxpts
!$OMP MASTER
      do iproc = 0, nproc - 1
        if( myid .eq. root ) then
          if( iproc .eq. myid ) then
!          write(6,*) startx_by_mpiID( iproc ), startx_by_mpiID( iproc ) + nxpts - 1
            re_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts - 1, ispn ) = re_l_bubble( : )
            im_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts - 1, ispn ) = im_l_bubble( : )
#ifdef MPI
          else
            call MPI_RECV( re_scratch( startx_by_mpiID( iproc ), ispn ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
                           iproc, iproc, comm, MPI_STATUS_IGNORE, ierr )
            call MPI_RECV( im_scratch( startx_by_mpiID( iproc ), ispn ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
                           iproc, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
#endif
          endif
#ifdef MPI
        elseif( iproc .eq. myid ) then
          call MPI_SEND( re_l_bubble, nxpts, MPI_DOUBLE_PRECISION, root, iproc, comm, ierr )
          call MPI_SEND( im_l_bubble, nxpts, MPI_DOUBLE_PRECISION, root, iproc+nproc, comm, ierr )
#endif
        endif
      enddo
! $OMP END MASTER

!! There is no implied barrier at the end of master !!
! $OMP BARRIER


      if( myid .eq. root ) then
        scratch(:) = cmplx( re_scratch(:, ispn), im_scratch(:,ispn), DP ) 
        call FFT_wrapper_single( scratch, OCEAN_BACKWARD, fo, .false. )
        scratch( : ) = scratch( : ) * bubble( : )
        call FFT_wrapper_single( scratch, OCEAN_FORWARD, fo, .false. )

!        if( ispn .eq. 1 ) then
          re_scratch(:,1) = real(scratch(:),DP) 
          im_scratch(:,1) = real(aimag(scratch(:)),DP)
!        else
!          re_scratch(:,1) = re_scratch(:,1) + real(scratch(:),DP)
!          im_scratch(:,1) = im_scratch(:,1) + real(aimag(scratch(:)),DP)
!        endif
      endif


!    enddo  ! ispn
    ! At this point re/im_scratch (on root) has the sum of the action of the exchange on 
    ! both the up-up and down-down electron-hole pairs



! $OMP MASTER
    do iproc = 0, nproc - 1
      if( myid .eq. root ) then
        if( iproc .eq. myid ) then
!          write(6,*) startx_by_mpiID( iproc ), startx_by_mpiID( iproc ) + nxpts - 1
          re_l_bubble( : ) = re_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts - 1, 1 )
          im_l_bubble( : ) = im_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts - 1, 1 )
#ifdef MPI
        else
          call MPI_SEND( re_scratch( startx_by_mpiID( iproc ), 1 ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
                         iproc, iproc, comm, ierr )
          call MPI_SEND( im_scratch( startx_by_mpiID( iproc ), 1 ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
                         iproc, iproc+nproc, comm, ierr )
#endif
        endif
#ifdef MPI
      elseif( iproc .eq. myid ) then
        call MPI_RECV( re_l_bubble, nxpts, MPI_DOUBLE_PRECISION, root, iproc, comm, MPI_STATUS_IGNORE, ierr )
        call MPI_RECV( im_l_bubble, nxpts, MPI_DOUBLE_PRECISION, root, iproc+nproc, comm, MPI_STATUS_IGNORE, ierr )
#endif
      endif
    enddo

!$OMP END MASTER 

!! There is no implied barrier at the end of master !!
!$OMP BARRIER

  do jspn = 1, sys%valence_ham_spin
    psi_spin2 = jspn * jspn  ! Either 1 or 4 
    dft_spin2 = min( jspn, sys%nspn )   ! The DFT states need not be different for spin up/down


!$OMP WORKSHARE
    re_amat = 0.0_dp
    im_amat = 0.0_dp
!$OMP END WORKSHARE

!$OMP DO COLLAPSE(3)
      do ik = 1, nkpts
        do bviter = 1, nbv
!        do ix = 1, nxpts, xwidth
!          xstop = min( nxpts, ix + xwidth - 1 )
!          do iix = 1, xstop
          do iix = 1, nxpts
              re_amat( iix, bviter, ik ) = re_l_bubble( iix ) * re_val( iix, bviter, ik, dft_spin2 )
              im_amat( iix, bviter, ik ) = im_l_bubble( iix ) * re_val( iix, bviter, ik, dft_spin2 )
!          enddo
!          do iix = 1, xstop
              re_amat( iix, bviter, ik ) = re_amat( iix, bviter, ik ) & 
                                         - im_l_bubble( iix ) * im_val( iix, bviter, ik, dft_spin2 )
              im_amat( iix, bviter, ik ) = im_amat( iix, bviter, ik ) &
                                         + re_l_bubble( iix ) * im_val( iix, bviter, ik, dft_spin2 )
          enddo
!          enddo
!        enddo
        enddo
      enddo
!$OMP END DO


    if( nxpts .gt. 0 ) then
!JTV add additional tiling
!$OMP DO 

      do ik = 1, nkpts
        call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, dft_spin2 ), nxpts_pad, & 
                    re_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%valr( 1, 1, ik, psi_spin2 ), psi_con_pad )
        call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, im_con( 1, 1, ik, dft_spin2 ), nxpts_pad, & 
                    im_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%valr( 1, 1, ik, psi_spin2 ), psi_con_pad )

        call DGEMM( 'T', 'N', nbc, nbv, nxpts, minus_spin_prefac, im_con( 1, 1, ik, dft_spin2 ), nxpts_pad, &
                    re_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%vali( 1, 1, ik, psi_spin2 ), psi_con_pad )
        call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, dft_spin2 ), nxpts_pad, &
                    im_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%vali( 1, 1, ik, psi_spin2 ), psi_con_pad )
      enddo
!$OMP END DO

    endif

  enddo ! jspn

  enddo ! ipsn

!$OMP END PARALLEL

    ! 
    !
    !
    deallocate( re_l_bubble, im_l_bubble, re_amat, im_amat )
111 continue
    if( ierr .ne. 0 ) then
      write(1000+myid,*) 'Error in AI_par_bubble.f90, quitting...'
    endif

  end subroutine AI_bubble_act




end module OCEAN_bubble
