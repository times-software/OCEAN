! This is the RPA bubble diagrams

module OCEAN_bubble

  use AI_kinds
  implicit none
  private
  save


  complex( DP ), allocatable :: scratch( : )
  real( DP ), allocatable  :: bubble( :, :, : )
  real(DP), allocatable :: re_scratch(:), im_scratch(:)

  integer, allocatable :: xstart_by_mpiID(:)

  integer*8        :: fplan
  integer*8        :: bplan

  public :: AI_bubble_prep

  contains

  subroutine AI_bubble_prep( sys, ierr )
    use OCEAN_mpi
    use OCEAN_val_states, only : nxpts_by_mpiID( : )
    !
    implicit none
    include 'fftw3.f03'
    !
    !
    integer :: i, j
!$  integer, external :: omp_get_max_threads
    real( DP ) :: length
    !
    if( myid .eq. root ) then
      allocate( bubble(  sys%nxpts ), &
                scratch( sys%nxpts ), re_scratch( sys%nxpts ), im_scratch( sys%nxpts ), &
                STAT=ierr )
      if( ierr .ne. 0 ) then
        write( 1000+myid ) 'Failed to allocate bubble and scratch'
        goto 111
      endif
      call dfftw_plan_with_nthreads( omp_get_max_threads() )
      call dfftw_plan_dft_3d( fplan, sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), &
                              scratch, scratch, FFTW_FORWARD, FFTW_PATIENT )
      call dfftw_plan_dft_3d( bplan, sys%xmesh(3), sys%xmesh(2), sys%xmesh(1), &
                              scratch, scratch, FFTW_BACKWARD, FFTW_PATIENT )
      call dfftw_plan_with_nthreads( 1 )
      call AI_bubble_wall_search( sys, length, ierr )
      if( ierr .ne. 0 ) goto 111
      call AI_bubble_populate( sys, length, ierr )
      if( ierr .ne. 0 ) goto 111

      allocate( xstart_by_mpiID( 0:nproc-1 ), STAT=ierr )
      if( ierr .ne. 0 ) goto 111
      j = 1
      do i = 0, nproc - 1
        xstart_by_mpiID( i ) = j
        j = j + nxpts_by_mpiID( i )
      enddo

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
    ystemype( o ), intent( in ) :: sys
    real( DP ), intent( out ) :: length
    integer( S_INT ), intent( inout ) :: ierr
    !
    integer :: ix, iy, iz, gvec( 3 )
    real( DP ) :: min_edge_length
    !
    real( DP ), external :: gvec_length
    !
    ! exclude anything greater than or equal to the min distance we missed.
    !  ie. if x goes from -5 to 5 than the |x| = 6 surface will tell us what we missed
    !  this may be slightly more accepting than taking everything equal or less than the |x|=5 surface
    !  since it may include some corners that would otherwise be thrown out.
    gvec = (/ floor( dble(sys%xmesh( 1 )) / 2d0 ) + 1, floor(dble(sys%xmesh( 2 )) / 2d0 )+ 1, &
              floor(dble(sys%xmesh( 3 )) / 2d0 ) + 1/)
    min_edge_length = gvec_length( gvec, sys%bvecs )
    do ix = ceiling(dble(-sys%xmesh(1))/2d0) , floor(dble(sys%xmesh(1))/2d0  )
      do iy = ceiling(dble(-sys%xmesh(2)) / 2d0) , floor(dble(sys%xmesh( 2 )) / 2d0 )
        gvec = (/ ix, iy, ceiling(dble(-sys%xmesh(3))/2d0)  /)
        length = gvec_length( gvec, sys%bvecs )
        if ( length .lt. min_edge_length ) min_edge_length = length
        gvec = (/ ix, iy, floor(dble(sys%xmesh(3))/2d0) /)
        length = gvec_length( gvec, sys%bvecs )
        if ( length .lt. min_edge_length ) min_edge_length = length
      enddo
      do iz = ceiling(dble(-sys%xmesh( 3 )) / 2d0) , floor(dble(sys%xmesh( 3 )) / 2d0  )
        gvec = (/ ix, ceiling(dble(-sys%xmesh(2))/2d0) , iz /)
        length = gvec_length( gvec, sys%bvecs )
        if ( length .lt. min_edge_length ) min_edge_length = length
        gvec = (/ ix, floor(dble(sys%xmesh(2))/2d0) , iz /)
        length = gvec_length( gvec, sys%bvecs )
        if ( length .lt. min_edge_length ) min_edge_length = length
      enddo
    enddo
    do iy = ceiling(dble(-sys%xmesh(2)) / 2d0) , floor(dble(sys%xmesh( 2 )) / 2d0 )
      do iz = ceiling(dble(-sys%xmesh( 3 )) / 2d0) , floor(dble(sys%xmesh( 3 )) / 2d0 )
        gvec = (/ ceiling(dble(-sys%xmesh(1))/2d0) , iy, iz /)
        length = gvec_length( gvec, sys%bvecs )
        if ( length .lt. min_edge_length ) min_edge_length = length
        gvec = (/ floor(dble(sys%xmesh(2))/2d0) , iy, iz /)
        length = gvec_length( gvec, sys%bvecs )
        if ( length .lt. min_edge_length ) min_edge_length = length
      enddo
    enddo
    length = min_edge_length
    write( 6, * ) 'Max gvec length: ', length
  end subroutine AI_bubble_wall_search
!
!
! 
  subroutine AI_bubble_populate( sys, length, ierr )
    use ai_kinds
    use m_ai_systems
    !
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    real( DP ), intent( in ) :: length
    integer( S_INT ), intent( inout ) :: ierr
    !
    integer :: ix, iy, iz, ij, ik, igvec, izz, iyy, ixx, temp_gvec( 3 ), iter
    real( DP ) :: mul, gsqd, qq( 3 ), bmet( 3, 3 ), pi
    real( DP ), external :: gvec_length
    ! 
    igvec = 0
    pi = 4.0d0 * atan( 1.0d0 )
    do ij = 1, 3
      do ik = 1, 3
        bmet(ij,ik) = dot_product( p_sys%bvecs(:,ij), p_sys%bvecs( :, ik ) )
      enddo
    enddo
    !
    iter = 0
    bubble( : ) = 0.d0
    write( 6, * ) 'qinb: ', probe%qinb( : )
    !
    do ix = ceiling( dble(-sys%xmesh(1))/2d0 ), floor( dble(sys%xmesh(1))/2d0 )
      temp_gvec( 1 ) = ix 
      do iy =  ceiling( dble(-sys%xmesh(2))/2d0 ), floor( dble(sys%xmesh(2))/2d0 )
        temp_gvec( 2 ) =  iy 
        do iz = ceiling( dble(-sys%xmesh(3))/2d0 ), floor( dble(sys%xmesh(3))/2d0 )
          temp_gvec( 3 ) = iz 
          qq( : ) = sys%qinunitsofbvectors( : ) + real( temp_gvec( : ) )
          gsqd = 0.d0
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
          if( gvec_length( temp_gvec, p_sys%bvecs ) .ge. length ) then
            mul = 0.d0
          else
            mul = 4.0d0 * pi * 27.2114d0 / ( sys%nkpts * sys%celvol * gsqd )
            igvec = igvec + 1
          endif
!          bubble( izz, iyy, ixx ) = mul
            bubble( iter ) = mul
        enddo
      enddo 
    enddo 
    igvec = igvec - 1
    bubble( 1 ) = 0.d0 ! skip G = 0
    write( 6, * )'Num gvecs retained: ', igvec, sys%nxpts
  end subroutine AI_bubble_populate
!
!
!
  subroutine AI_bubble_act( psi, psiout, ierr )
    use OCEAN_val_states
    use OCEAN_psi
    use OCEAN_mpi
    implicit none
    include 'fftw3.f03'
    !
    type( OCEAN_vector ), intent( in ) :: psi
    type( OCEAN_vector ), intent( inout ) :: psiout
    integer, intent( inout ) :: ierr
    !
    integer :: bciter, bviter, ik
    integer :: ix, iix, xstop, xwidth, nthreads
    integer, external :: omp_get_max_threads
    real( DP ), allocatable :: re_l_bubble( : ), im_l_bubble
    real( DP ), allocatable :: re_amat( :, : ), im_amat( :, : )
    !
    real( DP ), parameter :: one = 1.0_dp
    real( DP ), parameter :: zero = 0.0_dp
    real( DP ), parameter :: minusone = -1.0_dp
    !
    !
    allocate( re_l_bubble( nxpts ), im_l_bubble )
    re_l_bubble( : ) = 0.0_dp
    im_l_bubble( : ) = 0.0_dp


!$  nthreads = omp_get_max_threads()
    xwidth = max( (nxpts_pad/cache_double) * nkpts / nthreads, 1 )
    xwidth = xwidth * cache_double
    

    allocate( re_amat( nxpts_pad, nbv, nkpts ), im_amat( nxpts_pad, nbv, nkpts ) )

!$OMP PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( sys, val_ufunc, con_ufunc, psi, sum_l_bubble, xwidth, re_amat, im_amat ) &
!$OMP& PRIVATE( bciter, bviter, lcindx, xstop )



!$OMP DO COLLAPSE(2)
    do ik = 1, nkpts
      do ix = 1, nxpts, xwidth
        call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik), nxpts_pad, &
                    psi%valr( 1, 1, ik, 1 ), psi%cband, zero, re_amat( ix, 1, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', xwidth, nbv, nbc, minusone, re_con(ix, 1, ik), nxpts_pad, &
                    psi%valr( 1, 1, ik, 1 ), psi%cband, one, re_amat( ix, 1, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', xwdith, nbv, nbc, one, im_con(ix, 1, ik), nxpts_pad, &
                    psi%valr( 1, 1, ik, 1 ), psi%cband, zero, im_amat( ix, 1, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik), nxpts_pad, &
                    psi%vali( 1, 1, ik, 1 ), psi%cband, one, im_amat( ix, 1, ik ), nxpts_pad )
      enddo
    enddo
!$OMP END DO

!$OMP DO 
    do ix = 1, nxpts, xwidth
      xstop = min( nxpts, ix + xwidth - 1 )
      do bciter = 1, nbc
        do bviter = 1, nbv
          do ik = 1, nkpts
            do iix = ix, xstop
              re_l_bubble( iix ) = re_l_bubble( iix ) &
                                 + re_val( iix, bviter, ik, 1 ) * re_amat( iix, bviter, ik ) &
                                 + im_val( iix, bviter, ik, 1 ) * im_amat( iix, bciter, ik ) 
              im_l_bubble( iix ) = im_l_bubble( iix ) &
                                 + re_val( iix, bviter, ik, 1 ) * im_amat( iix, bviter, ik ) &
                                 - im_val( iix, bviter, ik, 1 ) * re_amat( iix, bciter, ik ) 
            enddo
          enddo
        enddo
      enddo
    enddo
!$OMP END DO


!$OMP MASTER
    do iproc = 0, nproc - 1
      if( myid .eq. root ) then
        if( iproc .eq. myid ) then
          re_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts ) = re_l_bubble( : )
          im_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts ) = im_l_bubble( : )
#ifdef MPI
        else
          call MPI_RECV( re_scratch( startx_by_mpiID( iproc ) ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
                         iproc, iproc, comm, MPI_STATUS_IGNORE, ierr )
          call MPI_RECV( im_scratch( startx_by_mpiID( iproc ) ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
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

    if( myid .eq. root )
      scratch(:) = cmplx( re_scratch(:), im_scratch(:) )
      call dfftw_execute_dft(fplan, scratch, scratch)
      scratch( : ) = scratch( : ) * bubble( : )
      call dfftw_execute_dft(bplan, scratch, scratch)
      re_scratch(:) = real(scratch(:))
      im_scratch(:) = aimag(scratch(:))
    endif



!$OMP MASTER
    do iproc = 0, nproc - 1
      if( myid .eq. root ) then
        if( iproc .eq. myid ) then
          re_l_bubble( : ) = re_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts )
          im_l_bubble( : ) = im_scratch( startx_by_mpiID( iproc ) : startx_by_mpiID( iproc ) + nxpts )
#ifdef MPI
        else
          call MPI_SEND( re_scratch( startx_by_mpiID( iproc ) ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
                         iproc, iproc, comm, ierr )
          call MPI_SEND( im_scratch( startx_by_mpiID( iproc ) ), nxpts_by_mpiID( iproc ), MPI_DOUBLE_PRECISION, &
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

!$OMP DO COLLAPSE(3)
    do ik = 1, nkpts
      do bviter = 1, nbv
        do ix = 1, nxpts, xwidth
          xstop = min( nxpts, ix + xwidth - 1 )
          do iix = 1, xstop
            re_amat( iix, bviter, ik ) = re_l_bubble( iix ) * re_val( iix, bviter, ik, 1 )
            im_amat( iix, bviter, ik ) = im_l_bubble( iix ) * re_val( iix, bviter, ik, 1 )
          enddo
          do iix = 1, xstop
            re_amat( iix, bviter, ik ) = re_amat( iix, bviter, ik ) - im_l_bubble( iix ) * im_val( iix, bviter, ik, 1 )
            im_amat( iix, bviter, ik ) = im_amat( iix, bviter, ik ) + re_l_bubble( iix ) * im_val( iix, bviter, ik, 1 )
          enddo
        enddo
      enddo
    enddo
!$OMP END DO

!JTV add additional tiling
!$OMP DO 
    do ik = 1, nkpts
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, one, re_con( 1, 1, ik, 1 ), nxpts_pad, re_amat( 1, 1, ik ), nxpts_pad, &
                  zero, psiout%rval( 1, 1, ik, 1 ), psi%cband )
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, one, im_con( 1, 1, ik, 1 ), nxpts_pad, im_amat( 1, 1, ik ), nxpts_pad, &
                  one, psiout%rval( 1, 1, ik, 1 ), psi%cband )
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, minusone, im_con( 1, 1, ik, 1 ), nxpts_pad, re_amat( 1, 1, ik ), nxpts_pad, &
                  zero, psiout%ival( 1, 1, ik, 1 ), psi%cband )
      call DGEMM( 'T', 'N', nbc, nbv, nxpts, one, re_con( 1, 1, ik, 1 ), nxpts_pad, im_amat( 1, 1, ik ), nxpts_pad, &
                  one, psiout%ival( 1, 1, ik, 1 ), psi%cband )
    enddo
!$OMP END DO


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
