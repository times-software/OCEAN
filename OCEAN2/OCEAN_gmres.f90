! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!> @brief This module contains the necessary pieces for carrying out GMRES
module OCEAN_gmres
  use OCEAN_psi, only : ocean_vector
  use AI_kinds

  implicit none
  private
  save
    
  complex(DP), allocatable :: c_matrix( : )

  type(ocean_vector), allocatable, target :: au_matrix(:)
  type(ocean_vector), allocatable, target :: u_matrix(:)

  real(DP) :: gmres_resolution !> imaginary part added to inverse, stored in Ha
  real(DP) :: gmres_preconditioner !> width of preconditioner 
  real(DP) :: gmres_convergence 
  integer :: gmres_depth ! max size of gmres Arnoldi? space before restart

  integer :: gmres_nsteps
  real(DP), allocatable :: gmres_energy_list(:)

  real(DP) :: interaction_scale = 1.0_DP
  logical :: echamp
  logical :: do_precondition = .true.
  logical :: allow_reuse_x = .true.
  logical :: force_reuse_x = .true.
  logical :: allow_recycle = .true.

  public :: OCEAN_gmres_do, OCEAN_gmres_setup, OCEAN_gmres_clean, OCEAN_gmres_do_recycle

  contains

  subroutine OCEAN_gmres_clean( ierr )
    use OCEAN_psi, only : OCEAN_psi_kill
    !
    integer, intent( inout ) :: ierr
    !
    integer :: i, n
    !
    if( allocated( c_matrix ) ) deallocate( c_matrix ) 
    if( allocated( gmres_energy_list ) ) deallocate( gmres_energy_list ) 

    if( allocated( au_matrix ) ) then
      n = size( au_matrix, 1 )
      do i = 1, n
        call OCEAN_psi_kill( au_matrix( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
      deallocate( au_matrix )
    endif

    if( allocated( u_matrix ) ) then
      n = size( u_matrix, 1 )
      do i = 1, n
        call OCEAN_psi_kill( u_matrix( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
      deallocate( u_matrix )
    endif

  end subroutine OCEAN_gmres_clean

  subroutine OCEAN_gmres_setup( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi!, only : myid, root, comm, MPI_BCAST, MPI_INTEGER, MPI_DOUBLE_PRECISION
    use OCEAN_corewidths, only : returnLifetime
    use OCEAN_constants, only : eV2Hartree
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    real(dp) :: e_start, e_stop, e_step
    character(len=4) :: inv_style
    integer :: ierr_, i
    !
    if( myid .eq. root ) then
!      open(unit=99,file='gmres.in',form='formatted',status='old',iostat=ierr,ERR=100)
      open(unit=99,file='bse.in',form='formatted',status='old',iostat=ierr,ERR=100)
      read(99,*)
      read(99,*)
      read(99,*)

      read(99,*) gmres_depth, gmres_resolution, gmres_preconditioner, gmres_convergence
      !
      ! if gres is negative fill it with core-hole lifetime broadening
      if( gmres_resolution .lt. 1.0d-20 ) then
        call returnLifetime( sys%cur_run%ZNL(1), sys%cur_run%ZNL(2), sys%cur_run%ZNL(3), gmres_resolution )
        if( gmres_resolution .le. 0.0_dp ) gmres_resolution = 0.1_dp
      endif
      gmres_resolution = gmres_resolution * eV2Hartree
      gmres_preconditioner = gmres_preconditioner * eV2Hartree
      !
      read(99,*) inv_style
      select case( inv_style )

        case('list')
          read(99,*) gmres_nsteps
          allocate( gmres_energy_list( gmres_nsteps ) )
          do i = 1, gmres_nsteps
            read(99,*) gmres_energy_list( i )
          enddo

        case('loop')
          read(99,*) e_start, e_stop, e_step
          gmres_nsteps = floor( ( e_stop - e_start + e_step*.9) / e_step )
          allocate( gmres_energy_list( gmres_nsteps ) )
          do i = 1, gmres_nsteps
            gmres_energy_list( i ) = e_start + ( i - 1 ) * e_step
          enddo

        case default
          write(6,*) 'Error in gmres.in! Unsupported energy point style', inv_style
          ierr = -1

      end select
      close( 99 )

      write(6,*) 'GMRES energy steps: ', gmres_nsteps

      inquire( file="echamp.inp", exist=echamp )
      if( echamp ) then
        open( unit=99, file="echamp.inp", form="formatted", status="old" )
        read( 99, * ) echamp
        close( 99 )
      endif
    endif
100 continue
    if( myid .eq. 0 .and. ierr .ne. 0 ) then
      write(6,*) 'Error reading gmres.in/bse.in!'
    endif
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    call MPI_BCAST( gmres_depth, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_resolution, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_convergence, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_preconditioner, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_nsteps, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( myid .ne. root ) allocate( gmres_energy_list( gmres_nsteps ) )
    call MPI_BCAST( gmres_energy_list, gmres_nsteps, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( echamp, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'Finished gmres setup'
    call initialize_gmres_storage( sys, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'Finished gmres setup'

  end subroutine OCEAN_gmres_setup

  subroutine initialize_gmres_storage( sys, ierr )
    use OCEAN_psi, only : OCEAN_psi_init, OCEAN_psi_new, OCEAN_psi_size_min
    use OCEAN_system, only : o_system
    use OCEAN_mpi, only : myid, root
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer :: i

    allocate( au_matrix( gmres_depth ) )
    allocate( u_matrix( gmres_depth ) )

    call OCEAN_psi_init( sys, ierr )
    if( ierr .ne. 0 ) return

    do i = 1, gmres_depth
      call OCEAN_psi_new( au_matrix( i ), ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_new( u_matrix( i ), ierr )
      if( ierr .ne. 0 ) return
    enddo

    if( myid .eq. root ) then
      i = OCEAN_psi_size_min(au_matrix(1))
      write(6,'(A,F8.2,A)') 'GMRES requires additional: ', gmres_depth*dble(i)*(32.0_dp/1048576.0_dp), ' MB'
    endif

    if( myid .eq. root ) then
      allocate( c_matrix( gmres_depth*(gmres_depth+1)/2 ) )
    else
      allocate( c_matrix( 1 ) )
    endif

  end subroutine initialize_gmres_storage

  subroutine make_recycle( current_iter, ener, prev_ener, gprc, hpsi1, psi_pcdiv, ierr )
    use OCEAN_psi, only : OCEAN_psi_axpy, OCEAN_psi_min_set_prec, OCEAN_psi_3element_mult
    integer, intent( in ) :: current_iter
    real(DP), intent( in ) :: ener, prev_ener, gprc
    type(OCEAN_vector), intent(inout) :: hpsi1
    type(OCEAN_vector), intent(inout) :: psi_pcdiv
    integer, intent( inout ) :: ierr
    !
    integer :: i
    real(DP) :: delta

    delta = ener - prev_ener

    do i = 1, current_iter
      call OCEAN_psi_axpy( delta, u_matrix( i ), au_matrix( i ), ierr )
      if( ierr .ne. 0 ) return
    enddo

    !JTv think this is not correct
!    call OCEAN_psi_min_set_prec( ener, gprc, hpsi1, psi_pcdiv, ierr, prev_ener )
!    do i = 1, current_iter
!      call OCEAN_psi_3element_mult( u_matrix(i), u_matrix(i), psi_pcdiv, ierr )
!      if( ierr .ne. 0 ) return
!      call OCEAN_psi_3element_mult( au_matrix(i), au_matrix(i), psi_pcdiv, ierr )
!      if( ierr .ne. 0 ) return
!    enddo
    

  end subroutine make_recycle

  subroutine update_gmres( current_iter, psi_g, psi_x, psi_ax, ierr)
    use OCEAN_mpi
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_dot, OCEAN_psi_axpy, OCEAN_psi_nrm
    use OCEAN_timekeeper, only : tk_inv, OCEAN_tk_start, OCEAN_tk_stop
    !
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent(inout ) :: psi_g, psi_x, psi_ax
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: c_temp( : ), coeff( : )
    real(DP), allocatable, dimension(:) :: re_coeff_vec, im_coeff_vec, re_rvec, im_rvec 
    integer, allocatable, dimension(:) :: re_coeff_request, im_coeff_request, re_rvec_request, im_rvec_request
    integer, allocatable :: ipiv( : )
    real(DP) :: tmp_r, tmp_i
    integer :: local_gmres_size, info, iter, c_size, iter_start
    !
    call OCEAN_tk_start( tk_inv )
    !
    allocate( re_coeff_vec( current_iter ), im_coeff_vec( current_iter ), & 
              re_rvec( current_iter ), im_rvec( current_iter ) )    

    allocate( re_coeff_request( current_iter ), im_coeff_request( current_iter ), &
              re_rvec_request( current_iter ), im_rvec_request( current_iter ) )
    !
    !
    ! For iter = current_iter we can call psi_nrm instead of psi_dot
    !  this will save a small amount of time, but need to correctly 
    !  set the imaginary value and request to 0/null 
    im_coeff_request( current_iter ) = MPI_REQUEST_NULL
    im_coeff_vec( current_iter ) = 0.0_DP
    !
    do iter = 1, current_iter - 1
      call OCEAN_psi_dot( au_matrix( current_iter ), au_matrix( iter ), &
                          re_coeff_request( iter ), re_coeff_vec( iter ), ierr, & 
                          im_coeff_request( iter ), im_coeff_vec( iter ) )

      call OCEAN_psi_dot( au_matrix( iter ), psi_g, &
                          re_rvec_request( iter ), re_rvec( iter ), ierr , &
                          im_rvec_request( iter ), im_rvec( iter ) )
    enddo

    call OCEAN_psi_nrm( re_coeff_vec( current_iter ), au_matrix( current_iter ), ierr, & 
                        re_coeff_request( current_iter ) )
    call OCEAN_psi_dot( au_matrix( current_iter ), psi_g, &
                        re_rvec_request( current_iter ), re_rvec( current_iter ), ierr , &
                        im_rvec_request( current_iter ), im_rvec( current_iter ) )

    call MPI_WAITALL( current_iter, re_coeff_request, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    call MPI_WAITALL( current_iter, im_coeff_request, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    !
    allocate( coeff( current_iter ) )
    if( myid .eq. root ) then
      iter_start = (current_iter * ( current_iter - 1 ) ) / 2
!      if( current_iter .eq. 1 ) then
!        write(6,*) iter_start, cmplx( re_coeff_vec( 1 ), im_coeff_vec( 1 ), DP )
!      endif
      do iter = 1, current_iter
        c_matrix( iter + iter_start ) = cmplx( re_coeff_vec( iter ), -im_coeff_vec( iter ), DP )
      enddo
      
      c_size = current_iter * (current_iter + 1 ) / 2 
      allocate( c_temp( c_size ), ipiv( current_iter ) )
      c_temp( : ) = c_matrix( 1 : c_size )
!      if( current_iter .le. 3 ) then
!        write(6,*) c_temp( : )
!      endif


!      if( current_iter .gt. 1 ) then
        call zhptrf( 'U', current_iter, c_temp, ipiv, info )
!      else 
!        info = 0
!      endif

      call MPI_WAITALL( current_iter, re_rvec_request, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( current_iter, im_rvec_request, MPI_STATUSES_IGNORE, ierr )

      ! Want to test after waitall so that all procs are on the same page if we abort
      if( info .ne. 0 ) then
        ierr = info
      else
        do iter = 1, current_iter
          coeff( iter ) = -cmplx( re_rvec( iter ), im_rvec( iter ), DP )
        enddo
!        if( current_iter .le. 3 ) then
!          write(6,*) coeff(:)
!        endif


!       if( current_iter .gt. 1 ) then
          call zhptrs( 'U', current_iter, 1, c_temp, ipiv, coeff, current_iter, info )
!        else
!          info = 0
!          coeff(1) = coeff(1)/c_temp(1)
!        endif
      endif
      deallocate( c_temp, ipiv )

    else
      call MPI_WAITALL( current_iter, re_rvec_request, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( current_iter, im_rvec_request, MPI_STATUSES_IGNORE, ierr )
    endif

    call MPI_BCAST( info, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    ierr = info
    if( ierr .ne. 0 ) return

    call MPI_BCAST( coeff, current_iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return

!    if( current_iter .le. 3 ) then
!      write(6,*) coeff(:)
!    endif


    ! update ax and x
    do iter = 1, current_iter
      tmp_r = real( coeff( iter ), DP )
      tmp_i = aimag( coeff( iter ) )
      call OCEAN_psi_axpy( tmp_r, u_matrix( iter ), psi_x, ierr, tmp_i )
      call OCEAN_psi_axpy( tmp_r, au_matrix( iter ), psi_ax, ierr, tmp_i )
    enddo

    deallocate( coeff, re_coeff_vec, im_coeff_vec, re_rvec, im_rvec, & 
                re_coeff_request, im_coeff_request, re_rvec_request, im_rvec_request )
    
    call OCEAN_tk_stop( tk_inv )
    !
  end subroutine update_gmres

  subroutine eigensystem_gmres( iter, keep, ierr )
!    use ocean_mpi, only : root, myid, MPI_REQUEST_NULL,  MPI_STATUSES_IGNORE
    use ocean_mpi
    use OCEAN_psi
    integer, intent( in ) :: iter
    integer, intent( inout ) :: keep
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: hess( :,: ), w(:), z(:,:), work(:), vl(:,:), vr(:,:)
    real(DP), allocatable :: re_hess(:,:), im_hess(:,:), re_coeff_vec(:), im_coeff_vec(:), rwork(:)
    integer, allocatable :: re_request(:,:), im_request(:,:), keep_list(:)
    real(DP) :: itmp, rtmp
    integer :: i, j, lwork, k, l, iter_start, ldvl
    type(ocean_vector), allocatable :: new_u(:)

    
    allocate( re_hess( iter, iter ), im_hess( iter, iter ),  &
              re_request( iter, iter ), im_request( iter, iter ), STAT=ierr )
    if( ierr .ne. 0 ) return

    re_hess = 0.0_DP
    im_hess = 0.0_DP
    re_request(:,:) = MPI_REQUEST_NULL
    im_request(:,:) = MPI_REQUEST_NULL
    
    do i = 1, iter
      do j = 1, min( i+1,iter)
        call OCEAN_psi_dot( u_matrix( j ), au_matrix( i ), &
                            re_request(j,i), re_hess( j,i ), ierr, &
                            im_request( j,i ), im_hess( j,i ) )
      enddo
    enddo

    call MPI_WAITALL( iter**2, re_request, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    call MPI_WAITALL( iter**2, im_request, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    
    allocate( w(iter), z(iter,iter) )
    if( myid .eq. root ) then
      lwork = iter * 12
      allocate( hess( iter, iter ), work( lwork ) )
      hess(:,:) = cmplx( re_hess(:,:), im_hess(:,:), DP )

      call zhseqr('S', 'I', iter, 1, iter, hess, iter, w, z, iter, work, lwork, ierr )
      if( ierr .ne. 0 ) then
        keep = 0 
      endif
      ldvl = iter
      allocate( vl(iter,iter), vr( iter, iter ), rwork( iter ) )
      call ztrevc( 'R', 'A', .true., iter, z, iter, vl, ldvl, vr, iter, iter, l, work, rwork, ierr)
      ierr = 0
      z(:,:) = vr(:,:)
      deallocate( hess, work, vl, vr )
    endif
    call MPI_BCAST( w, iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( z, iter*iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( keep, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    keep = min( keep, iter )
    allocate( keep_list( keep ) )
    write(6,*) myid, keep, iter
    keep_list(:) = 1
    do i = 2, keep
      k = i
      do j = 1, i-1
        if( abs(w(k)) .lt. abs(w(keep_list(j))) ) then
          l = keep_list(j)
          keep_list( j ) = k
          k = l
        endif
      enddo
      keep_list( i ) = k
    enddo

    do i = keep+1, iter
      k = i
      do j = 1, keep
        if( abs(w(k)) .lt. abs(w(keep_list(j))) ) then
          keep_list( j ) = k
          k = j
        endif
      enddo
    enddo


    do j = 1, keep
      i = keep_list(j)
      if( myid .eq. root ) write(6,*) 'EIG-S:', i, real(abs(w(i)))
       write(6,*) 'EIG-S:', myid, i, real(abs(w(i)))
    enddo

    call MPI_BARRIER( comm, ierr )

    allocate( new_u(keep) )
    do j = 1, keep
      call OCEAN_psi_new( new_u( j ), ierr )
      call OCEAN_psi_zero_min( new_u(j), ierr )
      if( ierr .ne. 0 ) return
      k = keep_list(j)
      do i = 1, iter
        if( k.lt.1 .or. k.gt.iter ) write(6,*) myid, j, k, keep_list(:)
        rtmp = real( z(i,k), DP )
        itmp = aimag( z(i,k) )
        call OCEAN_psi_axpy( rtmp, u_matrix(i), new_u(j), ierr, itmp )
        if( ierr .ne. 0 ) return
      enddo
    enddo

    ! update u and au 
    do j = 1, keep
      call OCEAN_psi_copy_min( u_matrix(j), new_u(j), ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_zero_min( au_matrix(j), ierr )
      if( ierr .ne. 0 ) return
      rtmp = real( w(keep_list(j) ), DP )
      itmp = aimag( w(keep_list(j) ) )
      call OCEAN_psi_axpy( rtmp, new_u(j), au_matrix(j), ierr, itmp )
      if( ierr .ne. 0 ) return
      
    enddo

    deallocate(keep_list, new_u)
    deallocate( re_hess, im_hess, re_request, im_request, w, z )

    keep = 1
    return
    allocate( re_coeff_vec(keep), im_coeff_vec(keep), re_request(keep,1), im_request(keep,1) )
    re_request(:,:) = MPI_REQUEST_NULL
    im_request(:,:) = MPI_REQUEST_NULL
    do j = 1, keep
      re_coeff_vec(:) = 0.0_DP
      im_coeff_vec(:) = 0.0_DP
      do i = 1, j-1
        call OCEAN_psi_dot( au_matrix( j ), au_matrix( i ), &
                          re_request( i,1 ), re_coeff_vec( i ), ierr, &
                          im_request( i,1 ), im_coeff_vec( i ) )
      enddo
      call OCEAN_psi_nrm( re_coeff_vec( j ), au_matrix( j ), ierr, &
                        re_request( j,1 ) )
      call MPI_WAITALL( j, re_request(1:j,1), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call MPI_WAITALL( j, im_request(1:j,1), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      iter_start = ( j * ( j - 1 ) ) / 2
      do i = 1, j
        c_matrix( i + iter_start ) = cmplx( re_coeff_vec( i ), -im_coeff_vec( i ), DP )
      enddo
    enddo

    deallocate( re_coeff_vec, im_coeff_vec, re_request, im_request )

  end subroutine eigensystem_gmres

  subroutine OCEAN_gmres_do_recycle( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_action, only : OCEAN_xact, OCEAN_action_h1
    use OCEAN_constants, only : eV2Hartree, Hartree2eV
    use OCEAN_mpi, only : myid, root, comm, MPI_STATUSES_IGNORE, MPI_REQUEST_NULL, MPI_INTEGER
    use OCEAN_energies, only : OCEAN_energies_allow
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
    type( ocean_vector ) :: psi_x, psi_ax, hpsi1, psi_g, psi_pcdiv
    type( ocean_vector), pointer :: psi_pg, psi_apg
    !
    real(DP) :: ener, rval, ival, gval, rel_error, fact, prev_energy
    real(dp), parameter :: one = 1.0_dp
    integer :: iter, step_iter, complete_iter, requests( 2 ), ierr_, prev_iter, global_iter, keep, i
    logical :: loud = .false.
    character( len=25 ) :: abs_filename
    integer, parameter :: abs_fh = 76
    type( ocean_vector ) :: hay_psi_x

    integer :: nhflag(6)
    nhflag(:) = sys%nhflag

    requests( : ) = MPI_REQUEST_NULL

    if( sys%cur_run%have_val ) then
      fact = hay_vec%kpref * 2.0_dp * sys%celvol
    else
      fact = hay_vec%kpref
    endif

    if( myid .eq. root ) then

      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref

      call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
      open( unit=abs_fh,file=abs_filename,form='formatted',status='unknown' )
      rewind( abs_fh )

    endif
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif


    ! Create all the psi vectors we need
    call OCEAN_psi_new( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( psi_g, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    !!!!


    ! Create the psi_min_arrays we will need
    ! AX etc

    ! Keep around H * (1,0) for making pcdiv for each loop
    call OCEAN_psi_one_full( psi_x, ierr )
    if( ierr .ne. 0 ) return

    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, psi_x, ierr )
      if( ierr .ne. 0 ) return
    endif

    call OCEAN_xact( sys, interaction_scale, psi_x, hpsi1, ierr,.false. , nhflag)
    if( ierr .ne. 0 ) return
    

    global_iter = 0
    prev_iter = 0
    prev_energy = -huge(prev_energy)

    do step_iter = 1, gmres_nsteps
      complete_iter = 0

!      ener = e_list( iter ) * eV2Hartree 
      ener = gmres_energy_list( step_iter ) * eV2Hartree
      if( myid .eq. root ) write(6,'(F20.8)') ener * Hartree2eV


      ! are we new or recycling?
      !   need to tune when we throw out because the preconditioner is no longer valid
      if( do_precondition .and. prev_iter .eq. 0 ) then 
        prev_energy = ener
        call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
        if( ierr .ne. 0 ) then
          if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
          return
        endif
      endif

      if( .not. allow_recycle ) then
        prev_iter = 0
      endif

      if( prev_iter .eq. 0 ) then
!        if( myid .eq. root ) write( 6, * ) 'no recycle: '
        call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, hay_vec, ierr , nhflag)
        if( ierr .ne. 0 ) then
          if( myid .eq. root ) write( 6,*) 'set_initial_vector failed'
          return
        endif
        if( loud ) then
          call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return

          call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
          if( ierr .ne. 0 ) return
          call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          if( myid .eq. 0 ) then
            write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, prev_iter, gmres_depth, &
                  gval, gmres_convergence, ener, (1.0_dp - rval), -ival
          endif
        endif
      else
        ! recycle au_matrix
        if( myid .eq. root ) write( 6, * ) 'recycle: ', prev_iter, (gmres_energy_list( step_iter ) - gmres_energy_list( step_iter -1))* eV2Hartree, (ener - prev_energy)
        call make_recycle( prev_iter, gmres_energy_list( step_iter ) * ev2Hartree, & 
                           gmres_energy_list( step_iter -1)* eV2Hartree, gmres_preconditioner, &
                           hpsi1, psi_pcdiv, ierr)
        if( ierr .ne. 0 ) return
!        call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
!        if( ierr .ne. 0 ) then
!          if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
!          return
!        endif
!        prev_energy = ener
        call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, hay_vec, ierr , nhflag)
        if( ierr .ne. 0 ) return
        do iter = 1, prev_iter
          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
        enddo
          !

        call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
        if( ierr .ne. 0 ) return
        if( loud ) then
          call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
          if( ierr .ne. 0 ) return
          call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          if( myid .eq. 0 ) then
            write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                  gval, gmres_convergence, ener, (1.0_dp - rval), -ival
          endif
        else
          if( myid .eq. 0 ) then
            write ( 6, '(1p,2x,3i5,3(1x,1e15.8),i8)' ) complete_iter, iter, gmres_depth, &
                  gval, gmres_convergence, ener, global_iter
          endif
        endif
!          if( myid.eq.root) write(*,*) 'conv check'
        if( gval .lt. gmres_convergence ) then
          goto 200 ! if convergered goto 200
        endif
        if( abs( prev_energy - ener ) > gmres_preconditioner/4.0_DP ) then
          if( myid .eq. root ) write(6,*) 'PRECONDITIONER OUT OF RANGE'
          prev_iter = 0
          prev_energy = ener
          call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
          if( ierr .ne. 0 ) then
            if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
            return
          endif
        endif
        !
      endif

      do ! outerloop for restarted GMRES
      
        do iter = prev_iter+1, gmres_depth
          complete_iter = complete_iter + 1
          global_iter = global_iter + 1
          psi_pg => u_matrix( iter )
          psi_apg => au_matrix( iter )
          ! pg = g * pcdiv
          if( do_precondition) then
            call OCEAN_psi_3element_mult( psi_pg, psi_g, psi_pcdiv, ierr )
          else
            call OCEAN_psi_copy_min( psi_pg, psi_g, ierr )
          endif
          if( ierr .ne. 0 ) then
            if( myid .eq. root ) write(6,*) 'do_precondition failed'
            return
          endif

          call OCEAN_psi_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return

          ! apg = H . pg
          call OCEAN_xact( sys, interaction_scale, psi_pg, psi_apg, ierr )
          if( ierr .ne. 0 ) return

          ! apg = (e+iG) * pg - apg
          call OCEAN_psi_axmy( psi_pg, psi_apg, ierr, ener, gmres_resolution )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_pg, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_apg, ierr )
          if( ierr .ne. 0 ) return


          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )

          if( iter .eq. gmres_depth ) then
            call OCEAN_psi_min2full( psi_x, ierr )
            if( ierr .ne. 0 ) return
            ! newi2loop section here
            call OCEAN_xact( sys, interaction_scale, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return

            ! psi_ax = (ener + iG ) * psi_x - psi_ax
            call OCEAN_psi_axmy( psi_x, psi_ax, ierr, ener, gmres_resolution )
            if( ierr .ne. 0 ) return
            !
            if( ener .ne. prev_energy ) then
              call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
              if( ierr .ne. 0 ) then
                if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
                return
              endif
              prev_energy = ener
            endif

            ! find approximate eigenvectors/values
!            keep = 5
!            call eigensystem_gmres( iter, keep, ierr )
!            prev_iter = keep
!            do i = 1, prev_iter
!              call update_gmres( i, psi_g, psi_x, psi_ax, ierr )
!              if( ierr .ne. 0 ) return
!              call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
!            enddo
!            if( myid .eq. 0 ) write( 6, * ) 'eigen'
            prev_iter = 0
          endif

          ! get new g
          ! g = ax - b   !!! y = x - z
          call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
          !

          call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return
          if( loud ) then
            call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
            if( ierr .ne. 0 ) return
            call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
            if( ierr .ne. 0 ) return
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener, (1.0_dp - rval), -ival
            endif
          else
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,3(1x,1e15.8),i8)' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener, global_iter
            endif
          endif
!          if( myid.eq.root) write(*,*) 'conv check'
          if( gval .lt. gmres_convergence ) then 
            prev_iter = iter
            if( iter .eq. gmres_depth ) prev_iter = 0 
            goto 200 ! if convergered goto 200
          endif
!          if(myid.eq.root) write(*,*) 'no conv: ', step_iter
        enddo


      enddo

!     Exit to here if we have converged
200   continue
!      prev_iter = iter

      call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
      if( ierr .ne. 0 ) return
      call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      if(myid.eq. root) then
        rel_error = -gval / ival
        write( abs_fh, '(1p,4(1e15.8,1x),1i6,x,i10)' ) ener * Hartree2eV, &
                    (1.0_dp - rval )*fact, -ival*fact, rel_error, complete_iter, global_iter
        flush(abs_fh)
      endif

      ! TODO: make this a separate flag
      if ( echamp ) then

        rval = 0.0_DP
        call OCEAN_psi_3element_mult( hay_psi_x, hay_vec, psi_x, ierr, rval, .true. )
        if( ierr .ne. 0 ) return
!JTV
!        call write_out_energy(sys, step_iter, hay_psi_x, nhflag, ierr)
!        if( ierr .ne. 0 ) return
      endif

      if( echamp ) then

        call write_out_echamp(sys, step_iter, psi_x, ierr)
        if( ierr .ne. 0 ) return
      endif
    enddo
!enddo 
   ! Create all the psi vectors we need

    call OCEAN_psi_kill( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_x, ierr )
       if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_g, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return


    if( myid .eq. root ) close( abs_fh )


  end subroutine OCEAN_gmres_do_recycle

  subroutine OCEAN_gmres_do( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_action, only : OCEAN_xact, OCEAN_action_h1
    use OCEAN_constants, only : eV2Hartree, Hartree2eV
    use OCEAN_mpi, only : myid, root, comm, MPI_STATUSES_IGNORE, MPI_REQUEST_NULL, MPI_INTEGER
    use OCEAN_energies, only : OCEAN_energies_allow
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
    type( ocean_vector ) :: psi_x, psi_ax, hpsi1, psi_g, psi_pcdiv
    !AK-NOTE:
    !! the strategy to output intermediate eigenvalues is:
    !! inefficient route: save as a new ocean_vector after the right combination
    !of hamiltonian elements have been applied and allow the minimization to
    !proceed as programmed.
    !! efficient route: use previously minimized eigenvalues as initial guesses
    !for the subsequent gmres of similar eigenvalues -- should provide a shorter
    !gmres time. 
    !! The elements of the BSE hamiltonian applied in OCEAN_gmres_do are: 
    type( ocean_vector), pointer :: psi_pg, psi_apg
    !
    real(DP) :: ener, rval, ival, gval, rel_error, fact
    real(dp), parameter :: one = 1.0_dp
    integer :: iter, step_iter, complete_iter, requests( 2 ), ierr_
    logical :: loud = .false.
    character( len=25 ) :: abs_filename
    integer, parameter :: abs_fh = 76
    type( ocean_vector ) :: hay_psi_x
    
    integer :: nhflag(6)
    nhflag(:) = sys%nhflag

    requests( : ) = MPI_REQUEST_NULL
    
    if( sys%cur_run%have_val ) then
      fact = hay_vec%kpref * 2.0_dp * sys%celvol
    else
      fact = hay_vec%kpref
    endif

    if( myid .eq. root ) then

      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref

      call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
      open( unit=abs_fh,file=abs_filename,form='formatted',status='unknown' )
      rewind( abs_fh )

    endif
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    
    ! Create all the psi vectors we need
    call OCEAN_psi_new( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_g, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    !!!!


    ! Create the psi_min_arrays we will need
    ! AX etc

    ! Keep around H * (1,0) for making pcdiv for each loop
    call OCEAN_psi_one_full( psi_x, ierr )
    if( ierr .ne. 0 ) return
    
    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, psi_x, ierr )
      if( ierr .ne. 0 ) return
    endif

    call OCEAN_xact( sys, interaction_scale, psi_x, hpsi1, ierr,.false. , nhflag)
    if( ierr .ne. 0 ) return




    do step_iter = 1, gmres_nsteps
      complete_iter = 0
      iter = 0

!      ener = e_list( iter ) * eV2Hartree
      ener = gmres_energy_list( step_iter ) * eV2Hartree
      if( myid .eq. root ) write(6,*) ener * Hartree2eV

      if( do_precondition ) then
        call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
        if( ierr .ne. 0 ) then
          if( myid .eq. root ) write(6,*) 'OCEAN_psi_min_set_prec failed'
          return
        endif
      endif
!      write(6,*) 'prec:', psi_pcdiv%min_r(1,1), psi_pcdiv%min_i(1,1)

!      zener = cmplx( ener, gres, DP )

      ! do we start with x = 0 or with x = previous or something else?
      call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, hay_vec, ierr , nhflag)
      if( ierr .ne. 0 ) then
        if( myid .eq. root ) write( 6,*) 'set_initial_vector failed'
        return
      endif

!      write(6,*) hay_vec%r(1,1,1), hay_vec%i(1,1,1)
!      write(6,*) psi_g%min_r(1,1), psi_g%min_i(1,1)

      if( loud ) then
        call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
        if( ierr .ne. 0 ) return

        call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
        if( ierr .ne. 0 ) return
        call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
        if( ierr .ne. 0 ) return
        if( myid .eq. 0 ) then
          write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                gval, gmres_convergence, ener, (1.0_dp - rval), -ival
        endif
      endif


      ! could have some out maximum here, like size_of_psi / gmres_depth 
      do ! outerloop for restarted GMRES
    ! if(myid .eq. root)  write(6, *) 'write out unmodified hay_vec with itter+200'
    !call write_out_echamp(sys, step_iter+200, hay_vec, ierr)
    !if( ierr .ne. 0 ) return
    
        do iter = 1, gmres_depth ! inner loop for restarted GMRES
          complete_iter = complete_iter + 1
          psi_pg => u_matrix( iter )
          psi_apg => au_matrix( iter )
          ! pg = g * pcdiv
          if( do_precondition) then
            call OCEAN_psi_3element_mult( psi_pg, psi_g, psi_pcdiv, ierr )
          else
            call OCEAN_psi_copy_min( psi_pg, psi_g, ierr )
          endif
          if( ierr .ne. 0 ) then
            if( myid .eq. root ) write(6,*) 'do_precondition failed'
            return
          endif

          call OCEAN_psi_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return
!          call OCEAN_psi_start_min2full( psi_pg, ierr )
!          if( ierr .ne. 0 ) return
!          call OCEAN_psi_finish_min2full( psi_pg, ierr )
!          if( ierr .ne. 0 ) return

          ! apg = H . pg
          call OCEAN_xact( sys, interaction_scale, psi_pg, psi_apg, ierr )
          if( ierr .ne. 0 ) return
        
          ! apg = (e+iG) * pg - apg
          call OCEAN_psi_axmy( psi_pg, psi_apg, ierr, ener, gmres_resolution )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_pg, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_free_fbe( psi_apg, ierr )
          if( ierr .ne. 0 ) return


          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )


          if( iter .eq. gmres_depth ) then
            call OCEAN_psi_min2full( psi_x, ierr )
            if( ierr .ne. 0 ) return
            ! newi2loop section here
            call OCEAN_xact( sys, interaction_scale, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return
            
            ! psi_ax = (ener + iG ) * psi_x - psi_ax
            call OCEAN_psi_axmy( psi_x, psi_ax, ierr, ener, gmres_resolution )
            if( ierr .ne. 0 ) return
            !
          endif

          ! get new g
          ! g = ax - b   !!! y = x - z
          call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
          !

          call OCEAN_psi_nrm( gval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return
          if( loud ) then
            call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
            if( ierr .ne. 0 ) return
            call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
            if( ierr .ne. 0 ) return
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener, (1.0_dp - rval), -ival
            endif
          else
            if( myid .eq. 0 ) then
              write ( 6, '(1p,2x,3i5,3(1x,1e15.8))' ) complete_iter, iter, gmres_depth, &
                    gval, gmres_convergence, ener
            endif
          endif
!          if( myid.eq.root) write(*,*) 'conv check'
          if( gval .lt. gmres_convergence ) goto 200 ! if convergered goto 200
!          if(myid.eq.root) write(*,*) 'no conv: ', step_iter
        enddo


      enddo

!     Exit to here if we have converged
200   continue

      call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
      if( ierr .ne. 0 ) return
      call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      if(myid.eq. root) then
        rel_error = -gval / ival
        write( abs_fh, '(1p,4(1e15.8,1x),1i6)' ) ener * Hartree2eV, & 
                    (1.0_dp - rval )*fact, -ival*fact, rel_error, complete_iter
        flush(abs_fh)
      endif

      ! TODO: make this a separate flag
      if ( echamp ) then
        ! Move create above
!        call OCEAN_psi_new( hay_psi_x, ierr, psi_x)
!        if( ierr .ne. 0 ) return
    
        rval = 0.0_DP
        call OCEAN_psi_3element_mult( hay_psi_x, hay_vec, psi_x, ierr, rval, .true. )
        if( ierr .ne. 0 ) return
!        call OCEAN_psi_dot_write(hay_vec,psi_x,hay_psi_x,requests(1),rval,ierr,requests(2),ival)
!        call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
!        if( ierr .ne. 0 ) return

!JTV
!        call write_out_energy(sys, step_iter, hay_psi_x, nhflag, ierr)
!        if( ierr .ne. 0 ) return
      endif

 
      if( echamp ) then
        
        call write_out_echamp(sys, step_iter, psi_x, ierr)
        if( ierr .ne. 0 ) return
      endif
    enddo    
!enddo 
   ! Create all the psi vectors we need
  
    call OCEAN_psi_kill( hay_psi_x, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_x, ierr )
       if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_g, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_kill( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return
    

    if( myid .eq. root ) close( abs_fh )
  
  end subroutine OCEAN_gmres_do

!> @brief Initializes X (and hence AX)
  subroutine set_initial_vector( sys, iter, psi_x, psi_g, psi_ax, hay_vec, ierr , nhflag)
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_zero_min, OCEAN_psi_axmz, OCEAN_psi_axmy, &
                          OCEAN_psi_min2full
    use OCEAN_action, only : OCEAN_xact
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : eV2Hartree
    use OCEAN_mpi, only : myid, root
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( in ) :: iter
    type( ocean_vector ), intent( inout ) :: psi_x, psi_g, psi_ax
    integer, intent( inout ) :: ierr
    !
    real(DP) :: ener
    character(len=4) :: initial = 'zero'
    integer, intent ( in ) :: nhflag(6)

    initial = 'zero'
    if( allow_reuse_x ) then
      ! option to reuse X
      if( iter .gt. 1 ) then
        if( ( gmres_energy_list( iter ) - gmres_energy_list( iter -1 ) ) * eV2Hartree & 
                .lt. 3.0_DP * gmres_resolution ) then
          initial = 'keep'
        endif
        if( force_reuse_x ) initial = 'keep'
      endif
    endif


    select case ( initial )

      case ( 'keep' )
        ! Keep X, generate psi_ax
        ! if( myid .eq. root ) write(6,*) '   Re-use previous X'
        call OCEAN_psi_min2full( psi_x, ierr )
        if( ierr .ne. 0 ) return
        call OCEAN_xact( sys, interaction_scale, psi_x, psi_ax, ierr, .false., nhflag )
        if( ierr .ne. 0 ) return
        ener = gmres_energy_list( iter ) * eV2Hartree
        call OCEAN_psi_axmy( psi_x, psi_ax, ierr, ener, gmres_resolution )

      case default 
        ! or 'zero'
        call OCEAN_psi_zero_min( psi_x, ierr )
        call OCEAN_psi_zero_min( psi_ax, ierr )

    end select
    !
    call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
    !
  end subroutine set_initial_vector

  subroutine write_out_energy( sys, current_iter, psi_x,hflag, ierr )
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_min2full, OCEAN_psi_vtor, OCEAN_psi_size_full
    use OCEAN_system, only : o_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_filenames, only : OCEAN_filenames_energy
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent( inout ) :: psi_x
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: out_vec( : )
    integer :: hflag(6)
    character(len=35) :: filename
    !
    call OCEAN_filenames_energy( sys, filename, current_iter,hflag, ierr )
    if( ierr .ne. 0 ) return
    !
   
    call OCEAN_psi_min2full( psi_x, ierr )
   !  write(6, *) 'FLAGAK 12 ', ierr
        if( ierr .ne. 0 ) return
    !
    if( myid .eq. root ) then
      allocate( out_vec( OCEAN_psi_size_full( psi_x ) ) )
      call OCEAN_psi_vtor( psi_x, out_vec ) 
      !
      open( unit=99, file=filename, form='formatted', status='unknown' )
      rewind( 99 )
      write(99, *) out_vec
      close( 99 )
      if( ierr .ne. 0 ) return
      deallocate( out_vec )
    endif
  end subroutine write_out_energy

  subroutine write_out_echamp( sys, current_iter, psi_x, ierr )
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_min2full, OCEAN_psi_vtor, OCEAN_psi_size_full
    use OCEAN_system, only : o_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_filenames, only : OCEAN_filenames_ehamp
    !
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent( inout ) :: psi_x
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: out_vec( : )
    character(len=28) :: filename
    !
    call OCEAN_filenames_ehamp( sys, filename, current_iter, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_min2full( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    if( myid .eq. root ) then
      allocate( out_vec( OCEAN_psi_size_full( psi_x ) ) )
      call OCEAN_psi_vtor( psi_x, out_vec ) 
      !
      open( unit=99, file=filename, form='unformatted', status='unknown' )
      rewind( 99 )
      write(99) out_vec
      close( 99 )
      deallocate( out_vec )
    endif
  end subroutine write_out_echamp

end module OCEAN_gmres
