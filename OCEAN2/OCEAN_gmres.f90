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

  public :: OCEAN_gmres_do, OCEAN_gmres_setup, OCEAN_gmres_clean

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
    call MPI_BCAST( gmres_preconditioner, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( gmres_nsteps, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( myid .ne. root ) allocate( gmres_energy_list( gmres_nsteps ) )
    call MPI_BCAST( gmres_energy_list, gmres_nsteps, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'Finished gmres setup'
    call initialize_gmres_storage( sys, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'Finished gmres setup'

  end subroutine OCEAN_gmres_setup

  subroutine initialize_gmres_storage( sys, ierr )
    use OCEAN_psi, only : OCEAN_psi_init, OCEAN_psi_new
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
      allocate( c_matrix( gmres_depth*(gmres_depth+1)/2 ) )
    else
      allocate( c_matrix( 1 ) )
    endif

  end subroutine


  subroutine update_gmres( current_iter, psi_g, psi_x, psi_ax, ierr)
    use OCEAN_mpi
    use OCEAN_psi, only : ocean_vector, OCEAN_psi_dot, OCEAN_psi_axpy, OCEAN_psi_nrm
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
    allocate( re_coeff_vec( current_iter ), im_coeff_vec( current_iter ), & 
              re_rvec( current_iter ), im_rvec( current_iter ) )    

    allocate( re_coeff_request( current_iter ), im_coeff_request( current_iter ), &
              re_rvec_request( current_iter ), im_rvec_request( current_iter ) )
    ! 
    im_coeff_request( current_iter ) = MPI_REQUEST_NULL
    im_coeff_vec( current_iter ) = 0.0_DP
    !
    do iter = 1, current_iter 
      call OCEAN_psi_dot( au_matrix( current_iter ), au_matrix( iter ), &
                          re_coeff_request( iter ), re_coeff_vec( iter ), ierr, & 
                          im_coeff_request( iter ), im_coeff_vec( iter ) )

      call OCEAN_psi_dot( au_matrix( iter ), psi_g, &
                          re_rvec_request( iter ), re_rvec( iter ), ierr , &
                          im_rvec_request( iter ), im_rvec( iter ) )
    enddo

!    call OCEAN_psi_nrm( re_coeff_vec( current_iter ), au_matrix( current_iter ), ierr, & 
!                        re_coeff_request( current_iter ) )
!    call OCEAN_psi_dot( au_matrix( current_iter ), psi_g, &
!                        re_rvec_request( current_iter ), re_rvec( current_iter ), ierr , &
!                        im_rvec_request( current_iter ), im_rvec( current_iter ) )

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
    

  end subroutine update_gmres


  subroutine OCEAN_gmres_do( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_action, only : OCEAN_xact
    use OCEAN_constants, only : eV2Hartree, Hartree2eV
    use OCEAN_mpi, only : myid, root, MPI_STATUSES_IGNORE, MPI_REQUEST_NULL
    use OCEAN_energies, only : OCEAN_energies_val_allow
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
    type( ocean_vector ) :: psi_x, psi_ax, hpsi1, psi_g, psi_pcdiv
    type( ocean_vector), pointer :: psi_pg, psi_apg
    !
    real(DP) :: ener, rval, ival, gval, rel_error, fact
    real(dp), parameter :: one = 1.0_dp
    integer :: iter, step_iter, complete_iter, requests( 2 )
    logical :: loud = .true.

    requests( : ) = MPI_REQUEST_NULL

    if( sys%cur_run%have_val ) then
      fact = hay_vec%kpref * 2.0_dp * sys%celvol
    else
      fact = hay_vec%kpref
    endif

    if( myid .eq. root ) then
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref
!      write(6,*) hay_vec%r(1,1,1), hay_vec%i(1,1,1)
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
    !!!!


    ! Create the psi_min_arrays we will need
    ! AX etc

    ! Keep around H * (1,0) for making pcdiv for each loop
    call OCEAN_psi_one_full( psi_x, ierr )
    if( ierr .ne. 0 ) return
    
    if( sys%cur_run%have_val ) then
      call OCEAN_energies_val_allow( sys, psi_x, ierr )
      if( ierr .ne. 0 ) return
    endif

    call OCEAN_xact( sys, one, psi_x, hpsi1, ierr )
    !!!!




    do step_iter = 1, gmres_nsteps
      complete_iter = 0

!      ener = e_list( iter ) * eV2Hartree
      ener = gmres_energy_list( step_iter ) * eV2Hartree
      if( myid .eq. root ) write(6,*) ener * Hartree2eV

      call OCEAN_psi_min_set_prec( ener, gmres_preconditioner, hpsi1, psi_pcdiv, ierr )
      if( ierr .ne. 0 ) return
!      write(6,*) 'prec:', psi_pcdiv%min_r(1,1), psi_pcdiv%min_i(1,1)

!      zener = cmplx( ener, gres, DP )

      ! do we start with x = 0 or with x = previous or something else?
      call set_initial_vector( sys, step_iter, psi_x, psi_g, psi_ax, hay_vec, ierr )
      if( ierr .ne. 0 ) return

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

        do iter = 1, gmres_depth ! inner loop for restarted GMRES
          complete_iter = complete_iter + 1
          psi_pg => u_matrix( iter )
          psi_apg => au_matrix( iter )
          ! pg = g * pcdiv
          call OCEAN_psi_element_mult( psi_pg, psi_g, psi_pcdiv, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_psi_prep_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_start_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_finish_min2full( psi_pg, ierr )
          if( ierr .ne. 0 ) return

          ! apg = H . pg
          call OCEAN_xact( sys, one, psi_pg, psi_apg, ierr )
          if( ierr .ne. 0 ) return
        
          ! apg = (e+iG) * pg - apg
          call OCEAN_psi_axmy( psi_pg, psi_apg, ierr, ener, gmres_resolution )
          if( ierr .ne. 0 ) return


          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )

          if( iter .eq. gmres_depth ) then
            ! newi2loop section here
            call OCEAN_xact( sys, one, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return
            !
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

          if( gval .lt. gmres_convergence ) goto 200 ! if convergered goto 200
  
        enddo


      enddo

!     Exit to here if we have converged
200   continue

      call OCEAN_psi_dot( hay_vec, psi_x, requests(1), rval, ierr, requests(2), ival )
      if( ierr .ne. 0 ) return
      call MPI_WAITALL( 2, requests, MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
  
      if( myid .eq. root ) then
        rel_error = -gval / ival
        write( 76, '(1p,1i5,4(1x,1e15.8))' ) complete_iter, ener * Hartree2eV, & 
                    (1.0_dp - rval )*fact, -ival*fact, rel_error
      endif
      

    enddo    

    ! Create all the psi vectors we need
    call OCEAN_psi_kill( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_kill( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_kill( hpsi1, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_kill( psi_g, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_kill( psi_pcdiv, ierr )
    if( ierr .ne. 0 ) return



  end subroutine OCEAN_gmres_do

!> @brief Initializes X (and hence AX)
  subroutine set_initial_vector( sys, iter, psi_x, psi_g, psi_ax, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    integer, intent( in ) :: iter
    type( ocean_vector ), intent( inout ) :: psi_x, psi_g, psi_ax
    integer, intent( inout ) :: ierr
    !
    if( .true. ) then
      call OCEAN_psi_zero_min( psi_x, ierr )
      call OCEAN_psi_zero_min( psi_ax, ierr )
    else
      ! Set psi_x to be something
!      call OCEAN_xact( sys, psi_x, psi_ax, ierr )
!      call OCEAN_psi_axmy( psi_x, psi_ax, ierr, ener, gmres_resolution )
    endif
    !
    call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
    !
  end subroutine set_initial_vector




end module OCEAN_gmres
