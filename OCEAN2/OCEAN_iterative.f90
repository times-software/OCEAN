module OCEAN_iterative

  use AI_kinds

  implicit none
  private
  save
    
  type(ocean_vector), allocatable :: au_matrix(:)
  type(ocean_vector), allocatable :: u_matrix(:)
  complex(DP) :: c_matrix( : )

!  integer :: 
  integer :: gmres_depth ! max size of gmres Arnoldi? space before restart

  public :: OCEAN_do_gmres

  contains

  subroutine initialize_gmres_storage( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
    integer :: local_gmres_size
    !
    local_gmres_size = OCEAN_psi_size_min( hay_vec )

    allocate( au_matrix( gmres_depth ) )
    allocate( u_matrix( gmres_depth ) )

    do iter = 1, gmres_depth
      call OCEAN_psi_new( au_matrix( iter ), ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_new( u_matrix( iter ), ierr )
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
    !
    integer, intent( in ) :: current_iter
    type( ocean_vector ), intent(inout ) :: psi_g, psi_x, psi_ax
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: c_temp( : ), coeff( : )
    integer, allocatable :: ipiv( : )
    integer :: local_gmres_size, info
    !
    
    allocate( re_coeff_request( current_iter ), im_coeff_request( current_iter ), &
              re_rvec_request( current_iter ), im_rvec_request( current_iter ) )
    !
    do iter = 1, current_iter 
      call OCEAN_psi_dot( au_matrix( current_iter ), au_matrix( iter ), &
                          re_coeff_request( iter ), re_coeff_vec( iter ), ierr, & 
                          im_coeff_request( iter ), im_coeff_vec( iter ) )

      call OCEAN_psi_dot( au_matrix( iter ), psi_g, &
                          re_rvec_request( iter ), re_rvec( iter ), ierr , &
                          im_rvec_request( iter ), im_rvec( iter ) )
    enddo

    call MPI_WAITALL( current_iter, re_coeff_request, MPI_STATUSES_IGNORE, ierr )
    call MPI_WAITALL( current_iter, im_coeff_request, MPI_STATUSES_IGNORE, ierr )
    !
    allocate( coeff( current_iter ) )
    if( myid .eq. root ) then
      iter_start = (current_iter * ( current_iter - 1 ) ) / 2
      do iter = 1, current_iter
        c_matrix( iter + iter_start ) = cmplx( re_coeff_vec( iter ), im_coeff_vec( iter ), DP )
      enddo
      
      c_size = current_iter * (current_iter + 1 ) / 2 ) 
      allocate( c_temp( c_size ), ipiv( current_iter ) )
      c_temp( : ) = c_matrix( 1 : c_size )

      call zhptrf( 'U', current_iter, c_temp, ipiv, info )

      call MPI_WAITALL( current_iter, re_rvec_request, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( current_iter, im_rvec_request, MPI_STATUSES_IGNORE, ierr )

      ! Want to test after waitall so that all procs are on the same page if we abort
      if( info .ne. 0 ) then
        ierr = info
      else
        do iter = 1, current_iter
          coeff( iter ) = -cmplx( re_coeff_vec( iter ), im_coeff_vec( iter ), DP )
        enddo

        call zhptrs( 'U', current_iter, 1, c_temp, ipiv, coeff, current_iter, info )
      endif
      deallocate( c_temp, ipiv )

    else
      call MPI_WAITALL( current_iter, re_rvec_request, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( current_iter, im_rvec_request, MPI_STATUSES_IGNORE, ierr )
    endif

    call BCAST( info, 1, MPI_INTEGER, root, comm, ierr )
    ierr = info
    if( ierr .ne. 0 ) return

    call MPI_BCAST( coeff, current_iter, MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. 0 ) return

    ! update ax and x
    do iter = 1, current_iter
      tmp_r = real( coeff( iter ), DP )
      tmp_i = aimag( coeff( iter ) )
      call OCEAN_psi_axpy( tmp_r, u_matrix( iter ), psi_x, ierr, tmp_i )
      call OCEAN_psi_axpy( tmp_r, au_matrix( iter ), psi_ax, ierr, tmp_i )
    enddo

    deallocate( coeff )
    

  end subroutine update_gmres


  subroutine OCEAN_do_gmres( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_haydock, only : OCEAN_xact
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
    type( ocean_vector ) :: psi_x, psi_ax, hpsi1, psi_g
    type( ocean_vector), pointer :: psi_pg, psi_apg
  

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

    call OCEAN_xact( sys, psi_x, hpsi1, ierr )
    !!!!



    do iter = 1, inv_loop

      ener = e_list( iter ) * eV2Hartree
      if( myid .eq. root ) write(6,*) ener * Hartree2eV

      call OCEAN_psi_min_set_prec( ener, gprc, hpsi1, psi_pcdiv, ierr )
      if( ierr .ne. 0 ) return

!      zener = cmplx( ener, gres, DP )

      ! do we start with x = 0 or with x = previous or something else?
      call set_initial_vector( sys, iter, psi_x, psi_g, psi_ax, hay_vec )


      ! could have some out maximum here, like size_of_psi / gmres_depth 
      do ! outerloop for restarted GMRES

        do iter = 1, gmres_depth ! inner loop for restarted GMRES
          psi_pg => u_matrix( iter )
          psi_apg => au_matrix( iter )
          ! pg = g * pcdiv
          call OCEAN_psi_element_mult( psi_pg, psi_g, psi_pcdiv, 0.0_dp, ierr )
          if( ierr .ne. 0 ) return

          ! apg = H . pg
          call OCEAN_xact( sys, psi_pg, psi_apg, ierr )
          if( ierr .ne. 0 ) return
        
          ! apg = (e-iG) * pg - apg
          call OCEAN_psi_axmy( ener, psi_pg, psi_apg, ierr, gres )
          if( ierr .ne. 0 ) return


          call update_gmres( iter, psi_g, psi_x, psi_ax, ierr )

          if( iter .eq. gmres_depth ) then
            ! newi2loop section here
            call OCEAN_xact( sys, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return
            !
            ! psi_ax = (ener + iG ) * psi_x - psi_ax
            call OCEAN_psi_axmy( ener, psi_x, psi_ax, ierr, gres )
            if( ierr .ne. 0 ) return
            !
          endif

          ! get new g
          ! g = ax - b   !!! y = x - z
          call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
          !

          call OCEAN_psi_nrm( rval, psi_g, ierr )  ! non-blocking wouldn't do any good
          if( ierr .ne. 0 ) return
          goto 200 if( rval .lt. ffff ) ! if convergered goto 200
  
        enddo


      enddo

!     Exit to here if we have converged
200   continue
      

    enddo    


  end subroutine OCEAN_do_gmres

!> @brief Initializes X (and hence AX)
  subroutine set_initial_vector( sys, iter, psi_x, psi_g, psi_ax, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( in ) :: iter
    type( ocean_vector ), intent( inout ) :: psi_x, psi_g, psi_ax
    integer, intent( inout ) :: ierr
    !
    if( .true. ) then
      call OCEAN_psi_zero_min( psi_x, ierr )
      call OCEAN_psi_zero_min( psi_ax, ierr )
    else
      ! Set psi_x to be something
      call OCEAN_xact( sys, psi_x, psi_ax, ierr )
      call OCEAN_psi_axmy( ener, psi_x, psi_ax, ierr, gres )
    endif
    !
    call OCEAN_psi_axmz( psi_ax, psi_g, hay_vec, ierr )
    !
  end subroutine set_initial_vector




end module OCEAN_iterative
