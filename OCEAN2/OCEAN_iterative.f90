module OCEAN_iterative

  use AI_kinds

  implicit none
  private
  save
    
  public :: OCEAN_do_gmres

  contains


  subroutine OCEAN_do_gmres( sys, hay_vec, ierr )
    use OCEAN_psi
    use OCEAN_system
    use OCEAN_haydock, only : OCEAN_xact
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr
    !
  

    ! Create all the psi vectors we need
    call OCEAN_psi_new( psi_x, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( psi_ax, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_new( hspi1, ierr )
    if( ierr .ne. 0 ) return
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

      zener = cmplx( ener, gres, DP )

      ! do we start with x = 0 or with x = previous or something else?
      call set_initial_vector( sys, iter, psi_x, psi_g, psi_ax, hay_vec )


      ! could have some out maximum here, like size_of_psi / nloop 
      do ! outerloop for restarted GMRES

        do iter = 1, nloop-1 ! inner loop for restarted GMRES
          ! pg = g * pcdiv
          call OCEAN_psi_element_multi( psi_pg, psi_g, psi_pcdiv, 0.0_dp, ierr )
          if( ierr .ne. 0 ) return

          ! apg = H . pg
          call OCEAN_xact( sys, psi_pg, psi_apg, ierr )
          if( ierr .ne. 0 ) return
        
          ! apg = (e-iG) * pg - apg
          call OCEAN_psi_zaxmy( zener, psi_pg, psi_apg, ierr )


          call update_gmres( )

          if( iter .eq. nloop ) then
            ! newi2loop section here
            call OCEAN_xact( sys, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return
            !
            call OCEAN_psi_zacmy( zener, psi_x, psi_ax, ierr )
            if( ierr .ne. 0 ) return
            !
          endif

          ! get new g
          ! g = ax - b
          call OCEAN_psi_axmy( psi_ax, hay_vec, psi_g )
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

end module OCEAN_iterative
