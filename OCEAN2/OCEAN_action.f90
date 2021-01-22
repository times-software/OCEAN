! Copyright (C) 2015 - 2017, 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! Was formerly part of OCEAN_haydock.f90
!
module OCEAN_action
  use AI_kinds
  use OCEAN_timekeeper

  implicit none
  private


  public :: OCEAN_xact, OCEAN_action_h1

  contains


subroutine OCEAN_action_h1(sys, inter_scale, psi, new_psi, ierr, backwards, hflag) 
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_bubble, only : AI_bubble_act
    use OCEAN_ladder, only : OCEAN_ladder_act
    use OCEAN_constants, only : Hartree2eV

    implicit none
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: inter_scale
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent(inout) :: new_psi
    integer, intent(inout) :: ierr
    logical, optional :: backwards
    !
    type(OCEAN_vector) :: psi_o, psi_i
    integer :: rrequest, irequest
    real(dp) :: rval, ival
    logical :: loud_valence = .true.
    logical :: back
    integer, intent( in ) :: hflag(6) !a number 1-7 that determines combination of
!hamiltonian terms 
    if( present( backwards ) ) then
      back = backwards
    else
      back = .false.
    endif
!sys%cur_run%have_val = .true.

    call OCEAN_psi_zero_full( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'Zero full'

    call OCEAN_psi_ready_buffer( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'Ready buffer'

!    call OCEAN_tk_stop( tk_psisum )

    if( sys%cur_run%have_core ) then
       
!      if( sys%e0 .and. myid .eq. 0) then

      if( sys%mult ) then
     
        call OCEAN_tk_start( tk_mult )
              

call OCEAN_mult_act( sys, inter_scale, psi, new_psi, back, hflag ) 
        call OCEAN_tk_stop( tk_mult )
        endif
      if( sys%long_range ) then
        call OCEAN_tk_start( tk_lr )
        if(hflag(1).eq.1)then
         call lr_act( sys, psi, new_psi, ierr )
        endif      


call OCEAN_tk_stop( tk_lr )
      endif

    endif  
    if( sys%cur_run%have_val ) then       
      if( loud_valence ) then
        !JTV
        ! Loud valence hasn't been tested in a while!!
        if( sys%cur_run%have_core ) then
        
          ierr = 1
          if( myid .eq. root ) write(6,*) 'This code pathway is currently disabled'
          return
        endif
        call OCEAN_psi_new( psi_o, ierr, psi )
        call OCEAN_psi_new( psi_i, ierr )

        if( sys%cur_run%bande ) then
          ! Might as well go ahead and put these on the correct new_psi min to start
          call  OCEAN_tk_start( tk_e0 )
          call ocean_energies_act( sys, psi, new_psi, back, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_tk_stop( tk_e0 )

          ! Then just to be sure, use the allow (should've been done before psi was passed in)
          call OCEAN_energies_allow( sys, new_psi, ierr )      
          if(ierr .ne. 0) return

!          call OCEAN_psi_dot( psi_o, psi_i, rrequest, rval, ierr, irequest, ival )
          call OCEAN_psi_dot( psi_o, new_psi, rrequest, rval, ierr, irequest, ival )
          call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
          call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
          if(myid.eq.root) write(6,'(A6,4X,E22.15,1X,E22.15)') 'one-el',rval*Hartree2eV, ival*Hartree2eV
          rval = 1.0_dp

        else
          call OCEAN_psi_zero_min( new_psi, ierr )
          if( ierr .ne. 0 ) return
        endif

        if( sys%cur_run%bflag ) then
          call OCEAN_psi_zero_full( psi_i, ierr )
          call OCEAN_psi_ready_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_zero_min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          if(hflag(3).eq.1) then
          call AI_bubble_act( sys, psi, psi_i, ierr )
          
          endif
        if( ierr .ne. 0 ) return

          call OCEAN_psi_send_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_buffer2min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_energies_allow( sys, psi_i, ierr )
          if( ierr .ne. 0 ) return
        
          call OCEAN_psi_dot( psi_o, psi_i, rrequest, rval, ierr, irequest, ival )
          call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
          call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
          if( myid .eq. root ) write(6,'(A6,4X,E22.15,1X,E22.15)') 'bubble', rval*Hartree2eV, ival*Hartree2eV
          rval = 1.0_dp
          call OCEAN_psi_axpy( rval, psi_i, new_psi, ierr )


        endif

        if( sys%cur_run%lflag ) then
                
          call OCEAN_psi_zero_full( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_ready_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_zero_min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          if(hflag(2).eq.1) then
          call OCEAN_ladder_act( sys, psi, psi_i, ierr )
              
        endif
        if( ierr .ne. 0 ) return

          call OCEAN_psi_send_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_buffer2min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_energies_allow( sys, psi_i, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_psi_dot( psi_o, psi_i, rrequest, rval, ierr, irequest, ival )
          call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
          call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
          if( myid .eq. root ) write(6,'(A6,4X,E22.15,1X,E22.15)') 'ladder', rval*Hartree2eV, ival*Hartree2eV
          rval = 1.0_dp
          call OCEAN_psi_axpy( rval, psi_i, new_psi, ierr )


        endif

        ! clean up aux psi vectors
        call OCEAN_psi_kill( psi_o, ierr )
        call OCEAN_psi_kill( psi_i, ierr )

        ! Need to do a better job figuring out what is shared between core/valence
        ! and what is then repeated between loud/quiet valence
        return

      else  ! loud_valence = false
        ! Option 2 doesn't give per-BSE hamiltonian values for E0, direct, and
        ! exchange. This should be faster because less communication needed.
        ! Only share the psi vectors at the end like in the core case.
    
        if( sys%cur_run%bflag ) then
         
          ! For now re-use mult timing for bubble
          call OCEAN_tk_start( tk_mult )
         ! call AI_bubble_act( sys, psi, new_psi, ierr )


          if(hflag(3).eq.1) then
          call AI_bubble_act( sys, psi, psi_i, ierr )
          endif

!          call OCEAN_energies_allow( sys, new_psi, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_tk_stop( tk_mult )
        endif

        if( sys%cur_run%lflag ) then
        
          ! For now re-use lr timing for ladder
          call OCEAN_tk_start( tk_lr )

          if(hflag(2).eq.1) then
          call OCEAN_ladder_act( sys, psi, new_psi, ierr )
          endif
!          call OCEAN_ladder_act( sys, psi, new_psi, ierr )
          if( ierr .ne. 0 ) return
!          call OCEAN_energies_allow( sys, new_psi, ierr )
          call OCEAN_tk_stop( tk_lr )
        endif

        ! This should be redundant
!        call OCEAN_energies_allow( sys, new_psi, ierr )
!        if( ierr .ne. 0 ) return
      
      endif
    
    endif ! sys%cur_run%have_val 


    call OCEAN_tk_start( tk_psisum )
    call OCEAN_psi_send_buffer( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    call OCEAN_tk_stop( tk_psisum )

    ! end

    !JTV future if we are doing multiplets as a two-step process then their
    !results get saved down to the local/min while _send_buffer is working
!      if( sys%mult .and. sys%cur_run%have_core ) then
!        call OCEAN_mult_finish
!      else
!    call OCEAN_psi_zero_min( new_psi, ierr )
!     endif
!    if( ierr .ne. 0 ) return



!   After this step the min storage component of new_psi contains the one-electron energy component.
!   The contents of new_psi's min storage are overwritten.
!   Or if one-e is disabled, the min storage is zeroed out.
    if( ( sys%e0 .and. sys%cur_run%have_core ) .or. ( sys%cur_run%bande .and. sys%cur_run%have_val ) ) then
        call OCEAN_tk_start( tk_e0 )
      call ocean_energies_act( sys, psi, new_psi, back, ierr )
      call OCEAN_tk_stop( tk_e0 )
    else

      call OCEAN_psi_zero_min( new_psi, ierr )
      if( ierr .ne. 0 ) return
    endif

!   The data sitting in new_psi's buffer are added to the min 
    call OCEAN_psi_buffer2min( new_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_tk_stop( tk_psisum )


!   The min storage is multiplied by the allow matrix (this enforces the Fermi-Dirac occpations at 0K)
    call OCEAN_energies_allow( sys, new_psi, ierr )
    if( ierr .ne. 0 ) return

end subroutine OCEAN_action_h1

! On entrance psi needs to be the same everywhere
! On exit new_psi is stored in min everywhere
  subroutine OCEAN_xact( sys, inter_scale, psi, new_psi, ierr, backwards )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_bubble, only : AI_bubble_act
    use OCEAN_ladder, only : OCEAN_ladder_act
    use OCEAN_constants, only : Hartree2eV

    implicit none
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: inter_scale
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent(inout) :: new_psi
    integer, intent(inout) :: ierr
    logical, optional :: backwards
    !
    type(OCEAN_vector) :: psi_o, psi_i
    integer :: rrequest, irequest
    real(dp) :: rval, ival
    logical :: loud_valence = .false.
    logical :: back 
    integer :: hflag(6)

    hflag(:) = sys%nhflag(:)

    if( present( backwards ) ) then
      back = backwards
    else
      back = .false.
    endif
    call OCEAN_psi_zero_full( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'Zero full'

    call OCEAN_psi_ready_buffer( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'Ready buffer'

!    call OCEAN_tk_stop( tk_psisum )

    if( sys%cur_run%have_core ) then
       
!      if( sys%e0 .and. myid .eq. 0) then

      if( sys%mult ) then
        call OCEAN_tk_start( tk_mult )
        call OCEAN_mult_act( sys, inter_scale, psi, new_psi, back, hflag )
        call OCEAN_tk_stop( tk_mult )
        endif

      if( sys%long_range .and. (hflag(1).eq.1) ) then
        
        call OCEAN_tk_start( tk_lr )
       
        call lr_act( sys, psi, new_psi, ierr )
        call OCEAN_tk_stop( tk_lr )
      endif

    endif  ! sys%cur_run%have_core
     !write OCEAN_vec after mult after all of the if statements to hopefully get
     !only the one final %r out.  
    if( sys%cur_run%have_val ) then
      if( loud_valence ) then
        !JTV
        ! Loud valence hasn't been tested in a while!!
        if( sys%cur_run%have_core ) then
          ierr = 1
          if( myid .eq. root ) write(6,*) 'This code pathway is currently disabled'
          return
        endif
  !      call OCEAN_energies_allow( sys, psi, ierr )
  !      if( ierr .ne. 0 ) return
        call OCEAN_psi_new( psi_o, ierr, psi )
        call OCEAN_psi_new( psi_i, ierr )

        if( sys%cur_run%bande ) then
          ! Might as well go ahead and put these on the correct new_psi min to start
          call  OCEAN_tk_start( tk_e0 )
          call ocean_energies_act( sys, psi, new_psi, back, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_tk_stop( tk_e0 )

          ! Then just to be sure, use the allow (should've been done before psi was passed in)
          call OCEAN_energies_allow( sys, new_psi, ierr )
          if( ierr .ne. 0 ) return

          
          call OCEAN_psi_dot( psi_o, new_psi, rrequest, rval, ierr, irequest, ival )
          call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
          call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
          if(myid.eq.root) write(6,'(A6,4X,E22.15,1X,E22.15)') 'one-el',rval*Hartree2eV, ival*Hartree2eV
        !  rval = 1.0_dp
        
         rval = 1.0_dp

!          call OCEAN_psi_axpy( rval, psi_i, new_psi, ierr )
        else
          call OCEAN_psi_zero_min( new_psi, ierr )
          if( ierr .ne. 0 ) return
        endif

        if( sys%cur_run%bflag ) then
          call OCEAN_psi_zero_full( psi_i, ierr )
          call OCEAN_psi_ready_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_zero_min( psi_i, ierr )
          if( ierr .ne. 0 ) return


          if( hflag(3) .eq. 1 ) then
            call AI_bubble_act( sys, psi, psi_i, ierr )
            if( ierr .ne. 0 ) return
          endif

          call OCEAN_psi_send_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_buffer2min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_energies_allow( sys, psi_i, ierr )
          if( ierr .ne. 0 ) return
        
          call OCEAN_psi_dot( psi_o, psi_i, rrequest, rval, ierr, irequest, ival )
          call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
          call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
          if( myid .eq. root ) write(6,'(A6,4X,E22.15,1X,E22.15)') 'bubble', rval*Hartree2eV, ival*Hartree2eV
          rval = 1.0_dp
          call OCEAN_psi_axpy( rval, psi_i, new_psi, ierr )


        endif

        if( sys%cur_run%lflag ) then
          call OCEAN_psi_zero_full( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_ready_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_zero_min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          if(hflag(2).eq.1) then
            call OCEAN_ladder_act( sys, psi, psi_i, ierr )
            if( ierr .ne. 0 ) return
          endif

          call OCEAN_psi_send_buffer( psi_i, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_psi_buffer2min( psi_i, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_energies_allow( sys, psi_i, ierr )
          if( ierr .ne. 0 ) return

          call OCEAN_psi_dot( psi_o, psi_i, rrequest, rval, ierr, irequest, ival )
          call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
          call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
          if( myid .eq. root ) write(6,'(A6,4X,E22.15,1X,E22.15)') 'ladder', rval*Hartree2eV, ival*Hartree2eV
          rval = 1.0_dp
          call OCEAN_psi_axpy( rval, psi_i, new_psi, ierr )


        endif

        ! clean up aux psi vectors
        call OCEAN_psi_kill( psi_o, ierr )
        call OCEAN_psi_kill( psi_i, ierr )

        ! Need to do a better job figuring out what is shared between core/valence
        ! and what is then repeated between loud/quiet valence
        return

      else  ! loud_valence = false
        ! Option 2 doesn't give per-BSE hamiltonian values for E0, direct, and
        ! exchange. This should be faster because less communication needed.
        ! Only share the psi vectors at the end like in the core case.
    
        if( sys%cur_run%bflag .and. (hflag(3).eq.1) ) then
          ! For now re-use mult timing for bubble
          call OCEAN_tk_start( tk_mult )
          call AI_bubble_act( sys, psi, new_psi, ierr )
!          call OCEAN_energies_allow( sys, new_psi, ierr )
          if( ierr .ne. 0 ) return
          call OCEAN_tk_stop( tk_mult )
        endif

        if( sys%cur_run%lflag .and. (flag(2).eq.1)) then
          ! For now re-use lr timing for ladder
          call OCEAN_tk_start( tk_lr )
          call OCEAN_ladder_act( sys, psi, new_psi, ierr )
          if( ierr .ne. 0 ) return
!          call OCEAN_energies_allow( sys, new_psi, ierr )
          call OCEAN_tk_stop( tk_lr )
        endif

        ! This should be redundant
!        call OCEAN_energies_allow( sys, new_psi, ierr )
!        if( ierr .ne. 0 ) return
      
      endif
    
    endif ! sys%cur_run%have_val 


    call OCEAN_tk_start( tk_psisum )
    call OCEAN_psi_send_buffer( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    call OCEAN_tk_stop( tk_psisum )

    ! end

    !JTV future if we are doing multiplets as a two-step process then their
    !results get saved down to the local/min while _send_buffer is working
!      if( sys%mult .and. sys%cur_run%have_core ) then
!        call OCEAN_mult_finish
!      else
!    call OCEAN_psi_zero_min( new_psi, ierr )
!     endif
!    if( ierr .ne. 0 ) return



!   After this step the min storage component of new_psi contains the one-electron energy component.
!   The contents of new_psi's min storage are overwritten.
!   Or if one-e is disabled, the min storage is zeroed out.
    if( ( sys%e0 .and. sys%cur_run%have_core ) .or. ( sys%cur_run%bande .and. sys%cur_run%have_val ) ) then
     if((sys%e0 .and. sys%cur_run%have_core))then
        endif
     if((sys%cur_run%bande .and. sys%cur_run%have_val))then
        
        endif
        call OCEAN_tk_start( tk_e0 )
      call ocean_energies_act( sys, psi, new_psi, back, ierr )
      call OCEAN_tk_stop( tk_e0 )
    else

      call OCEAN_psi_zero_min( new_psi, ierr )
      if( ierr .ne. 0 ) return
    endif

!   The data sitting in new_psi's buffer are added to the min 
    call OCEAN_psi_buffer2min( new_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_tk_stop( tk_psisum )


!   The min storage is multiplied by the allow matrix (this enforces the Fermi-Dirac occpations at 0K)
    call OCEAN_energies_allow( sys, new_psi, ierr )
    if( ierr .ne. 0 ) return

  end subroutine OCEAN_xact

end module OCEAN_action
