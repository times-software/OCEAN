! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine OCEAN_load_data( sys, hay_vec, ierr )
  use OCEAN_system
  use OCEAN_psi
  use OCEAN_energies
  use OCEAN_mpi
  use OCEAN_multiplet
  use OCEAN_long_range
  use OCEAN_val_states, only : OCEAN_val_states_load, OCEAN_val_states_init
  use OCEAN_bubble, only : AI_bubble_prep
  use OCEAN_ladder, only : OCEAN_ladder_init, OCEAN_ladder_new
  use OCEAN_fxc, only : OCEAN_fxc_prep

  implicit none
  integer, intent( inout ) :: ierr
  type( o_system ), intent( inout ) :: sys
  type(ocean_vector), intent( inout ) :: hay_vec

  logical :: complex_bse

  complex_bse = .false.

#ifdef MPI
!    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
#endif

  if( myid .eq. root ) write(6,*) 'Calc Type = ', sys%cur_run%calc_type
  if( myid .eq. root ) write(6,*) 'Init matrix elements 1'

  call ocean_psi_init( sys, ierr )
  if( ierr .ne. 0 ) return

  if( myid .eq. root ) write(6,*) 'Init matrix elements 2'
  call ocean_psi_new( hay_vec, ierr )
  if( ierr .ne. 0 ) return

  if( myid .eq. root ) then 
    if( sys%cur_run%have_val )write(6,*) 'Load matrix elements for valence calculation'
    if( sys%cur_run%have_core )write(6,*) 'Load matrix elements for core calculation'
  endif
  call ocean_psi_load( sys, hay_vec, ierr )
  if( ierr .ne. 0 ) return

  if( myid .eq. root ) write(6,*) 'Matrix elements loaded'
  if( sys%write_rhs ) then
    call ocean_psi_write( sys, hay_vec, 'rhs_', .true., ierr )
    if( ierr .ne. 0 ) return
  endif


  if( sys%cur_run%have_val) then
    ! Energy is how we construct the allow and sfact arrays
    !  these are needed even if (for some reason) we did not want to run energies
    if( myid .eq. root ) write(6,*) 'Init energies'
    call OCEAN_energies_val_load( sys, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) write(6,*) 'Trim & scale matrix elements'
    ! Now trim the hay_vec by the allow array 
    !  This 1) cuts off over-lapped states valence above Fermi/conduction below
    !       2) Uniform energy cutoff for upper bands
    call OCEAN_energies_allow( sys, hay_vec, ierr, sfact=sys%bwflg )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_pnorm( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) return
    if( myid .eq. root ) write(6,*) 'Trim & scale complete'


    if( sys%cur_run%bflag .or. sys%cur_run%lflag .or. sys%cur_run%aldaf ) then
      if( myid .eq. root ) write(6,*) 'Init val states'
      call OCEAN_val_states_init( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Load val states'
      call OCEAN_val_states_load( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Val states loaded'

    endif

    if( sys%cur_run%bflag ) then
      if( myid .eq. root ) write(6,*) 'Init bubble'
      call AI_bubble_prep( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Bubble prepped'
    endif

    if( sys%cur_run%lflag ) then
      if( myid .eq. root ) write(6,*) 'Init ladder'
      call OCEAN_ladder_init( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Load ladder'
      call OCEAN_ladder_new( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Ladder loaded'
    endif

    if( sys%cur_run%aldaf ) then
      if( myid .eq. root ) write(6,*) 'Init ALDA'
      call OCEAN_fxc_prep( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'ALDA initialized'
    endif

  endif 
  

  if( sys%cur_run%have_core ) then

!    if( sys%e0 ) then
      if( myid .eq. root ) write(6,*) 'Init energies'
      call ocean_energies_init( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Load energies'
      call ocean_energies_load( sys, complex_bse, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Energies loaded'
 
    if( myid .eq. root ) write(6,*) 'Trim & scale matrix elements'
    call OCEAN_energies_initial_allow( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_pnorm( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) return
!    endif
!    call OCEAN_energies_resetAllow( ierr )
!    if( ierr .ne. 0 ) return

    if( sys%mult ) then
      if( myid .eq. root ) write(6,*) 'Init mult'
      call OCEAN_create_central( sys, complex_bse, ierr )
      if( myid .eq. root ) write(6,*) 'Mult loaded'
      
    endif

    if( sys%long_range ) then
      if( myid .eq. root ) write(6,*) 'Init long_range'
      call lr_init( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Long_range loaded'
    endif

  endif

  if( sys%obf ) then
    if( myid .eq. root ) write(6,*) 'Init OBFs'
    if( myid .eq. root ) write(6,*) 'Load OBFs'
    if( myid .eq. root ) write(6,*) 'OBFs loaded'
  endif


  sys%cur_run%complex_bse = complex_bse

  if( myid .eq. root ) write(6,*) 'Initialization complete', ierr
  if( sys%cur_run%complex_bse ) then
    if( myid .eq. root ) write(6,*) '  Will use bi-lanczos for non-Hermitian Hamiltonian'
  endif
#ifdef MPI
!    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
#endif


end subroutine OCEAN_load_data

subroutine OCEAN_reload_val( sys, hay_vec, ierr )
  use OCEAN_system
  use OCEAN_psi
  use OCEAN_energies
  use OCEAN_mpi
  use OCEAN_multiplet
  use OCEAN_long_range
  use OCEAN_val_states, only : OCEAN_val_states_load, OCEAN_val_states_init
  use OCEAN_bubble, only : AI_bubble_prep
  use OCEAN_ladder, only : OCEAN_ladder_init, OCEAN_ladder_new, OCEAN_ladder_kill

  implicit none
  integer, intent( inout ) :: ierr
  type( o_system ), intent( inout ) :: sys
  type(ocean_vector), intent( inout ) :: hay_vec

  if( sys%cur_run%have_val) then

#ifdef MPI
!    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
#endif

    call ocean_psi_load( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) write(6,*) 'Trim & scale matrix elements'
    ! Now trim the hay_vec by the allow array 
    !  This 1) cuts off over-lapped states valence above Fermi/conduction below
    !       2) Uniform energy cutoff for upper bands
    call OCEAN_energies_allow( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_pnorm( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) return
    if( myid .eq. root ) write(6,*) 'Trim & scale complete'
#ifdef MPI
!    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
#endif
    
    if( sys%cur_run%lflag ) then
      call OCEAN_ladder_kill()
      if( myid .eq. root ) write(6,*) 'Init ladder'
      call OCEAN_ladder_init( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Load ladder'
      call OCEAN_ladder_new( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Ladder loaded'

    endif


  else
    return
  endif

end subroutine OCEAN_reload_val
