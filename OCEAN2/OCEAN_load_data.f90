! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine OCEAN_load_data( sys, hay_vec, lr, ierr )
  use OCEAN_system
  use OCEAN_psi
  use OCEAN_energies
  use OCEAN_mpi
  use OCEAN_multiplet
  use OCEAN_long_range

  implicit none
  integer, intent( inout ) :: ierr
  type( o_system ), intent( in ) :: sys
  type(ocean_vector), intent( inout ) :: hay_vec
  type( long_range ), intent(out) :: lr

#ifdef MPI
!    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
#endif


  if( myid .eq. root ) write(6,*) 'Calc Type = ', sys%calc_type
  if( myid .eq. root ) write(6,*) 'Init matrix elements 1'

  call ocean_psi_init( sys, ierr )
  if( ierr .ne. 0 ) return

  if( myid .eq. root ) write(6,*) 'Init matrix elements 2'
  call ocean_psi_new( hay_vec, ierr )
  if( ierr .ne. 0 ) return

  if( myid .eq. root ) write(6,*) 'Load matrix elements'
  call ocean_psi_load( sys, hay_vec, ierr )
  if( ierr .ne. 0 ) return

  if( myid .eq. root ) write(6,*) 'Matrix elements loaded'
  call ocean_psi_write( sys, hay_vec, ierr )
  if( ierr .ne. 0 ) return


  if( sys%cur_run%have_val) then
    if( sys%cur_run%bande ) then
      if( myid .eq. root ) write(6,*) 'Init energies'
      call OCEAN_energies_val_load( sys, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( sys%cur_run%bflag ) then
      if( myid .eq. root ) write(6,*) 'Init bubble'
!      call AI_bubble_prep( sys, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Bubble prepped'
    endif

    if( sys%cur_run%lflag ) then
!      if( myid .eq. root ) write(6,*) 'Init ladder'
!      call

    endif

  endif 
  

  if( sys%cur_run%have_core ) then

    if( sys%e0 ) then
      if( myid .eq. root ) write(6,*) 'Init energies'
      call ocean_energies_init( sys, ierr )
      if( myid .eq. root ) write(6,*) 'Load energies'
      call ocean_energies_load( sys, ierr )
      if( myid .eq. root ) write(6,*) 'Energies loaded'
    endif

    if( sys%mult ) then
      if( myid .eq. root ) write(6,*) 'Init mult'
      call OCEAN_create_central( sys, ierr )
      if( myid .eq. root ) write(6,*) 'Mult loaded'
      
    endif

    if( sys%long_range ) then
      if( myid .eq. root ) write(6,*) 'Init long_range'
      call lr_init( sys, lr, ierr )
      if( ierr .ne. 0 ) return
      if( myid .eq. root ) write(6,*) 'Long_range loaded'
    endif

  endif

  if( sys%obf ) then
    if( myid .eq. root ) write(6,*) 'Init OBFs'
    if( myid .eq. root ) write(6,*) 'Load OBFs'
    if( myid .eq. root ) write(6,*) 'OBFs loaded'
  endif



  if( myid .eq. root ) write(6,*) 'Initialization complete'


end subroutine OCEAN_load_data
