subroutine OCEAN_load_data( sys, ierr )
  use OCEAN_system
  use OCEAN_mpi

  implicit none
  integer, intent( inout ) :: ierr
  type( ocean_system ), intent( in ) :: sys


  if( myid .eq. root ) write(6,*) 'Calc Type = ', sys%calc_type

  if( myid .eq. root ) write(6,*) 'Init matrix elements'
  call ocean_psi_init( sys, ierr )
  if( myid .eq. root ) write(6,*) 'Load matrix elements'
  call ocean_psi_load( sys, ierr )
  if( myid .eq. root ) write(6,*) 'Matrix elements loaded'

  if( sys%e0 ) then
    if( myid .eq. root ) write(6,*) 'Init energies'
    call ocean_energies_init( sys, ierr )
    if( myid .eq. root ) write(6,*) 'Load energies'
    call ocean_energies_load( sys, ierr )
    if( myid .eq. root ) write(6,*) 'Energies loaded'
  endif

  if( sys%mult ) then
    if( myid .eq. root ) write(6,*) 'Init mult'
    if( myid .eq. root ) write(6,*) 'Load mult'
    if( myid .eq. root ) write(6,*) 'Mult loaded'
    
  endif

  if( sys%long_range ) then
    if( myid .eq. root ) write(6,*) 'Init long_range'
    if( myid .eq. root ) write(6,*) 'Load long_range'
    if( myid .eq. root ) write(6,*) 'Long_range loaded'
  endif

  if( sys%obf ) then
    if( myid .eq. root ) write(6,*) 'Init OBFs'
    if( myid .eq. root ) write(6,*) 'Load OBFs'
    if( myid .eq. root ) write(6,*) 'OBFs loaded'
  endif


  if( myid .eq. root ) write(6,*) 'Initialization complete'


end subroutine OCEAN_load_data
