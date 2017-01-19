! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine qe_grabwf(ikpt, isppol, nsppol, maxband, maxnpw, kg, eigen, occ, cg, npw, ierr )
  use iotk_module !, only : 
  implicit none
  integer, intent( in ) :: ikpt, isppol, nsppol, maxband, maxnpw
  integer, intent( out ) :: npw
  integer, intent( inout ) :: ierr
  real(kind=kind(1.0d0)), intent( out ) :: eigen(maxband), occ(maxband), cg(maxband,2*maxnpw)
  integer, intent( out ) :: kg( 3, maxnpw )
  character(len=22) :: dirname
  character(len=16) :: prefix = 'Out/system.save/' 
  character(len=128) :: filename
  integer :: i, j
  complex(kind=kind(1.0d0)), allocatable :: tbuffer(:)

  write( dirname, '(a16,a1,i5.5)') prefix, 'K', ikpt




  ! Open eigval.xml
  if( nsppol .eq. 1) then
    write( filename, '(a,a,a)' ) trim( dirname ), '/', 'eigenval.xml'
  else
    if( isppol .eq. 1 ) then
      write( filename, '(a,a,a)' ) trim( dirname ), '/', 'eigenval1.xml'
    else
      write( filename, '(a,a,a)' ) trim( dirname ), '/', 'eigenval2.xml'
    endif
  endif
  call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) 'Failed to open ', trim(filename), ierr
    return
  endif
  call iotk_scan_dat( 99, "EIGENVALUES", eigen, IERR=ierr )
  if( ierr .ne. 0 ) then
     write(6,*) 'Failed to read eigenvalues ', ierr
    return
  endif
  call iotk_scan_dat( 99, "OCCUPATIONS", occ, IERR=ierr )
  if( ierr .ne. 0 ) then
     write(6,*) 'Failed to read occupations ', ierr
    return
  endif
  call iotk_close_read( 99 )



  ! gkvectors.dat
  write( filename, '(a,a,a)' ) trim( dirname ), '/', 'gkvectors.dat'
  call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr )
  if( ierr .ne. 0 ) then 
    write(6,*) 'Failed to open ', trim(filename), ierr
    return
  endif
  call iotk_scan_dat( 99, "NUMBER_OF_GK-VECTORS", npw, IERR=ierr )
  if( ierr .ne. 0 ) then
     write(6,*) 'Failed to read number of gk-vectors ', ierr
    return
  endif
  call iotk_scan_dat( 99, "GRID", kg(1:3,1:npw), IERR=ierr )
  if( ierr .ne. 0 ) then
     write(6,*) 'Failed to read grid ', ierr
    return
  endif
  call iotk_close_read( 99 )  



  allocate( tbuffer( npw ) )
  if( nsppol .eq. 1) then
    write( filename, '(a,a,a)' ) trim( dirname ), '/', 'evc.dat'
  else 
    if( isppol .eq. 1 ) then
      write( filename, '(a,a,a)' ) trim( dirname ), '/', 'evc1.dat'
    else
      write( filename, '(a,a,a)' ) trim( dirname ), '/', 'evc2.dat'
    endif
  endif
  call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr ) 
  if( ierr .ne. 0 ) then 
    write(6,*) 'Failed to open ', trim(filename), ierr
    return
  endif
  do i = 1, maxband
    call iotk_scan_dat( 99, "evc" // trim(iotk_index( i ) ), tbuffer(1:npw), IERR=ierr )
    if( ierr .ne. 0 ) return
    do j = 1, npw
      cg( i, 2*j-1 ) = real( tbuffer( j ) )
      cg( i, 2*j ) = aimag( tbuffer( j ) )
    enddo
  enddo
  call iotk_close_read( 99 )
  deallocate( tbuffer )



end subroutine
  
