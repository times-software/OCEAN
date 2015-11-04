! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine qe_grabwf(ikpt, isppol, nsppol, maxband, maxnpw, kg_unshift, kg_shift, eigen_un, eigen_sh, occ_un, &
                     occ_sh, cg_un, cg_sh, occ_max, unocc_max, nband, un_npw, sh_npw, noshift, ierr )
  use iotk_module
  implicit none
  integer :: ikpt, isppol, nsppol, maxband, maxnpw, nband(2), iband,ii,un_npw,sh_npw,nspinor, occ_max, unocc_max
  double precision :: eigen_un(maxband), eigen_sh(maxband),         &
    occ_un(maxband), occ_sh(maxband), cg_un(maxband,2*maxnpw),      &
    cg_sh(maxband,2*maxnpw)
  integer :: kg_unshift(3,maxnpw), kg_shift(3,maxnpw)
  logical :: noshift
  character(len=22) :: dirname
  character(len=16) :: prefix = 'Out/system.save/'
  character(len=128) :: filename
  integer :: ierr, npw, i, j, ib, ig
  complex(kind=kind(1.0d0)), allocatable :: tbuffer(:), cg_un_tmp(:,:)

  write( dirname, '(a16,a1,i5.5)') prefix, 'K', ikpt




  ! Open eigval.xml
  if( nsppol .eq. 1) then
    write( filename, '(a,a,a)' ) trim( dirname ), '/', 'eigenval.xml'
  else
    if( isppol .eq. 1 ) then
      write( filename, '(a,a,a)' ) trim( dirname ), '/', 'eigenval1.dat'
    else
      write( filename, '(a,a,a)' ) trim( dirname ), '/', 'eigenval2.dat'
    endif
  endif
  call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) 'Failed to open ', trim(filename), ierr
    return
  endif
  call iotk_scan_dat( 99, "EIGENVALUES", eigen_un, IERR=ierr )
  if( ierr .ne. 0 ) then
     write(6,*) 'Failed to read eigenvalues ', ierr
    return
  endif
  call iotk_scan_dat( 99, "OCCUPATIONS", occ_un, IERR=ierr )
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
  call iotk_scan_dat( 99, "GRID", kg_unshift(1:3,1:npw), IERR=ierr )
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
      cg_un( i, 2*j-1 ) = real( tbuffer( j ) )
      cg_un( i, 2*j ) = aimag( tbuffer( j ) )
      cg_sh( i, 2*j-1 ) = real( tbuffer( j ) )
      cg_sh( i, 2*j ) = aimag( tbuffer( j ) )
    enddo
  enddo
  call iotk_close_read( 99 )
  deallocate( tbuffer )


  un_npw = npw
  sh_npw = npw

end 
  
