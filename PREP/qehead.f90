! Uses iotk to read data-file.xml to determine the needed information about the wavefunctions
! open up datafile
! 
! return the other ints
!
! maxband gives the largest number of bands any kpt might have
! maxnpw  "" for planewaves
! nspin is 1 or 2 for spin polarized calcs
! nspinor is for spinors which are *not* supported
! nkpt gives the number of kpts
!
!
subroutine qehead( datafile, maxband, maxnpw, nspin, nspinor, nkpt, ierr )
!
!
  use iotk_module
  implicit none
  !
  character(len=128), intent(in) :: datafile
  integer, intent( out ) :: maxband
  integer, intent( out ) :: maxnpw
  integer, intent( out ) :: nspin
  integer, intent( out ) :: nspinor
  integer, intent( out ) :: nkpt
  integer, intent(inout) :: ierr

  logical :: non_colinear

  
  call iotk_open_read ( 99, FILE = trim(datafile), IERR=ierr )
  if( ierr .ne. 0 ) return

  ! PLANE_WAVES
  call iotk_scan_begin( 99, "PLANE_WAVES", IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) 'Failed at begin PLANE_WAVES', ierr
    return
  endif

  call iotk_scan_dat( 99, "MAX_NUMBER_OF_GK-VECTORS", maxnpw, IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) "Failed to grab gvects inside plane_waves", ierr
    return
  endif

  call iotk_scan_end( 99, "PLANE_WAVES", IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) 'Failed to end PLANE_WAVES', ierr
    return
  endif
  ! END PLANE_WAVES

  ! BAND_STRUCTURE_INFO
  call iotk_scan_begin( 99, "BAND_STRUCTURE_INFO", IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) 'Failed at begin BAND_STRUCTURE_INFO', ierr
    return
  endif

  call iotk_scan_dat( 99, "NUMBER_OF_K-POINTS", nkpt, IERR=ierr )
  if( ierr .ne. 0 ) return
  
  call iotk_scan_dat( 99, "NUMBER_OF_SPIN_COMPONENTS", nspin, IERR=ierr )
  if( ierr .ne. 0 ) return
  
  call iotk_scan_dat( 99, "NON-COLINEAR_CALCULATION", non_colinear, IERR=ierr )
  if( ierr .ne. 0 ) return

  if( non_colinear ) then
    nspinor = 2
  else
    nspinor = 1
  endif


  if( nspin .eq. 1 ) then
    call iotk_scan_dat( 99, "NUMBER_OF_BANDS", maxband, IERR=ierr )
    if( ierr .ne. 0 ) return
  elseif( nspin .eq. 2 ) then
    ! Not sure yet
  else
    ierr = -1
    return
  endif

  call iotk_scan_end( 99, "BAND_STRUCTURE_INFO", IERR=ierr )
  if( ierr .ne. 0 ) then
    write(6,*) 'Failed to end BAND_STRUCTURE_INFO', ierr
    return
  endif
  ! END BAND

  call iotk_close_read( 99 )

end subroutine qehead

