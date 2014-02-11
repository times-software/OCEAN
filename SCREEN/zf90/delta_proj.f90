program delta_proj
  implicit none

  integer, parameter :: DP = kind(1.d0)  

  integer :: ZNL(3), lmin, lmax, nr, nrtot
  integer, allocatable :: nprj( : )

  real(DP) :: rmax


  character(len=11) :: prj_filename, rad_filename


  integer :: il


  open( unit=99, file='ZNL', status='old', form='formatted' )
  rewind( 99 )
  read( 99, * ) ZNL(:)
  close(99)

  write(prj_filename,'(a8,i3.3)') 'prjfilez', ZNL(1)
  open( unit=99, file=prj_filename, status='old', form='formatted' )
  rewind(99)
  read(99,*) lmin, lmax 
  allocate( nprj( lmin : lmax ) )
  do il = lmin, lmax
    read( 99, * ) nprj( il )
  enddo
  close( 99 )



  write(rad_filename,'(a8,i3.3)') 'radfilez', ZNL(1)
  open( unit=99, file=rad_filename, status='old', form='formatted' )
  rewind(99)
  read(99,*) rmax, nr, nrtot
  


  write(6,*) lmin, lmax
  write(6,*) nprj(:)
  write(6,*) rmax, nr, nrtot



  do il = lmin, lmax
    ! Read in projectors
    if( il .eq. lmin ) then
      allocate( ae_prj( 0 : nproj( il ) ), ps_prj( 0 : nproj( il ) ) )
    else

    endif

    do iproj = 1, nprj( il )


  








  deallocate( nprj )

end program delta_proj
