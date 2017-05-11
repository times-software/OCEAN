! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine qe_grabwf_noiotk(ikpt, isppol, nsppol, maxband, maxnpw, kg, eigen, cg, npw, ierr )
  implicit none
  integer, intent( in ) :: ikpt, isppol, nsppol, maxband, maxnpw
  integer, intent( out ) :: npw
  integer, intent( inout ) :: ierr
  real(kind=kind(1.0d0)), intent( out ) :: eigen(maxband), cg(maxband,2*maxnpw)
  integer, intent( out ) :: kg( 3, maxnpw )
  !
  character(len=22) :: dirname
  character(len=16) :: prefix = 'Out/system.save/' 
  character(len=128) :: filename
  integer :: i, j, crap, max_ngvecs
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
  open( unit=99, file=filename, form='formatted', status='old', IOSTAT=ierr )
  do i = 1, 9
    read(99,*) 
  enddo

  do i = 1, maxband
    read(99,*) eigen( i )
  enddo

  close( 99 )
  



  ! gkvectors.dat
  write( filename, '(a,a,a)' ) trim( dirname ), '/', 'gkvectors.dat'
  open( unit=99, file=filename, form='unformatted', status='old' )
  do i = 1, 12
    read( 99 )
  enddo
  read(99) crap, npw
  do i = 1, 4
    read(99)
  enddo
  read(99) crap, max_ngvecs
  if( max_ngvecs .gt. maxnpw ) then
    write(6,*) 'K-point specific max_ngevcs is larger than expected'
    write(6,*) '  ', max_ngvecs, maxnpw
  endif

  do i = 1, 19
    read( 99 )
  enddo

  read( 99 ) crap, kg( :, 1:npw )
  close( 99 )

  


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

  open( unit=99, file = filename, form='unformatted', status='old' )
  do i = 1, 12
    read( 99 )
  enddo

  do i = 1, maxband
    read(99)
    read(99)
!    call clear_tag( 99 )
    read(99) crap, tbuffer
    do j = 1, npw
      cg( i, 2*j-1 ) = real( tbuffer( j ) )
      cg( i, 2*j ) = aimag( tbuffer( j ) )
    enddo
    read(99)
    read(99)
  enddo
  close(99)
  deallocate( tbuffer )



end subroutine
  

subroutine clear_tag( fh )
  implicit none
  integer, intent( in ) :: fh
  !
  integer :: counter
  integer(kind=4) :: line_count, crap
  character( len=1 ), allocatable :: io_holder( : )

  read( fh ) crap
  counter = crap / 256
  allocate( io_holder( counter ) )
  read( fh ) line_count, io_holder
!  write( 6, * )crap,  line_count
  line_count = line_count - 128
  line_count = line_count / 256
!  write(6,*) counter, line_count
  write(6,*) io_holder(2:counter-1)
  deallocate( io_holder )

end subroutine clear_tag

