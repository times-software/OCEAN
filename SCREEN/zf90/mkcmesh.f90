! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine mkcmesh( nang_in, nr, rad, drad, element, indx, posn, wpt, drel, avec, ierr )
  implicit none
!
  integer, intent( in ) :: nang_in, nr, indx
  real( kind = kind( 1.0d0 )), intent( in ) :: avec( 3, 3 )
  real( kind = kind( 1.0d0 )), intent( in ), dimension( nr ) :: rad, drad
  real( kind = kind( 1.0d0 )), intent( out ) :: posn( 3, nang_in*nr ), wpt( nang_in*nr ), &
                                                drel( nang_in*nr )
  character(len=2), intent( in ) :: element
  integer, intent( inout ) :: ierr
!
  integer :: i, j, ii, nang
  real( kind = kind( 1.0d0 ) ) :: pi, su, tmp, alpha(3), tau(3), dr
  real( kind = kind( 1.0d0 ) ), allocatable :: x( :, : ), wang( : ), wr( : )
!
  write ( 6, * ) ' we are in mkcmesh'
  pi = 4.0d0 * atan( 1.0d0 )
  allocate( x( 3, nang_in ), wang( nang_in ), wr( nr ) )
  open( unit=99, file='specpnt', form='formatted', status='old')
  rewind 99
  read ( 99, * ) nang
  if( nang .ne. nang_in ) then
    ierr = -1
    return
  endif
  su = 0
  do j = 1, nang
    read ( 99, * ) x( :, j ), wang( j )
    su = su + wang( j )
    tmp = dot_product( x( :, j ), x( :, j ) )
    x( :, j ) = x( :, j ) / sqrt( tmp )
  end do
  close( unit=99 )
  wang( : ) = wang( : ) * ( 4.0d0 * pi ) / su
  do i = 1, nr
    wr( i ) = drad( i ) * rad( i ) ** 2
  end do
  call snatch( element, indx, alpha )
  tau = 0
  do i = 1, 3
    tau( : ) = tau( : ) + alpha( i ) * avec( :, i )
  end do
  write(6,*) alpha(:)
  write(6,*) tau(:)
  ii = 0
  do i = 1, nr
    do j = 1, nang
      ii = ii + 1
      drel( ii ) = rad( i )
      wpt( ii ) = wr( i ) * wang( j )
      posn( :, ii ) = tau( : ) + rad( i ) * x( :, j )
    end do
  end do
  deallocate( x, wang, wr )
!
  return
end subroutine mkcmesh
