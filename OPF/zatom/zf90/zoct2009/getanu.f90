subroutine getanu( l, nr, r, rcut, dl, phi )
  implicit none
  integer :: l, nr, i, j, nproj, irc, lmin, lmax
  double precision :: dl, rcut, phi( nr ), r( nr ), su
  double precision, allocatable :: psproj( :, : ), f( : ), a( : )
  character * 3 :: name
  open( unit=99, file='prjfile', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lmin, lmax
  do i = lmin, l
     read ( 99, * ) nproj
  end do
  close( unit=99 )
  do i = 1, nr, 4
     if ( r( i ) .lt. rcut ) irc = i
  end do
  allocate( f( nr ), psproj( nr, nproj ), a( nproj ) )
  write ( unit=name, fmt='(1a2,1i1)' ) 'ps', l ! mknam
  open( unit=99, file=name, form='formatted', status='unknown' )
  do i = 1, nr
     read ( 99, * ) r( i ), ( psproj( i, j ), j = 1, nproj )
  end do
  close( unit=99 )
  do i = 1, nproj
     f( : ) = psproj( :, i ) * phi( : ) / r( : )
     call bintegrate( nr, r, dl, f, a( i ), irc )
  end do
  write ( unit=name, fmt='(1a2,1i1)' ) 'as', l ! mknam
  open( unit=99, file=name, form='formatted', status='unknown' )
  do i = 1, nproj
     write ( 99, '(2x,1e15.8)' ) a( i ), 0
  end do
  close( unit=99 )
  write ( unit=name, fmt='(1a2,1i1)' ) 'vu', l ! mknam
  open( unit=99, file=name, form='formatted', status='unknown' )
  do i = 1, irc
     su =0.d0
     do j = 1, nproj
        su = su + a( j ) * psproj( i, j )
     end do
     write ( 99, '(3(2x,1e15.8))' ) r( i ), phi( i ) / r( i ), su
  end do
  close( unit=99 )
  deallocate( psproj, f, a )
  return
end subroutine getanu
