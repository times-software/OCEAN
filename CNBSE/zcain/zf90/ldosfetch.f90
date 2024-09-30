subroutine ldosfetch( spcttype, lmin, lmax, npmax, nproj, rmax, ldosmel )
  implicit none
  !
  character(len=10), intent( in ) :: spcttype
  integer, intent( in ) :: lmin, lmax, npmax, nproj( lmin : lmax )
  real( kind = kind( 1.0d0 ) ), intent( in ) :: rmax
  real( kind = kind( 1.0d0 ) ), intent( out ) :: ldosmel( npmax  )
  !
  real( kind = kind( 1.0d0 ) ) :: dum, wv( npmax ), rmax2, dl, su
  real( kind = kind( 1.0d0 ) ), allocatable :: phps (:, : ), r(:)
  integer :: targetL, zz, nn, ll, nr, dumi, i, j, nr2, ii
  character( len=8 ) :: s8
  character( len=11 ) :: s11


  ldosmel( : ) = 0.0d0

  select case( spcttype(5:5) )
    case ('0')
      targetL = 0
    case ('1')
      targetL = 1
    case ('2')
      targetL = 2
    case ('3')
      targetL = 3
  end select

  if( targetL .lt. lmin .or. targetL .gt. lmax) then
    write(6,*) 'WARNING: ldos out of range of OPFS', targetL, lmin, lmax
    return
  endif
    

  open( unit=99, file='ZNL', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) zz, nn, ll
  close( unit=99 )

  write( s11, '(1a8,1i3.3)') 'radfilez', zz
  open(unit=99, file=s11, form='formatted', status='old' )
  rewind 99
  read(99,*) rmax2, dumi, nr
  close(99)

  if( rmax2 .lt. rmax ) then
    write(6,*) 'Using smaller maximum r for l-dos than requested:', rmax2, rmax
  else
    rmax2 = rmax
  endif

  allocate( phps( nr, nproj(targetL) ), r( nr ) )

  write ( s8, '(1a2,1i1,1a1,1i3.3)' ) 'ps', targetL, 'z', zz
  write(6 ,* ) '!! ', s8
  open( unit=99, file=s8, form='formatted', status='old' )
  rewind 99
  do j = 1, nr
    read ( 99, * ) r(j), wv( 1 : nproj( targetL ) )
    do i = 1, nproj( targetL )
       phps( j, i ) = wv( i )
    end do
  end do
  close( unit=99 )

  nr2 = nr
  do j = 1, nr
    if( r( j ) .gt. rmax2 ) then
      nr2 = j-1
      exit
    endif
  enddo

  ii = nr / 2
  dl = log( r( 1 + ii ) / r( 1 ) ) / dble( ii )
  do i = 1, nproj( targetL )
    su = 0.0d0
    do j = 1, nr2
      su = su + dl * r( j ) * (phps( j, i )*r(j))**2
    enddo
    ldosmel( i ) = su
    write(6,*) 'LDOSMEL:', su
  enddo

  deallocate( r, phps )

end subroutine
