subroutine sdump( l, psflag, nel, nl, nr, r, phe )
  implicit none
  !
  integer :: l, nel, nr
  logical :: psflag
  integer :: nl( nel )
  double precision :: r( nr )
  double precision :: phe( nr, nel )
  !
  integer ::  i, j, k, nco
  integer, allocatable :: ilp( : )
  !
  character * 2 :: str
  character * 3 :: name
  !
  allocate( ilp( nel ) )
  !
  if ( psflag ) then
     nco = 0
     str = 'ps'
  else
     open(unit=99,file='corcon', form='formatted',status='unknown')
     rewind 99
     read ( 99, * ) nco
     close( unit=99 )
     str = 'ae'
  end if
  k = 0
  do i = nco + 1, nco + nel
     if ( nl( i ) .eq. l ) then
        k = k + 1 
        ilp( k ) = i
     end if
  end do
  write ( 6, * ) k, nco, nel
  write ( unit=name, fmt='(1a2,1i1)' ) str, l
  open(unit=99,file=name, form='formatted',status='unknown')
  rewind 99
  do i=1,nr
     write (99,'(1x,20f25.15)') r(i),(phe(i,ilp(j))/r(i),j=1,k)
  end do
  close(unit=99)
  !
  deallocate( ilp )
  !
  return
end subroutine sdump
