subroutine amult( l, norb, nr )
  implicit none
  !
  ! This subr. multiplies the A operator by the basis vector of
  ! partial, pseudo wavefunctions (psorb).  The output are the projector
  ! functions (proj).  The A operator is read in from file A and
  ! the wavefunctions are read in from pseudo. 
  ! 
  integer :: i, j, k, l, nr, norb
  real( kind = kind( 1.0d0 ) ) :: su
  real( kind = kind( 1.0d0 ) ), allocatable :: psorb(:,:), A(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: proj(:,:), r(:)
  character * 3 :: name
  !
  allocate( a( norb, norb) )
  write ( unit=name, fmt='(1a2,1i1)' ) 'am', l
  open( unit=99, file=name, form='formatted', status='unknown' )
  rewind 99
  do i = 1, norb
     do j = 1, norb
        read ( 99, * ) a( i, j )
     end do
  end do
  close( unit=99 )
  !
  allocate( psorb( nr, norb ), r( nr ), proj( nr, norb ) )
  write ( unit=name, fmt='(1a2,1i1)' ) 'ps', l
  open( unit=99, file=name, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     read ( 99, * ) r(i), ( psorb( i, k ), k = 1, norb )
  end do
  close ( unit=99 )
  !
  ! ** Projector functions are formed by multiplying A matrix by
  ! ** pseudo basis functions
  !  
  do i = 1, norb
     do j = 1, nr
        su = 0
        do k = 1, norb
           su = su + a( i, k ) * psorb( j, k )
        end do
        proj( j, i ) = su
     end do
  end do
  !
  write ( unit=name, fmt='(1a2,1i1)' ) 'pr', l
  open( unit=99, file=name, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     write ( 99, '(8f22.11)' ) r(i), ( proj( i, j ), j = 1, norb )
  end do
  close ( unit=99 )
  !
  deallocate( psorb, a, proj, r )
  !
  return
end subroutine amult
