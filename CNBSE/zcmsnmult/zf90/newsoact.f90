subroutine newsoact( nc, somel, n, v, hv )
  implicit none
  !
  integer :: nc, n
  real( kind = kind( 1.0d0 ) ) :: v( n, nc, 2 ), hv( n, nc, 2 ), somel( nc, nc, 2 )
  !
  integer :: ic, jc
  !
! Limit this to only the 1 loop; do not unroll
! Really want to have some way to explouit cache lines, doesn't matter that much, doesn't take long
!$OMP PARALLEL DO &
!$OMP DEFAULT( NONE ) &
!$OMP PRIVATE( ic, jc ) &
!$OMP SHARED( hv, v, nc, somel )
  do ic = 1, nc
! Could manually unroll this to groups of 4
     do jc = 1, nc
        hv( :, ic, 1 ) = hv( :, ic, 1 ) + somel( ic, jc, 1 ) * v( :, jc, 1 ) - somel( ic, jc, 2 ) * v( :, jc, 2 )
        hv( :, ic, 2 ) = hv( :, ic, 2 ) + somel( ic, jc, 1 ) * v( :, jc, 2 ) + somel( ic, jc, 2 ) * v( :, jc, 1 )
     end do
  end do
!$OMP END PARALLEL DO
  !
  return
end subroutine newsoact
