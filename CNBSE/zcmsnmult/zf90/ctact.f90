subroutine ctact( nc, n, nq, nbd, nspn, celvol, inter, v, hv, lmin, lmax, npmax, nproj, mpcr, mpci, mpm )
  implicit none
  !
  integer :: nc, n, nq, nbd, nspn, lmin, lmax, npmax
  real( kind = kind( 1.0d0 ) ) :: celvol, inter 
  integer :: nproj( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: mpcr( nq * nbd, npmax, -lmax : lmax, lmin : lmax, nspn )
  real( kind = kind( 1.0d0 ) ) :: mpci( nq * nbd, npmax, -lmax : lmax, lmin : lmax, nspn )
  real( kind = kind( 1.0d0 ) ) :: mpm( npmax, npmax, lmin : lmax )
  real( kind = kind( 1.0d0 ) ), dimension( n, nc, 2 ) :: v, hv
  !
  integer :: ic, l, m, i, nu, ispn
  real( kind = kind( 1.0d0 ) ) :: mul
  real( kind = kind( 1.0d0 ) ), dimension( npmax ) :: ampr, ampi, hampr, hampi
  !
  mul = inter / ( dble( nq ) * celvol )
!$OMP PARALLEL DO &
!$OMP DEFAULT( NONE ) &
!$OMP PRIVATE( ic, l, i, nu, m, ampr, ampi, hampr, hampi, ispn ) &
!$OMP SHARED( mul, nc, lmin, lmax, n, nproj, nspn ) &
!$OMP SHARED( hv, v, mpcr, mpci, mpm )
  do ic = 1, nc
    ispn = 2 - mod( ic, 2 )
     if( nspn .eq. 1 ) then
       ispn = 1
     endif
     do l = lmin, lmax
        do m = -l, l
           ampr( : ) = 0
           ampi( : ) = 0
           do i = 1, n
              ampr( : ) = ampr( : ) + v( i, ic, 1 ) * mpcr( i, :, m, l, ispn ) - v( i, ic, 2 ) * mpci( i, :, m, l, ispn )
              ampi( : ) = ampi( : ) + v( i, ic, 1 ) * mpci( i, :, m, l, ispn ) + v( i, ic, 2 ) * mpcr( i, :, m, l, ispn )
           end do
           hampr( : ) = 0
           hampi( : ) = 0
           do nu = 1, nproj( l )
              hampr( : ) = hampr( : ) - mpm( :, nu, l ) * ampr( nu )
              hampi( : ) = hampi( : ) - mpm( :, nu, l ) * ampi( nu )
           end do
           hampr( : ) = hampr( : ) * mul
           hampi( : ) = hampi( : ) * mul
           do nu = 1, nproj( l )
              hv( :, ic, 1 ) = hv( :, ic, 1 ) &
                             + mpcr( :, nu, m, l, ispn ) * hampr( nu ) + mpci( :, nu, m, l, ispn ) * hampi( nu )
              hv( :, ic, 2 ) = hv( :, ic, 2 ) &
                             + mpcr( :, nu, m, l, ispn ) * hampi( nu ) - mpci( :, nu, m, l, ispn ) * hampr( nu )
           end do
        end do
     end do
  end do
!$OMP END PARALLEL DO
  !
  return
end subroutine ctact 
