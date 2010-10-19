subroutine ctact( nc, n, nq, nbd, celvol, inter, v, hv, lmin, lmax, npmax, nproj, mpcr, mpci, mpm )
  implicit none
  !
  integer :: nc, n, nq, nbd, lmin, lmax, npmax
  real( kind = kind( 1.0d0 ) ) :: celvol, inter 
  integer :: nproj( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: mpcr( nq * nbd, npmax, -lmax : lmax, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: mpci( nq * nbd, npmax, -lmax : lmax, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: mpm( npmax, npmax, lmin : lmax )
  real( kind = kind( 1.0d0 ) ), dimension( n, nc, 2 ) :: v, hv
  !
  integer :: ic, l, m, i, nu 
  real( kind = kind( 1.0d0 ) ) :: mul
  real( kind = kind( 1.0d0 ) ), dimension( npmax ) :: ampr, ampi, hampr, hampi
  !
  mul = inter / ( dble( nq ) * celvol )
  do ic = 1, nc
     do l = lmin, lmax
        do m = -l, l
           ampr( : ) = 0
           ampi( : ) = 0
           do i = 1, n
              ampr( : ) = ampr( : ) + v( i, ic, 1 ) * mpcr( i, :, m, l ) - v( i, ic, 2 ) * mpci( i, :, m, l )
              ampi( : ) = ampi( : ) + v( i, ic, 1 ) * mpci( i, :, m, l ) + v( i, ic, 2 ) * mpcr( i, :, m, l )
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
              hv( :, ic, 1 ) = hv( :, ic, 1 ) + mpcr( :, nu, m, l ) * hampr( nu ) + mpci( :, nu, m, l ) * hampi( nu )
              hv( :, ic, 2 ) = hv( :, ic, 2 ) + mpcr( :, nu, m, l ) * hampi( nu ) - mpci( :, nu, m, l ) * hampr( nu )
           end do
        end do
     end do
  end do
  !
  return
end subroutine ctact 
