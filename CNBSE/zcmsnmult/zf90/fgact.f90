subroutine fgact( lvl, lvh, lmin, lmax, lc, nc, npmax, nproj, jbeg, mham, nq, itot, jtot, &
     mhr, mhi, n, v, hv, mpcr, mpci, inter, celvol )
  implicit none
  !
  integer :: lvl, lvh, lmin, lmax, lc, nc, npmax, nq, itot, jtot, n
  integer, dimension( lvl : lvh ) :: nproj, jbeg, mham
  real( kind = kind( 1.0d0 ) ) :: inter, celvol
  real( kind = kind( 1.0d0 ) ), dimension( n, nc, 2 ) :: v, hv
  real( kind = kind( 1.0d0 ) ), dimension( n, npmax, -lmax : lmax, lmin : lmax ) :: mpcr, mpci
  real( kind = kind( 1.0d0 ) ), dimension( jtot ) :: mhr, mhi
  !
  integer :: lv, ii, ic, ivml, nu, j1, jj
  real( kind = kind( 1.0d0 ) ) :: mul
  real( kind = kind( 1.0d0 ) ), dimension( itot ) :: pwr, pwi, hpwr, hpwi
  !
  mul = inter / ( dble( nq ) * celvol )
  do lv = lvl, lvh
     pwr( : ) = 0
     pwi( : ) = 0
     ii = 0
     do ic = 1, nc
        do ivml = -lv, lv
           do nu = 1, nproj( lv )
              ii = ii + 1
              pwr( ii ) = sum( v( :, ic, 1 ) * mpcr( :, nu, ivml, lv ) - v( :, ic, 2 ) * mpci( :, nu, ivml, lv ) )
              pwi( ii ) = sum( v( :, ic, 1 ) * mpci( :, nu, ivml, lv ) + v( :, ic, 2 ) * mpcr( :, nu, ivml, lv ) )
           end do
        end do
     end do
     hpwr( : ) = 0
     hpwi( : ) = 0
     j1 = jbeg( lv )
     do jj = 1, mham( lv )
        do ii = 1, mham( lv )
           hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
           hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
           j1 = j1 + 1
        end do
     end do
     hpwr = hpwr * mul
     hpwi = hpwi * mul
     ii = 0
     do ic = 1, nc
        do ivml = -lv, lv
           do nu = 1, nproj( lv )
              ii = ii + 1
              hv( :, ic, 1 ) = hv( :, ic, 1 ) + hpwr( ii ) * mpcr( :, nu, ivml, lv ) + hpwi( ii ) * mpci( :, nu, ivml, lv )
              hv( :, ic, 2 ) = hv( :, ic, 2 ) + hpwi( ii ) * mpcr( :, nu, ivml, lv ) - hpwr( ii ) * mpci( :, nu, ivml, lv )
           end do
        end do
     end do
  end do
  !
  return
end subroutine fgact
