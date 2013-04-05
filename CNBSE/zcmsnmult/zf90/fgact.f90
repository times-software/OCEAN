subroutine fgact( lvl, lvh, lmin, lmax, lc, nc, npmax, nproj, jbeg, mham, nq, itot, jtot, &
     mhr, mhi, n, nspn, v, hv, mpcr, mpci, inter, celvol )
  implicit none
  !
  integer :: lvl, lvh, lmin, lmax, lc, nc, npmax, nq, itot, jtot, n, nspn
  integer, dimension( lvl : lvh ) :: nproj, jbeg, mham
  real( kind = kind( 1.0d0 ) ) :: inter, celvol
  real( kind = kind( 1.0d0 ) ), dimension( n, nc, 2 ) :: v, hv
  real( kind = kind( 1.0d0 ) ), dimension( n, npmax, -lmax : lmax, lmin : lmax, nspn ) :: mpcr, mpci
  real( kind = kind( 1.0d0 ) ), dimension( jtot ) :: mhr, mhi
  !
  integer :: lv, ii, ic, ivml, nu, j1, jj, ispn
  real( kind = kind( 1.0d0 ) ) :: mul
  real( kind = kind( 1.0d0 ) ), dimension( itot ) :: pwr, pwi, hpwr, hpwi
  !
  ! lvl, lvh is not enough for omp
  ! should probably pull thread spawning out of the do loop though

  mul = inter / ( dble( nq ) * celvol )
  do lv = lvl, lvh
!     pwr( : ) = 0
!     pwi( : ) = 0
     hpwr( : ) = 0
     hpwi( : ) = 0
!$OMP PARALLEL &
!$OMP PRIVATE( ic, ivml, nu, ii, jj, j1, ispn ) &
!$OMP SHARED( mul, nproj, nc, nspn, lv, mham, jbeg, pwr, pwi, mhr, mhi, mpcr, mpci, v, hv, hpwr, hpwi ) &
!$OMP DEFAULT( NONE )

!$OMP DO 
     do ic = 1, nc
       ispn = 2 - mod( ic, 2 )
!       ispn = 1 + mod( ic, 2 )
       if( nspn .eq. 1 ) then
         ispn = 1
       endif
        do ivml = -lv, lv
           do nu = 1, nproj( lv )
              ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
              pwr( ii ) = sum( v( :, ic, 1 ) * mpcr( :, nu, ivml, lv, ispn ) - v( :, ic, 2 ) * mpci( :, nu, ivml, lv, ispn ) )
              pwi( ii ) = sum( v( :, ic, 1 ) * mpci( :, nu, ivml, lv, ispn ) + v( :, ic, 2 ) * mpcr( :, nu, ivml, lv, ispn ) )
           end do
        end do
     end do
!$OMP END DO


!$OMP DO
     do ii = 1, mham( lv )
        do jj = 1, mham( lv )
           j1 = jbeg( lv ) + ( jj -1 ) + ( ii - 1 ) * mham( lv )
           hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
           hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
        end do
        hpwr( ii ) = hpwr( ii ) * mul
        hpwi( ii ) = hpwi( ii ) * mul
     end do
!$OMP END DO


!$OMP DO
     do ic = 1, nc
       ispn = 2 - mod( ic, 2 )
!       ispn = 1 + mod( ic, 2 )
       if( nspn .eq. 1 ) then
         ispn = 1
       endif
        do ivml = -lv, lv
           do nu = 1, nproj( lv )
              ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
              hv( :, ic, 1 ) = hv( :, ic, 1 ) + hpwr( ii ) * mpcr( :, nu, ivml, lv, ispn ) + hpwi( ii ) * mpci( :, nu, ivml, lv, ispn )
              hv( :, ic, 2 ) = hv( :, ic, 2 ) + hpwi( ii ) * mpcr( :, nu, ivml, lv, ispn ) - hpwr( ii ) * mpci( :, nu, ivml, lv, ispn )
           end do
        end do
     end do
!$OMP END DO
!$OMP END PARALLEL 
  end do
  !
  return
end subroutine fgact
