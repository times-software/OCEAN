subroutine cainxact( n, nc, lc, nq, nbd, e0, v, hv, pref, &
     iq, qphys, zn, inter, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, &
     vms, cml, cms, xi, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, mhr, mhi, somel )
  implicit none
  !
  integer :: n, nc, lc, nq, nbd, nx, ny, nz, lmin, lmax, lvl, lvh, itot, jtot, zn( 3 ), npmax
  real( kind = kind( 1.0d0 ) ) :: pref, bvec( 3, 3 ), celvol, amet( 3, 3 ), tau( 3 ), rcut, rzero, ptab( 100 ), inter, renorm, xi
  integer :: iq( 3, nq ), nproj( lmin : lmax ), jbeg( lvl : lvh ), mham( lvl : lvh )
  real( kind = kind( 1.0d0 ) ) :: vms( nc ), cml( nc ), cms( nc ), mhr( jtot ), mhi( jtot )
  real( kind = kind( 1.0d0 ) ) :: qphys( 3, nq ), e0( n ), v( n, nc, 2 ), hv( n, nc, 2 )
  real( kind = kind( 1.0d0 ) ) :: ur( nx * ny * nz, nbd, nq ), ui( nx * ny * nz, nbd, nq )
  real( kind = kind( 1.0d0 ) ), dimension( n, npmax, -lmax : lmax, lmin : lmax ) :: mpcr, mpci
  real( kind = kind( 1.0d0 ) ) :: mpm( npmax, npmax, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: somel( nc, nc, 2 )
  !
  integer :: ic
  !
  ! band + monopole attraction
  do ic = 1, nc
     call cainchact( v( :, ic, 1 ), v( :, ic, 2 ), n, e0, hv( :, ic, 1 ), hv( :, ic, 2 ), nq, nbd, zn, inter, amet, nx, ny, nz, &
          ur, ui, tau, rcut, rzero, ptab )
  end do
  !
  ! spin-orbit
  call newsoact( nc, somel, n, v, hv )
  !
  if ( inter .gt. 0.01d0 ) then
     ! 
     ! central term
     call ctact( nc, n, nq, nbd, celvol, inter, v, hv, lmin, lmax, npmax, nproj, mpcr, mpci, mpm )
     !
     ! multiplet
     call fgact( lvl, lvh, lmin, lmax, lc, nc, npmax, nproj( lvl ), jbeg, mham, nq, itot, jtot, &
          mhr, mhi, n, v, hv, mpcr, mpci, inter, celvol )
  end if
  !
  return
end subroutine cainxact
