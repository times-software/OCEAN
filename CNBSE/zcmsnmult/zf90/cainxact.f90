! Copyright (C) 2013 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cainxact( n, nspn, nc, lc, nq, nbd, e0, self_energy, v, hv, pref, &
     iq, qphys, zn, inter, bvec, celvol, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, renorm, ptab, &
     vms, cml, cms, xi, lmin, lmax, npmax, nproj, mpcr, mpci, mpm, lvl, lvh, jbeg, mham, itot, jtot, mhr, mhi, somel, core_offset )
  implicit none
  !
  integer :: n, nc, lc, nq, nbd, nx, ny, nz, lmin, lmax, lvl, lvh, itot, jtot, zn( 3 ), npmax, spin_test, nspn
  real( kind = kind( 1.0d0 ) ) :: pref, bvec( 3, 3 ), celvol, amet( 3, 3 ), tau( 3 ), rcut, rzero, ptab( 100 ), inter, renorm, xi
  integer :: iq( 3, nq ), nproj( lmin : lmax ), jbeg( lvl : lvh ), mham( lvl : lvh )
  real( kind = kind( 1.0d0 ) ) :: vms( nc ), cml( nc ), cms( nc ), mhr( jtot ), mhi( jtot )
  real( kind = kind( 1.0d0 ) ) :: qphys( 3, nq ), e0( n, nspn ), v( n, nc, 2 ), hv( n, nc, 2 )
  complex( kind = kind( 1.0d0 ) ) :: self_energy( n, nspn )
  real( kind = kind( 1.0d0 ) ) :: ur( nx * ny * nz, nbd, nq, nspn ), ui( nx * ny * nz, nbd, nq, nspn )
  real( kind = kind( 1.0d0 ) ), dimension( n, npmax, -lmax : lmax, lmin : lmax ) :: mpcr, mpci
  real( kind = kind( 1.0d0 ) ) :: mpm( npmax, npmax, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: somel( nc, nc, 2 )
  real( kind = kind( 1.0d0 ) ) :: core_offset
  !
  integer :: ic
  real( kind = kind( 1.0d0 ) ) :: eps

  !
  ! band + monopole attraction
  open( unit=99, file='epsilon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) eps
  close( unit=99 )
!
! Need to switch from icms, icml ordering to test wether it is spin up or down from the wfns
!
!$OMP PARALLEL DO  & 
!$OMP& SCHEDULE( STATIC  ) &
!$OMP& PRIVATE( ic, spin_test ) &
!$OMP& SHARED( nc, v, n, e0, hv, nq, nbd, zn, inter, amet, nx, ny, nz, ur, ui, tau, rcut, rzero, ptab, eps, nspn, self_energy, core_offset) &
!$OMP& DEFAULT ( NONE )
  do ic = 1, nc
     ! if ic is odd, spin_test is 1, if even it is 2
     spin_test = 2 - mod( ic, 2 )
     if( nspn .eq. 1 ) then
       spin_test = 1
     endif
     call cainchact( v( :, ic, 1 ), v( :, ic, 2 ), n, e0( :, spin_test ), self_energy( :, spin_test ), &
                     hv( :, ic, 1 ), hv( :, ic, 2 ), nq, nbd, zn, inter, amet, nx, ny, nz, &
                     ur( :, :, :, spin_test), ui( :, :, :, spin_test ), tau, rcut, rzero, ptab, eps, core_offset )
  end do
!$OMP END PARALLEL DO
  !
  ! spin-orbit
  call newsoact( nc, somel, n, v, hv )
  !
  if ( inter .gt. 0.01d0 ) then
     ! 
     ! central term
     call ctact( nc, n, nq, nbd, nspn, celvol, inter, v, hv, lmin, lmax, npmax, nproj, mpcr, mpci, mpm )
     !
     ! multiplet
     call fgact( lvl, lvh, lmin, lmax, lc, nc, npmax, nproj( lvl ), jbeg, mham, nq, itot, jtot, &
          mhr, mhi, n, nspn, v, hv, mpcr, mpci, inter, celvol )
  end if
  !
  return
end subroutine cainxact
