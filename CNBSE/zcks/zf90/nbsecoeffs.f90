!
! ng = number of Gvecs for the writing of each wavevector
! kvc = inverse of the g list
! bmet 
! bvec
! ck = complex wavefunction
! ibl = start band
! ibh = stop band
subroutine nbsecoeffs( ng, kvc, bmet, bvec, ck, ibl, ibh, q, ntau, tau, lmin, lmax, nproj, npmax, nqproj, dqproj, fttab, coeff, &
     prefs )
  implicit none
  !
  integer :: ng, ibl, ibh, ntau 
  integer :: lmin, lmax, npmax, nqproj
  integer :: kvc( 3, ng ), nproj( lmin : lmax )
  double precision :: dqproj
  double precision :: bvec( 3, 3 ), bmet( 3, 3 ), q( 3 ), prefs( 0 : 1000 )
  double precision :: tau( 3, ntau )
  double precision :: fttab( nqproj, npmax, lmin : lmax )
  double complex :: coeff( -lmax:lmax, ibl:ibh, npmax, lmin:lmax, ntau )
  !
  integer itau, l
  complex( kind = kind( 1.0d0 ) ) :: ck( ng, ibl : ibh )
  double complex, allocatable :: tauphs( : , : )
  double complex, allocatable :: ylmfac( : , : )
  double precision, allocatable :: sfq( : , : , : )
  !
  allocate( tauphs( ng, ntau ) )
  allocate( sfq( ng, npmax, lmin : lmax ) )
  allocate( ylmfac( ng * ( 2 * lmax + 1 ), lmin : lmax ) )
  !
  do itau = 1, ntau
     call getphase( ng, kvc, q, tau( 1, itau ), tauphs( 1, itau ) )
  end do
  do l = lmin, lmax
     call seanitup( ng, q, kvc, bmet, l, nproj( l ), sfq( 1, 1, l ), npmax, nqproj, dqproj, fttab( 1, 1, l ) )
     call getylmfac( ng, kvc, q, bvec, l, ylmfac( 1, l ), prefs )
  end do
  do itau = 1, ntau
     do l = lmin, lmax
        call getcoeff( l, ng, 1 + ibh - ibl, ck, tauphs( 1, itau ), ylmfac( 1, l ), coeff( -lmax, ibl, 1, l, itau ), lmax, &
             npmax, sfq( 1, 1, l ), nproj( l ) )
     end do
  end do
  !
  deallocate( tauphs, ylmfac, sfq )
  !
  return
end subroutine nbsecoeffs
