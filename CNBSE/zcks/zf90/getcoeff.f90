subroutine getcoeff( l, ng, nbd, ck, tauphs, ylmfac, coeff, lmax, npmax, seanfq, nproj )
  implicit none
  !
  integer :: l, ng, nbd, nproj, lmax, npmax
  double precision :: seanfq( ng, nproj )
  double complex :: ck( ng, nbd )
  double complex :: tauphs( ng ), ylmfac( -l : l, ng )
  double complex :: coeff( -lmax : lmax, nbd, npmax )
  !
  integer :: m, ig, ibd, iproj
  double complex :: su
  !
  do ibd = 1, nbd
     do iproj = 1, nproj
        do m = -l, l
           su = 0
           do ig = 1, ng
              su = su + ylmfac( m, ig ) * tauphs( ig ) * ck( ig, ibd ) * seanfq( ig, iproj )
           end do
           coeff( m, ibd, iproj ) = su
        end do
     end do
  end do
  !
  return
end subroutine getcoeff
