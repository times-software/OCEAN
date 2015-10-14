! Copyright (C) 2011 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getcoeff( l, ng, nbd, ck, tauphs, ylmfac, coeff, lmax, npmax, seanfq, nproj, temperature, efermi, ww )
  implicit none
  !
  integer :: l, ng, nbd, nproj, lmax, npmax
  double precision :: seanfq( ng, nproj )
  double complex :: ck( ng, nbd )
  double complex :: tauphs( ng ), ylmfac( -l : l, ng )
  double complex :: coeff( -lmax : lmax, nbd, npmax )
  double precision :: efermi, temperature, ww( nbd )
  !
  integer :: m, ig, ibd, iproj
  double complex :: su
  double precision :: ffactor
  !
  do ibd = 1, nbd
     if( temperature .gt. 0.000001) then
       ffactor = 1.d0 - 1.d0 / ( exp( ( ww( ibd ) - efermi ) / temperature ) + 1 )
       write(6,*) ffactor, ww( ibd )
     else
       ffactor = 1
     endif
     do iproj = 1, nproj
        do m = -l, l
           su = 0
           do ig = 1, ng
              su = su + ylmfac( m, ig ) * tauphs( ig ) * ck( ig, ibd ) * seanfq( ig, iproj )
           end do
           coeff( m, ibd, iproj ) = su * ffactor
        end do
     end do
  end do
  !
  return
end subroutine getcoeff
