! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine seanitup( ng, q, gvec, bmet, l, nproj, &
     &               seanfq, npmax, nq, dq, fttab )
  implicit none
  !
  integer :: ng, l, nproj, nq, npmax
  integer :: gvec( 3, ng )
  !
  double precision :: dq, q( 3 ), bmet( 3, 3 )
  double precision :: seanfq( ng, nproj )
  double precision :: fttab( nq, nproj )
  !
  integer :: i, j, i1, i2, ii
  double precision :: gc( 3 ), qbar, qii, x
  double precision :: pm1, pz0, pp1, pp2
  double precision :: wm1, wz0, wp1, wp2
  !
  pm1 = - 1.d0 / 6.d0
  pz0 = 1.d0 / 2.d0
  pp1 = - 1.d0 / 2.d0
  pp2 = 1.d0 / 6.d0
  do j = 1, ng
     gc( : ) = q( : ) + gvec( :, j )
     qbar = 0.d0
     do i1 = 1, 3
        do i2 = 1, 3
           qbar = qbar + gc( i1 ) * gc( i2 ) * bmet( i1, i2 )
        end do
     end do
     qbar = sqrt( qbar )
     ii = 1 + qbar / dq
     if ( ii .lt. 2 ) ii = 2
     if ( ii .gt. nq - 2 ) ii = nq - 2
     qii = dq * ( ii - 1 )
     x = ( qbar - qii ) / dq
     wm1 = pm1 * x * ( x - 1.d0 ) * ( x - 2.d0 )
     wz0 = pz0 * ( x + 1.d0 ) * ( x - 1.d0 ) * ( x - 2.d0 ) 
     wp1 = pp1 * ( x + 1.d0 ) * x * ( x - 2.d0 )
     wp2 = pp2 * ( x + 1.d0 ) * x * ( x - 1.d0 )
     do i = 1, nproj
        seanfq( j, i ) = fttab( ii - 1, i ) * wm1 + &
             &           fttab( ii, i ) * wz0 + &
             &           fttab( ii + 1, i ) * wp1 + &
             &           fttab( ii + 2, i ) * wp2
     end do
  end do
  !
  return
end subroutine seanitup
