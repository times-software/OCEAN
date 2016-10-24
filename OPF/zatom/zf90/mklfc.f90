! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine mklfc( lfc, ratio, nr, r, dl, cq, nel, np, occ, phe )
  implicit none
  !
  integer :: lfc, nr, nel, np
  real( kind = kind( 1.0d0 ) ) :: dl, ratio, occ( nel ), phe( nr, nel ), r( nr ), cq( nr ) 
  !
  integer :: i, im
  real( kind = kind( 1.0d0 ) ) :: pi, vq( nr )
  real( kind = kind( 1.0d0 ) ) :: cor0, corp, corm, f, fp, rhs, xl, xh, br, diff, aa, bb, tmp
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  cq = 0
  vq = 0
  if ( lfc .ne. 0 ) then
     do i = 1, np - 1
        cq( : ) = cq( : ) + occ( i ) * phe( :, i ) ** 2
     end do
     do i = np, nel
        vq( : ) = vq( : ) + occ( i ) * phe( :, i ) ** 2
     end do
     if ( ratio .ge. 0 ) then
        im = 0
        do i = nr, 1, -1
           if ( ( im .eq. 0 ) .and. ( ( cq( i ) * ratio ) .gt. vq( i ) ) ) im = i
        end do
     else
        im = log( -ratio / r( 1 ) ) / dl
     end if
     write ( 6, '(1i5,1e15.8)' ) im, r( im )
     cor0 = cq( im ) / r( im )
     corp = cq( im + 1 ) / r( im + 1 )
     corm = cq( im - 1 ) / r( im - 1 )
     f = cor0
     fp = ( corp - corm ) / ( 2 * dl * r( im ) )
     rhs = r( im ) * fp / f
     if ( rhs .gt. 0 ) then
        xl =0
     else
        xl = pi / 2
     end if
     xh = xl + pi / 2
     do
        br = ( xl + xh ) / 2
        diff = tan( br ) - ( br / rhs )
        write ( 6, '(1x,2f20.10)' ) br, diff
        if ( diff .ge. 0 ) xh = br
        if ( diff .le. 0 ) xl = br
        if ( abs( xh - xl ) .lt. 0.0000001d0 ) exit
     end do
     bb = br / r( im )
     aa = f / sin( bb * r( im ) )
     write ( 6, '(1x,1a6,1f10.6,2x,1a2,1f10.6)' ) 'lfc a=', aa, 'b=', bb
     open( unit=99, file='corchg', form='formatted', status='unknown' )
     rewind 99
     do i = 1, nr
       tmp = cq( i )
       if ( i .lt. im ) cq( i ) = r( i ) * aa * sin( bb * r( i ) )
       write ( 99, '(1x,3f20.10)' ) r( i ), tmp, cq( i )
     end do
     close( unit=99 )
  end if
  !
  return
end subroutine mklfc
