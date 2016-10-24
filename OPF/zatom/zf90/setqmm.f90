! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine setqmm( i, orb, l, xj, idoflag, v, zef, zorig, rel, nr, r, r2, dl, xm1, xm2, njrc, vi, psflag )
! subroutine setqmm( orb, l, xj, idoflag, v, zef, zorig, rel, nr, r, r2, dl, xm1, xm2, njrc, vi, psflag ) 
  implicit none
  integer i, l, idoflag, nr
  real( kind = kind( 1.0d0 ) ) xj, zef, zorig, rel, dl
  integer njrc( 4 )
  real( kind = kind( 1.0d0 ) ) v( nr ), r( nr ), r2( nr ), orb( nr )
  real( kind = kind( 1.0d0 ) ) xm1( nr ), xm2( nr ), vi( nr, 7 )
  integer j, lp, lpx, lp2, l1, l2, lu, ij
  real( kind = kind( 1.0d0 ) ) alpha, aa, a2, zaa, za2, d1, d2, w1, w2
  real( kind = kind( 1.0d0 ) ) dvdl, ddvdll, dvdr, ddvdrr, c
  logical psflag
  !
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: vu, o1, o2, o3, vief
  !
  ! setting up constants which might be used
  include 'alfinv.h'
  alpha = rel / c
  aa = alpha * alpha
  a2 = 0.5d0 * aa
  lp = l + 1
  lpx = lp
  if ( lp .gt. 4 ) lpx = 4
  lp2 = l + l + 1
  if ( lp2 .gt. 7 ) lp2 = 7
  zef = zorig
  if ( ( njrc( lpx ) .ne. 0 ) .and. ( psflag ) ) zef = 0
  zaa = zef * aa
  za2 = zef * a2
  ! we carry on only if idoflag is not zero
  if ( idoflag .eq. 0 ) return
  d1 = 0.5d0 / dl
  d2 = 1.d0 / ( dl * dl )
  o1( : ) = 1.0d0 / r( : ); o2( : ) = o1( : ) ** 2; o3( : ) = o1( : ) * o2( : )
  ! below, first fork = full potential, second fork = pseudopotential
  if ( ( njrc( lpx ) .eq. 0 ) .or. ( .not. psflag ) ) then
     if ( idoflag .eq. 1 ) v( : ) = - zef * o1( : ) + orb( : )
     vu( : ) = orb( : )
  else
     if ( idoflag .eq. 1 ) then
        lu = 0
        ij = 2.1 * ( abs( xj ) - dble( l ) )
        if ( l .eq. 0 ) then
           lu = 1
        else
           if ( ij .lt. 0 ) lu = 2 * l
           if ( ij .gt. 0 ) lu = 2 * l + 1
        end if
        if ( lu .gt. 0 ) then
           vief( : ) = vi( :, lu )
!          write ( 6, * ) 'SQMM', lu, vief( 5 )
        else
           l1 = 2 * l
           l2 = 2 * l + 1
           w1 = dble( l ) / dble( 2 * l + 1 )
           w2 = dble( l + 1 ) / dble( 2 * l + 1 )
           vief( : ) = w1 * vi( :, l1 ) + w2 * vi( :, l2 )
!          write ( 6, * ) 'SQMM', l1, l2, w1, w2, vi( 5, l1 ), vi( 5, l2 ), vief( 5 )
        end if
        v( : ) = vief( : ) + orb( : )
     end if
     vu( : ) = v( : )
  end if
  ! following indept of full potential versus pseudopential
  do j = 3, nr - 2
     dvdl = ( 8.d0 * ( vu( j + 1 ) - vu( j - 1 ) ) - ( vu( j + 2 ) - vu( j - 2 ) ) ) / ( 12.d0 * dl )
     ddvdll = ( 16.d0 * ( vu( j + 1 ) + vu( j - 1 ) ) - ( vu( j + 2 ) + vu( j - 2 ) ) - 30.d0 * vu( j ) ) / ( 12.d0 * dl * dl )
     dvdr = dvdl * o1( j )
     ddvdrr = ( ddvdll - dvdl ) * o2( j )
     xm1( j ) = -a2 * dvdr - za2 * o2( j )
     xm2( j ) = -a2 * ddvdrr + zaa * o3( j )
  end do
  ! inner end point
  xm1( 1 ) = xm1( 3 ) + za2 * ( o2( 3 ) - o2( 1 ) )
  xm2( 1 ) = xm2( 3 ) - zaa * ( o3( 3 ) - o3( 1 ) )
  xm1( 2 ) = xm1( 3 ) + za2 * ( o2( 3 ) - o2( 2 ) )
  xm2( 2 ) = xm2( 3 ) - zaa * ( o3( 3 ) - o3( 2 ) )
  ! outer end point
  xm1( nr - 1 ) = xm1( nr - 2 )
  xm2( nr - 1 ) = xm2( nr - 2 )
  xm1( nr ) = xm1( nr - 1 )
  xm2( nr ) = xm2( nr - 1 )
  !
  return
end subroutine setqmm
