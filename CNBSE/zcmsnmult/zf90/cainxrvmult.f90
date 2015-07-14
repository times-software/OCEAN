! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cainxrvmult( nx, ny, nz, zn, tau, amet, epsi, xwrkr, xwrki, ptab )
  implicit none
  !
  integer :: nx, ny, nz, zn( 3 )
  real( kind = kind( 1.0d0 ) ) :: epsi, tau( 3 ), amet( 3, 3 ), ptab( 100 )
  real( kind = kind( 1.0d0 ) ), dimension( zn( 1 ) * zn( 2 ) * zn( 3 ), nx * ny * nz ) :: xwrkr, xwrki
  !
  integer i, ix, iy, iz, j, k1, k2, k3, kk1, kk2, kk3, ii
  real( kind = kind( 1.0d0 ) ) :: fr( 3 ), xk( 3 ), alf( 3 ), r, frac, potn
  !
  i = 0
  do iz = 1, nz
     fr( 3 ) = dble( iz - 1 ) / dble( nz )
     do iy = 1, ny
        fr( 2 ) = dble( iy - 1 ) / dble( ny )
        do ix = 1, nx
           fr( 1 ) = dble( ix - 1 ) / dble( nx )
           i = i + 1
           j = 0
           do k3 = 1, zn( 3 )
              kk3 = k3 - 1
              if ( kk3 .ge. zn( 3 ) / 2 ) kk3 = kk3 - zn( 3 )
              xk( 3 ) = kk3
              do k2 = 1, zn( 2 )
                 kk2 = k2 - 1
                 if ( kk2 .ge. zn( 2 ) / 2 ) kk2 = kk2 - zn( 2 )
                 xk( 2 ) = kk2
                 do k1 = 1, zn( 1 )
                    kk1 = k1 - 1
                    if ( kk1 .ge. zn( 1 ) / 2 ) kk1 = kk1 - zn( 1 )
                    xk( 1 ) = kk1
                    j = j + 1
                    alf( : ) = xk( : ) + fr( : ) - tau( : )
                    r = sqrt( dot_product( alf, matmul( amet, alf ) ) )
                    if ( r .ge. 9.9d0 ) then
                       potn = epsi / r
                    else
                       ii = 1.0d0 + 10.0d0 * r
                       frac = 10.d0 * ( r - 0.1d0 * dble( ii - 1 ) )
                       potn = ptab( ii ) + frac * ( ptab( ii + 1 ) - ptab( ii ) )
                    end if
                    xwrkr( j, i ) = xwrkr( j, i ) * potn
                    xwrki( j, i ) = xwrki( j, i ) * potn
                 end do
              end do
           end do
        end do
     end do
  end do
  !
  return
end subroutine cainxrvmult
