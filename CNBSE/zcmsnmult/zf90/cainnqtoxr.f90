subroutine cainnqtoxr( nx, ny, nz, nfft, nq, nbd, vecr, veci, xwrkr, xwrki, tau, nn1, zn, wrk, jfft, ur, ui )
  implicit none
  !
  integer :: nx, ny, nz, nfft, nq, nbd, nn1, zn( 3 ), jfft
  real( kind = kind( 1.0d0 ) ) :: tau( 3 ), wrk( jfft )
  double precision, dimension( nbd, nq ) :: vecr, veci
  double precision, dimension( nfft, nx, ny, nz ) :: xwrkr, xwrki
  double precision, dimension( nx, ny, nz, nbd, nq ) :: ur, ui
  !
  integer ix, iy, iz, iq1, iq2, iq3, iq, ibd, ii
  !
  xwrkr( :, :, :, : ) = 0
  xwrki( :, :, :, : ) = 0
  iq = 0
  do iq1 = 1, zn( 1 )
     do iq2 = 1, zn( 2 )
        do iq3 = 1, zn( 3 )
           iq = iq + 1
           ii = 1 + ( iq1 - 1 ) + zn( 1 ) * ( ( iq2 - 1 ) + zn( 2 ) * ( iq3 - 1 ) )
           do ibd = 1, nbd
              xwrkr( ii, :, :, : ) = xwrkr( ii, :, :, : ) + vecr( ibd, iq ) * ur( :, :, :, ibd, iq ) &
                   - veci( ibd, iq ) * ui( :, :, :, ibd, iq )
              xwrki( ii, :, :, : ) = xwrki( ii, :, :, : ) + vecr( ibd, iq ) * ui( :, :, :, ibd, iq ) &
                   + veci( ibd, iq ) * ur( :, :, :, ibd, iq )
           end do
        end do
     end do
  end do
  !
  do iz = 1, nz
     do iy = 1, ny
        do ix = 1, nx
           call cfft( xwrkr( 1, ix, iy, iz ), xwrki( 1, ix, iy, iz ), nn1, zn( 1 ), zn( 2 ), zn( 3 ), -1, wrk, jfft )
        end do
     end do
  end do
  !
  return
end subroutine cainnqtoxr
