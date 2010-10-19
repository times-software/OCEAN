subroutine loadux( nx, ny, nz, nbd, nq, zn, ur, ui )
  implicit none
  !
  integer, parameter :: u2dat = 35
  !
  integer :: nx, ny, nz, nbd, nq, zn( 3 )
  real( kind = kind( 1.0d0 ) ), dimension( nx, ny, nz, nbd, nq ) :: ur, ui
  !
  integer :: iq, ibd, ig, idum( 3 ), ix, iy, iz, ivl, ivh, icl, ich
  integer :: iq1, iq2, iq3, dumint, ivh2
  real( kind = kind( 1.0d0 ) ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
  real( kind = kind( 1.0d0 ) ) :: su, sul, suh
  logical :: metal
  !
  sul = 1.0d0
  suh = 1.0d0
  open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ivl, ivh, icl, ich
  close( unit=99 )
!  ivh2 = ivh
!  open( unit=99, file='metal', form='formatted', status='old')
!  read( 99, * ) metal
!  close( 99 )
!  if( metal ) then
!    open( unit=36, file='ibeg.h', form='formatted', status='old' )
!  endif
  if ( nbd .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'
  open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
  rewind u2dat
  do iq = 1, nq
     open( unit=99, file='gumatprog', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(2i8)' ) iq, nq
     close( unit=99 )
!     if( metal ) then
!       read( 36, * ) dumint, ivh2
!       ivh2 = ivh2 - 1
!     endif
     do ibd = ivl, ivh
        do ig = 1, nx * ny * nz
           read ( u2dat )
        end do
     end do
     do ibd = 1, nbd 
        do ix = 1, nx
           do iy = 1, ny
              do iz = 1, nz
                 read ( u2dat ) idum( 1 : 3 ), ur( ix, iy, iz, ibd, iq ), ui( ix, iy, iz, ibd, iq ) 
              end do
           end do
        end do
        su = sum( ur( :, :, :, ibd, iq ) ** 2 + ui( :, :, :, ibd, iq ) ** 2 )
        sul = min( su, sul )
        suh = max( su, suh )
     end do
  end do
!  if( metal ) then
!    close( 36 )
!  endif
  close( unit=u2dat )
  write ( 6, '(1a16,2f20.15)' ) 'norm bounds ... ', sul, suh
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  iq = 0
  do iq1 = 1, zn( 1 )
     do iq2 = 1, zn( 2 )
        do iq3 = 1, zn( 3 )
           iq = iq + 1
           do iz = 1, nz
              do iy = 1, ny
                 do ix = 1, nx
                    phsx = 2.0d0 * pi * dble( ( ix - 1 ) * ( iq1 - 1 ) ) / dble( nx * zn( 1 ) )
                    phsy = 2.0d0 * pi * dble( ( iy - 1 ) * ( iq2 - 1 ) ) / dble( ny * zn( 2 ) )
                    phsz = 2.0d0 * pi * dble( ( iz - 1 ) * ( iq3 - 1 ) ) / dble( nz * zn( 3 ) )
                    cphs = cos( phsx + phsy + phsz )
                    sphs = sin( phsx + phsy + phsz )
                    do ibd = 1, nbd
                       psir = cphs * ur( ix, iy, iz, ibd, iq ) - sphs * ui( ix, iy, iz, ibd, iq )
                       psii = cphs * ui( ix, iy, iz, ibd, iq ) + sphs * ur( ix, iy, iz, ibd, iq )
                       ur( ix, iy, iz, ibd, iq ) = psir
                       ui( ix, iy, iz, ibd, iq ) = psii
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  !
  return
end subroutine loadux
