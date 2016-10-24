subroutine escanner( ii, psiat, nel, orb, nl, xnj, zorig, rel, nr, r, r2, dl, njrc, vi, psflag, de, eh, rcl, rch )
  implicit none
  !
  integer :: ii, njrc( 4 ), nel, nr
  real( kind = kind( 1.0d0 ) ) :: de, eh, rcl, rch, rel
  integer :: nl( nel )
  real( kind = kind( 1.0d0 ) ) :: orb( nr, nel )
  real( kind = kind( 1.0d0 ) ) :: vi( nr, 7 ), zorig
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr ), psiat( nr ), dl
  real( kind = kind( 1.0d0 ) ) :: xnj( nel ) 
  logical :: psflag
  !
  integer :: istop, i, n, nn, ief
  real( kind = kind( 1.0d0 ) ) :: zeff,kappa, plead, ee, k, psi, psip, xdummy, mel, norm
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: xm1, xm2, v, phi, f
  !
  istop = 1
  do while ( rch .gt. r( istop ) )
     istop = istop + 1
  end do
  n = 2 * nr
  nn = 2 * nr
  call setqmm(ii,orb(1,ii),nl(ii),xnj(ii),1,v,zeff, zorig,rel,nr,r,r2,dl,xm1,xm2,njrc,vi,psflag)
  call getkap( xnj( ii ), nl( ii ), kappa )
  call getplead( nl( ii ), xnj( ii ), rel, kappa, plead, zeff )
  ee = 0.5d0 * dble( de )
  open( unit=98, file='edep', form='formatted', status='unknown' )  
  rewind 98
  do while ( ee .lt. eh )
     k = sqrt( 2.0d0 * ee )
     phi( : ) = 0.0d0
     call integ( ee, nl( ii ), kappa, n, nn, istop, ief, xdummy, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
     open( unit=99, file='nck', form='formatted', status='unknown' )  
     rewind 99
     do i = 3, istop - 2
        if ( ( r( i - 2 ) .gt.  rcl ) .and. ( r( i + 2 ) .lt. rch ) ) then
           psi = phi( i )
           psip = ( 8.0d0 * ( phi( i + 1 ) - phi( i - 1 ) ) - ( phi( i + 2 ) - phi( i - 2 ) ) ) / ( 12.0d0 * dl * r( i ) )
           norm = psi ** 2 + ( psip / k ) ** 2
           write ( 99, '(5(1x,1e15.8))' ) r( i ), psi, psip, k, norm
        end if
     end do
     close( unit=99 )
     phi = phi / sqrt( norm )
     do i = 1, nr
        f( i ) = psiat( i ) * phi( i ) / r( i ) ! r for matrix element, divide by r squared
     end do
     call bintegrate( nr, r, dl, f, mel, nr - 5 )
     write ( 98, * ) ee, mel
     ee = ee + de
  end do
  close( unit=98 )
  ! 
  return
end subroutine escanner
