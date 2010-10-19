program condens
  implicit none
  !
  integer :: ng, i, nx( 3 ), nxprod
  real( kind = kind( 1.0d0 ) ) :: vol
  integer, allocatable :: gvec( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: rhogr( : ), rhogi( : ), rhocell( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: comprho( : )
  !
  real( kind = kind( 1.0d0 ) ) :: n, su, ux1, ux2, uc1, uc2, kxc, dn
  real( kind = kind( 1.0d0 ) ) :: valp2, valp1, valz0, valm1, valm2, np2, np1, nz0, nm1, nm2
  real( kind = kind( 1.0d0 ) ) :: exp2, exp1, exz0, exm1, exm2, ecp2, ecp1, ecz0, ecm1, ecm2
  !
  read ( 5, * ) ng
  allocate( gvec( ng, 3 ), rhogr( ng ), rhogi( ng ) )
  do i = 1, ng
     read ( 5, * ) gvec( i, : ), rhogr( i ), rhogi( i )
  end do
  open( unit=99, file='omega.h', form='formatted', status='unknown' )
  rewind 99
  call rgetval( vol, 'volbo' )
  close( unit=99 )  
  rhogr( : ) = rhogr( : ) / vol
  rhogi( : ) = rhogi( : ) / vol
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  rewind 99
  call igetval( nx( 1 ), 'ngx  ' )
  call igetval( nx( 2 ), 'ngy  ' )
  call igetval( nx( 3 ), 'ngz  ' )
  close( unit=99 )
  write ( 6, '(1a4,1i6)' ) 'ng= ', ng
  write ( 6, '(1a6,3i10)' ) 'nx(:)=', nx(:)
  write ( 6, '(1a15,3i6)' ) 'gvec( 1, : ) = ', gvec( 1, : )
  write ( 6, '(10x,1a12,2(1x,1e15.8))' ) 'rhog( 1 ) = ', rhogr( 1 ), rhogi( 1 )
  write ( 6, '(1a15,3i6)' ) 'gvec( 2, : ) = ', gvec( 2, : )
  write ( 6, '(10x,1a12,2(1x,1e15.8))' ) 'rhog( 2 ) = ', rhogr( 2 ), rhogi( 2 )
  open( unit=99, file='rho.xpts', form='unformatted', status='unknown' )
  rewind 99
  call gentoreal( nx, 1, rhogr, rhogi, ng, gvec, 99, .true. )
  deallocate( gvec, rhogr, rhogi )
  nxprod = nx( 1 ) * nx( 2 ) * nx( 3 )
  allocate( rhocell( nxprod ), comprho( nxprod ) )
  rewind 99
  read ( 99 ) comprho( : )
  close( unit=99 )
  do i = 1, nxprod
     rhocell( i ) = comprho( i )
  end do
  deallocate( comprho )
  call brcapper( nxprod, rhocell )
  !
  su = 0
  open( unit=99, file='rhokxc.txt', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nxprod
     n = rhocell( i )
     su = su + n
     dn = 0.01d0 * n
     np2 = n + 2 * dn
     np1 = n + dn
     nz0 = n
     nm1 = n - dn
     nm2 = n - 2 * dn
     call cacorr( np2, exp2, ecp2, ux1, ux2, uc1, uc2 )
     call cacorr( np1, exp1, ecp1, ux1, ux2, uc1, uc2 )
     call cacorr( nz0, exz0, ecz0, ux1, ux2, uc1, uc2 )
     call cacorr( nm1, exm1, ecm1, ux1, ux2, uc1, uc2 )
     call cacorr( nm2, exm2, ecm2, ux1, ux2, uc1, uc2 )
     valp2 = np2 * ( exp2 + ecp2 )
     valp1 = np1 * ( exp1 + ecp1 )
     valz0 = nz0 * ( exz0 + ecz0 )
     valm1 = nm1 * ( exm1 + ecm1 )
     valm2 = nm2 * ( exm2 + ecm2 )
     kxc = ( 16.0d0 * ( valm1 + valp1 ) - ( valm2 + valp2 ) - 30.0d0 * valz0 ) / ( 12.0d0 * dn ** 2 )
     write ( 99, '(1p,2(1x,1e15.8))' ) n, kxc
  end do
  close( unit=99 )
  su = su / dble( nxprod ) * vol
  write ( 6, '(1a11,1x,1f20.6)' ) 'checksum = ', su
  deallocate( rhocell )
  !
end program condens
