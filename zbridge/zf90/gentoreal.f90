subroutine gentoreal( nx, nfcn, fcnr, fcni, ng, gvec, iu, loud )
  implicit none
  !
  integer :: nx( 3 ), nfcn, ng, iu
  integer :: gvec( ng, 3 )
  real( kind = kind( 1.0d0 ) ) :: fcnr( ng, nfcn ), fcni( ng, nfcn )
  logical :: loud
  !
  integer :: ix, iy, iz, i1, i2, i3, nfft( 3 ), idwrk, igl, igh, nmin, ii
  integer :: fac( 3 ), j, ig, i, locap( 3 ), hicap( 3 ), toreal, torecp, nftot
  real( kind = kind( 1.0d0 ) ) :: normreal, normrecp
  complex( kind = kind( 1.0d0 ) ) :: rm1
  character * 80 :: fstr
  !
  integer, parameter :: nfac = 3
  integer, allocatable :: ilist( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: zr( :, :, : ), zi( :, :, : ), wrk( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: cres( :, :, : )
  integer, external :: optim
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  fac( 1 ) = 2; fac( 2 ) = 3; fac( 3 ) = 5
  hicap( 1 ) = 20; hicap( 2 ) = 3; hicap( 3 ) = 1
  !
  allocate( ilist( max( nx( 1 ), nx( 2 ), nx( 3 ) ), 3 ) )
  do j = 1, 3
     igl = gvec( 1, j )
     igh = gvec( 1, j )
     do ig = 2, ng
        igl = min( igl, gvec( ig, j ) )
        igh = max( igh, gvec( ig, j ) )
     end do
     nmin = 1 + igh - igl
     call facpowfind( nx( j ), nfac, fac, locap ) 
     nfft( j ) = optim( nmin, nfac, fac, locap, hicap )
     fstr = '(1i1,1a1,5x,1a4,1i4,5x,1a4,1i4,5x,1a3,1i4,5x,1a5,1i4)'
     if ( loud ) write ( 6, fstr ) j, ':', 'igl=', igl, 'igh=', igh, 'nx=', nx( j ), 'nfft=', nfft( j )
     if ( igl * igh .ge. 0 ) stop 'zero not included!'
     do i = 1, nx( j )
        ilist( i, j ) = 1 + ( i - 1 ) * nfft( j ) / nx( j )
     end do
     if ( loud ) write ( 6, '(15i5)' ) ilist( 1 : nx( j ), j )
  end do
  ! 
  i = max( nfft( 1 ), nfft( 2 ), nfft( 3 ) )
  idwrk = 2 * i * ( i + 1 )
  allocate( zr( nfft( 1 ), nfft( 2 ), nfft( 3 ) ) )
  allocate( zi( nfft( 1 ), nfft( 2 ), nfft( 3 ) ) )
  allocate( wrk( idwrk ) )
  nftot = nfft( 1 ) * nfft( 2 ) * nfft( 3 ) 
  !
  ! the indices only of the output are reversed here.
  allocate( cres( nx( 3 ), nx( 2 ), nx( 1 ) ) )
  call chkfftreal( toreal, normreal, loud )
  call chkfftrecp( torecp, normrecp, loud )
  !
  do i = 1, nfcn
     zr( :, :, : ) = 0.0d0; zi( :, :, : ) = 0.0d0
     fstr = '(3(1a1,2i8,2x),1a5,2i8)'
     do ig = 1, ng
        i1 = 1 + gvec( ig, 1 )
        if ( i1 .le. 0 ) i1 = i1 + nfft( 1 )
        i2 = 1 + gvec( ig, 2 )
        if ( i2 .le. 0 ) i2 = i2 + nfft( 2 )
        i3 = 1 + gvec( ig, 3 )
        if ( i3 .le. 0 ) i3 = i3 + nfft( 3 )
        zr( i1, i2, i3 ) = fcnr( ig, i )
        zi( i1, i2, i3 ) = fcni( ig, i )
        if ( loud .and. ( ig .le. 10 ) .and. ( i .le. 3 ) ) then
           if ( loud ) write ( 6, fstr ) 'x', gvec( ig, 1 ), i1, 'y', gvec( ig, 2 ), i2, 'z', gvec( ig, 3 ), i3, 'ig, i', ig, i
        end if
     end do
     call cfft( zr, zi, nfft( 1 ), nfft( 1 ), nfft( 2 ), nfft( 3 ), toreal, wrk, idwrk )
     zr = zr / dble( nftot ) ** normreal
     zi = zi / dble( nftot ) ** normreal
     ii = 0
     fstr = '(2(1a9,3i5,5x),1a9,2(1x,1e15.8))'
     do iz = 1, nx( 3 )
        i3 = ilist( iz, 3 )
        do iy = 1, nx( 2 )
           i2 = ilist( iy, 2 )
           do ix = 1, nx( 1 )
              i1 = ilist( ix, 1 )
              ii = ii + 1
              if ( loud .and. ( ii .le. 10 ) .and. ( i .le. 3 ) ) then
                 write ( 6, fstr ) 'mesh ind.', i1, i2, i3, 'cell ind.', ix, iy, iz, 'value = ', zr( i1, i2, i3 ), zi( i1, i2, i3 )
              end if
              cres( iz, iy, ix ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
           end do
        end do
     end do 
     write ( iu ) cres
  end do
  deallocate( zr, zi, wrk, cres, ilist )
  !
  return
end subroutine gentoreal
