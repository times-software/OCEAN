program formreytab
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  !
  integer :: i, npt, nptvu, indx
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: xpt, ypt, zpt, wpt
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: xptvu, yptvu, zptvu, wptvu
  !
  integer :: ltab( 9 ), mtab( 9 )
  complex( kind = kind( 1.0d0 ) ) :: rm1, c( 9, 9 )
  !
  real( kind = kind( 1.0d0 ) ) :: rey, imy, rpos, rneg
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  call getprefs( prefs )
  !
  open( unit=99, file='specpnt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) npt
  allocate( xpt( npt ), ypt( npt ), zpt( npt ), wpt( npt ) )
  do i = 1, npt
     read ( 99, * ) xpt( i ), ypt( i ), zpt( i ), wpt( i )
  end do
  close( unit=99 )
  !
  open( unit=99, file='specvu', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nptvu
  allocate( xptvu( nptvu ), yptvu( nptvu ), zptvu( nptvu ), wptvu( nptvu ) )
  do i = 1, nptvu
     read ( 99, * ) xptvu( i ), yptvu( i ), zptvu( i ), wptvu( i )
  end do
  close( unit=99 )
  !
  ltab( 1 ) = 0; ltab( 2 : 4 ) = 1; ltab( 5 : 9 ) = 2
  mtab( 1 ) = 0; mtab( 2 ) = -1; mtab( 3 ) = 0; mtab( 4 ) = +1
  mtab( 5 ) = -2; mtab( 6 ) = -1; mtab( 7 ) = 0; mtab( 8 ) = +1; mtab( 9 ) = +2
  !
  c( :, : ) =  0.0d0
  c( 1, 1 ) = +1.0d0
  c( 3, 2 ) = +1.0d0
  c( 2, 3 ) = +1.0d0 / sqrt( 2.0d0 )
  c( 4, 3 ) = -1.0d0 / sqrt( 2.0d0 )
  c( 2, 4 ) = +rm1 / sqrt( 2.0d0 )
  c( 4, 4 ) = +rm1 / sqrt( 2.0d0 )
  c( 7, 5 ) = +1.0d0
  c( 5, 6 ) = +1.0d0 / sqrt( 2.0d0 )
  c( 9, 6 ) = +1.0d0 / sqrt( 2.0d0 )
  c( 5, 7 ) = +rm1 / sqrt( 2.0d0 )
  c( 9, 7 ) = -rm1 / sqrt( 2.0d0 )
  c( 6, 8 ) = +1.0d0 / sqrt( 2.0d0 )
  c( 8, 8 ) = -1.0d0 / sqrt( 2.0d0 )
  c( 6, 9 ) = +rm1 / sqrt( 2.0d0 )
  c( 8, 9 ) = +rm1 / sqrt( 2.0d0 )
  !
  do indx = 1, 9
     do i = 1, npt
        call getrey( rey, imy, c, indx, ltab, mtab, xpt( i ), ypt( i ), zpt( i ), prefs )
        write ( 6, '(2(1x,1e15.8))' ) rey, imy
     end do
     open( unit=20+indx, form='formatted', status='unknown' )
     rewind 20 + indx
     do i = 1, nptvu
        call getrey( rey, imy, c, indx, ltab, mtab, xptvu( i ), yptvu( i ), zptvu( i ), prefs )
        rpos = max( 0.0d0, rey ); rneg = abs( min( 0.0d0, rey ) )
        write ( 20+indx, '(6(1x,1f10.5))' ) xptvu( i ), yptvu( i ), zptvu( i ), rpos, rneg, imy
     end do
     close( unit=20+indx )
  end do
  !
end program formreytab
  !
subroutine getrey( rey, imy, c, indx, ltab, mtab, x, y, z, prefs )
  implicit none
  !
  integer :: indx, ltab( 9 ), mtab( 9 )
  real( kind = kind( 1.0d0 ) ) :: rey, imy, x, y, z, prefs( 0 : 1000 )
  complex( kind = kind( 1.0d0 ) ) :: c( 9, 9 )
  !
  integer :: i
  complex( kind = kind( 1.0d0 ) ) :: rm1, ylm
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  rey = 0.0d0
  imy = 0.0d0
  do i = 1, 9
     call ylmeval( ltab( i ), mtab( i ), x, y, z, ylm, prefs )
     rey = rey + c( i, indx ) * ylm
     write ( 40 + indx, '(8(1x,1f10.5))' ) x, y, z, ylm * c( i, indx ), ylm, c( i, indx )
     imy = imy - rm1 * c( i, indx ) * ylm
  end do
  !
  return
end subroutine getrey
