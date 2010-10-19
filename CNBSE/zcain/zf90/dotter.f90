program dotter
  implicit none
  !
  integer :: nptot, ntot, mc, i, lc, idum
  real( kind = kind( 1.0d0 ) ) :: rr, ri, ir, ii, br, bi, tau( 3 )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: mer, mei
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: pcr, pci
  character * 6 :: str
  character * 7 :: str2
  character * 9 :: filnam
  !
  read ( 5, * ) filnam
  open( unit=99, file=filnam, form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) nptot, ntot
  read ( 99 ) tau( : )
  allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
  read ( 99 ) pcr
  read ( 99 ) pci
  close( unit=99 )
  !
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) idum, idum, lc
  close( unit=99 )
  !
  allocate( mer( nptot ), mei( nptot ) )
  open( unit=99, file='mels', form='formatted', status='unknown' )
  rewind 99
  do mc = -lc, lc
     do i = 1, nptot
        read ( 99, * ) mer( i ), mei( i )
     end do
     write ( str, '(1a4,1i2.2)' ) 'beff', 1 + mc + lc
     write ( str2, '(1a5,1i2.2)' ) 'abeff', 1 + mc + lc
     open( unit=98, file=str, form='unformatted', status='unknown' )
     open( unit=97, file=str2, form='formatted', status='unknown' )
     rewind 98
     write ( 98 ) tau( : )
     rewind 97
     write ( 97, '(3(1x,1e15.8))' ) tau( : )
     do i = 1, ntot
        rr = dot_product( mer( : ), pcr( :, i ) )
        ri = dot_product( mer( : ), pci( :, i ) )
        ir = dot_product( mei( : ), pcr( :, i ) )
        ii = dot_product( mei( : ), pci( :, i ) )
        br = rr - ii
        bi = -ri - ir
        write ( 98 ) br, bi
        write ( 97, '(2(1x,1e15.8))' ) br, bi
     end do
  end do
  !
end program dotter
