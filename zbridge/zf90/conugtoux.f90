 program conugtoux
   implicit none
   !
   integer, parameter :: u1dat = 22, iwf = 23
   integer :: nk, nb, nx( 3 ), ng, i, nwfile, iwfile
   integer, allocatable :: g( :, : ), ikl( :, : )
   real( kind = kind( 1.0d0 ) ), allocatable :: zr( :, : ), zi( :, : )
   character * 20, allocatable :: wfnam( : )
   !
   open( unit=99, file='kandb.h', form='formatted', status='unknown' )
   call igetval( nk, 'nk   ' )
   call igetval( nb, 'nb   '  )
   call igetval( nx( 1 ), 'ngx  ' )
   call igetval( nx( 2 ), 'ngy  ' )
   call igetval( nx( 3 ), 'ngz  ' )
   close( unit=99 )
   !
   open( unit=99, file='masterwfile', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) nwfile
   close( unit=99 )
   allocate( wfnam( nwfile ), ikl( 2, nwfile ) )
   open( unit=99, file='listwfile', form='formatted', status='unknown' )
   rewind 99
   do i = 1, nwfile
      read ( 99, * ) ikl( 1, i ), wfnam( i )
      if ( i .gt. 1 ) ikl( 2, i - 1 ) = ikl( 1, i ) - 1
   end do
   close( unit=99 )
   ikl( 2, nwfile ) = nk
   !
   open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown' )
   rewind u1dat
   iwfile = 1
   do i = 1, nk
      open( unit=99, file='prog', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(2x,1a6,3x,2i5)' ) 'conv03g', i,nk
      close( unit=99 )
      if ( i .eq. ikl( 1, iwfile ) ) then
         open ( unit=iwf, file=wfnam( iwfile ), form='unformatted',status='unknown' )
         rewind iwf
         read ( iwf ) ng
         allocate( g( ng, 3 ), zr( ng, nb ), zi( ng, nb ) )
         read ( iwf ) g
      end if
      read ( iwf ) zr
      read ( iwf ) zi
      call gentoreal( nx, nb, zr, zi, ng, g, u1dat, i .eq. 1 )
      if ( i .eq. ikl( 2, iwfile ) ) then
         deallocate( g, zr, zi )
         iwfile = iwfile + 1
      end if
   end do
   close( unit=u1dat )
   !
   deallocate( wfnam, ikl )
   !
 end program conugtoux
