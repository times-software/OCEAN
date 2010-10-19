subroutine vcbder( ibl, ibh, vlev, vhev, clev, chev )
  implicit none
  !
  integer :: ibl, ibh
  real( kind = kind( 1.0d0 ) ) :: vlev, vhev, clev, chev
  !
  integer :: i, j, ii
  real( kind = kind( 1.0d0 ) ) :: efryd, efev, dum, ee( 2 )
  logical :: vhit, chit
  !
  open( unit=99, file='fermi.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) efryd
  read ( 99, * ) efev
  close( unit=99 )
  if ( abs( efryd * 13.6057d0 - efev ) .gt. 0.01d0 ) stop 'bad match vcbder'
  vhit = .false.
  chit = .false.
  open( unit=99, file='bbds', form='formatted', status='unknown' )
  rewind 99
  do i = 1, ibh
     read ( 99, * ) ii, dum, dum, ee( : )
     do j = 1, 2
        if ( ( ii .ge. ibl ) .and. ( ii .le. ibh ) ) then
           if ( ee( j ) .lt. efev ) then
              if ( vhit ) then
                 vlev = min( vlev, ee( j ) )
                 vhev = max( vhev, ee( j ) )
              else
                 vlev = ee( j )
                 vhev = ee( j )
              end if
              vhit = .true.
           else
              if ( chit ) then
                 clev = min( clev, ee( j ) )
                 chev = max( chev, ee( j ) )
              else
                 clev = ee( j )
                 chev = ee( j )
              end if
              chit = .true.
           end if
        end if
     end do
  end do
  close( unit=99 )
  !
  return
end subroutine vcbder
