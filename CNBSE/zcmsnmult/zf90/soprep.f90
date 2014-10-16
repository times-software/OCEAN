subroutine soprep( nc, vms, cml, cms, lc, xi, somel )
  implicit none
  !
  integer :: nc, lc
  real( kind = kind( 1.0d0 ) ) :: xi
  real( kind = kind( 1.0d0 ) ) :: cms( nc ), cml( nc ), vms( nc )
  real( kind = kind( 1.0d0 ) ) :: somel( nc, nc, 2 )
  !
  integer :: ic, jc, i
  complex( kind = kind( 1.0d0 ) ) :: ctmp, rm1, vrslt( 3 )
  complex( kind = kind( 1.0d0 ) ), external :: jimel
  real( kind = kind( 1.d0 ) )  :: life_time( 2 ), l_alpha, l_beta, delta_so( 2 )
  logical :: broaden_exist
  !
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 ) 
  include 'sphsetnx.h.f90'
  !
  include 'sphsetx.h.f90'
  call newgetprefs( prefs, lc, nsphpt, wsph, xsph, ysph, zsph )
  !
  ! 
  inquire( file='core_broaden.ipt', exist= broaden_exist )
  if( broaden_exist ) then
    write(6,*) 'spin-orbit dependent broadening'
    open( unit=99, file='core_broaden.ipt', form='formatted', status='old' )
    read( 99, * ) life_time( : )
    close( 99 )
    life_time( : ) = life_time( : ) / 27.2114d0
 
    delta_so( 1 ) = 0.5d0 * ( ( real( lc ) + 0.5 ) * ( real( lc ) + 1.5 ) - real( lc ) * real( lc + 1 ) - 0.75d0 )
    delta_so( 2 ) = 0.5d0 * ( ( real( lc ) - 0.5 ) * ( real( lc ) + 0.5 ) - real( lc ) * real( lc + 1 ) - 0.75d0 )
    !
    l_alpha = ( life_time( 2 ) - life_time( 1 ) ) / ( - xi * ( delta_so( 2 ) - delta_so( 1 ) ) )
    l_beta = ( life_time( 2 ) * delta_so( 1 ) - life_time( 1 ) * delta_so( 2 ) ) / ( delta_so( 1 ) - delta_so( 2 ) )
    write(6,*) l_alpha, l_beta
  else
    l_alpha = 0.d0
    l_beta  = 0.d0
  endif

  somel( :, :, : ) = 0.0d0
  rm1 = -1; rm1 = sqrt( rm1 )
  do ic = 1, nc
     do jc = 1, nc
        if ( vms( ic ) .eq. vms( jc ) ) then
           call limel( lc, nint( cml( jc ) ), nint( cml( ic ) ), vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
           ctmp = 0
           do i = 1, 3
              ctmp = ctmp + vrslt( i ) * jimel( 0.5d0, cms( jc ), cms( ic ), i )
           end do
           write(20,*) ic, jc, real( ctmp )
!           write(21,*) 
           ctmp = -xi * ctmp
!           somel( ic, jc, 1 ) = ctmp 
           somel( ic, jc, 1 ) = real( ctmp ) - aimag( ctmp ) * l_alpha
           somel( ic, jc, 2 ) = aimag( ctmp ) + real( ctmp ) * l_alpha 
!           somel( ic, jc, 2 ) = -rm1 * ctmp + real( ctmp ) * l_alpha 
           if( ic .eq. jc ) then
             somel( ic, jc, 2 ) = somel( ic, jc, 2 ) + l_beta !* real( ctmp )
!             write( 6, * ) ctmp, l_alpha, l_beta
!             somel( ic, jc, 1 )  = somel( ic, jc, 1 ) 
           endif
        end if
     end do
  end do
  !
  return
end subroutine soprep
