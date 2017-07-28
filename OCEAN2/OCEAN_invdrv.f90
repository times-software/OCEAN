! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine OCEAN_invdrv( x, b, n, i1, i2, n2, need, iwrk, w, bv, av, bs, as, req, ct, ev, f )
  use AI_kinds, only : DP
  implicit none
  !
  integer :: n, i1, i2, n2, need, iwrk
  complex(DP) :: x( n ), b( n ), bv( n ), av( n )
  character(len=3) :: req, bs, as
  character(len=5) :: ev
  character(len=9) :: ct
  complex(DP) :: w( 0 : iwrk - 1 )
  real(DP) :: f( 2 )
  !
  integer :: ax, g, pg, apg, u, au, c
  !
  ax =   0
  g =  ax + n
  pg =   g + n
  apg =  pg + n
  u = apg + n
  au =   u + n * n2
  c =  au + n * n2
  need =   c + n2 * n2
  !
  if ( iwrk .lt. need ) then
     if ( ct .eq. 'beginning' ) then
        req = 'all'
     else
        stop 'invdrv: late adjustment of iwrk sensed'
     end if
     return
  end if
  call backend( av, as, bv, bs, x, b, w( ax ), w( g ), w( pg ), w( apg ), w( u ), w( au ), & 
                w( c ), n, i1, i2, n2, ev, ct, req, f )
  !
  return
end subroutine OCEAN_invdrv
!
subroutine backend( aftvec, as, befvec, bs, x, b, ax, g, pg, apg, u, au, c, n, i1, i2, n2, &
                    ev, ct, req, f )
  use AI_kinds, only : DP
  implicit none
  !
  integer n, i1, i2, n2
  real(DP) :: f( 2 )
  character(len=3) req, as, bs
  character(len=5) ev
  character(len=9) ct
  complex(DP) :: x( n ), b( n )
  complex(DP) :: aftvec( n ), befvec( n )
  complex(DP) :: ax( n ), g( n ), pg( n ), apg( n )
  complex(DP) :: u( n, n2 ), au( n, n2 ), c( n2, n2 )
  !
  complex(DP), external :: hilbcd
  !
  select case( as )
  case ( 'ax_' )
     ax = aftvec
  case ( 'pg_' )
     pg = aftvec
  case ( 'apg' )
     apg = aftvec
  end select
  req = '---'
  do while ( req .eq. '---' )
     select case( ct )
     case ( 'beginning' )
        i1 = 0
        select case ( ev )
        case ( 'zerox' )
           x = 0
        case ( 'loadx' )
           call loadx( n, x )
        case ( 'havex' )
           write ( 6, * ) 'x is assumed'
        end select
        ct = 'newi2loop'
     case ( 'newi2loop' )
        i2 = 0
        req= 'act'
        bs = 'x__'
        as = 'ax_'
        ct = 'checkconv'
     case ( 'checkconv' )
        g = ax - b
        f( 2 ) = real( hilbcd( n, g, g ), DP )
!        write( 6, * ) f( 2 )
        if ( f( 2 ) .lt. f( 1 ) ) then
           req = 'end'
        else
           req = 'prc'
           bs = 'g__'
           as = 'pg_'
           ct = 'havenewpg'
        end if
     case ( 'havenewpg' )
        req = 'act'
        bs = 'pg_'
        as = 'apg'
        ct = 'runupdate' 
     case ( 'runupdate' )
        call update( x, ax, g, pg, apg, u, au, c, n, i1, i2, n2 )
!        call savex( n, x )
        if ( i2 .eq. n2 ) then
           ct = 'newi2loop'
        else
           ct = 'checkconv'
        end if
     end select
  end do
  if ( ( req .eq. 'prc' ) .or. ( req .eq. 'act' ) ) then
     select case( bs )
     case ( 'x__' )
        befvec = x
     case ( 'g__' )
        befvec = g
     case ( 'pg_' )
        befvec = pg
     end select
  end if
  !
  return
end subroutine backend
!
! here we try to do some GMRES-like miminization
subroutine update( x, ax, g, pg, apg, u, au, c, n, i1, i2, n2 )
  use AI_kinds, only : DP
  implicit none
  !
  integer :: i1, i2, n2, n
  complex(DP) :: x( n ), ax( n ), g( n ), pg( n ), apg( n )
  complex(DP) :: u( n, n2 ), au( n, n2 )
  complex(DP) :: c( n2, n2 )
  !
  integer :: i, ip, nstart, nstop, info, i3
  integer, allocatable :: ipiv( : )
  complex(DP) :: coeff(i2+1)
  complex(DP), allocatable :: c2( : , : ), cinv( : , : ), r( : ), cwork( : )
  complex(DP) :: dumc(i2+1)
  complex(DP), external :: hilbcd
  complex(DP), parameter :: cm1 = -1.0_dp
  complex(DP), parameter :: cp1 = 1.0_dp
  complex(DP), parameter :: c0  = 0.0_dp
  !
!  call sizereport( 16 * n2 ** 2, 'c2........' ); allocate( c2( n2, n2 ) )
!  call sizereport( 16 * n2 ** 2, 'cinv......' ); allocate( cinv( n2, n2 ) )
!  call sizereport( 16 * n2, 'r.........' ); allocate( r( n2 ) )
  allocate( c2( n2, n2 ), cinv( n2, n2 ), r( n2 ), ipiv( n2 ) )
  i1 = i1 + 1
  i2 = i2 + 1
  u( : , i2 ) = pg
  au( : , i2 ) = apg
if( 0 ) then
  r = 0
  dumc=0

!!$OMP PARALLEL DO &
!!$OMP DEFAULT( NONE ) &
!!$OMP SCHEDULE( STATIC ) &
!!$OMP FIRSTPRIVATE( g, i2, n ) &
!!$OMP PRIVATE( ip, nstart, nstop ) &
!!$OMP SHARED( au, c ) &
!!$OMP REDUCTION(+:r,dumc)
do nstart = 1, n, 256
  nstop = min( nstart+255,n)

  do ip = 1, i2
    r( ip ) = r( ip ) - dot_product( au( nstart:nstop, ip ), g(nstart:nstop) )
    dumc(ip) = dumc(ip) + dot_product( au(nstart:nstop, i2), au( nstart:nstop, ip ) )
  enddo
enddo
!!$OMP END PARALLEL DO

  c( 1:i2, i2 ) = conjg( dumc(:) )
  c( i2, 1:i2 ) = dumc(:)

  call invert( i2, n2, c, c2, cinv )

  coeff(1:i2) = matmul( cinv(1:i2,1:i2), r(1:i2) )


!!$OMP PARALLEL DO &
!!$OMP DEFAULT( NONE ) &
!!$OMP SCHEDULE( STATIC ) &
!!$OMP FIRSTPRIVATE( coeff, n, i2 ) &
!!$OMP PRIVATE( nstart, nstop, i ) &
!!$OMP SHARED( x, ax, u, au )
  do nstart = 1, n, 256
    nstop = min( nstart+255,n)
    do i = 1, i2
       x(nstart:nstop) = x(nstart:nstop) + coeff(i) * u( nstart:nstop , i )
       ax(nstart:nstop) = ax(nstart:nstop) + coeff(i) * au( nstart:nstop , i )
    end do
  enddo
!!$OMP END PARALLEL DO


else
  call ZGEMV( 'C', n, i2, cm1, au, n, g, 1, c0, r, 1 )


  if( 0 ) then
    call ZGEMV( 'C', n, i2, cp1, au, n, au(:,i2), 1, c0, dumc, 1 )
    c( 1:i2, i2 ) = dumc(:)
    c( i2, 1:i2 ) = conjg( dumc(:) )
    call invert( i2, n2, c, c2, cinv )
    coeff(1:i2) = matmul( cinv(1:i2,1:i2), r(1:i2) )
  else
!    call ZGEMM( 'C', 'N', i2, i2, n, cp1, au, n, au, n, c0, c, n2 )
    call ZGEMV( 'C', n, i2, cp1, au, n, au(:,i2), 1, c0, dumc, 1 )
    c( 1:i2, i2 ) = dumc(:)
    c( i2, 1:i2 ) = conjg( dumc(:) )
    c2( 1:i2, 1:i2 ) = c( 1:i2, 1:i2 )

    if( 0 ) then
      call ZGETRF( i2, i2, c, n2, ipiv, info )
      call ZGETRS( 'N', i2, 1, c, n2, ipiv, r, n2, info )
    else
      i3 = 64*i2
      allocate( cwork( i3 ) )
      call ZHETRF( 'L', i2, c, n2, ipiv, cwork, i3, info )
      call ZHETRS( 'L', i2, 1, c, n2, ipiv, r, n2, info )
      deallocate( cwork )
    endif
    coeff(1:i2) = r(1:i2)
    c( 1:i2, 1:i2 ) = c2( 1:i2, 1:i2 )
  endif

  call ZGEMV( 'N', n, i2, cp1, u, n, coeff, 1, cp1, x, 1 )
  call ZGEMV( 'N', n, i2, cp1, au, n, coeff, 1, cp1, ax, 1 )

endif





!  call sizereport( 0, 'c2........' ); deallocate( c2 ) 
!  call sizereport( 0, 'cinv......' ); deallocate( cinv ) 
!  call sizereport( 0, 'r.........' ); deallocate( r ) 
  deallocate( c2, cinv, r, ipiv )
  return
end subroutine update
!
function hilbcd( n, l, r )
  use AI_kinds, only : DP
  implicit none
  !
  integer, intent(in) :: n 
!, i
  complex(DP), intent(in) ::  l( n ), r( n )
  complex(DP) :: hilbcd
  !
!  hilbcd = 0
!  do i = 1, n
!     hilbcd = hilbcd + conjg( l( i ) ) * r( i )
!  end do
  hilbcd = dot_product( l, r )
  !
  return
end function hilbcd
!
subroutine loadx( n, x )
  use AI_kinds, only : DP
  implicit none
  !
  integer n
  complex(DP) :: x( n )
  !
  open( unit=99, file='xfile', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) x
  close( unit=99 )
  !
  return
end subroutine loadx
!
subroutine savex( n, x )
  use AI_kinds, only : DP
  implicit none
  !
  integer n
  complex(DP) :: x( n )
  !
  open( unit=99, file='xfile', form='unformatted', status='unknown' )
  rewind 99
  write ( 99 ) x
  close( unit=99 )
  !
  return
end subroutine savex
