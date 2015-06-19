subroutine AI_kmapr( kmesh, xk, kk, ladcap )
  implicit none
  !
  integer, intent( in ) :: kmesh(3)
  integer, intent( out ) :: kk( kmesh(3), kmesh(2), kmesh(1), 3 )
  integer, intent( in ) ::  ladcap(2,3)
  real(kind=kind(1.d0)), intent( out ) :: xk( kmesh(3), kmesh(2), kmesh(1), 3 )
  !
  integer :: ikx, iky, ikz, j
  real(kind=kind(1.d0)) :: c( 3 ), pi
  !
  pi = 4.0d0 * datan( 1.0d0 )
  !
  c( : ) = 2.0d0 * pi / dble( kmesh( : ) )
  do ikx = 1, kmesh(1)
    do iky = 1, kmesh(2)
      do ikz = 1, kmesh(3)
        xk( ikz, iky, ikx, : ) = c( : ) * dble( (/ ikx-1, iky-1, ikz-1 /) ) 
        kk( ikz, iky, ikx, : ) = (/ ikx-1, iky-1, ikz-1 /)
        ! mapping k points onto (approximately) zero-centered grid
        do j = 1, 3
          if( kk( ikz, iky, ikx, j ) .gt. ladcap( 2, j ) ) &
              kk( ikz, iky, ikx, j ) = kk( ikz, iky, ikx, j ) - kmesh( j )
        enddo
!        write(74,*) xk( ikz, iky, ikx, : )
      enddo
    enddo
  enddo
  !

end subroutine AI_kmapr

!  open( unit=99, file='ladcap', form='formatted', status='unknown' )
!  rewind 99
!  read ( 99, * ) kl( : ), kh( : )
!  close ( unit=99 )
!  open( unit=99, file='kpts.dat', form='formatted', status='unknown' )
!  rewind 99
!  nkv( 1 ) = nkx; nkv( 2 ) = nky; nkv( 3 ) = nkz
!  c( : ) = 2.0d0 * pi / dble( nkv( : ) )
!  ik = 0
!  do ikx = 1, nkx
!     do iky = 1, nky
!        do ikz = 1, nkz
!           ikv( 1 ) = ikx - 1; ikv( 2 ) = iky - 1; ikv( 3 ) = ikz - 1
!           ik = ik + 1
!           xk( ik, : ) = c( : ) * dble( ikv( : ) )
!           kk( ik, : ) = ikv( : )
!           read ( 99, * ) ick, xck( : )
!           if ( ick .ne. ik ) stop 'bad ick'
           ! mapping k points onto (approximately) zero-centered grid
!           do j = 1, 3
!              if ( abs( xck( j ) - xk( ik, j ) ) .gt. 0.0001d0 ) stop 'bad kpt'
!              if ( ikv( j ) .gt. ladcap( 2, j ) ) kk( ik, j ) = ikv( j ) - nkv( j )
!           end do
!        end do
!     end do
!  end do
!  close( unit=99 )
!  !
!  return
!end subroutine AI_kmapr
