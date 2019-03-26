subroutine makeProjectorRadialIntegral( npmax, lmin, lmax, nproj, atno, powmax, radialPrj, ierr )
  implicit none
  !
  integer,  parameter :: DP = kind(1.0d0)
  integer, intent( in ) :: npmax, lmin, lmax, nproj( lmin : lmax ), atno, powmax
  integer, intent( inout ) :: ierr
  real(DP), intent( out ) :: radialPrj( npmax, lmin:lmax, npmax, lmin:lmax, 0 : powmax )
  !
  integer :: dumi, nrad, ll, ll1, ll2, nu1, nu2, ip, ir
  real(DP) :: dumf, area, dl, mel
  character( len=11 ) :: filnam
  character( len=17 ) :: cornam

  real(DP), allocatable :: projectors( :, :, : ), radialGrid(:), projLine(:), f(:), core(:)
  !
  write( filnam, '(A,I3.3)' ) 'radfilez', atno
  open( unit=99, file=filnam, form='formatted', status='old' )
  read( 99, * ) dumf, dumi, nrad
  close( 99 )

  allocate( projectors( nrad, npmax, lmin:lmax ), radialGrid( nrad ) )
  
  do ll = lmin, lmax
    write( filnam, '(A,I1.1,A,I3.3)' ) 'ae', ll, 'z', atno
  
    allocate( projLine( nproj( ll ) ) )
    open(unit=99,file=filnam, form='formatted', status='old' )
    do ir = 1, nrad
      read(99,*) radialGrid(ir), projLine(:)
      projectors( ir, 1:nproj(ll), ll ) = projLine( : ) 
    enddo
    deallocate( projLine )

    close( 99 )
  enddo

  write(6,*) 'Projectors read in'


  radialPrj( :, :, :, :, : ) = 0.0_dp
  allocate( f( nrad ) )
  dl = 0.01d0 * dlog(radialGrid(101)/radialGrid(1))
  do ip = 0, powmax

    do ll1 = lmin, lmax
      do nu1 = 1, nproj( ll1 )

        do ll2 = lmin, lmax
          do nu2 = 1, nproj( ll2 )
            f(:) = projectors( :, nu2, ll2 ) * projectors( :, nu1, ll1 ) * radialGrid(:)**ip
            call bintegrate( nrad, radialGrid, dl, f, area )
            call rpower( nrad, radialGrid, dl, projectors( :, nu2, ll2 ), projectors( :, nu1, ll1 ), ip, mel )
!            if( ll2 .eq. ll1 ) then
              write(6,'(5I8,2(X,E20.12))') ip, ll1, nu1, ll2, nu2, area, mel
!            endif
            radialPrj( nu2, ll2, nu1, ll1, ip ) = area
          enddo
        enddo
      
      enddo
    enddo
  enddo
    
  allocate( core( nrad ) )
  write(cornam, '(A,I3.3,A)' ) 'coreorbz', atno, 'n01l00'
  open( unit=99, file=cornam, form='formatted', status='old' )
    read( 99, * ) 
  do ir = 1, nrad
    read( 99, * ) dumf, core( ir )
  enddo
  close(99)

  do ll1 = lmin, lmax
    do ip = 0, powmax
      do nu1 = 1, nproj( ll1 )
        f(:) = core(:) * projectors( :, nu1, ll1 ) * radialGrid(:)**(ip-1)
        call bintegrate( nrad, radialGrid, dl, f, area )
        call rpower( nrad, radialGrid, dl, core, projectors( :, nu1, ll1 ), ip-1, mel )
        write(6,'(3I8,2(X,E20.12))') ll1, ip, nu1, area, mel
      enddo
    enddo
  enddo

end subroutine

subroutine bintegrate( nr, r, dl, f, area )
  implicit none
  !
  integer :: nr,  j, k, m
  real( kind = kind( 1.0d0 ) ) :: dl, area, su, f( nr ), r( nr )
  !  
  real( kind = kind( 1.0d0 ) ), allocatable :: rrr( : )
  !
  ! f is multiplied by r^3,
  ! where one r comes from each of the radial wf's,
  ! and one r comes from using a log mesh.
  ! 
  !
  su = 0.0d0
  k = nr
  do while ( k .gt. 0 )
     if ( ( k .eq. nr ) .or. ( k .lt. 5 ) ) then
        su = su + 14.0d0 * f( k ) * r( k )**3 !r( k ) ** 3
     else
        su = su + 28.0d0 * f( k ) * r( k )**3 !r( k ) ** 3
     end if
     k = k - 4
  end do
  k = k + 4
  j = nr - 1
  do while ( j .gt. k )
     su = su + 64.0d0 * f( j ) * r( j )**3 !r( j ) ** 3
     j = j - 2
  end do
  m = nr - 2
  do while ( m .gt. k )
     su = su + 24.0d0 * f( m ) * r( m )**3 !r( m ) ** 3
     m = m - 4
  end do
  area = su * dl / 45.0d0
  !
  return
end subroutine

subroutine rpower( nr, r, dl, ph1, ph2, npowr, mel )
  implicit none
  !
  integer :: nr, npowr
  real( kind = kind( 1.0d0 ) ) :: dl, mel, r( nr ), ph1( nr ), ph2( nr )
  !
  integer :: i, j, ii
  real( kind = kind( 1.0d0 ) ) :: f( 0 : 4 )
  !
  mel = 0.d0

  do i = 1, nr - 4, 4
    do j = 0, 4
      ii = i + j
      f( j ) = ph1( ii ) * ph2( ii ) * r( ii ) ** ( npowr + 2 )
    end do
    mel = mel + 14.0d0 / 45.0d0 * dl * r( i ) * ( f( 0 ) + f( 4 ) )
    mel = mel + 64.0d0 / 45.0d0 * dl * r( i ) * ( f( 1 ) + f( 3 ) )
    mel = mel + 24.0d0 / 45.0d0 * dl * r( i ) * f( 2 )
  end do
  !
  return
end subroutine rpower
