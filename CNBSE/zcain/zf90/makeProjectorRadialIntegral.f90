subroutine makeProjectorRadialIntegral( npmax, lmin, lmax, nproj, atno, powmax, radialPrj, ierr )
  implicit none
  !
  integer,  parameter :: DP = kind(1.0d0)
  integer, intent( in ) :: npmax, lmin, lmax, nproj( lmin : lmax ), atno, powmax
  integer, intent( inout ) :: ierr
  real(DP), intent( out ) :: radialPrj( npmax, lmin:lmax, npmax, lmin:lmax, 0 : powmax )
  !
  integer :: dumi, nrad, ll, ll1, ll2, nu1, nu2, ip
  real(DP) :: dumf, area, dl
  character( len=11 ) :: filnam

  real(DP), allocatable :: projectors( :, :, : ), radialGrid(:), projLine(:), f(:)
  !
  write( filnam, '(A,I3.3)' ) 'radfilez', atno
  open( unit=99, file=filnam, form='formatted', status='old' )
  read( 99, * ) dumf, dumi, nrad
  close( 99 )

  allocate( projectors( nrad, npmax, lmin:lmax ), radialGrid( nrad ) )
  
  do ll = lmin, lmax
    write( filnam, '(A,I1.1,A,I3.3)' ) 'ae', ll, 'z', atno
  
    allocate( projLine( nproj( ll )
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
  dl = 
  do ip = 0, powmax

    do ll1 = lmin, lmax
      do nu1 = 1, proj( ll1 )

        do ll2 = lmin, lmax
          do nu2 = 1, proj( ll2 )
            f(:) = projectors( :, nu2, ll2 ) * projectors( :, nu1, ll1 )
            call bintegrate( nrad, radialGrid, dl, f, area )
            if( ll2 .eq. ll1 ) then
              write(6,'(4I8,X,20.12E)') ip, ll1, nu1, nu2, area
            endif
            radialPrj( nu2, ll2, nu1, ll1, ip ) = area
          enddo
        enddo
      
      enddo
    enddo
  enddo
    

end subroutine

subroutine bintegrate( nr, r, dl, f, area, irc )
  implicit none
  !
  integer :: nr, irc, j, k, m
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
  k = irc
  do while ( k .gt. 0 )
     if ( ( k .eq. irc ) .or. ( k .lt. 5 ) ) then
        su = su + 14.0d0 * f( k ) * r( k ) !r( k ) ** 3
     else
        su = su + 28.0d0 * f( k ) * r( k ) !r( k ) ** 3
     end if
     k = k - 4
  end do
  k = k + 4
  j = irc - 1
  do while ( j .gt. k )
     su = su + 64.0d0 * f( j ) * r( j ) !r( j ) ** 3
     j = j - 2
  end do
  m = irc - 2
  do while ( m .gt. k )
     su = su + 24.0d0 * f( m ) * r( m ) !r( m ) ** 3
     m = m - 4
  end do
  area = su * dl / 45.0d0
  !
  return
end subroutine
