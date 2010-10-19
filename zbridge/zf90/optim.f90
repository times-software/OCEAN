! the function optim returns the length [optim] for a fft grid, 
! based on a mininum length [nmin].
! the length can have [nfac] factors stored in aray [fac], 
! such that each factor is within the [locap] and [hicap] bounds.
! having [locap] can help ensure divisibility by desired factors.
!
! optim is the smallest value possible.  
!
!---------------------------------------------------------
!
function optim( nmin, nfac, fac, locap, hicap )
  implicit none
  !
  integer :: optim, nmin, nfac
  integer :: fac( nfac ), locap( nfac ), hicap( nfac )
  !
  integer :: i, j, n, nbase, ntest
  logical :: valid
  integer, allocatable :: pow( : ), powmin( : ), powmax( : )
  !
  allocate( pow( nfac ), powmin( nfac ), powmax( nfac ) )
  powmin( : ) = locap( : )
  nbase = 1
  do j = 1, nfac
     do i = 1, powmin( j )
        nbase = nbase * fac( j )
     end do
  end do
  powmax( : ) = powmin( : ) 
  do j = 1, nfac
     i = nbase
     do while ( i .lt. nmin )
        powmax( j ) = powmax( j ) + 1
        i = i * fac( j )
     end do
     powmax( j ) = min( powmax( j ), hicap( j ) )
     if ( powmax( j ) .lt. powmin( j ) ) stop 'impossibility in optim!'
  end do
  n = 1
  do j = 1, nfac
     n = n * fac( j ) ** powmax( j )
  end do
  call reginit( nfac, pow, powmin )
  valid = .true.
  do while( valid )
     ntest = 1
     do j = 1, nfac
        ntest = ntest * fac( j ) ** pow( j )
     end do
     if ( ( ntest .gt. nmin ) .and. ( ntest .le. n ) ) n = ntest
     call reginc( nfac, pow, powmin, powmax, valid )
  end do
  deallocate( pow, powmin, powmax )
  optim = n
  !
  return
end function optim
!
!-----------------------------------------------------------
!
subroutine reginit( nreg, reg, regmin )
  implicit none
  !
  integer :: nreg
  integer :: reg( nreg ), regmin( nreg )
  !
  integer :: i
  !
  do i = 1, nreg
     reg( i ) = regmin( i )
  end do
  !
  return
end subroutine reginit
!
!-----------------------------------------------------------
!
subroutine reginc( nreg, reg, regmin, regmax, valid )
  implicit none
  !
  integer :: nreg
  integer :: reg( nreg ), regmin( nreg ), regmax( nreg )
  logical :: valid
  !
  integer :: foc
  logical :: done
  !
  done = .false.
  foc = 1
  do while ( ( .not. done ) .and. ( foc .le. nreg ) )
     reg( foc ) = reg( foc ) + 1
     if ( reg( foc ) .le. regmax( foc ) ) then
        done = .true.
     else
        reg( foc ) = regmin( foc )
        foc = foc + 1
        done = .false.
     end if
  end do
  valid = ( foc .le. nreg )
  !
  return
end subroutine reginc
