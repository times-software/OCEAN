! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine orthred( nr, irc, ntot, nnew, dl, r, pheps, pheae, pspr, aepr, prec )
  implicit none
  !
  integer :: nr, irc, ntot, nnew
  real( kind = kind( 1.0d0 ) ) :: dl, r( irc ), pheps( nr, ntot ), pheae( nr, ntot )
  real( kind = kind( 1.0d0 ) ) :: pspr( irc, ntot ), aepr( irc, ntot )
  real( kind = kind( 1.0d0 ) ) :: prec
  !
  integer :: i, j, ierr
  real( kind = kind( 1.0d0 ) ) :: tmp, tmp1, tmp2, err1, err2, err3
  real( kind = kind( 1.0d0 ) ) :: w( ntot ), fv1( ntot ), fv2( ntot ), fm1( 2 * ntot )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: ar, ai, zr, zi
  !
  allocate( ar( ntot, ntot ), ai( ntot, ntot ), zr( ntot, ntot ), zi( ntot, ntot ) )
  !
  ai = 0
  err1 = 0.0d0; err2 = 0.0d0
  do i = 1, ntot
     do j = 1, ntot
        call radint( irc, r, dl, pheps( :, i ), pheps( :, j ), tmp )
        if ( i .eq. j ) then 
           err1 = max( err1, abs( tmp - 1.0d0 ) )
        else
           err2 = max( err2, abs( tmp ) )
        end if
        ar( i, j ) = -tmp
     end do
  end do
  write ( 6, '(1a10,2(1x,1e15.8))' ) 'norm err =', err1, err2
  call elsch( ntot, ntot, ar, ai, w, 1, zr, zi, fv1, fv2, fm1, ierr )
  !
  nnew = 0
  do i = 1, ntot
     if ( abs( w( i ) ) .gt. prec ) nnew = nnew + 1
  end do
  write ( 6, '(6(1x,1e15.8))' ) w( 1 : nnew )
  write ( 6, '(2x,1a9,1f10.4)' ) 'runsu = ', sum( w( : ) )
  !
  err1 = 0.0d0; err2 = 0.d0; err3 = 0.0d0
  do i = 1, ntot
     err1 = max( err1, abs( sum( zr( :, i ) ** 2 ) - 1.0d0 ) )
     err2 = max( err2, sum( zi( :, i ) ** 2 ) )
     do j = 1, ntot
        if ( j .ne. i ) err3 = max( err3, dot_product( zr( :, i ), zr( :, j ) ) )
     end do
  end do
  write ( 6, '(1a7,5x,3(1x,1e15.8))' ) 'errs = ', err1, err2, err3
  !
  do i = 1, ntot
     zr( :, i ) = zr( :, i ) / sqrt( abs( w( i ) ) )
  end do
  pspr( :, : ) = 0.d0; aepr( :, : ) = 0.0d0
  do i = 1, nnew
     do j = 1, ntot
        pspr( 1 : irc, i ) = pspr( 1 : irc, i ) + zr( j, i ) * pheps( 1 : irc, j )
        aepr( 1 : irc, i ) = aepr( 1 : irc, i ) + zr( j, i ) * pheae( 1 : irc, j )
     end do
  end do
  err1 = 0.0d0
  err2 = 0.0d0
  do i = 1, nnew
     do j = 1, nnew
        call radint( irc, r, dl, pspr( :, i ), pspr( :, j ), tmp1 )
        call radint( irc, r, dl, aepr( :, i ), aepr( :, j ), tmp2 )
        tmp1 = abs( tmp1 ); tmp2 = abs( tmp2 ); 
        if ( j .eq. i ) then
           tmp1 = abs( tmp1 - 1.0d0 ); tmp2 = abs( tmp2 - 1.0d0 ); 
        end if
        err1 = max( err1, tmp1 ); err2 = max( err2, tmp2 )
     end do
  end do
  write ( 6, '(1a20,2(1x,1e15.8))' ) 'new overlap error = ', err1, err2
  !
  return
end subroutine orthred
