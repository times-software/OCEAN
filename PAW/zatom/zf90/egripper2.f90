! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine egripper2( nr, nelmax, r, dr, r2, dl, nl, pheps, pheae, njrc, vi, zorig, iuflag, orb, &
     l, ntest, irc, frac, rel, eint, emin, emax, north, prec, prec2 )
  implicit none
  !
  integer :: nr, nelmax, njrc( 4 )
  integer :: l, iuflag, ntest, irc, north
  !
  real( kind = kind( 1.0d0 ) ) :: dl, zorig, emin, emax, prec2
  real( kind = kind( 1.0d0 ) ) :: rel
  real( kind = kind( 1.0d0 ) ) :: r( nr ),dr( nr ),r2( nr )
  real( kind = kind( 1.0d0 ) ) :: pheps( nr, nelmax )
  real( kind = kind( 1.0d0 ) ) :: pheae( nr, nelmax )
  real( kind = kind( 1.0d0 ) ) :: orb( nr, nelmax )
  integer :: nl( nelmax )
  real( kind = kind( 1.0d0 ) ) :: vi( nr, 7 )
  !
  real( kind = kind( 1.0d0 ) ) :: frac( 0 : ntest )
  real( kind = kind( 1.0d0 ) ) :: eint( 0 : ntest )
  !
  integer :: i, j, k, nco, jj, ierr
  real( kind = kind( 1.0d0 ) ) :: f, runf
  real( kind = kind( 1.0d0 ) ), allocatable :: enew( : ), dif( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: map( : , : )
  real( kind = kind( 1.0d0 ) ), external :: efrac
  real( kind = kind( 1.0d0 ) ) :: prec, runsu, su, nrm
  real( kind = kind( 1.0d0 ) ) :: fac1, fac2
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: ff( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: orthps( : , : ), workps( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: orthae( : , : ), workae( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ar(:,:), ai(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: zr(:,:), zi(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: w(:),fv1(:),fv2(:),fm1(:)
  !
  integer :: nnew
  real( kind = kind( 1.0d0 ) ) :: sm1, sm2
  real( kind = kind( 1.0d0 ) ), allocatable :: newps(:,:),newae(:,:)
  !
  open(unit=99,file='corcon', form='formatted',status='unknown')
  rewind 99
  read ( 99, * ) nco
  close( unit=99 )
  !
  allocate( enew( 0 : ntest ), dif( 0 : ntest ) )
  allocate( map( 0 : ntest, 0 : ntest ) )
  do i = 0, ntest
     eint( i ) = emin + ( emax - emin ) * dble( i ) / dble( ntest )
  end do
  write ( 77, * ) 'calling grips up front'
  call grips( ntest + 1, l, nr, nelmax, r, dr, r2, dl, nl, pheps, njrc, vi, zorig, eint, orb, iuflag, &
       irc, map, dif, 0.d0, .true. )
  call optdif( ntest, eint, 'dif1', l, dif )
  call optmap( ntest, 'map1', l, map )
  !
  su = 0
  do i = 1, ntest - 1, 2
     dif( i ) = abs( dif( i ) )
     su = su + sqrt( dif( i ) )
  end do
  !
  frac( 0 ) =0.d0
  frac( ntest ) = 1.d0
  runf = 0.d0
  do i = 1, ntest - 3, 2
     runf = runf + sqrt( dif( i ) )
     frac( i + 1 ) = runf / su
  end do
  do i = 0, ntest
  end do
  !
  do i = 0, ntest
     f = dble( i ) / dble( ntest )
     enew( i ) = efrac( ntest, f, frac, eint )
  end do
  write ( 77, * ) 'calling grips below 1'
  call grips( ntest + 1, l, nr, nelmax, r, dr, r2, dl, nl, pheps, njrc, vi, zorig, enew, orb, iuflag, &
       irc, map, dif, 0.d0, .true. )
  call optdif( ntest, enew, 'dif2', l, dif )
  call optmap( ntest, 'map2', l, map )
  write ( 77, * ) 'calling grips below 2'
  call grips( ntest + 1, l, nr, nelmax, r, dr, r2, dl, nl, pheae, njrc, vi, zorig, enew, orb, iuflag, &
       irc, map, dif, rel, .false. )
  write ( 77, * ) 'back from grips'
  !
  allocate( orthps( irc, ntest + 2 ), orthae( irc, ntest + 2 ) )
  allocate( workps( irc ), workae( irc ) )
  allocate( ff( irc ) )
  north = 0
  do i = 1, ntest + 2
     do j = 1, irc
        workps( j ) = pheps( j, i )
        workae( j ) = pheae( j, i + nco )
     end do
     nrm = 1.d0
     if ( north .gt. 0 ) then
        do k = 1, north
           do j = 1, irc
              ff( j ) = workps( j ) * orthps( j, k ) / r( j ) ** 2
           end do
           call bintegrate( irc, r, dl, ff, su, irc )
           do j = 1, irc
              workps( j ) = workps( j ) - su * orthps( j, k )
           end do
           nrm = nrm - su ** 2
           do j = 1, irc
              ff( j ) = workae( j ) * orthae( j, k ) / r( j ) ** 2
           end do
           call bintegrate( irc, r, dl, ff, su, irc )
           do j = 1, irc
              workae( j ) = workae( j ) - su * orthae( j, k )
           end do
        end do
     end if
     if ( nrm .gt. prec ) then
        north = north + 1
        do j = 1, irc
           ff( j ) = ( workps( j ) / r( j ) ) ** 2
        end do
        call bintegrate( irc, r, dl, ff, su, irc )
        su = 1.d0 / sqrt( su )
        do j = 1, irc
           orthps( j, north ) = workps( j ) * su
        end do
        do j = 1, irc
           ff( j ) = ( workae( j ) / r( j ) ) ** 2
        end do
        call bintegrate( irc, r, dl, ff, su, irc )
        su = 1.d0 / sqrt( su )
        do j = 1, irc
           orthae( j, north ) = workae( j ) * su
        end do
     else
     end if
  end do
  !
  runsu = 0.d0
  do i = 1, ntest + 2
     do k = 1, north
        do j = 1, irc
           ff( j ) = orthps( j, k ) * pheps( j, i ) / r( j ) ** 2
        end do
        call bintegrate( irc, r, dl, ff, su, irc )
        runsu = runsu + su ** 2
     end do
  end do
  write ( 6, '(2x,1a9,1f10.4)' ) 'runsu = ', runsu
  !
  allocate( ar( north, north ) )
  allocate( ai( north, north ) )
  allocate( zr( north, north ) )
  allocate( zi( north, north ) )
  allocate( fv1( north ), fv2( north ), fm1( 2 * north ) )
  allocate( w( north ) )
  do i = 1, north
     do j = 1, north
        runsu = 0.d0
        do k = 1, ntest + 2
           do jj = 1, irc
              ff( jj ) = orthps( jj, i ) * pheps( jj, k ) / r( jj ) ** 2
           end do
           call bintegrate( irc, r, dl, ff, fac1, irc )
           do jj = 1, irc
              ff( jj ) = orthps( jj, j ) * pheps( jj, k ) / r( jj ) ** 2
           end do
           call bintegrate( irc, r, dl, ff, fac2, irc )
           runsu = runsu + fac1 * fac2
        end do
        ar( i, j ) = - runsu
     end do
  end do
  ai = 0
  call elsch( north, north, ar, ai, w, 1, zr, zi, fv1, fv2, fm1, ierr )
  write ( 6, '(2x,1f20.10)' ) w
  do i = 1, north
     sm1 = 0.d0
     sm2 = 0.d0
     do j = 1, north
        sm1 = sm1 + zr( j, i ) ** 2
        sm2 = sm2 + zi( j, i ) ** 2
     end do
     if ( sm1 .gt. sm2 ) then
        su = sqrt( ( sm1 + sm2 ) / sm1 )
        do j = 1, north
           zr( j, i ) = zr( j, i ) * su
        end do
     else
        su = sqrt( ( sm1 + sm2 ) / sm2 )
        do j = 1, north
           zr( j, i ) = zi( j, i ) * su
        end do
     end if
  end do
  !
  allocate( newps( irc, north ), newae( irc, north ) )
  nnew = 0
  do i = 1, north
     if ( abs( w( i ) ) .gt. prec2 ) then
        nnew = nnew + 1
        do j = 1, irc
           newps( j, nnew ) = 0.d0
           newae( j, nnew ) = 0.d0
           do k = 1, north
              newps( j, nnew ) = newps( j, nnew ) + zr( k, i ) * orthps( j, k )
              newae( j, nnew ) = newae( j, nnew ) + zr( k, i ) * orthae( j, k )
           end do
        end do
     end if
  end do
  !
  deallocate( ar, ai, zr, zi, fv1, fv2, fm1, w )
  !
  do k = 1, nnew  
     pheps( :, k ) = 0
     pheae( :, k + nco ) = 0
     do j = 1, irc
        pheps( j, k ) = newps( j, k )
        pheae( j, k + nco ) = newae( j, k )
     end do
  end do
  !
  do i = 1, nnew
     do j = 1, nnew
        do jj = 1, irc
           ff( jj ) = pheps( jj, i ) * pheps( jj, j ) / r( jj ) ** 2
        end do
        call bintegrate( irc, r, dl, ff, su, irc )
        write ( 6, '(2x,1a6,2i5,1x,1f10.6)' ) 'ovlp: ', i, j, su
     end do
  end do
  !
  north = nnew
  !
  deallocate( dif, map, enew, orthps, orthae, workps, workae, ff )
  deallocate( newps, newae )
  !
  return
end subroutine egripper2
!
subroutine optmap( ntest, name, l, map )
  implicit none
  !
  integer :: ntest, l
  character * 4 :: name
  character * 5 :: name2
  real( kind = kind( 1.0d0 ) ) :: map( 0 : ntest, 0 : ntest ) 
  !
  integer :: i, j
  !
  write ( unit=name2, fmt='(1a4,1i1)' ) name, l ! mknam
  open( unit=99, file=name2, form='formatted', status='unknown' )
  rewind 99
  do i = 0, ntest
     do j = 0, ntest
        write ( 99, '(2x,2i5,1(2x,1e15.8))' ) i, j, map( i, j )
     end do
     write ( 99, * )
  end do
  close( unit=99 )
  !
  return
end subroutine optmap
!
subroutine optdif( ntest, e, name, l, dif )
  implicit none
  !
  integer :: ntest, l
  character * 4 :: name
  character * 5 :: name2
  real( kind = kind( 1.0d0 ) ) :: e( 0 : ntest ), dif( 0 : ntest )
  !    
  integer :: i
  !
  write ( unit=name2, fmt='(1a4,1i1)' ) name, l
  open( unit=99, file=name2, form='formatted', status='unknown' )
  rewind 99
  do i = 1, ntest - 1, 2
     write ( 99, '(2x,1i5,2(2x,1e15.8))' ) i, e( i ), dif( i )
  end do
  close( unit=99 )
  !
  return
end subroutine optdif
!
function efrac( ntest, f, frac, eint )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: efrac
  !
  integer :: ntest
  real( kind = kind( 1.0d0 ) ) :: f
  real( kind = kind( 1.0d0 ) ) :: frac( 0 : ntest ), eint( 0 : ntest )
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: el, eh, fl, fh
  !
  if ( ( f .ge. frac( 2 ) ) .and. ( f .le. frac( ntest - 2 ) ) ) then
     i = 0
     do while ( f .ge. frac( i + 2 ) )
        i = i + 2
     end do
  else
     if ( f .le. frac( 2 ) ) then
        i = 0
     else
        i = ntest - 2
     end if
  end if
  el = eint( i )
  eh = eint( i + 2 )
  fl = frac( i )
  fh = frac( i + 2 )
  efrac = el + ( eh - el ) * ( f - fl ) / ( fh - fl )
  !
  return
end function efrac
!
subroutine grips( ne, l, nr, nelmax, r, dr, r2, dl, nl, phe, njrc, vi, zorig, eint, orb, iuflag, &
     irc, smtab, diff, rel, mode )
  implicit none
  !
  integer :: ne, l, nr, nelmax, njrc( 4 ), iuflag, irc
  !
  logical :: mode
  !
  real( kind = kind( 1.0d0 ) ) :: dl, zorig
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  real( kind = kind( 1.0d0 ) ) :: phe( nr, nelmax ), vi( nr, 7 )
  real( kind = kind( 1.0d0 ) ) :: orb( nr, nelmax )
  !
  real( kind = kind( 1.0d0 ) ) :: smtab( ne, ne ), diff( ne )
  real( kind = kind( 1.0d0 ) ) :: eint( ne )
  integer :: nl( nelmax )
  !
  integer :: nel
  real( kind = kind( 1.0d0 ) ) :: etot, rel, alfa, xntot
  !
  integer, allocatable :: no( : ), nm( : ), is( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xnj( : ), occ( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ev( : ), ek( : ), cq( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: f( : )
  !
  integer :: nval, ncor, ntot
  real( kind = kind( 1.0d0 ) ) :: xj
  integer, allocatable :: ntmp( : ), ltmp( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: occtmp( : )
  !
  integer :: i, j, i1, i2
  real( kind = kind( 1.0d0 ) ) :: su
  !
  allocate( cq( nr ) )
  allocate( f( irc ) )
  !
  alfa = 0.d0
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ncor
  close( unit=99 )
  !
  open( unit=99, file='valcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nval
  allocate ( ntmp( nval ), ltmp( nval ), occtmp( nval ) )
  do i = 1, nval
     read ( 99, * ) ntmp( i ), ltmp( i ), occtmp( i )
  end do
  close( unit=99 )
  !
  ntot = ne + ncor + nval
  allocate( no( ntot ), nm( ntot ), is( ntot ), xnj( ntot ) )
  allocate( occ( ntot ), ev( ntot ), ek( ntot ) )
  !
  open( unit=99, file='config', form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) ne + nval, 0.4d0, 0.0001d0, 100.d0, 0, 2
  xj = - l
! if ( l .eq. 0 ) xj = -0.5d0
  do i = 1, nval
     if ( ltmp( i ) .eq. l ) then
        write ( 99, * ) ntmp( i ), l, 0, xj, 1, occtmp( i )
     end if
  end do
  do i = 1, ne
     write ( 99, * ) -1, l, 0, xj, 1, 0.d0
     write ( 99, * ) eint( i )
  end do
  do i = 1, nval
     if ( ltmp( i ) .ne. l ) then
        xj = - ltmp( i )
!       if ( ltmp( i ) .eq. 0 ) xj = -0.5d0
        write ( 99, * ) ntmp( i ), ltmp( i ), 0, xj, 1, occtmp( i )
     end if
  end do
  close( unit=99)
  deallocate( ntmp, ltmp, occtmp )
  !
  call abinitio(etot,rel,alfa,nr,r,dr,r2,dl, phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj, &
       ev,occ,is,ek,orb,iuflag,cq,mode,nelmax)
  !
  do i = 1, ne
     do j = 1, irc
        f( j ) = phe( j, i + 1 ) ** 2 / r( j ) ** 2
     end do
     call bintegrate( irc, r, dl, f, su, irc )
     su = 1.d0 / sqrt( su )
     do j = 1, nr
        phe( j, i + 1 ) = phe( j , i + 1 ) * su
!       write ( 123, '(2(1x,1e15.8),1x,1a20)' ) r( j ), phe( j, i + 1 ), mode
     end do
  end do
  do i1 = 1, ne
     do i2 = 1, ne
        smtab( i1, i2 ) = 0.d0
        if ( iabs( i1 - i2 ) .lt. 2 ) then
           do j = 1, irc
              f( j ) = phe( j, i1 + 1 ) * phe( j, i2 + 1 ) / r( j ) ** 2
           end do
           call bintegrate( irc, r, dl, f, smtab( i1, i2 ), irc )
        end if
     end do
  end do
  diff( 1 ) = 0.d0
  diff( ne ) = 0.d0
  do i = 2, ne - 1
     diff( i ) = 2.d0 * smtab( i, i ) - smtab( i, i + 1 ) - smtab( i, i - 1 )
  end do
  !
  deallocate( no, nm, is, xnj, occ, ev, ek, cq, f )
  !
  return
end subroutine grips
