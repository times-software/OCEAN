! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine mwreduce( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rcocc, lmin, lmax ) 
  implicit none
  !
  double precision :: etot, rel, alfa
  double precision :: dl, zorig, xntot, rcocc( 0 : 3 )
  integer :: nr, njrc( 4 ), nel, isuse, lmin, lmax, skips( 0 : 3 )
  integer :: iuflag
  double precision :: r( nr ),dr( nr ),r2( nr )
  double precision :: vi( nr, 7 ), cq( nr ), etmin( lmin : lmax )
  !
  integer, allocatable :: no( : ), nl( : ), nm( : )
  integer, allocatable :: is( : )
  double precision, allocatable :: ev( : ), occ( : ), ek( : )
  double precision, allocatable :: xnj( : )
  !
  integer :: ncore, lcore, icore, nq, nco, n2, is2, nuse, luse, iuse
  integer :: l, ne, lcmin, lcmax, k
  integer :: i, j, lev, ilev, jlev, irc
  integer :: tfac, nt, i1, i2, j1, j2, lc
  double precision :: emin, emax, emid, dq, qmax, prec, prec2, rc, ctol, condnum, su, dummy
  double precision, external :: efrac
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: c( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: meltmp( : )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: orb, proj, pheae, pheps, aesav
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: meltrue, melfake, melval
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: etrue, efake, exch
  !
  integer, parameter :: prjfile = 93, melfile = 94
  !
  integer  :: ntest, nelmax
  double precision, allocatable :: eint( : ), frac( : )
  double precision, dimension( nr ) :: ff, copy, reconned, phespc
  !
  integer :: nval, npower
  double precision :: xj
  integer, allocatable, dimension( : ) :: ntmp, ltmp
  double precision, allocatable :: occtmp( : )
  !
  real( kind = kind( 1.0d0 ) ) :: ru, au, nulocc( 0 : 3 ), hptab( nr, 2 )
  !
  character * 1 :: ltr
  character * 3 :: nam3
  character * 4 :: nam4
  character * 6 :: s6
  character * 10 :: addend, adduse
  character * 11 :: filnam11
  character * 16 :: s16
  !
  read ( 5, * ) ncore, lcore
  write ( addend, '(1a1,1i3.3,2(1a1,1i2.2))' ) 'z', nint( zorig ), 'n', ncore, 'l', lcore
  call geticore( nco, icore, ncore, lcore )
  write ( 6, * ) ' all stuffs = ', icore, ncore, lcore
  !
  open( unit=99, file='skip', form='formatted', status='unknown' )
  read ( 99, * )
  do l = 0, 3
     read ( 99, * ) i, skips( l )
  end do
  !
  ntest = 80
  allocate( eint( 0 : ntest ), frac( 0 : ntest ) )
  !
  n2 = 1 + lmax - lmin
  nelmax = 1 + ntest + nco + ( 1 + lmax - lmin ) + 10
  allocate( no( nelmax ), nl( nelmax ) )
  allocate( nm( nelmax ), is( nelmax ) )
  allocate( ev( nelmax ), ek( nelmax ) )
  allocate( occ( nelmax ), xnj( nelmax ) )
  allocate( pheps( nr, nelmax ), pheae( nr, nelmax ) )
  allocate( orb( nr, nelmax ) )
  !
  call freshen( lmin, lmax, rcocc, skips )
  call abinitio(etot,rel,alfa,nr,r,dr,r2,dl,pheae,njrc,vi,zorig,xntot,nel,no, &
                nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.true.,nelmax)
  write ( 6, * ) nl( 1 : 3 ), ev( 1 : 3 )
  do i = 1, n2
     etmin( nl( i ) ) = ev( i ) - 0.15d0
  end do
  do l = lmin, lmax
     write ( 6, '(1i3,5x,1a8,1f10.5)' ) l, 'etmin = ', etmin( l )
  end do
  !
  read ( 5, * ) npower
  read ( 5, * ) emin, emax, prec, prec2, tfac
  read ( 5, * ) rc, ctol
  read ( 5, * ) dq, qmax
  nq = 1 + qmax / dq
  !
  ! reset cutoff radius to radius value on the radial grid
  irc = 0
  do i = 1, nr, 4
     if ( r( i ) .lt. rc ) irc = i
  end do
  write ( filnam11, '(1a8,1i3.3)' ) 'radfilez', nint( zorig )
  open( unit=99, file=filnam11, form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) r( irc ), nr, irc
  close( unit=99 )
  !
  write ( filnam11, '(1a8,1i3.3)' ) 'prjfilez', nint( zorig )
  open( unit=prjfile, file=filnam11, form='formatted', status='unknown' )
  write ( filnam11, '(1a8,1i3.3)' ) 'melfilez', nint( zorig )
  open( unit=melfile, file=filnam11, form='formatted', status='unknown' )
  rewind prjfile 
  rewind melfile
  write ( prjfile, '(3i5,2x,1e15.8)' ) lmin, lmax, nq, dq
  !
  ! loop over values of valence angular momentum ... 
  do l = lmin, lmax
     !
     ! determine suitable project functions
     write ( 6, * ) 'calling egripper'
     call egripper2( nr, nelmax, r, dr, r2, dl, nl, &
          &                pheps, pheae, njrc, vi, zorig, iuflag, orb, &
          &                l, ntest, irc, frac, rel, &
          &                eint, emin, emax, ne, prec, prec2 ) 
     write ( 6, * ) 'back from egripper'
     write ( 6, *)  'project functions found'
     write ( 6, * ) 'l, ne ... ', l, ne
     !
     ! output projector functions and their fourier transforms 
     nl( 1 : ne ) = l 
     call sdump( l, .true., ne, nl, nr, r, pheps, nelmax, zorig )
     nl( nco + 1 : nco + ne ) = l
     call sdump( l, .false., ne, nl, nr, r, pheae, nelmax, zorig )
     !
     allocate( aesav( irc, ne ) )
     aesav( 1 : irc, 1 : ne ) = pheae( 1 : irc, 1 + nco : ne + nco )
     call breader( l, irc, condnum, nr, dl, ne, zorig )
     !
     lcmin = iabs( lcore - l )
     lcmax = lcore + l 
     write ( 6, '(1a14,2i5)' ) 'lcmin, lcmax: ', lcmin, lcmax
     write ( prjfile, '(1i5)' ) ne
     allocate( exch( lcmin : lcmax, ne, ne ) )
     phespc = pheae( :, icore )
     write ( 6, * ) 'icore = ', icore, 'ncore = ', ncore, 'lcore = ', lcore
     nt = tfac * ne
     allocate( melval( 0 : npower, ne ), meltrue( 0 : npower, nt ), melfake( 0 : npower, nt ), meltmp( 0 : npower ) )
     write ( 6, * ) ' for melfile ...', ne
     do i = 1, ne
        lev = i + nco
        call rpower( nr, r, dl, phespc, pheae( 1, lev ), npower, meltmp )
        melval( :, i ) = meltmp( : )
     end do
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
     open( unit=99, file='config', form='formatted', status='unknown' )
     rewind 99
     write ( 99, * ) nt + nval, 0.4d0, 0.0001d0, 100.d0, 0, 0
     do i = 1, 2 * nt, 2
        write ( 99, * ) -1, l, l, -dble( l ), 1, 0.d0
        emid = etmin( l ) + ( emax - etmin( l ) ) * dble( i ) / dble( 2 * nt )
        write ( 99, * ) emid
     end do
     do i = 1, nval
        xj = -ltmp( i )
        if ( ltmp( i ) .eq. 0 ) xj = -0.5d0
        write ( 99, * ) ntmp( i ), ltmp( i ), 0, xj, 1, occtmp( i )
     end do
     close( unit=99)
     deallocate( ntmp, ltmp, occtmp )
     !
     call abinitio(etot,rel,alfa,nr,r,dr,r2,dl,pheae,njrc,vi,zorig,xntot,nel,no, &
                   nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     !
     allocate( etrue( lcmin : lcmax, nt, nt ) )
     do i = 1, nt
        ilev = i + ( nel - nval - nt )
        do j = 1, nt
           jlev = j + ( nel - nval - nt )
           call getexch( nr, r, dl, lcmin, lcmax, phespc, pheae( 1, ilev ), pheae( 1, jlev ), etrue( lcmin, i, j ) )
        end do
     end do
     do i = 1, nt
        lev = i + ( nel - nval - nt )
        call rpower( nr, r, dl, phespc, pheae( 1, lev ), npower, meltmp )
        meltrue( :, i ) = meltmp( : )
     end do
     !
     call amult( l, ne, nr, zorig )
     allocate( proj( nr, ne ), c( ne, nt ) )
     write ( unit=nam3, fmt='(1a2,1i1)' ) 'pr', l
     open( unit=99, file=nam3, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nr
        read ( 99, * ) dummy, ( proj( i, j ), j = 1, ne )
     end do
     close( unit=99 )
     !
     call transp( l, ne, nr, qmax, dq, nq, dl, irc, zorig )
     !
     do j = 0, npower 
        write ( melfile, '(8(2x,1e15.8))' ) melval( j, : )
     end do
     !
     call abinitio(etot,rel,alfa,nr,r,dr,r2,dl,pheps,njrc,vi,zorig,xntot,nel,no, &
                   nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.true.,nelmax)
     !
     ! reconstruct wave functions.  insodoing, output diagnostic file showing actual and reconstructed wave functions
     write ( unit=nam3, fmt='(1a2,i1)' ) 'di', l
     open( unit=99, file=nam3, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nt
        do k = 1, irc
           copy( k ) = pheps( k, i ) / r( k )
           reconned( k ) = 0.d0
        end do
        do j = 1, ne
           do k = 1, irc
              ff( k ) = proj( k, j ) * pheps( k, i ) / r( k )
           end do
           call bintegrate( nr, r, dl, ff, su, irc )
           do k = 1, irc
              copy( k ) = copy( k ) - su * proj( k, j )
              reconned( k ) = reconned( k ) + su * aesav( k, j ) / r( k )
           end do
        end do
        call recon( nr, r, dl, irc, pheps( 1, i ), ne, proj, c( 1, i ) )
        do j = 1, irc
           write ( 99, '(2x,2i5,5(2x,1e15.8))' ) is2, j, r( j ), copy( j ), pheps( j, i ),  reconned( j ), pheae( j, i + nco )
        end do
        write ( 99, * )
     end do
     !
     ! from reconstruction coefficients, calculate trumped-up transition matrix elements
     do i = 1, nt
        do j = 0, npower
           su = 0.d0
           do k = 1, ne
              su = su + c( k, i ) * melval( j, k )
           end do
           melfake( j, i ) = su
        end do
     end do
     close( unit=99 )
     !
     ! output diagnostic file to show quality of reconstructed transition matrix elements
     do j = 0, npower
        write ( unit=nam4, fmt='(1a2,2i1)' ) 'mt', j, l
        open( unit=99, file=nam4, form='formatted', status='unknown' )
        rewind 99
        do i = 1, nt
           write ( 99, '(2x,1i5,3(2x,1e15.8))' ) i, ev( i ), meltrue( j, i ), melfake( j, i )
        end do
        close( unit=99 )
     end do
     !
     ! from reconstruction coefficients, calculate trumped-up exchange integrals
     allocate( efake( lcmin : lcmax, nt, nt ) )
     do i1 = 1, nt
        do i2 = 1, nt
           do lc = lcmin, lcmax, 2
              su = 0.d0
              do j1 = 1, ne
                 do j2 = 1, ne
                    su = su + exch( lc, j1, j2 ) * c( j1, i1 ) * c( j2, i2 )
                 end do
              end do
              efake( lc, i1, i2 ) = su
           end do
        end do
     end do
     !
     ! output diagnostic file to show quality of reconstructed exchange integrals
     write ( unit=nam3, fmt='(1a2,1i1)' ) 'ex', l
     open( unit=99, file=nam3, form='formatted', status='unknown' )
     rewind 99
     do i1 = 1, nt
        do i2 = 1, nt
           do lc = lcmin, lcmax, 2
              write ( 99, '(2x,2i5,4(2x,1e15.8))' ) i1, i2, ev( i1 ), ev( i2 ), efake( lc, i1, i2 ), etrue( lc, i1, i2 )
           end do
        end do
        write ( 99, * )
     end do
     close( unit=99 )
     !
     deallocate( melval, meltrue, melfake, meltmp, proj, c, exch, etrue, efake, aesav )
     !
  end do
  close( unit=prjfile )
  close( unit=melfile )
  !
  ru = 1.0d0
  au = 0.0d0
  nulocc( : ) = 0.0d0
! call freshen( lmin, lmax, nulocc, skips )
  !
  ! loop over all core levels
! do istub = 1, nco
!    iuse = nco + 1 - istub
  do iuse = 1, nco
     !
     ! acquiesce N and L
     open( unit=99, file='corcon', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * )
     do j = 1, iuse
        read ( 99, * ) nuse, luse
     end do
     close( unit=99 )
     !
     ! fix file addend
     write ( adduse, '(1a1,1i3.3,2(1a1,1i2.2))' ) 'z', nint( zorig ), 'n', nuse, 'l', luse
     !
     ! make guidefile
     select case( luse )
     case( 0 )
        ltr = 's'
     case( 1 )
        ltr = 'p'
     case( 2 )
        ltr = 'd'
     case( 3 )
        ltr = 'f'
     end select
     write ( s6, '(1a4,1i1,1a1)' ) 'znl.', nuse, ltr
     open( unit=99, file=s6, form='formatted', status='unknown' )
     write ( 99, '(3i5)' ) nint( zorig ), nuse, luse
     close( unit=99 )
     !
     ! output deflin...
     write ( s16, '(1a6,1a10)' ) 'deflin', adduse
     open( unit=99, file=s16, form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(3i5)' ) luse, lmin, lmax
     close( unit=99 )
     !
     ! output real input files
     write ( 6, * ) 'valence null occ for core part begin'
     call freshen( lmin, lmax, nulocc, skips )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,pheae,njrc,vi,zorig,xntot,nel,no, &
                   nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     call hpload( hptab, nr, 1 )
     write ( 6, '(2i5,1x,1f10.5,1x,1a4)' ) iuse, nco, hptab( nr - 10, 1 ) * r( nr - 10 ), 'tail'
     call chgocc( iuse, -1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,pheae,njrc,vi,zorig,xntot,nel,no, &
                   nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     call hpload( hptab, nr, 2 )
     call chgocc( iuse, +1.0d0 )
     write ( 6, '(2i5,1x,1f10.5,1x,1a4)' ) iuse, nco, hptab( nr - 10, 2 ) * r( nr - 10 ), 'tail'
     call optript( nr, r, hptab, adduse, 'vcxxxxx' )
     write ( 6, * ) 'valence null occ for core part end'
     !
     ! output Slater Fk and Gk integrals
     call freshen( lmin, lmax, rcocc, skips )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,pheae,njrc,vi,zorig,xntot,nel,no, &
                   nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     call getfgnew( nr, dl, r, nel, nl, pheae, zorig, iuse, nuse, luse )
     !
     ! output radial melfiles
     call getmeznl( nr, dl, r, nel, nl, pheae, zorig, iuse, nuse, luse, npower )
     !
  end do
  !
  deallocate( pheps, pheae, orb, no, nl, nm, is, ev, ek, occ, xnj )
  !
  return
end subroutine mwreduce
