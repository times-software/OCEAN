subroutine melwriter( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse )
  implicit none
  !
  double precision :: etot, rel, alfa
  double precision :: dl, zorig, xntot
  integer :: nr, njrc( 4 ), nel, isuse
  integer :: iuflag
  double precision :: r( nr ),dr( nr ),r2( nr )
  double precision :: vi( nr, 7 ), cq( nr )
  !
  integer, allocatable :: no( : ), nl( : ), nm( : )
  integer, allocatable :: is( : )
  double precision, allocatable :: ev( : ), occ( : ), ek( : )
  double precision, allocatable :: xnj( : )
  !
  character * 3 :: nam3
  character * 4 :: nam4
  character * 6 :: meltyp
  integer :: icore, lcore, nq, nco, n2
  integer :: lmin, lmax
  integer :: l, ne, lcmin, lcmax, k, is2
  integer :: i, j, lev, ilev, jlev, irc
  integer :: tfac, nt, i1, i2, j1, j2, lc
  integer :: idum
  double precision, allocatable :: etmin( : )
  double precision :: emin, emax, emid
  double precision :: dq, qmax, prec, prec2
  double precision :: rc, ctol, condnum
  double precision :: su, dummy
  double precision, external :: efrac
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: c( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: meltmp( : )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: orb, proj, pheae, pheps, aesav
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: meltrue, melfake, melval
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: etrue, efake, exch
  !
  integer, parameter :: prjfile = 93, melfile = 94, excfile = 95
  !
  integer  :: ntest, nelmax
  double precision, allocatable :: eint( : ), frac( : )
  double precision, allocatable :: ff( : ), copy( : ), reconned( : )
  !
  double precision, allocatable :: phespc( : )
  logical :: special, okay
  !
  integer :: nval, llo, lhi, numl
  double precision :: xj, q
  integer, allocatable, dimension( : ) :: ntmp, ltmp
  double precision, allocatable :: occtmp( : )
  !
  allocate( ff( nr ), copy( nr ), reconned( nr ) )
  !
  read ( 5, '(1a6)' ) meltyp
  okay = .false.
  if ( meltyp .eq. 'kshell' ) meltyp = 'ssubsh'
  if ( meltyp .eq. 'psubsh' ) meltyp = 'dipole'
  if ( meltyp .eq. 'dipole' ) okay = .true.
  if ( meltyp .eq. 'ssubsh' ) okay = .true.
  if ( meltyp .eq. 'longit' ) then
     okay = .true.
     read ( 5, * ) q
  end if
  if ( .not. okay ) stop 'bad mel type'
  !
  allocate( phespc( nr ) )
  read ( 5, * ) icore
  special = ( icore .lt. 0 )
  if ( special ) then
     icore = abs( icore )
     open ( unit=99, file='spec', form='unformatted', status='unknown' )
     rewind 99
     read ( 99 ) phespc
     close( unit=99 )
     write ( 6, * ) 'we are special!!!'
  end if
  read ( 5, * ) lmin, lmax
  allocate( etmin( lmin : lmax ) )
  read ( 5, * ) etmin
  read ( 5, * ) emin, emax, prec, prec2, tfac
  read ( 5, * ) rc, ctol
  irc = 0
  do i = 1, nr, 4
     if ( r( i ) .lt. rc ) irc = i
  end do
  open( unit=99, file='radfile', form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) r( irc ), nr
  close( unit=99 )
  read ( 5, * ) dq, qmax
  !
  open( unit=prjfile, file='prjfile', form='formatted', status='unknown' )
  open( unit=melfile, file='melfile', form='formatted', status='unknown' )
  open( unit=excfile, file='excfile', form='formatted', status='unknown' )
  rewind prjfile 
  rewind melfile
  rewind excfile
  !
  ntest = 80
  allocate( eint( 0 : ntest ), frac( 0 : ntest ) )
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  do i = 1, icore
     read ( 99, * ) idum, lcore 
  end do
  write ( 6, '(1a8,2x,1i4)' ) 'lcore = ', lcore
  close( unit=99 )
  open( unit=99, file='deflin', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(3i5)' ) lcore, lmin, lmax
  close( unit=99 )
  open( unit=99, file='config', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) n2
  close( unit=99 )
  nelmax = 1 + ntest + nco + n2 + 10
  allocate( no( nelmax ), nl( nelmax ) )
  allocate( nm( nelmax ), is( nelmax ) )
  allocate( ev( nelmax ), ek( nelmax ) )
  allocate( occ( nelmax ), xnj( nelmax ) )
  allocate( pheps( nr, nelmax ), pheae( nr, nelmax ) )
  allocate( orb( nr, nelmax ) )
  !
  nq = 1 + qmax / dq
  write ( prjfile, '(3i5,2x,1e15.8)' ) lmin, lmax, nq, dq
  !
  ! loop over values of valence angular momentum ... 
  do l = lmin, lmax
     !
     ! detmerine suitable project functions
     call egripper2( nr, nelmax, r, dr, r2, dl, nl, &
          &                pheps, pheae, njrc, vi, zorig, iuflag, orb, &
          &                l, ntest, irc, frac, rel, &
          &                eint, emin, emax, ne, prec, prec2 ) 
     write ( 6, *)  'project functions found'
     write ( 6, * ) 'l, ne ... ', l, ne
     !
     ! output projector functions and their fourier transforms 
     nl( 1 : ne ) = l 
     call sdump( l, .true., ne, nl, nr, r, pheps, nelmax )
     nl( nco + 1 : nco + ne ) = l
     call sdump( l, .false., ne, nl, nr, r, pheae, nelmax )
     !
     allocate( aesav( irc, ne ) )
     do i = 1, ne
        do j = 1, irc
           aesav( j, i ) = pheae( j, i + nco )
        end do
     end do
     call breader( l, irc, condnum, nr, dl, ne )
     !
     lcmin = iabs( lcore - l )
     lcmax = lcore + l 
     write ( 6, '(1a14,2i5)' ) 'lcmin, lcmax: ', lcmin, lcmax
     write ( prjfile, '(1i5)' ) ne
     write ( excfile, '(2i5)' ) lcmin, lcmax
     allocate( exch( lcmin : lcmax, ne, ne ) )
     if ( .not. special ) phespc = pheae( : , icore )
     do i = 1, ne
        ilev = i + nco
        do j = 1, ne
           jlev = j + nco
           call getexch( nr, r, dl, lcmin, lcmax, &
                &                    phespc, &
                &                    pheae( 1, ilev ), pheae( 1, jlev ), &
                &                    exch( lcmin, i, j ) )
           write ( excfile, '(2i3,7(2x,1f12.8))' ) &
                &        i, j, ( exch( lc, i, j ), lc = lcmin, lcmax, 2 )
        end do
     end do
     if ( meltyp .eq. 'longit' ) then
        llo = abs( l - lcore )
        lhi = l + lcore
        numl = 1 + ( lhi - llo ) / 2
     else
        numl = 1
     end if
     nt = tfac * ne
     allocate( melval( numl, ne ), meltrue( numl, nt ), melfake( numl, nt ), meltmp( numl ) )
     do i = 1, ne
        lev = i + nco
        select case( meltyp )
        case( 'dipole' )
           call dipole( nr, r, dl, phespc, pheae( 1, lev ), meltmp )
        case( 'ssubsh' )
           call ssubsh( nr, r, dl, phespc, pheae( 1, lev ), meltmp, l )
        case( 'longit' )
           do j = 1, numl
              call longit( nr, r, dl, phespc, pheae( 1, lev ), meltmp( j ), llo + 2 * ( j - 1 ), q )
           end do
        end select
        do j = 1, numl
           write ( 6, '(2x,2i5,2x,1e15.8,2x,1a6)' ) i, j, meltmp( j ), 'melval'
        end do
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
        xj = - ltmp( i )
        if ( ltmp( i ) .eq. 0 ) xj = -0.5d0
        write ( 99, * ) ntmp( i ), ltmp( i ), 0, xj, 1, occtmp( i )
     end do
     close( unit=99)
     deallocate( ntmp, ltmp, occtmp )
     !
     call abinitio(etot,rel,alfa,nr,r,dr,r2,dl, &
          &              pheae,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj, &
          &              ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     !
     allocate( etrue( lcmin : lcmax, nt, nt ) )
     do i = 1, nt
        ilev = i + ( nel - nval - nt )
        do j = 1, nt
           jlev = j + ( nel - nval - nt )
           call getexch( nr, r, dl, lcmin, lcmax, &
                &                    phespc, &
                &                    pheae( 1, ilev ), pheae( 1, jlev ), &
                &                    etrue( lcmin, i, j ) )
        end do
     end do
     do i = 1, nt
        lev = i + ( nel - nval - nt )
        select case( meltyp )
        case( 'dipole' )
           call dipole( nr, r, dl, phespc, pheae( 1, lev ), meltmp )
        case( 'ssubsh' )
           call ssubsh( nr, r, dl, phespc, pheae( 1, lev ), meltmp, l )
        case( 'longit' )
           do j = 1, numl
              call longit( nr, r, dl, phespc, pheae( 1, lev ), meltmp( j ), llo + 2 * ( j - 1 ), q )
           end do
        end select
        meltrue( :, i ) = meltmp( : )
     end do
     !
     call amult( l, ne, nr )
     allocate( proj( nr, ne ), c( ne, nt ) )
     write ( unit=nam3, fmt='(1a2,1i1)' ) 'pr', l
     open( unit=99, file=nam3, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nr
        read ( 99, * ) dummy, ( proj( i, j ), j = 1, ne )
     end do
     close( unit=99 )
     !
     call transp( l, ne, nr, qmax, dq, nq, dl, irc )
     !
     select case( meltyp )
     case ( 'dipole' )
        if ( ( l .eq. lcore + 1 ) .or. ( l .eq. lcore - 1 ) ) then
           write ( melfile, '(8(2x,1e15.8))' ) melval( :, : )
        else
           write ( melfile, '(8(2x,1e15.8))' ) ( 0.d0, i = 1, ne )
        end if
     case( 'ssubsh' )
        if ( ( l .eq. 1 ) .or. ( l.eq. 2 ) ) then
           write ( melfile, '(8(2x,1e15.8))' ) melval( :, : )
        else
           write ( melfile, '(8(2x,1e15.8))' ) ( 0.d0, i = 1, ne )
        end if
     case( 'longit' )
        do j = 1, numl 
           write ( melfile, '(8(2x,1e15.8))' ) melval( j, : )
        end do
     end select
     !
     call abinitio(etot,rel,alfa,nr,r,dr,r2,dl,pheps,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.true.,nelmax)
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
        do j = 1, numl
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
     do j = 1, numl
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
  close( unit=excfile )
  close( unit=melfile )
  !
  deallocate( pheps, pheae, orb, no, nl, nm, is, ev, ek, occ, xnj )
  !
  return
end subroutine melwriter
