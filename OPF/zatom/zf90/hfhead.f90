! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program hartfock
  implicit none
  !
  integer, allocatable :: no(:),nl(:),nm(:),is(:)
  integer, allocatable :: ilp(:)
!
  double precision, allocatable :: r(:),dr(:),r2(:)
  double precision, allocatable :: ev(:),occ(:),xnj(:)
  double precision, allocatable :: ek(:),phe(:,:),orb(:,:)
  double precision, allocatable :: vi(:,:),cq(:),vctab(:,:)
  double precision, allocatable :: wgts(:)
!
  integer :: ii, ic
  real( kind = kind( 1.0d0 ) ) :: de, eh, rcl, rch
!
  integer i,j,nel,nr,nst,iuflag,vtry,isuse,nwgt, l, njrc( 4 )
  integer oldnel, newnel, nelmax
  character( len = 8 ) :: nam
!
  integer :: lmin, lmax, skips( 0 : 3 ), idum
  double precision rel,alfa,etot,dl,zorig,xntot,rmin,rmax, zc, mul, pertstr
  real( kind = kind( 1.0d0 ) ) :: rlast, potn, rad, vpert, dum, tmp, rcocc( 0 : 3 ), rsocc( 0 : 3 )
!
  logical :: done, psflag, rset
!
  character( len = 20 ) :: cmd
  character( len = 3 ) :: mode
 !
!     logical found
!     double precision rcut
! integer :: ntest, nsmax, irc
! real( kind = kind( 1.0d0 ) ) :: cush, emax, rrmd, prec
!
  integer, parameter :: nwgmx = 10, stdin = 5
!
!
      allocate( wgts( nwgmx ) )
!
      rset = .false.
!
      open( unit=77, file='tmp', form='formatted', status='unknown' )
      rewind 77
      rel = 0.d0
      nst = 2
      
      read ( stdin, * ) cmd
!
!====================================
      do while ( cmd .ne. 'quit' )
!====================================
!
      select case( cmd )
      case ( '#' )
        ! this is the comment option'
      case ( 'output_psi' )
        read ( stdin, * ) i, nam, mul
        open( unit=99, file=nam, form='formatted', status='unknown' )
        rewind 99
        do j = 1, nr
           tmp = phe( j, i )
           if ( abs( tmp ) .lt. 1.0d-30 ) tmp = 0.0d0
           write ( 99,'(3(1e15.8,1x))' ) r( j ), tmp * mul, ev( i )
        end do
        close( unit=99 )
      case ( 'b' )
        call bachelet(vi,r,njrc,nr)
      case ( 'X' )
        read (5,'(2x,1i4,2x,1a3)') nwgt,mode
        if (nwgt.gt.nwgmx) then
          write (6,*) 'bad nwgt'
          stop
        end if
        read (5,*) (wgts(i),i=1,nwgt)
        close(unit=77)
        call grandop(nwgt,wgts,mode)
        open(unit=77,file='tmp',form='formatted',status='unknown')
        rewind 77
      case ( 'd' )
        read (5,*) rel
      case ( 'j' )
        call besft( nr, r, dr, nel, nl, phe )
      case ( 'x' )
        read  (5,*) alfa
      case ( 'Q' )
        call alloy( vi, nr, r )
      case ( '!' ) ! another comment option
      case ( 'K' )
        call config
!     case ( 'A' )
!       read ( stdin, * ) l
!       open( unit=99, file='radfile', form='formatted', status='unknown' )
!       rewind 99
!       read ( 99, * ) rcut
!       rcut = rcut + 0.000001d0
!       close( unit=99 )
!       found = .false.
!       i = 0
!       do while ( .not. found )
!          i = i + 1
!          if ( l .eq. nl( i ) ) found = .true.
!       end do
!       call getanu( nl( i ), nr, r, rcut, dl, phe( 1, i ) )
      case ( 'a' )
        write ( 6, * ) 'nelmax = ', nelmax
        call abinitio(etot,rel,alfa,nr,r, dr,r2,dl,phe,njrc,vi,zorig,xntot,nel, &
             no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq, psflag, nelmax )
      case ( 'initgrid' )
         read ( stdin, * ) zorig, nr, nelmax
         if ( rset ) then
            deallocate( r, dr, r2, phe, orb, vi, cq, vctab )
            deallocate( no, nl, nm, is, ilp, ev, occ, xnj, ek )
         else
            rset = .true.
         end if
         allocate( r( nr ),dr( nr ), r2( nr ) )
         allocate( phe( nr, nelmax ), orb( nr, nelmax ) )
         allocate( vi( nr, 7 ), cq( nr ), vctab( nr, 0 : 3 ) )
         allocate( no( nelmax ), nl( nelmax ) )
         allocate( nm( nelmax ), is( nelmax ) )
         allocate( ilp( nelmax ) )
         allocate( ev( nelmax ), occ( nelmax ) )
         allocate( xnj( nelmax ), ek( nelmax ) )
         call initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,xntot,nel)
         vctab( :, : ) = 0.0d0
      case ( 'e' )
        read ( stdin, * ) j
        open( unit=99, file='spec', form='unformatted', status='unknown' )
        rewind 99
        write ( 99 ) phe( :, j )
        close( unit=99 )
! defunct defunct defunct
!     case ( 's' )
!       read ( stdin, * ) l
!       if ( 1 .lt. 2 ) stop 'sdump interface must be modified.'
!       call sdump( l, psflag, nel, nl, nr, r, phe )
! defunct defunct defunct
      case ( 'S' )
        read ( stdin, * ) ii, ic, de, eh, rcl, rch
        call escanner( ii, phe( :, ic ), nel, orb, nl, xnj, zorig, rel, nr, r, r2, dl, njrc, vi, psflag, de, eh, rcl, rch )
!     case ( 'W' )
!       call melwriter( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse )
      case ( 'u' )
        write (6,*) 'please enter iuflag. (0=u, 1=su, 2=r).'
        read ( stdin, * ) iuflag
      case ( 'P' )
        read ( stdin, * ) psflag
      case ( 'p' )
        write ( 6, * ) 'calling lead, before'
        call leadbehv(nel,nl,nr,phe,r,dl,.false.)
        oldnel = nel
        write ( 6, * ) 'calling pseudo'
        call pseudo(etot,rel,alfa,nr,rmin, rmax,r,dr,r2,dl, phe,orb,njrc,vi,cq,zorig,xntot,nel, &
                    no,nl,nm,xnj,ev,occ,is,ek,iuflag,vctab, vtry,isuse)
        read ( 5, * ) nel
        newnel = nel
!       call leadbehv(nel,nl,nr,phe,r,dl,.true.)
!       call mkkxfile( oldnel, newnel )
        rcocc( : ) = 0.0d0
        lmin = nl( 1 ); lmax = nl( 1 )
        do i = 1, nel
           rcocc( nl( i ) ) = occ( i )
           lmin = min( lmin, nl( i ) ) 
           lmax = max( lmin, nl( i ) ) 
        end do
        write ( 6, '(1a12,4f10.4)' ) 'refconocc = ', rcocc( : )
        write ( 6, '(1a7,1i5)' ) 'lmin = ', lmin
        write ( 6, '(1a7,1i5)' ) 'lmax = ', lmax
!     case( 'fillinpaw' )
!       call fillinpaw( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rcocc, lmin, lmax )
      case( 'continuum' )
         call continuum( nel, no, nl, rel, nr, r, r2, dr, dl, phe, orb( :, 1 ), zorig, njrc, vi, psflag )
      case( 'projaug' )
         call projaug( nr, r, dl )
      case( 'trck' )
         call trck( nel, nl, phe, nr, r, dl )
      case( 'spartanfip' )
         call spartanfip( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rcocc, lmin, lmax )
      case( 'screencore' )
         call screencore( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rsocc, lmin, lmax )
      case( 'ppdump' )
         call ppdump( njrc, vi, nr, r )
      case( 'ppload' )
         call ppload( njrc, vi, nr, r, zc )
         write ( 6, '(4i8,5x,1a4)' ) njrc( : ), 'njrc'
         do l = 0, 3
            if ( njrc( l + 1 ) .ne. 0 ) lmax = l
         end do
         do l = 3, 0, -1
            if ( njrc( l + 1 ) .ne. 0 ) lmin = l
         end do
         write ( 6, '(1a7,1i5)' ) 'lmin = ', lmin
         write ( 6, '(1a7,1i5)' ) 'lmax = ', lmax
      case( 'calcso' )
        call mkcorcon( alfa, rel, zorig, zc, rcocc, rsocc, .true. )
        open( unit=99, file='skip', form='formatted', status='unknown' )
        rewind 99
        read ( 99, * )
        do l = 0, 3
           read ( 99, * ) idum, skips( l )
        end do
        close( unit=99 )
        call freshen( lmin, lmax, rcocc, skips, 1, -0.01d0 )
        call abinitio(etot,rel,alfa,nr,r, dr,r2,dl,phe,njrc,vi,zorig,xntot,nel, &
             no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
!     case ( 'fipfront' )
!        read ( 5, * ) l, cush, emax, rcut, nsmax, rrmd, prec, ntest
!        call pwavesetter( l, ntest, cush, emax, rcut, nsmax, rrmd, prec, irc, nr, r, dr, r2, dl, vi, iuflag, cq, isuse, njrc )
!        call fipfront( lmin, lmax, nr, r, dr, r2, dl, zorig, vi, iuflag, cq, njrc )
      case( 'mkcorcon' )
        call mkcorcon( alfa, rel, zorig, zc, rcocc, rsocc, .false. )  
        write ( 6, '(1a12,4f10.4)' ) 'refconocc = ', rcocc( : )
      case ( 'C' )
        call vpscpp( nr, r, r2, vi )
      case ( 'v' )
        call realspace(nr,r,vi)
      case ( 'V' )
        call fourier(nr,r,dr,r2,vi)
      case ( 'g' )
        call ppopt(nr,r,vi,cq,nel,phe,occ)
        read ( stdin, * ) rad
        write ( 6, * ) 'rad = ', rad
        open( unit=99, file='precpw', form='formatted', status='unknown' )
        rewind 99
        write ( 99, * ) nel, rad
        write ( 6, * ) nel, rad
        do i = 1, nel
          write ( 6, * ) i, nel
          write ( 99, * ) nl( i )
          rlast = r( 1 )
          j = 2         
          done = .false.
          do while ( .not. done )
            if ( r( j ) .gt. rlast + 0.01d0 ) then
              potn = vi( j, 2 * nl( i ) + 1 ) + orb( j, i )
              write ( 99, '(3(2x,1e15.8))' ) r( j ), phe( j, i ), potn
              rlast = r( j )
              done = ( rlast .gt. rad )
            end if
            j = j + 1
          end do
        end do
      case ( 'k' )
        call mkkbfile(nel,nl,xnj,ev,dl,nr,r,dr,r2,cq, vi,orb,zorig,rel,njrc,psflag)
      case ( 'F' )
        call sigfit(zorig,nr,r,dr,nel,nl,ev,phe)
      case ( 'f' )
        call sigfit2(nel,nl,ev,nr,dr,r,phe)
!     case ( 'G' )
!       call getfg( nr, dl, r, nel, nl, phe, zorig )
      case ( 'c' )
        call mkvctab( nr, vctab, r, r2 )
      case ( 'z' )
        call corepot( nr, vi )
      case ( 'Y' )
        read ( 5, * ) pertstr
        open( unit=99, file='Real.Input', form='formatted', status='unknown' )
        rewind 99
        do i = 1, nr
           read ( 99, * ) dum, vpert 
           if ( abs( dum - r( i ) ) .gt. 0.0001d0 ) stop 'bad r'
           vi( i, : ) = vi( i, : ) + vpert * pertstr
        end do
        close( unit=99 )
      case( 'moment' )
         call orbmom( nr, r, dl, nel, phe )
      case ( 'Z' )
        do i = 1, nr
           do j = 1, 7 
              vi( i, j ) = vi( i, j ) - 1 / r( i )
           end do
        end do
      end select
!
!====================================
        read ( stdin, * ) cmd
      end do
      write ( 6, * ) 'hartfock terminus achieved'
      stop
!====================================
!
      end
