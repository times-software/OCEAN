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
      integer :: ii
      real( kind = kind( 1.0d0 ) ) :: de, eh, rcl, rch
!
      integer i,j,nel,nr,nst,iuflag,vtry,isuse,nwgt
      integer oldnel, newnel, nelmax
      integer l
      integer njrc( 4 )
      character * 8 :: nam
!
      double precision rel,alfa,etot,dl,zorig,xntot,rmin,rmax
      double precision rlast, potn, rad, rcut, vpert, dum, tmp
!
      logical done, psflag, rset, found
!
      character * 20 cmd
      character * 3 mode
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
      case ( 'o' )
        read ( stdin, * ) i, nam
        open( unit=99, file=nam, form='formatted', status='unknown' )
        rewind 99
        do j = 1, nr
           tmp = phe( j, i )
           if ( abs( tmp ) .lt. 1.0d-30 ) tmp = 0.0d0
           write ( 99,'(2(1e15.8,1x))' ) r( j ), tmp
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
      case ( 'ppdump' )
      case ( 'ppload' )
      case ( 'Q' )
        call alloy( vi, nr, r )
      case ( '!' ) ! another comment option
      case ( 'K' )
        call config
      case ( 'A' )
        read ( stdin, * ) l
        open( unit=99, file='radfile',
     &        form='formatted', status='unknown' )
        rewind 99
        read ( 99, * ) rcut
        rcut = rcut + 0.000001d0
        close( unit=99 )
        found = .false.
        i = 0
        do while ( .not. found )
           i = i + 1
           if ( l .eq. nl( i ) ) found = .true.
        end do
        call getanu( nl( i ), nr, r, rcut, dl, phe( 1, i ) )
      case ( 'a' )
        call abinitio(etot,rel,alfa,nr,r,
     &                dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,
     &                no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,
     &                psflag, nelmax )
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
        open( unit=99, file='spec',
     &        form='unformatted', status='unknown' )
        rewind 99
        write ( 99 ) phe( :, j )
        close( unit=99 )
      case ( 's' )
        read ( stdin, * ) l
        if ( 1 .lt. 2 ) stop 'sdump interface must be modified.'
        call sdump( l, psflag, nel, nl, nr, r, phe )
      case ( 'S' )
        read ( stdin, * ) ii, de, eh, rcl, rch
        call escanner( ii, nel, orb, nl, xnj, zorig, rel, nr, r, r2, 
     &        dl, njrc, vi, psflag, de, eh, rcl, rch )
      case ( 'W' )
        call melwriter( etot, rel, alfa, nr, r, dr, r2, dl,
     &                  njrc, vi, zorig, xntot, nel,
     &                  iuflag, cq, isuse )
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
        call pseudo(etot,rel,alfa,nr,rmin,
     &              rmax,r,dr,r2,dl,
     &              phe,orb,njrc,vi,cq,zorig,xntot,nel,
     &              no,nl,nm,xnj,ev,occ,is,ek,iuflag,vctab,
     &              vtry,isuse)
        newnel = nel
        write ( 6, * ) 'calling lead, after'
        call leadbehv(nel,nl,nr,phe,r,dl,.true.)
        write ( 6, * ) 'calling mkkxfile'
        call mkkxfile( oldnel, newnel )
        write ( 6, * ) 'done ... '
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
        open( unit=99, file='precpw',
     &        form='formatted', status='unknown' )
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
        call mkkbfile(nel,nl,xnj,ev,dl,nr,r,dr,r2,cq,
     &                vi,orb,zorig,rel,njrc,psflag)
      case ( 'F' )
        call sigfit(zorig,nr,r,dr,nel,nl,ev,phe)
      case ( 'f' )
        call sigfit2(nel,nl,ev,nr,dr,r,phe)
      case ( 'G' )
        call getfg( nr, dl, r, nel, nl, phe )
      case ( 'c' )
        call mkvctab( nr, vctab, r, r2 )
      case ( 'z' )
        call corepot( nr, vi )
      case ( 'Y' )
        open( unit=99, file='Real.Input', 
     &        form='formatted', status='unknown' )
        rewind 99
        do i = 1, nr
           read ( 99, * ) dum, vpert 
           if ( abs( dum - r( i ) ) .gt. 0.0001d0 ) stop 'bad r'
           do j = 1, 7 
              vi( i, j ) = vi( i, j ) + vpert
           end do
        end do
        close( unit=99 )
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
c-------------------------------------------------------------------
      subroutine sigfit( z, nr, r, dr, nlev, llev, eps, phi )
      implicit none
c
      integer nr, nlev
      double precision z
      integer llev( nlev )
      double precision r( nr ), dr( nr ), eps( nlev )
      double precision phi( nr, nlev )
c
      integer i, j, k, n
      double precision corpol, expt, rl, rh, sd, sl, sh, rc, sc
c
      double precision, external :: sig
c
      k = z + 0.1d0
      read  (5,*) corpol,n
      do i=1,n
        read  (5,*) j,expt,rl,rh
        sd=dabs(dabs(expt)-dabs(eps(j)))
        sl=sig(corpol,rl,nr,r,dr,phi(1,j))
        sh=sig(corpol,rh,nr,r,dr,phi(1,j))
        if (((sd-sh)*(sd-sl).ge.0.d0).and.
     &      (abs(rl-rh).gt.0.0001d0)      ) then
          rc=-1.d0
        else
  324     rc=rl+(rh-rl)/2.d0
          sc=sig(corpol,rc,nr,r,dr,phi(1,j))
          write (6,'(1x,8f9.4)') rl,sl,rh,sh,rc,sc,sd,eps(j)-sc
          if (sc.gt.sd) then
            rl=rc
            sl=sc
          else
            rh=rc
            sh=sc
          end if
          if ((abs(sc-sd).gt.0.000001d0).and.
     &        (abs(rl-rh).gt.0.0001d0)       ) go to 324
        end if
        write (6,'(1x,2i4,1f10.4)') k,llev(j),rc
        write (16,'(1x,2i4,1f10.4)') k,llev(j),rc
      end do
c
      return
      end
c-----------------------------------------------------------------
      subroutine sigfit2( nlev, llev, eps, nr, dr, r, phi )
      implicit none
      integer nlev, nr
      integer llev( nlev )
      double precision eps( nlev ), dr( nr ), r( nr )
      double precision phi( nr, nlev )
      integer i, k, l, n
      double precision lim
      double precision tmp, eav, e1, e2, corpol
      double precision sc, sd, sl, sh, rc, rl, rh
      character * 10 limunit, eunit
      read  (5,*) corpol
      read  (5,*) limunit, lim
      read  (5,*) i, n
      l = llev( i )
      if (n.eq.1) then
        read (5,*) eunit, eav
      else
        read (5,*) eunit, e1, e2
        eav=(e1*dble(l)+e2*dble(l+1))/dble(2*l+1)
      end if
      call convert( limunit, lim )
      call convert( eunit, eav )
      eav=dabs(lim)-dabs(eav)
      sd=dabs(dabs(eav)-dabs(eps(i)))
      rl= 0.d0
      rh=10.d0
      sl= 0.d0
      sh= 0.d0
      sc = sd + 1
      do while ( abs( sc - sd ) .gt. 0.000001d0 )
         if (sl*sh.le.0.00000001d0) then
            rc=rl+(rh-rl)/2.d0
         else
            rc=rl+(rh-rl)*(sd-sl)/(sh-sl)
         end if
         sc=0.d0
         do k=1,nr
            tmp=(1.d0-dexp(-(r(k)/rc)**2))**2
            tmp=corpol/(2.d0*r(k)**4)*tmp*tmp
            sc=sc+dr(k)*phi(k,i)*phi(k,i)*tmp
         end do
         if (sc.gt.sd) then
            rl=rc
            sl=sc
         else
            rh=rc
            sh=sc
         end if
      end do
      write (6,'(1x,3f14.4)') rc,sc,sc*27.2114d0
      return
      end
c------------------------------------------------------------
      subroutine convert( unit, value )
      implicit none
      character * 10 unit
      double precision value, mult
      mult = -1.d0
      if ( unit .eq. 'hartree' ) mult = 1.d0
      if ( unit .eq. 'Ryd' ) mult = 0.5d0
      if ( unit .eq. 'eV' ) mult = 1.d0 / 27.2114d0
      if ( unit .eq. 'invcm' ) mult = 0.000123985d0 / 27.2114d0
      if ( mult .lt. 0.d0 ) stop 'unsupported unit'
      value = value * mult
      unit = 'hartree' 
      return
      end
c------------------------------------------------------------
      subroutine mkvctab( nr, vctab, r, rsqd )
      implicit none
      integer nr
      double precision vctab( nr, 0 : 3 ), r( nr ), rsqd( nr )
      integer k
      double precision corpol, rs, rp, rd
      double precision pref, rssqd, rpsqd, rdsqd, fssqd, fpsqd, fdsqd
      read (5,*) corpol,rs,rp,rd
      rssqd = rs ** 2
      rpsqd = rp ** 2
      rdsqd = rd ** 2
      do k=1,nr
        pref = - 0.5d0 * corpol / r( k ) ** 4
        fssqd=(1.d0-dexp(-rsqd(k)/rssqd))**4
        fpsqd=(1.d0-dexp(-rsqd(k)/rpsqd))**4
        fdsqd=(1.d0-dexp(-rsqd(k)/rdsqd))**4
        vctab(k,0)= pref * fssqd
        vctab(k,1)= pref * fpsqd
        vctab(k,2)= pref * fdsqd
        vctab(k,3)= 0.d0
      end do
      return
      end
