  inquire(file='core_offset', exist=exst)
  if( exst ) then
    open(unit=99,file='core_offset',  form='formatted', status='old' )
    rewind 99
    read( 99, * ) core_offset
    core_offset = core_offset / 27.2114d0
    close( 99 )
  else
    core_offset = 0
  endif
  pi = 4.0d0 * atan( 1.0d0 )
  open( unit=99, file='epsilon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) eps
  close( unit=99 )
  xi = xi / 27.2114d0
  nc = 4 * ( 2 * lc + 1 )
  write ( s11, '(1a7,1a4)' ) 'prjfile', add04
  open( unit=99, file=s11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lmin, lmax
  allocate( nproj( lmin : lmax ) )
  read ( 99, * ) nproj
  close( unit=99 )
  if ( ( lvl .lt. lmin ) .or. ( lvh .gt. lmax ) ) stop 'bad lvl lvh'
  !
  ! generalize to multiplet hamiltonians for several values of lv
  allocate( ibeg( lvl : lvh ), jbeg( lvl : lvh ), mham( lvl : lvh ) )
  itot = 0
  jtot = 0
  do lv = lvl, lvh
     ibeg( lv ) = itot + 1
     jbeg( lv ) = jtot + 1
     mham( lv ) = 4 * ( 2 * lc + 1 ) * ( 2 * lv + 1 ) * nproj( lv )
     itot = itot + mham( lv )
     jtot = jtot + mham( lv ) ** 2
  end do
  allocate( cml( nc ), cms( nc ), vms( nc ) )
  allocate( hcml( itot ), hcms( itot ) )
  allocate( hvml( itot ), hvms( itot ) )
  allocate( hvnu( itot ) )
  allocate( mhr( jtot ), mhi( jtot ) )
  allocate( pwr( itot ), pwi( itot ) )
  allocate( hpwr( itot ), hpwi( itot ) )
  ! 
  open( unit=99, file='channelmap', form='formatted', status='unknown' )
  rewind 99
  ic = 0
  do icms = -1, 1, 2
     do icml = -lc, lc
        do ivms = -1, 1, 2
           ic = ic + 1 
           cms( ic ) = 0.5d0 * icms
           cml( ic ) = icml
           vms( ic ) = 0.5d0 * ivms
           write ( 99, '(3f8.2,2i5,1a11)' ) cms( ic ), cml( ic ), vms( ic ), ic, nc, 'cms cml vms'
        end do
     end do
  end do
  !
  ii = 0 
  do lv = lvl, lvh
     do ic = 1, 4 * ( 2 * lc + 1 )
        do ivml = -lv, lv
           do nu = 1, nproj( lv )
              ii = ii + 1
              hcms( ii ) = cms( ic )
              hcml( ii ) = cml( ic )
              hvms( ii ) = vms( ic )
              hvml( ii ) = ivml
              hvnu( ii ) = nu
           end do
        end do
     end do
  end do
  !
  call nbsemkcmel( add04 )
  do lv = lvl, lvh
     ii = ibeg( lv )
     jj = jbeg( lv )
     call nbsemhsetup( lc, lv, nproj( lv ), mham( lv ), hcms( ii ) , hcml( ii ), hvms( ii ), hvml( ii ), hvnu( ii ), &
          mhr( jj ), mhi( jj ), add10 )
  end do
  mhr = mhr / 27.2114d0
  mhi = mhi / 27.2114d0
  write ( 6, * ) 'multiplet hamiltonian set up'
  !
  ! !!!!!!!!!
  ! In time this will be deprecated by making sure all the helper files 
  !   always say 1 or 2 spins. For now, check here.
  inquire(file='nspn.ipt', exist = nspn_exist )
  if( nspn_exist ) then
    open(unit=99,file='nspn.ipt',form='formatted',status='old')
    read(99,*) nspn
    close(99)
  else
    nspn = 1
  endif
  write(6,*) 'Number of spins = ', nspn
  !
  open(unit=99, file='cks.normal', form='formatted', status='old')
  rewind 99
  read(99, * ) conduct
  close(99)
  ! read in generic information
  if (conduct) then
    infoname = 'wvfcninfo'
  else
    infoname = 'wvfvainfo'
  endif
  open( unit=99, file=infoname, form='unformatted', status='unknown' )
  rewind 99
  if( nspn .eq. 2 ) then
    read( 99 ) nbd, nq, nspn
  else
    read ( 99 ) nbd, nq
  endif
  n = nq * nbd
  write ( 6, * ) 'n, nc, nspn', n, nc, nspn
  allocate( iq( 3, nq ), e0( n, nspn ), v( n, nc, 2 ), hv( n, nc, 2 ) )
!  if (.not. allocated (cor) ) stop
!  cor = 0; coi = 0
  allocate( qphys( 3, nq ) )
  v( :, :, : )  = 0
  read ( 99 ) e0( :, : )
  close( unit=99 )
  !
  ic = 0
  if( nspn .eq. 1 ) then
    do icms = -1, 1, 2
      do icml = -lc, lc
        write ( str, '(1a4,1i2.2)' ) 'beff', 1 + icml + lc
        open( unit=99, file=str, form='unformatted', status='unknown' )
        rewind 99
        read ( 99 ) tau( : )
        do ivms = -1, 1, 2
           ic = ic + 1
           write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
           if ( icms .eq. ivms ) then
              do j = 1, n
                 read ( 99 ) v( j, ic, 1 ), v( j, ic, 2 )
              end do 
           end if
        end do
        close( unit=99 )
      end do
    end do
  else
    do icms = 1, 2
      do icml = -lc, lc
        do ivms = 1, 2
           write ( str, '(1a4,1i2.2,1a1,1i1.1)' ) 'beff', 1 + icml + lc, '.', ivms
           open( unit=99, file=str, form='unformatted', status='unknown' )
           rewind 99
           read ( 99 ) tau( : )
           ic = ic + 1
           write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
           if ( icms .eq. ivms ) then
              do j = 1, n
                 read ( 99 ) v( j, ic, 1 ), v( j, ic, 2 )
              end do
           end if
           close( unit=99 )
        end do
      end do
    end do
  endif
  !
  write (6,*) 'band states have been read in'
  npmax = maxval( nproj( lmin : lmax ) )
! v( :, :, : ) = 0.0d0
! call vfilter( nq * nbd, nc, v, muroot, lmin, lmax, nproj, npmax, lc )
  !
  allocate( mpcr( nq * nbd, npmax, -lmax : lmax, lmin : lmax, nspn ) )
  allocate( mpci( nq * nbd, npmax, -lmax : lmax, lmin : lmax, nspn ) )
  allocate( ampr( npmax ), ampi( npmax ), hampr( npmax ), hampi( npmax ) )
  !
  ip = 0
  do l = lmin, lmax
     do m = -l, l
        do nu = 1, nproj( l )
           ip = ip + 1
        end do
     end do
  end do
  nptot = ip
! Spin needed here too!
  allocate( pcr( nptot, n, nspn ), pci( nptot, n, nspn ) )
  open( unit=99, file='ufmi', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 )
  read ( 99 )
  read ( 99 ) pcr
  read ( 99 ) pci
  do ispn = 1, nspn
    do i = 1, nq * nbd
      ip = 0
      do l = lmin, lmax
        do m = -l, l
           do nu = 1, nproj( l )
              ip = ip + 1
              mpcr( i, nu, m, l, ispn ) = pcr( ip, i, ispn )
              mpci( i, nu, m, l, ispn ) = pci( ip, i, ispn )
           end do
        end do
      end do
    end do
  enddo
  close( unit=99 )
  write ( 6, * ) 'projector coefficients have been read in'
  !
  allocate( list( npmax * npmax ), mpm( npmax, npmax, lmin : lmax ) )
  mpm( :, :, : ) = 0
  do l = lmin, lmax
     if ( l .gt. 9 ) stop 'l exceeds 9 in multip' 
     write ( s5, '(1a4,1i1)' ) 'cmel', l
     open( unit=99, file=s5, form='formatted', status='unknown' )
     rewind 99
     read( 99, * ) list( 1: nproj( l ) * nproj( l ) )
     i = 0
     do nu = 1, nproj( l )
        mpm( 1 : nproj( l ), nu, l ) = list( i + 1 : i + nproj( l ) )
        i = i + nproj( l )
     end do
     close( unit=99 )
  end do
  mpm( :, :, : ) = mpm( :, :, : ) / 27.2114d0
  write ( 6, * ) 'central projector matrix elements have been read in'
  !
  allocate( ur( nx * ny * nz, nbd, nq, nspn ), ui( nx * ny * nz, nbd, nq, nspn ) )
  call loaduxzerotau( nx, ny, nz, nbd, nq, nspn, zn, tau, ur, ui )
  write ( 6, * ) 'ur has been read in'
  !
  open(unit=99,file='rpottrim',form='formatted',status='unknown')
  rewind 99
  do i = 1, 100
     read( 99, * ) ptab( i )
  end do
  close(unit=99)
  !
  do i = 1, 3
     do j = 1, 3
        amet( i , j ) = dot_product( avec( :, i ), avec( :, j ) )
     end do
  end do
  !
  call getomega( avec, celvol )
  pref = 1.0d0 / ( dble( nq ) * celvol )
  rzero = ( 3.0d0 / ( 4.0d0 * pi ) * celvol / dble( nx * ny * nz ) ) ** ( 1.0d0 / 3.0d0 )
  rcut = ( 0.2d0 * rzero ) ** 2
  !
  allocate( somel( nc, nc, 2 ) )
