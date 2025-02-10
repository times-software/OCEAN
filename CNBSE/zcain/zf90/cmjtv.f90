! Copyright (C) 2010, 2024 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cmjtv( nsphpt, xsph, ysph, zsph, wsph, prefs )
  implicit none
  !
  integer :: nsphpt
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  !
  integer :: l, m, iproj, mc, npmax, nqproj, atno, nc, lc, lmin, lmax, ip, powmax, &
             i, nsemi, lsemi, ltest, ntest, j, stat
  real( kind = kind( 1.0d0 ) ) :: dqproj, qhat( 3 ), ehat( 3 ), ephotev, lam, pi, q, dummy
  integer, allocatable :: nproj( : )
  character(len=2) :: stringNL
  character(len=4) :: add04
  character(len=15) :: spcttype
  character(len=17) :: f17
  character(len=13) :: f13
  character(len=12) :: f12
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: ifcn
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: semifcn
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: c2c
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: atomEnergy
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: semiReducedEnergy
  complex( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: nbsemel
  integer, allocatable, dimension( :, : ) :: atomNL
  logical :: ex
  integer, parameter :: atomNLMax = 20
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  ! get core level info
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) atno, nc, lc
  close( unit=99 )
  write ( add04, '(1a1,1i3.3)' ) 'z', atno
  !
  ! input information about the projectors 
  call nbseprjprep( lmin, lmax, npmax, dqproj, nqproj, add04 )
  allocate( nproj( lmin : lmax ) )
  call nbseprjnproj( lmin, lmax, nproj, add04 )
  allocate( nbsemel( npmax, -lmax : lmax, lmin : lmax, -lc : lc ) )
  !
  ! input information about radial part of matrix elements
  write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', atno, 'n', nc, 'l', lc
  open( unit=99, file=f17, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) powmax   
  allocate( ifcn( npmax, 0 : powmax, lmin : lmax ) ) 
  ifcn( :, :, : ) = 0.0d0 
  do l = lmin, lmax
     do ip = 0, powmax
        read ( 99, * ) ifcn( 1 : nproj( l ), ip, l )
     end do
  end do 
  close( unit=99 )
  !
  ! input information about the probe
  open( unit=99, file='spectfile', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) spcttype
  call fancyvector( ehat, dummy, 99 )
  call fancyvector( qhat, q, 99 )
  read ( 99, * ) ephotev
  close( unit=99 )
  ehat( : ) = ehat( : ) / sqrt( sum( ehat( : ) ** 2 ) )
  qhat( : ) = qhat( : ) / sqrt( sum( qhat( : ) ** 2 ) )
  write ( 6, '(3(1f10.5))' ) ehat( : ), qhat( : )
  select case( spcttype )
  case( 'dipole', 's-dipole', 'p-dipole', 'd-dipole' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'quad', 's-quad', 'p-quad', 'd-quad' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'quadalone', 's-quadalone', 'p-quadalone', 'd-quadalone' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'NRIXS' )
      write(6,*) 'NRIXS is deprecated in favor of qRaman. Consider changing your inputs'
     ! nothing needs to be done!!!
  case( 'ldos0', 'ldos1', 'ldos2', 'ldos3', 'ldos00', 'ldos1-1', 'ldos10', 'ldos11', &
        'ldos2-2', 'ldos2-1', 'ldos20', 'ldos21', 'ldos22' )
     q = ephotev
  case( 'qRaman', 'qRs', 'qRp', 'qRd', 'qRf', 'tp', 'tpq', 'tp1', 'tp2' )
     ! nothing needs to be done!!!

  case default
     write(6,*) "Photon type was not recognized. Must be one of the following:"
     write(6,*) "      dipole, quad, quadalone, NRIXS, qRaman: ", spcttype
     stop

  end select
  write ( 6, '(1a4,1f10.5)' ) 'q = ', q
  !
  if( spcttype(1:2) .eq. 'tp' ) then
    write(f12, '(1A9I3.3)') 'corezetaz', atno
    open( unit=98, file=f12, form='formatted', status='old' )
    allocate( atomNL( 2, atomNLMax ), atomEnergy( atomNLMax ) )
    atomNL(:,:) = -1
    atomEnergy(:) = 0.0d0
    do i = 1, atomNLMax
      read(98,*,iostat=stat) dummy, dummy, stringNL, atomEnergy( i )
      if( stat == 0 ) then
        call decodeStringNL( stringNL, atomNL(:,i) )
      else
        exit
      endif
    enddo
    close( 98 )
        
    ! twp-photon radial info
    nsemi = 0
    lsemi = 1
    do i = 1,5
      write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', atno, 'n', i, 'l', lsemi
      inquire( file=f17, exist=ex )
      if( .not. ex ) then
        write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melsemi', 'z', atno, 'n', i, 'l', lsemi
        inquire( file=f17, exist=ex )
      endif
      if( ex ) nsemi = nsemi + 1
    enddo
    write(6,*) 'Found semi core options: ', nsemi
    allocate( semifcn( npmax, 0 : powmax, lmin : lmax, nsemi ), c2c( 0:powmax, nsemi) )
    allocate( semiReducedEnergy( nsemi ) )
    nsemi = 0
    do i = 1,5
      write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', atno, 'n', i, 'l', lsemi
      inquire( file=f17, exist=ex )
      if( .not. ex ) then
        write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melsemi', 'z', atno, 'n', i, 'l', lsemi
        inquire( file=f17, exist=ex )
      endif
      if( ex ) then
        nsemi = nsemi + 1
        open( unit=99, file=f17, form='formatted', status='old' )
        rewind 99
        read ( 99, * ) powmax
        do l = lmin, lmax
           do ip = 0, powmax
              read ( 99, * ) semifcn( 1 : nproj( l ), ip, l, nsemi )
           end do
        end do
        close( unit=99 )

        write( f13, '(1a3,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'c2c', 'z', atno, 'n', i, 'l', lsemi
        inquire( file=f13, exist=ex)
        if( .not. ex ) write( f13, '(1a3,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 's2c', 'z', atno, 'n', &
                                                                           i, 'l', lsemi
        open( unit=99, file=f13, form='formatted', status='old' )
        rewind(99)
        read(99,*) powmax
        do
          read(99,*) c2c(:,nsemi), ntest, ltest
          if( ntest .eq. nc .and. ltest .eq. lc ) exit
        enddo
        close(99)
        semiReducedEnergy( nsemi ) = 1.0d0
        do j = 1, atomNLMax
          if( atomNL(1,j) .eq. i .and. atomNL(2,j) .eq. lsemi ) then
            semiReducedEnergy( nsemi ) = -atomEnergy( j ) /  ( ephotev + atomEnergy( j ) )
            write( 6, * ) 'SRE', i, lsemi, semiReducedEnergy(nsemi), -atomEnergy( j )
            exit
          endif
        enddo
      endif
    enddo
  endif
  ! calculate full matrix elements
  call jtvsub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, ifcn, spcttype, ehat, qhat, q, &
       nsphpt, xsph, ysph, zsph, wsph, prefs, nsemi, semifcn, semiReducedEnergy, c2c )
  !
  ! output matrix elements in desired order
  open( unit=99, file='mels', form='formatted', status='unknown' )
  rewind 99
  do mc = -lc, lc
     do l = lmin, lmax
        do m = -l, l
           do iproj = 1, nproj( l )
              write ( 99, '(2(1x,1e15.8),5i4)' ) nbsemel( iproj, m, l, mc ), nproj( l ), iproj, &
                                                 l, m, mc
           end do
        end do
     end do
  end do
  close( unit=99 )
  !

!  deallocate( ifcn, semifcn, c2c, atomEnergy, semiReducedEnergy, nbsemel, nproj, atomNL )
  return

  contains

    subroutine decodeStringNL( st, nl )
      character( len=2 ), intent( in ) :: st
      integer, intent( out ) :: nl(2)
      !
      select case( st(2:2) )
        case( 's' )
          nl( 2 ) = 0
        case( 'p' )
          nl( 2 ) = 1
        case( 'd' )
          nl( 2 ) = 2
        case( 'f' )
          nl( 2 ) = 3
        case default
          nl( 2 ) = 4
      end select

      read(st(1:1),*) nl(1)
    end subroutine decodeStringNL

end subroutine cmjtv
