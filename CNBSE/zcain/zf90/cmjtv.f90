! Copyright (C) 2010 OCEAN collaboration
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
  integer :: l, m, iproj, mc, npmax, nqproj, atno, nc, lc, lmin, lmax, ip, powmax
  real( kind = kind( 1.0d0 ) ) :: dqproj, qhat( 3 ), ehatr( 3 ), ehati( 3 )
  real( kind = kind( 1.0d0 ) ) :: ephotev, lam, pi, q, dummy, enorm
  complex( kind = kind( 1.0d0 ) ) :: rm1, ehat( 3 )
  integer, allocatable :: nproj( : )
  character(len=4) :: add04
  character(len=15) :: spcttype
  character(len=17) :: f17
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: ifcn
  complex( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: nbsemel
  !
  pi = 4.0d0 * atan( 1.0d0 )
  rm1 = -1.0d0
  rm1 = sqrt( rm1 )
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
  call fancyvector( ehatr, dummy, 99 )
  call fancyvector( ehati, dummy, 99 )
  call fancyvector( qhat, q, 99 )
  read ( 99, * ) ephotev
  close( unit=99 )
  ehat( 1 ) = ehatr( 1 ) + rm1 * ehati( 1 )
  ehat( 2 ) = ehatr( 2 ) + rm1 * ehati( 2 )
  ehat( 3 ) = ehatr( 3 ) + rm1 * ehati( 3 )
  enorm = ehat(1)*conjg(ehat(1)) + ehat(2)*conjg(ehat(2)) + ehat(3)*conjg(ehat(3))
  ehat(:) = ehat(:) / sqrt( enorm )
!  ehat( : ) = ehat( : ) / sqrt( sum( ehat( : ) ** 2 ) )
  qhat( : ) = qhat( : ) / sqrt( sum( qhat( : ) ** 2 ) )
  write ( 6, '(3(1f10.5))' ) ehat( : ), qhat( : )
  select case( spcttype )
  case( 'dipole' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'quad' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'quadalone' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'NRIXS' )
      write(6,*) 'NRIXS is deprecated in favor of qRaman. Consider changing your inputs'
     ! nothing needs to be done!!!
  case( 'ldos0' )
     ! nothing needs to be done!!!
  case( 'ldos1' )
     ! nothing needs to be done!!!
  case( 'ldos2' )
     ! nothing needs to be done!!!
  case( 'qRaman' )
     ! nothing needs to be done!!!

  case default
     write(6,*) "Photon type was not recognized. Must be one of the following:"
     write(6,*) "      dipole, quad, quadalone, NRIXS, qRaman"
     stop

  end select
  write ( 6, '(1a4,1f10.5)' ) 'q = ', q
  !
  ! calculate full matrix elements
  call jtvsub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, ifcn, spcttype, ehat, qhat, q, &
       nsphpt, xsph, ysph, zsph, wsph, prefs )
  !
  ! output matrix elements in desired order
  open( unit=99, file='mels', form='formatted', status='unknown' )
  rewind 99
  do mc = -lc, lc
     do l = lmin, lmax
        do m = -l, l
           do iproj = 1, nproj( l )
              write ( 99, '(2(1x,1e15.8),5i4)' ) nbsemel( iproj, m, l, mc ), nproj( l ), iproj, l, m, mc
           end do
        end do
     end do
  end do
  close( unit=99 )
  !
  return
end subroutine cmjtv
