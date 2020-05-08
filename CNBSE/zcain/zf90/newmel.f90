! Copyright (C) 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program newmel
  implicit none
  !
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: powmax = 2
  integer, parameter :: ylmax = 5
  !
  character(len=4) :: add04
  character(len=15) :: spcttype
  character(len=17) :: f17
  integer :: atno, nc, lc, lmin, lmax, nqproj, npmax, ierr, nsphpt, isphpt, ll, mc, iprj, ll2, mc2, iprj2, ip
  real(DP) :: dqproj, qhat( 3 ), ehat( 3 ), ephotev, lam, pi, q, dummy, sphsu, prefs( 0 : 1000 )
  integer, allocatable :: nproj( : )
  real(DP), allocatable :: radialPrj( :, :, :, :, : ), ifcn( :, :, : ), xsph(:), ysph(:), zsph(:), wsph(:)
  complex(DP), allocatable :: tempMels( :, :, :, : ), mels( :,:,:,:,:,:)

  open( unit=99, file='sphpts', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nsphpt
  allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
  do isphpt = 1, nsphpt
     read ( 99, * ) xsph( isphpt ), ysph( isphpt ), zsph( isphpt ), wsph( isphpt )
  end do
  close( unit=99 )
  sphsu = sum( wsph( : ) )
  wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / sphsu )
  write ( 6, * ) nsphpt, ' points with weights summing to four pi '

  call getprefs( prefs, ylmax, nsphpt, wsph, xsph, ysph, zsph )

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
  ! Need the radial overlaps of everybody
  allocate( radialPrj( npmax, lmin:lmax, npmax, lmin:lmax, 0: powmax ) )

  
  call makeProjectorRadialIntegral( npmax, lmin, lmax, nproj, atno, powmax, radialPrj, ierr )

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
  case( 'dipole' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'quad' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )
  case( 'quadalone' )
     lam = 10000.0d0 / ( ephotev * 0.806554445d0 )
     q = 0.529177d0 * ( 2.0d0 * pi / lam )

  case default
     write(6,*) "Only dipole, quad, and quadalone are implemented in newmel"

  end select
  write ( 6, '(1a4,1f10.5)' ) 'q = ', q
  
  allocate( ifcn( npmax, 0:powmax, lmin:lmax ) )

  write( f17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', atno, 'n', nc, 'l', lc
  open( unit=99, file=f17, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ll2
  do ll = lmin, lmax
     do ip = 0, powmax
        read ( 99, * ) ifcn( 1 : nproj( ll ), ip, ll )
     end do
  end do
  close( unit=99 )


  allocate( mels( npmax, -lmax : lmax, lmin:lmax, npmax, -lmax : lmax, lmin:lmax ) )
  mels = 0.0_dp
  do ll = lmin, lmax
    allocate( tempMels( npmax, -lmax : lmax, lmin:lmax, -ll : ll ) )
    do iprj = 1, nproj( ll )
      tempMels = 0.0_dp
      do ll2 = lmin, lmax
        do ip = 0, powmax
          ifcn( :, ip, ll2 ) = radialPrj( :, ll2, iprj, ll, ip )
        enddo
      enddo
      call jtvsub( lmin, lmax, nproj, npmax, ll, tempMels, powmax, ifcn, spcttype, ehat, qhat, q, &
                   nsphpt, xsph, ysph, zsph, wsph, prefs )
      do mc = -ll, ll
        mels( :, :, :, iprj, mc, ll ) = tempMels( :, :, :, mc )
      enddo
    enddo
    deallocate( tempMels )
  enddo


  open( unit=99, file='newmels', form='formatted', status='unknown' )
  rewind 99

  do ll = lmin, lmax
    do mc = -ll, ll
      do iprj = 1, nproj(ll)

        do ll2 = lmin, lmax
          do mc2 = -ll2, ll2
            do iprj2 = 1, nproj(ll2) 

              write(99, '(2(1x,1e26.18),6i4)' ) mels( iprj2, mc2, ll2, iprj, mc, ll ), iprj2, mc2, ll2, iprj, mc, ll
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo

  



  deallocate( radialPrj, ifcn, mels )
end program
