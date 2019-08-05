! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program vhommod
  implicit none
  integer, parameter :: stdin = 5, stdout = 6
  integer, parameter :: dp = kind( 1.0d0 )
  !
  integer :: i, j, k, nr, nq, idum, nshells, nsites, isite, ierr, xmesh( 3 )
  integer, allocatable :: siteIndex( : )
  real( kind=dp ) :: epsq0, qf, wpsqd, wf, lam, pi, eps
  real( kind=dp ) :: qq, q, dq, qmax, r1, r2, dr, rmax
  real( kind=dp ) :: jqr1, jqr2, tj1, tj2, th1, th2, j1, s
  real( kind=dp ) :: nav, ds, ss, dn, n, avec( 3, 3 ), omega, nel
  !
  real( kind=dp ), allocatable :: v( :, : ), vh( :, : ), dv( : )
  real( kind=dp ), allocatable :: tab1( :, : ), bq( : ), shells( : )
  !
  character(len=2), allocatable :: siteSymbol( : )
  character(len=4) :: smode
  character(len=80) :: avgName, reoptName
  !
  real( kind=dp ), external :: levlou, sphj0, sphj1
  !
  logical :: havenav, valenceGrid
  character(len=80) :: dummy
  !
  pi = 4.d0 * datan( 1.d0 )
  !
!  read ( stdin, * ) r2 !, dr, nr
!  read ( stdin, * ) dq, qmax   ! dq and qmax in units of qfermi
  !
  open( unit=99, file='shells', form='formatted', status='old' )
  rewind 99
  read( 99, * ) nshells
  allocate( shells( nshells ) )
  do i = 1, nshells
    read( 99, * ) shells( i )
  enddo
  close( 99 )

  inquire(file='screen.mode', exist=valenceGrid )
  if( valenceGrid ) then
    open( unit=99, file='screen.mode', form='formatted', status='old' )
    read( 99, * ) smode
    close( 99 )
    if( smode .eq. 'grid' ) then
      valenceGrid = .true.
    else
      valenceGrid = .false.
    endif
  endif

  if( valenceGrid ) then
    open( unit=99, file='xmesh.ipt', form='formatted', status='old' )
    read( 99, * ) xmesh(:)
    close( 99 )
    nsites = product( xmesh( : ) )
    allocate( siteSymbol( 0 ), siteIndex( 0 ) )
  else  !core-level
    open( unit=99, file='sitelist', form='formatted', status='old' )
    read( 99, * ) nsites
    allocate( siteSymbol( nsites ), siteIndex( nsites ) )
    do i = 1, nsites
      read( 99, * ) siteSymbol( i ), idum, siteIndex( i )
    enddo
    close( 99 )
  endif

  open( unit=99, file='screen.final.rmax', form='formatted', status='old' )
  rewind 99
  read( 99, * ) rmax
  close( 99 )
  open( unit=99, file='screen.final.dr', form='formatted', status='old' )
  read( 99, * ) dr
  close( 99 )
  nr = rmax / dr

  open( unit=99, file='screen.model.dq', form='formatted', status='old' )
  rewind 99
  read( 99, * ) dq
  close( 99 )
  open( unit=99, file='screen.model.qmax', form='formatted', status='old' )
  rewind 99
  read( 99, * ) qmax
  close( 99 )
  open( unit=99, file='epsilon', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) epsq0
  close( 99 )
  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) avec( :, : )
  close( unit=99 )
  call getomega( avec, omega )
  open ( unit=99, file='rhoofg', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) idum
  read ( 99, * ) idum, idum, idum, nel
  close( unit=99 )
  inquire( file='fake_nav.ipt', exist=havenav)
  if( havenav ) then
    open( unit=99, file='fake_nav.ipt', form='formatted', status='old' )
    rewind 99
    read( 99, * ) nav
  else
    nav = dble( nel ) / omega
  endif
  !
!  read ( stdin, * ) epsq0 
  allocate( v( nr, nshells ), vh( nr, nshells ), dv( nr ) )
  !
  qf = ( 3 * pi ** 2 * nav ) ** ( 1.d0 / 3.d0 )
  wpsqd = 4 * pi * nav
  wf = qf ** 2 / 2
  lam = dsqrt( wpsqd / ( wf ** 2 * ( epsq0 - 1 ) ) )
  !
  nq = qmax / dq
  dq = qf * dq
  allocate( tab1( nr, nq ), bq( nq ) )
  !
  v( :, : ) = 0
  do k = 1, nshells
    r2 = shells( k )
    do i = 1, nq
       q = dq * ( i - 0.5d0 )
       jqr2 = sphj0( q * r2 )
       qq = q / qf
       eps = 1 / levlou( qq, qf, lam )
       bq( i ) = ( 1 - 1 / eps ) / ( 4 * pi )
       r1 = dr / 2
       do j = 1, nr
          jqr1 = sphj0( q * r1 )
          v( j, k ) = v( j, k ) + 8 * dq * bq( i ) * jqr1 * jqr2
          r1 = r1 + dr
       end do
    end do
  end do
  vh( :, : ) = v( :, : )
  !
  do i = 1, nq 
     q = dq * ( i - 0.5d0 )
     r1 = dr / 2
     do j = 1, nr
        tab1( j, i ) = sphj0( q * r1 )
        r1 = r1 + dr
     end do
  end do

  do isite = 1, nsites
    s = 0.00001d0
    ds = 0.10d0

    if( isite .gt. 1 ) then
      v( :, : ) = vh( :, : )
    endif

    if( valenceGrid ) then
      write( avgName, '(A3,I6.6)' ) 'avg', isite
    else
      write( avgName, '(A3,A2,I4.4)' ) 'avg', siteSymbol( isite ), siteIndex( isite )
    endif
    write( 6, * ) trim( avgName )

    open( unit=99, file=trim(avgName), form='formatted', status='old' )
    rewind 99

    do 
      read( 99, *, iostat=ierr ) dummy, dummy, ss, n
      if( ierr .lt. 0 ) then
        goto 100
      endif
      dn = n - nav
      if ( abs( s - ss ) .gt. 0.000001d0 ) then
        write(6,*) s, ss
        stop 'jive!'
      endif

      do k = 1, nshells
        r2 = shells( k )
        th2 = 0
        if ( s .gt. r2 ) th2 = 1
        dv = 0
        do i = 1, nq
          q = dq * ( i - 0.5d0 )
          j1 = sphj1( q * s )
          tj2 = j1 * th2
          jqr2 = sphj0( q * r2 )
          r1 = dr / 2
          do j = 1, nr
            th1 = 0
            if ( s .gt. r1 ) th1 = 1
            tj1 = j1 * th1
            jqr1 = tab1( j, i )
            dv( j ) = dv( j ) + q * bq( i ) * ( tj1 * jqr2 + tj2 * jqr1 )
            r1 = r1 + dr
          end do
        end do
        v( :, k ) = v( :, k ) + dv( : ) * ds * dn * 4 * dq / nav
      enddo

      s = s + ds
    end do
100 continue
    close( 99 )

    do k = 1, nshells
      if( valenceGrid ) then
        write( reoptName, '(A5,I6.6,A2,F4.2)' ) 'reopt', isite, &
                                                   '.R', shells( k )
      else
        write( reoptName, '(A5,A2,I4.4,A2,F4.2)' ) 'reopt', siteSymbol( isite ), siteIndex( isite ), &
                                                   '.R', shells( k )
      endif

      open( unit=99, file=reoptName, form='formatted', status='unknown' )
      rewind( 99 )
      r1 = dr / 2
      do j = 1, nr
        write ( 99, '(3(1x,1e15.8))' ) r1, vh( j, k ), v( j, k )
        r1 = r1 + dr
      end do
      close( unit=99 )
    enddo

  enddo

#if 0
!  open( unit=99, file='avden', form='formatted', status='unknown' )
!  rewind 99
  do while ( s .lt. 40 )
     read ( 99, * ) dummy, dummy, ss, n
     dn = n - nav
     if ( abs( s - ss ) .gt. 0.000001d0 ) stop 'jive!'
     th2 = 0
     if ( s .gt. r2 ) th2 = 1
     dv = 0
     do i = 1, nq
        q = dq * ( i - 0.5d0 )
        j1 = sphj1( q * s )
        tj2 = j1 * th2
        jqr2 = sphj0( q * r2 )
        r1 = dr / 2
        do j = 1, nr
           th1 = 0
           if ( s .gt. r1 ) th1 = 1
           tj1 = j1 * th1
           jqr1 = tab1( j, i )
           dv( j ) = dv( j ) + q * bq( i ) * ( tj1 * jqr2 + tj2 * jqr1 )
           r1 = r1 + dr
        end do
     end do
     v( : ) = v( : ) + dv( : ) * ds * dn * 4 * dq / nav
     s = s + ds
  end do
  close( unit=99 )
  !
  open( unit=99, file='reopt', form='formatted', status='unknown' )
  rewind 99
  r1 = dr / 2
  do j = 1, nr
     write ( 99, '(3(1x,1e15.8))' ) r1, vh( j ), v( j )
     r1 = r1 + dr
  end do
  close( unit=99 )
#endif
  !
  deallocate( v, vh, tab1, shells, siteSymbol, siteIndex )
  !
end program vhommod
