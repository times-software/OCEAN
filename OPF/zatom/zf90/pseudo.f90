! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine pseudo( etot, rel, alfa, nr, rmin, rmax, r, dr, r2, dl, phe, orb, njrc, vi, cq, zorig, xntot, nel, &
     no, nl, nm, xnj, ev, occ, is, ek, iuflag, vctab, vtry, isuse )
  implicit none
  !
  integer :: nr, nel, iuflag, vtry, isuse, njrc( 4 )
  integer, dimension( nel ) :: no, nl, nm, is
  real( kind = kind( 1.0d0 ) ) :: etot, rel, alfa, rmin, rmax, dl
  real( kind = kind( 1.0d0 ) ) :: zorig, xntot
  real( kind = kind( 1.0d0 ) ) :: vi( nr, 7 ), vctab( nr, 0 : 3 )
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: r, dr, r2, cq
  real( kind = kind( 1.0d0 ) ), dimension( nr, nel ) :: phe, orb
  real( kind = kind( 1.0d0 ) ), dimension( nel ) :: xnj, ek, ev, occ
  !
  integer :: i, i1, np, lmax, lfc, lu, ido, mjrc( 4 ), csh( 0 : 3 ), ntar
  double precision :: pi, corpol, rnorm, ratio, zuse, zeff, rcut, rtest, factor
  double precision, allocatable :: xm1( : ), xm2( : )
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  do i = 1, 4
     if ( njrc( i ) .gt. 0 ) stop 'cannot repseudize as of now'
     mjrc( i ) = 0
  end do
  !
  read ( 5, * ) np, corpol, rnorm
  read ( 5, * ) lfc,ratio
  zuse = zorig - sum( occ( 1 : np - 1 ) )
  write ( 6, '(1a14,1f10.5)' ) 'core charge = ', zuse
  call mklfc( lfc, ratio, nr, r, dl, cq, nel, np, occ, phe )
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) np - 1
  csh( : ) = 0
  do i = 1, np - 1
     csh( nl( i ) ) = csh( nl( i ) ) + nint( occ( i ) )
     write ( 99, * ) no( i ), nl( i ), nm( i ), xnj( i ), is( i ), occ( i )
  end do
  close( unit=99 )
  do i = 0, 3
     csh( i ) = csh( i ) / ( 2 * ( 2 * i + 1 ) )
  end do
  write ( 6, '(1a13,4i3)' ) 'core skips = ', csh( : )
  !
  open( unit=99, file='valcon', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(2x,1i5)' ) 1 + nel - np
  lmax = 0
  do i = np, nel
     if ( nl( i ) .gt. lmax ) lmax = nl( i )
     write ( 99, '(2x,2i5,1f10.6)' ) no( i ), nl( i ), occ( i )
  end do
  close( unit=99 )
  xntot = sum( occ( np : nel ) )
  !
  open( unit=99, file='skip', form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) lmax
  do i = 0, lmax
     write ( 99, '(2i5)' ) i, csh( i )
  end do
  close( unit=99)
  !
  allocate( xm1( nr ), xm2( nr ) )
  !
  ! conduct pseudization
  do i = np, nel
     read ( 5, * ) ido, rcut, rtest, factor
     if ( ido .gt. 0 ) then
        write ( 6, '(1a11,5x,2i3,1f6.1,3f12.8)' ) 'pseudizing ', no( i ), nl( i ), xnj( i ), rcut, rtest, factor
        lu = 2 * nl( i ) + 1
        if ( abs( xnj( i ) ) + 0.25d0 .lt. dble( nl( i ) ) ) lu = 2 * nl( i )
        orb( :,i ) = orb( :, i ) + vctab( :, nl( i ) )
        call setqmm( i, orb( 1, i ), nl( i ), xnj( i ), 1, vi( 1, lu ), zeff, zorig, rel, nr, r, r2, dl, xm1, xm2, &
             mjrc, vi, .false. )
        orb( :, i ) = 0
        ntar = no( i ) - csh( nl( i ) ) - nl( i ) - 1
        write ( 6, '(10x,1a7,1i5)' ) 'ntar = ', ntar
        call pseudize( i, orb( 1, i ), ev( i ), nl( i ), xnj( i ), no( i ), njrc, zeff, vi( 1, lu ), xm1, xm2, &
             nr, rmin, rmax, r, dr, r2, dl, rel, rcut, rtest, factor, ntar )
     end if
  end do
  !
  ! step down shell structure
  do i = np, nel
     lu = 2 * nl( i ) + 1
     write ( 6, '(3i5,5x,1a12,2i5)' ) i, np, nel, 'old n and l ', no( i ), nl( i )
     no( i ) = no( i ) - csh( nl( i ) )
     write ( 6, '(3i5,5x,1a12,2i5)' ) i, np, nel, 'new n and l ', no( i ), nl( i )
     call elsolve( i, no( i ), nl( i ), -1.0d0, xnj( i ), zorig, zeff, ev( i ), phe( 1, i ), vi( 1, lu ), &
          xm1, xm2, nr, r, dr, r2, dl, 0.0d0, 1, 0 )
     write ( 6, '(1x,1a13,1i5,3x,1f4.1,3x,1f14.6)' ) 'solution ... ', nl( i ), xnj( i ), ev( i )
     i1 = i - ( np - 1 )
     phe( :, i1 ) = phe( :, i )
     no( i1 ) = no( i )
     nl( i1 ) = nl( i )
     nm( i1 ) = nm( i )
     xnj( i1 ) = xnj( i )
     is( i1 ) = is( i )
     ev( i1 ) = ev( i )
     occ( i1 ) = occ( i )
  end do
  nel = 1 + nel - np
  deallocate( xm1, xm2 )
  !
  call unscreen( nr, r, dr, r2, dl, etot, alfa, xntot, nel, phe, orb, occ, is, nl, nm, no, xnj, iuflag, cq, ev, vi, zuse, corpol, &
       rnorm )
  !
  return
end subroutine pseudo
