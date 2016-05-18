! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine escanprep( nr, r, dr, r2, dl, vi, iu, cq, njrc, irc, lmin, lmax, sig, emin, cush, kappa, &
     aexm1, aexm2, aepot, psxm1, psxm2, pspot, aepl, pspl )
  implicit none
  !
  integer :: nr, iu, njrc( 4 ), irc, lmin, lmax
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr ), dl, vi( nr, 7 ), cq( nr )
  real( kind = kind( 1.0d0 ) ) :: aexm1( nr, lmin : lmax ), aexm2( nr, lmin : lmax ), aepot( nr, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: psxm1( nr, lmin : lmax ), psxm2( nr, lmin : lmax ), pspot( nr, lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: aepl( lmin : lmax ), pspl( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: sig( lmin : lmax ), emin( lmin : lmax ), cush, kappa( lmin : lmax )
  !
  integer :: i, ll, nv, nco, skips( 0 : 3 )
  real( kind = kind( 1.0d0 ) ) :: etot, rel, alfa, zorig, zc, rcocc( 0 : 3 ), rsocc( 0 : 3 ), ntot
  integer, allocatable :: no( : ), nl( : ), nm( : ), is ( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ev( : ), occ( : ), ek( : ), xnj( : ), phe( :, : ), pot( :, : )
  logical :: dirac
  real( kind = kind( 1.0d0 ) ), parameter :: zpseu = 0.0d0, nrl = 0.0d0
  logical, parameter :: AE = .false., PS = .true.
  !
  call mkcorcon( alfa, rel, zorig, zc, rcocc, rsocc, dirac )
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  close( unit=99 )
  open( unit=99, file='skip', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * )
  do ll = 0, 3
     read ( 99, * ) i, skips( ll )
  end do
  close( unit=99 )
  write ( 6, '(1a9,4i5)' ) 'skip   = ', skips( : )
  !
  nv = nco + ( 1 + lmax - lmin )
  allocate( no( nv ), nl( nv ), nm( nv ), is( nv ), ev( nv ), ek( nv ), occ( nv ), xnj( nv ), phe( nr, nv ), pot( nr, nv ) )
  call freshen( lmin, lmax, rcocc, skips, 1, -0.01d0 )
  call abinitio(etot,rel,alfa,nr,r,dr,r2,dl,phe,njrc,vi,zorig,ntot,nv,no,nl,nm,xnj,ev,occ,is,ek,pot,iu,cq,AE,nv)
  call optradf( nr, irc, zorig, nco, no, nl, r, phe, lmin, lmax )
  do i = nco + 1, nco + 1 + lmax - lmin
     emin( nl( i ) ) = ev( i ) - cush
     call setqmm(i,pot(:,i),nl(i),xnj(i),1,aepot(:,nl(i)),zorig,zorig,rel,nr,r,r2,dl,aexm1(:,nl(i)),aexm2(:,nl(i)),njrc,vi,AE)
     call getkap( xnj( i ), nl( i ), kappa( nl( i ) ) )
     call getplead( nl( i ), xnj( i ), rel, kappa( nl( i ) ), aepl( nl( i ) ), zorig )
     sig( nl( i ) ) = ( -1 ) ** skips( nl( i ) )
  end do
  !
  call freshen( lmin, lmax, rcocc, skips, 1, -0.01d0 )
  call abinitio(etot,nrl,alfa,nr,r,dr,r2,dl,phe,njrc,vi,zorig,ntot,nv,no,nl,nm,xnj,ev,occ,is,ek,pot,iu,cq,PS,nv)
  do i = 1, 1 + lmax - lmin
     call setqmm(i,pot(:,i),nl(i),xnj(i),1,pspot(:,nl(i)),zpseu,zorig,nrl,nr,r,r2,dl,psxm1(:,nl(i)),psxm2(:,nl(i)),njrc,vi,PS)
     call getkap( xnj( i ), nl( i ), kappa( nl( i ) ) )
     call getplead( nl( i ), xnj( i ), nrl, kappa( nl( i ) ), pspl( nl( i ) ), zpseu )
  end do
  !
  return
end subroutine escanprep
