! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine pseff( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rsocc, lmin, lmax ) 
  implicit none
  !
  double precision :: etot, rel, alfa
  double precision :: dl, zorig, xntot, rsocc( 0 : 3 )
  integer :: nr, njrc( 4 ), nel, isuse, lmin, lmax, skips( 0 : 3 ), noskip( 0 : 3 )
  integer :: iuflag
  double precision :: r( nr ),dr( nr ),r2( nr )
  double precision :: vi( nr, 7 ), cq( nr )
  integer, allocatable :: no( : ), nl( : ), nm( : )
  integer, allocatable :: is( : )
  double precision, allocatable :: ev( : ), occ( : ), ek( : )
  double precision, allocatable :: xnj( : )
  integer :: nco, nuse, luse, iuse, i, j, l, nelmax
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: orb, phe
  real( kind = kind( 1.0d0 ) ) :: ru, au, nulocc( 0 : 3 ), hptab( nr, 2 )
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: vcg, vcx, vvg, vvx, vorb
  real( kind = kind( 1.0d0 ) ), dimension( nr, 7 ) :: visav
  character * 10 :: adduse
  !
  open( unit=99, file='skip', form='formatted', status='unknown' )
  read ( 99, * )
  do l = 0, 3
     read ( 99, * ) i, skips( l )
  end do
  noskip( : ) = 0.0d0
  write ( 6, * ) 'skip = ', skips( : )
  write ( 6, * ) 'noskip = ', noskip( : )
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  close( unit=99 )
  write ( 6, * ) 'nco = ', nco 
  !
  nelmax = 1 + nco + ( 1 + lmax - lmin ) + 10
  allocate( no( nelmax ), nl( nelmax ) )
  allocate( nm( nelmax ), is( nelmax ) )
  allocate( ev( nelmax ), ek( nelmax ) )
  allocate( occ( nelmax ), xnj( nelmax ) )
  allocate( phe( nr, nelmax ), orb( nr, nelmax ) )
  !
  ru = 1.0d0
  au = 0.0d0
  !
  nulocc( : ) = 0.0d0
  do iuse = 1, nco
     !
     ! acquiesce N and L
     open( unit=99, file='corcon', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * )
     do j = 1, iuse
        read ( 99, * ) nuse, luse
     end do
     close( unit=99 )
     !
     ! fix file addend
     write ( adduse, '(1a1,1i3.3,2(1a1,1i2.2))' ) 'z', nint( zorig ), 'n', nuse, 'l', luse
     !
     ! output real input files
     call freshen( lmin, lmax, nulocc, skips )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     vcg( : ) = 0.0d0
     do i = 1, nco
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vcg( : ) = vcg( : ) + occ( i ) * vorb( : )
     end do
     call chgocc( iuse, -1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     vcx( : ) = 0.0d0
     do i = 1, nco
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vcx( : ) = vcx( : ) + occ( i ) * vorb( : )
     end do
     call chgocc( iuse, +1.0d0 )
     !
     hptab( :, 1 ) = vcg( : )
     hptab( :, 2 ) = vcx( : )
     call optript( nr, r, hptab, adduse, 'vc_bare' )
     !
     ! output files that have valence electrons
     call freshen( lmin, lmax, rsocc, skips )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     vcg( : ) = 0.0d0
     do i = 1, nco
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vcg( : ) = vcg( : ) + occ( i ) * vorb( : )
     end do
     vvg( : ) = 0.0d0
     do l = lmin, lmax
        i = nco + 1 + ( l - lmin )   
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vvg( : ) = vvg( : ) + occ( i ) * vorb( : )
     end do
     call chgocc( iuse, -1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.false.,nelmax)
     vcx( : ) = 0.0d0
     do i = 1, nco
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vcx( : ) = vcx( : ) + occ( i ) * vorb( : )
     end do
     vvx( : ) = 0.0d0  
     do l = lmin, lmax
        i = nco + 1 + ( l - lmin )   
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vvx( : ) = vvx( : ) + occ( i ) * vorb( : )
     end do
     call chgocc( iuse, +1.0d0 )
     ! 
     hptab( :, 1 ) = vcg( : )
     hptab( :, 2 ) = vcx( : )
     call optript( nr, r, hptab, adduse, 'vcallel' )
     ! 
     hptab( :, 1 ) = vvg( : )
     hptab( :, 2 ) = vvx( : )
     call optript( nr, r, hptab, adduse, 'vvallel' )
     !
     call freshen( lmin, lmax, rsocc, skips )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.true.,nelmax)
     vvg( : ) = 0.0d0
     do l = lmin, lmax
        i = 1 + ( l - lmin )
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vvg( : ) = vvg( : ) + occ( i ) * vorb( : )
     end do
     !
     visav( :, : ) = vi( :, : )
     call addpot( zorig, nuse, luse, 'vc_bare', nr, vi )
     !
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,.true.,nelmax)
     vvx( : ) = 0.0d0
     do l = lmin, lmax
        i = 1 + ( l - lmin )
        call orbcont( nr, r, dl, phe( 1, i ), vorb )
        vvx( : ) = vvx( : ) + occ( i ) * vorb( : )
     end do
     !
     hptab( :, 1 ) = vvg( : )
     hptab( :, 2 ) = vvx( : )
     call optript( nr, r, hptab, adduse, 'vvpseud' )
     !
     vi( :, : ) = visav( :, : )
     !
  end do
  !
  deallocate( phe, orb, no, nl, nm, is, ev, ek, occ, xnj )
  !
  return
end subroutine pseff
