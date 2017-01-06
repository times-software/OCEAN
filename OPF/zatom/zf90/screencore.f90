! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine screencore( etot, rel, alfa, nr, r, dr, r2, dl, njrc, vi, zorig, xntot, nel, iuflag, cq, isuse, rsocc, lmin, lmax ) 
  implicit none
  !
  integer :: iuflag, nr, njrc( 4 ), nel, isuse, lmin, lmax, skips( 0 : 3 ), noskip( 0 : 3 )
  real( kind = kind( 1.0d0 ) ) :: etot, rel, alfa, dl, zorig, xntot, rsocc( 0 : 3 )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr ), vi( nr, 7 ), cq( nr )
  !
  integer :: nco, nuse, luse, iuse, i, j, l, nelmax, nclist, nvlist
  real( kind = kind( 1.0d0 ) ) :: ru, au, nulocc( 0 : 3 )
  real( kind = kind( 1.0d0 ) ) :: vcg( nr ), vcx( nr ), vvg( nr ), vvx( nr ), visav( nr, 7 ), gtot( nr ), xtot( nr )
  real( kind = kind( 1.0d0 ) ), allocatable :: ev( : ), occ( : ), ek( : ), xnj( : ), orb( :, : ), phe( :, : )
  integer, allocatable :: no( : ), nl( : ), nm( : ), is( : ), vlist( : ), clist( : ), vplist( : )
  character * 7 :: cp7, vp7
  character * 10 :: adduse, vs10
  logical, parameter :: AE = .false., PS=.true.
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
  allocate( no( nelmax ), nl( nelmax ), nm( nelmax ), is( nelmax ), vlist( nelmax ), clist( nelmax ), vplist( nelmax ) )
  allocate( ev( nelmax ), ek( nelmax ), occ( nelmax ), xnj( nelmax ), phe( nr, nelmax ), orb( nr, nelmax ) )
  !
  ru = rel
  au = alfa
  visav( :, : ) = vi( :, : )
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
     ! output core hole potl files
     nclist = nco
     do i = 1, nclist
        clist( i ) = i
     end do
     nvlist = 1 + lmax - lmin
     do i = 1, nvlist
        vlist( i ) = nco + i
        vplist( i ) = i
     end do
     ! 
     ! output files that have no valence electrons
     cp7 = 'vc_bare'
     write ( 6, * ); write ( 6, '(1a5,2(1x,1a7))' ) 'begin', cp7;
     call freshen( lmin, lmax, nulocc, skips, 1, +0.01d0 )
     write ( 6, '(1a14)' ) 'Excited state:'
     call chgocc( iuse, -1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,AE,nelmax)
     call potfigure( nclist, clist, nr, r, dl, phe, nelmax, occ, vcx )
     write ( 6, '(1a14)' ) 'Ground state: '
     call chgocc( iuse, +1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,AE,nelmax)
     call potfigure( nclist, clist, nr, r, dl, phe, nelmax, occ, vcg )
     call optript( nr, r, vcg, vcx, adduse, cp7 )
     write ( 6, '(1a5,2(1x,1a7))' ) 'end  ', cp7; write ( 6, * ); 
     !
     ! output files that have valence electrons
     cp7 = 'vcallel'; vp7 = 'vvallel'
     write ( 6, * ); write ( 6, '(1a5,2(1x,1a7))' ) 'begin', cp7, vp7;
     call freshen( lmin, lmax, rsocc, skips, 1, +0.01d0 )
     write ( 6, '(1a14)' ) 'Excited state:'
     call chgocc( iuse, -1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,AE,nelmax)
     call potfigure( nclist, clist, nr, r, dl, phe, nelmax, occ, vcx )
     call potfigure( nvlist, vlist, nr, r, dl, phe, nelmax, occ, vvx )
     write ( 6, '(1a14)' ) 'Ground state: '
     call chgocc( iuse, +1.0d0 )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,AE,nelmax)
     call potfigure( nclist, clist, nr, r, dl, phe, nelmax, occ, vcg )
     call potfigure( nvlist, vlist, nr, r, dl, phe, nelmax, occ, vvg )
     call optript( nr, r, vcg, vcx, adduse, cp7 )
     call optript( nr, r, vvg, vvx, adduse, vp7 )
     write ( 6, '(1a5,2(1x,1a7))' ) 'end  ', cp7, vp7; write ( 6, * );
     ! write sum of all electron results; this is helpful for core screening in atom prog from 
     ! valence but frozen out bands in the solid, like Ti(3s,3p), Mg(2s,2p), etc.
     cp7 = 'aetotal'
     write ( 6, * ); write ( 6, '(1a5,2(1x,1a7))' ) 'begin', cp7;
     gtot = vcg + vvg
     xtot = vcx + vvx
     call optript( nr, r, gtot, xtot, adduse, cp7 )
     write ( 6, '(1a5,2(1x,1a7))' ) 'end  ', cp7; write ( 6, * );
     !
     ! output files for changes in val occs
     cp7 = 'valence'
     do l = 0, 3
        if ( rsocc( l ) .gt. 1.0d0 ) then
           write ( vs10, '(1a1,1i3.3,1a5,1i1.1)' ) 'z', nint( zorig ), 'semil', l
           write ( 6, * ); write ( 6, '(1a5,1x,1a7,1x,1a10)' ) 'begin', cp7, vs10;
           write ( 6, '(1a14)' ) 'Excited state:'
           rsocc( l ) = rsocc( l ) - 1.0d0
           call freshen( lmin, lmax, rsocc, skips, 1, +0.010d0 )
           call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,AE,nelmax)
           call potfigure( nclist, clist, nr, r, dl, phe, nelmax, occ, vcx )
           call potfigure( nvlist, vlist, nr, r, dl, phe, nelmax, occ, vvx )
           write ( 6, '(1a14)' ) 'Ground state: '
           rsocc( l ) = rsocc( l ) + 1.0d0
           call freshen( lmin, lmax, rsocc, skips, 1, +0.010d0 )
           call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,AE,nelmax)
           call potfigure( nclist, clist, nr, r, dl, phe, nelmax, occ, vcg )
           call potfigure( nvlist, vlist, nr, r, dl, phe, nelmax, occ, vvg )
           gtot = vcg + vvg
           xtot = vcx + vvx
           call optript( nr, r, gtot, xtot, vs10, cp7 )
           write ( 6, * ); write ( 6, '(1a5,1x,1a7,1x,1a10)' ) 'end  ', cp7, vs10;
        end if
     end do
     !
     ! output file with pseudo val response to bare core screening
     cp7 = 'vc_bare'; vp7 = 'vpseud0'
     write ( 6, * ); write ( 6, '(1a5,2(1x,1a7))' ) 'begin', cp7, vp7;
     call freshen( lmin, lmax, rsocc, skips, 1, +0.01d0 )
     write ( 6, '(1a14)' ) 'Excited state:'
     ru = 0; au = alfa
     write ( 6, * ) 'xx 1'
     call addpot( zorig, nuse, luse, cp7, nr, vi )
     write ( 6, * ) 'xx 2'
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,PS,nelmax)
     write ( 6, * ) 'xx 3'
     call potfigure( nvlist, vplist, nr, r, dl, phe, nelmax, occ, vvx )
     write ( 6, * ) 'xx 4'
     write ( 6, '(1a14)' ) 'Ground state: '
     ru = 0; au = alfa
     vi( :, : ) = visav( :, : )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,PS,nelmax)
     call potfigure( nvlist, vplist, nr, r, dl, phe, nelmax, occ, vvg )
     call optript( nr, r, vvg, vvx, adduse, vp7 )
     write ( 6, '(1a5,2(1x,1a7))' ) 'end  ', cp7, vp7; write ( 6, * );
     !
     ! output file with pseudo val response to dressed core screening
     cp7 = 'vcallel'; vp7 = 'vpseud1'
     write ( 6, * ); write ( 6, '(1a5,2(1x,1a7))' ) 'begin', cp7, vp7;
     call freshen( lmin, lmax, rsocc, skips, 1, +0.01d0 )
     write ( 6, '(1a14)' ) 'Excited state:'
     ru = 0; au = alfa
     call addpot( zorig, nuse, luse, cp7, nr, vi )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,PS,nelmax)
     call potfigure( nvlist, vplist, nr, r, dl, phe, nelmax, occ, vvx )
     write ( 6, '(1a14)' ) 'Ground state: '
     ru = 0; au = alfa
     vi( :, : ) = visav( :, : )
     call abinitio(etot,ru,au,nr,r,dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,orb,iuflag,cq,PS,nelmax)
     call potfigure( nvlist, vplist, nr, r, dl, phe, nelmax, occ, vvg )
     call optript( nr, r, vvg, vvx, adduse, vp7 )
     write ( 6, '(1a5,2(1x,1a7))' ) 'end  ', cp7, vp7; write ( 6, * );
     !
  end do
  !
  return
end subroutine screencore
