! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine fsassem( nv, avec, bvec, bmet, zn, cfr, cfi, ng, kvc )
  implicit none
  !
  integer, parameter :: stdin = 5
  !
  integer :: nv, zn( 3 ), ng
  integer :: kvc( 3, ng )
  real( kind = kind( 1.0d0 ) ) :: avec(3,3),bvec(3,3),bmet(3,3)
  real( kind = kind( 1.0d0 ) ) :: cfr(nv,ng),cfi(nv,ng)
  !
  !
  integer :: i, j, ib, lc, mc, nbv, nbc, nq, ig
  character * 5 :: fnrootv, fnrootc
  character * 6 :: fnv, fnc
  character * 11 :: qfil
  real( kind = kind( 1.0d0 ) ) :: tauv( 3 ), tauc( 3 ), tmp( 3 ), dum, umklapp( 3 )
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
  integer, dimension( ng, 3 ) :: kvcopt
  real( kind = kind( 1.0d0 ) ), dimension( nv, nv ) :: zr, zi
  real( kind = kind( 1.0d0 ) ), dimension( ng, nv ) :: redump, imdump
  complex( kind = kind( 1.0d0 ) ), dimension( ng, nv ) :: ck
  real( kind = kind( 1.0d0 ) ), allocatable :: ev( :, : ), ec( :, : )
  real, allocatable :: vcf( :, :, :, : ), ccf( :, :, : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: conmel1( :, :, :, :, :, : ), conmel( :, :, :, : ), valmel( :, :, : , : )
  integer, parameter :: uv = 98, uc = 97, enkfile = 96, wvu = 95, wvup = 94, qsft=93
  !
! integer :: ibv, ibc
! real( kind = kind( 1.0d0 ) ) :: sur, sui
  integer :: ibc, iq
  !
  rm1 = -1; rm1 = sqrt( rm1 )
  ! 
  read ( stdin, '(1a5)' ) fnrootv
  read ( stdin, '(1a5)' ) fnrootc
  write ( fnv, '(1a5,1a1)' ) fnrootv, '1'
  write ( fnc, '(1a5,1a1)' ) fnrootc, '1'
  !
  open( unit=uv, file=fnv, form='unformatted', status='unknown' )
  rewind uv
  read ( uv ) tauv
  write ( 6, '(1a7,3f10.5)' ) 'vtau = ', tauv
  !
  open( unit=uc, file=fnc, form='unformatted', status='unknown' )
  rewind uc
  read ( uc ) tauc
  write ( 6, '(1a7,3f10.5)' ) 'ctau = ', tauc
  !
  tmp = tauv - tauc
  if ( dot_product( tmp, tmp ) .gt. 0.000001d0 ) stop 'tau mismatch'
  !
  read ( uv ) nv, nbv, nq
  read ( uc ) nv, nbc, nq
  !
  allocate( ev( nbv, nq ), ec( nbc, nq ) )
  open( unit=99, file='vraw', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) ev
  close( unit=99 )
  open( unit=99, file='craw', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) ec
  close( unit=99 )
  !
  allocate( ccf( nbc, nq, nv, 2 ), vcf( nbv, nq, nv, 2 ) )
  !
  read ( uv )
  read ( uv )
  read ( uv )
  read ( uv )
  read ( uv )
  ! 
  read ( uc )
  read ( uc )
  read ( uc )
  read ( uc )
  read ( uc )
  !
  open( unit=enkfile, file='enkfile', form='formatted', status='unknown' )
  rewind enkfile
  do i = 1, nq
     write ( enkfile, '(4f20.12,1x)' ) ev( 1 : nbv, i )
     write ( enkfile, '(4f20.12,1x)' ) ec( 1 : nbc, i )
  end do
  close( unit=enkfile )
  !
  do i = 1, nv
     read ( uv ) vcf( :, :, i, 1 )
     read ( uv ) vcf( :, :, i, 2 )
     read ( uc ) ccf( :, :, i, 1 )
     read ( uc ) ccf( :, :, i, 2 )
  end do
  !
  close( unit=uv )
  close( unit=uc )
  !
  ! chug-a-lug the wave functions out onto disk
  open( unit=wvu, file='masterwfile', form='formatted', status='unknown' )
  rewind wvu
  write ( wvu, '(1i8)' ) nq
  close( unit=wvu )
  open( unit=qsft, file='shifter', form='formatted', status='unknown' )
  rewind qsft
  open( unit=wvu, file='listwfile', form='formatted', status='unknown' )
  rewind wvu
  do i = 1, nq
     write ( qfil, '(1a5,1i6.6)' ) '.rixs', i
     write ( wvu, '(1i8,1x,1a11)' ) i, qfil
     do j = 1, 3
        read ( qsft, * ) dum, dum, umklapp( j )
     end do
     open( unit=wvup, file=qfil, form='unformatted', status='unknown' )
     rewind wvup
     write ( wvup ) ng
!    write ( 6, * ) ' ... ', qfil, ng
     do j = 1, 3
        do ig = 1, ng
           kvcopt( ig, j ) = kvc( j, ig ) + umklapp( j )
        end do
     end do
     write ( wvup ) kvcopt
     do ib = 1, nbv
        zr( :, ib ) = vcf( ib, i, :, 1 )
        zi( :, ib ) = vcf( ib, i, :, 2 )
     end do
     do ib = 1, nbc
        zr( :, ib + nbv ) = ccf( ib, i, :, 1 )
        zi( :, ib + nbv ) = ccf( ib, i, :, 2 )
     end do
     call getck( nv, ng, zr, zi, cfr, cfi, 1, nbv + nbc, ck )
     redump( :, : ) = ck( :, : )
     imdump( :, : ) = -rm1 * ck( :, : )
     write ( wvup ) redump( :, : )
     write ( wvup ) imdump( :, : )
     close( unit=wvup )
  end do
  close( unit=wvu )
  close( unit=qsft )
  !
  ! this section is for preparation of the file, tmels ...
  read ( stdin, * ) lc
  allocate( valmel( nbv, nq, 2 * lc + 1, 2 ) )
  allocate( conmel( nbc, nq, 2 * lc + 1, 2 ) )
  allocate( conmel1( nbc, nq, 2, 2 * lc + 1, 2, 2 ) )
  !
  ! read in cond. electron-core hole amplitudes
  open( unit=99, file='echamp', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) conmel1( :, :, :, :, :, : )
  close( unit=99 )
  open( unit=99, file='ascchamp', form='formatted', status='unknown' )
  rewind 99
  do ibc = 1, nbc
     do iq = 1, nq
        write ( 99, '(3(1x,1e15.8))' ) ec( ibc, iq ), conmel1( ibc, iq, 1, 1, 1, : )
     end do
  end do
  close( unit=99 )
  !
  !
  ! read in val. electron-core hole transition matrix elements
  do mc = -lc, lc
     write ( fnv, '(1a5,1i1)' ) fnrootv, 1 + mc + lc
     open( unit=99, file=fnv, form='unformatted', status='unknown' )
     rewind 99
     read ( 99 )
     read ( 99 )
     read ( 99 )
     read ( 99 )
     read ( 99 ) valmel( :, :, 1 + mc + lc, 1 ) 
     read ( 99 ) valmel( :, :, 1 + mc + lc, 2 ) 
     close( unit=99 )
  end do
  !
  ! this convertion projects out the total S=0, Sz=0 part for spin.
  conmel( :, :, :, : ) = sqrt( 0.5d0 ) * ( conmel1( :, :, 1, :, 1, : ) + conmel1( :, :, 2, :, 2, : ) ) 
  call mktmel( 'sing1', nq, nbc, nbv, lc, valmel, conmel )
  !
  ! this convertion projects out the total S=1, Sz=0 part for spin.
  conmel( :, :, :, : ) = sqrt( 0.5d0 ) * ( conmel1( :, :, 1, :, 1, : ) - conmel1( :, :, 2, :, 2, : ) ) 
  call mktmel( 'trip1', nq, nbc, nbv, lc, valmel, conmel )
  !
  ! this convertion projects out the total S=1, Sz=0 part for spin.
  conmel( :, :, :, : ) = conmel1( :, :, 1, :, 2, : )
  call mktmel( 'trip2', nq, nbc, nbv, lc, valmel, conmel )
  !
  ! this convertion projects out the total S=1, Sz=0 part for spin.
  conmel( :, :, :, : ) = conmel1( :, :, 2, :, 1, : )
  call mktmel( 'trip3', nq, nbc, nbv, lc, valmel, conmel )
  !
  !
  return
end subroutine fsassem
