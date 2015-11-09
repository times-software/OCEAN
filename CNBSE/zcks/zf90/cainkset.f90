! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cainkset( avec, bvec, bmet, prefs )
  implicit none
  integer, parameter :: stdin = 5, stdout = 6, mubase = 80
  character * 7, parameter :: u7 = 'unknown'
  character * 9, parameter :: f9 = 'formatted'
  character * 11, parameter :: u11 = 'unformatted'
  !
  real( kind=kind( 1.0d0 ) ) :: bmet( 3, 3 ), bvec( 3, 3 ), avec( 3, 3 ), prefs( 0 : 1000 )
  !
  integer :: zn( 3 )
  integer :: ntau, itau
  real( kind = kind( 1.0d0 ) ), allocatable :: tau( :, : )
  !
  integer :: i, ii, ik1, ik2, ik3, nq, nktot, ntot, ivl, ivh, icl, ich, nbtot, ng
  integer :: nbd, ibd, atno, nc, lc, ibeg, nspin, ispin
  real( kind = kind( 1.0d0 ) ) :: qbase( 3 )
  real( kind = kind( 1.0d0 ) ) :: dbeta( 3 ), qraw( 3 ), dq( 3 )
  real( kind = kind( 1.0d0 ) ) :: eshift, edge, sc, pi, celvol
  !
  integer, allocatable :: iq( :, : ), g( :, : ), kvc( :, : ), flip_g( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: qphys( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: e0( : ), eraw( : )
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: zzr( :, : ), zzi( :, : ), ww( :, :, : ), w( : )
  !
  integer :: iproj, indx
  logical :: metal, conduct
  real( kind = kind( 1.0d0 ) ) :: efermi, temperature
  !
  character * 2 :: element
  character * 4 :: add04
  character * 5 :: fnroot
  character * 9, allocatable, dimension( : ) :: fntau
  character * 10 :: add10, infoname
  character * 12, allocatable, dimension( :, : ) :: wnam
  !
  real( kind = kind( 1.0d0 ) ) :: dqproj
  integer, allocatable :: nproj( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: fttab( :, :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: coeff( :, :, :, :, : ), cwgt( :, :, :, : ), ck( :, : )
  !
  integer :: nprj, ip
  complex( kind = kind( 1.0d0 ) ), allocatable :: prj( : , : , : )
  !
  integer :: l, m, lmin, lmax, npmax, nqproj, nptot, nlm, ilm
  complex( kind = kind( 1.0d0 ) ) :: rm1
  integer, allocatable :: lml( : ), lmm( : ), ibeg_array( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: pcoefr( :, :, : ), pcoefi( :, :, : )  
  !
  integer :: j
  real( kind = kind( 1.0d0 ) ) :: su, qsqd, betot( 3 ), nrm
  integer :: OMP_GET_THREAD_NUM, fh
  logical :: temp_exist, is_jdftx
  character(len=3) :: DFT
  !
  write ( stdout, * ) 'warning: this assumes one-component system'
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  pi = 4.d0 * datan( 1.d0 )
  call getcel( celvol, avec )
  !
  !
! read ( stdin, * ) qbase( : )
  open( unit=99, file='scaledkzero.ipt', form=f9, status='old' )
  rewind 99
  read ( 99, * ) qbase( : )
  close( unit=99 )
  ! 
! read ( stdin, * ) nbd
  open( unit=99, file='nbuse.ipt', form=f9, status='old' )
  rewind 99
  read ( 99, * ) nbd
  close( unit=99 )
  ! 
! read ( stdin, * ) cs
  open( unit=99, file='eshift.ipt', form=f9, status='old' )
  rewind 99
  read ( 99, * ) eshift
  close( unit=99 )
  ! 
!  read ( stdin, * ) metal ! new
  open( unit=99, file='metal', form=f9, status='old')
  rewind 99
  read(99, * ) metal
  close( 99 )
  !
  open(unit=99, file='cks.normal', form=f9, status='old')
  rewind 99
  read(99, * ) conduct
  close(99)
!  conduct = .true.
  !
  open( unit=99, file='nspin', form=f9, status='old' )
  rewind 99
  read( 99, * ) nspin
  close( unit=99 )
  !
  open( unit=99, file='kmesh.ipt', form=f9, status='old' )
  rewind 99
  read ( 99, * ) zn( : )
  close( unit=99 )
  nktot = product( zn( : ) )
  ntot = nbd * nktot * nspin
  !
  ! set energy zero to Fermi level or conduction band minimum
  open( unit=99, file='efermiinrydberg.ipt', form=f9, status='old' )
  rewind 99
  read ( 99, * ) efermi
  close( unit=99 )
  !
  inquire( file='temperature.ipt', exist = temp_exist )
  if( temp_exist ) then
    open( unit=99, file='temperature.ipt', form=f9, status='old' )
    read( 99, * ) temperature
    close( 99 )
  else
    temperature = 0.d0
  endif
  !
  open( unit=99, file='brange.ipt', form=f9, status='old' )
  rewind 99
  read ( 99, * ) ivl, ivh, icl, ich
  close( unit=99 )
  nbtot = 1 + ivh - ivl + 1 + ich - icl
  allocate( ww( nbtot, nktot, nspin ) )
  do ispin = 1, nspin
    do i = 1, nktot
       call enkread( 99, i, .true., 1, nbtot, ww( 1, i, ispin ) ) 
    end do
  enddo
  close( unit=99 )
  write ( stdout, * ) maxval( ww( ivh, :, : ) )
  write ( stdout, * ) minval( ww( icl, :, : ) )
  write ( stdout, '(1a5,1f10.3)' ) 'cs = ', eshift
  !
  ! energy shifted and stretched
!  read ( stdin, * ) edge, sc
  open( unit=99, file='cksshift', form=f9, status='old')
  read( 99, * ) edge
  close( 99 )
  !
  open( unit=99, file='cksstretch', form=f9, status='old')
  read( 99, * ) sc
  close( 99 )
  !
!  read ( stdin, * ) dq( : )
  open( unit=99, file='cksdq', form=f9, status='old')
  read( 99, * ) dq( : )
  close( 99 )
  if (conduct) then
    do i = 1, 3
      dbeta( i ) = sum( dq( : ) * avec( :, i ) ) / ( 2.0d0 * pi )
    end do
  else
    dbeta(:) = 0.0d0
  endif
  !
  allocate( iq( 3, nktot ), qphys( 3, nktot ), e0( ntot ), eraw( ntot ) )
  !
  ! get core level info
  open( unit=99, file='ZNL', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) atno, nc, lc
  close( unit=99 )
  write ( add04, '(1a1,1i3.3)' ) 'z', atno
  write ( add10, '(1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'z', atno, 'n', nc, 'l', lc
  !
  ! get number of sites
!  open(unit=99, file='sitelist', form='formatted', status='old')
!  read( 99, * ) ntau
  read ( stdin, * ) ntau
  allocate( tau( 3, ntau ), fntau( ntau ) )
  itau = 0
  if (conduct) then
    fnroot = 'cksc.'
  else
    fnroot = 'cksv.'
  endif
  do itau = 1, ntau
     read ( stdin, * ) element, indx, fnroot
!     read( 99, * ) element, indx
     call snatch( element, indx, tau( 1, itau ) )
     write ( fntau( itau ), '(1a5,1i4.4)' ) fnroot, itau
  end do
!  close(99)
  !
  call nbseprjprep( lmin, lmax, npmax, dqproj, nqproj, add04 )
  allocate( nproj( lmin : lmax ) )
  call nbseprjnproj( lmin, lmax, nproj, add04 )
  nprj = 0
  do l = lmin, lmax
     nprj = nprj + nproj( l ) * ( 2 * l + 1 )
  end do
  allocate( prj( ntot, nprj, ntau ) )
!  allocate( coeff( -lmax : lmax, 1 : nbd, npmax, lmin : lmax, ntau ) )
  allocate( cwgt( -lmax : lmax, npmax, lmin : lmax, ntau ) )
  allocate( fttab( nqproj, npmax, lmin : lmax ) )
  call nbseftread( nqproj, npmax, nproj, lmin, lmax, fttab, add04 )
  nlm = nint( ( 1 + lmin + lmax ) * dble( 1 + lmax - lmin ) )
  allocate( lml( nlm ), lmm( nlm ) )
  call setlm( nlm, lml, lmm, lmin, lmax, nproj, nptot )
  allocate( pcoefr( nptot, ntot, ntau ), pcoefi( nptot, ntot, ntau ) )
  allocate( wnam( nktot, nspin ) ) !, w( nbd ) )
  open( unit=99, file='listwfile', form=f9, status='old' )
  rewind 99
  do ispin = 1, nspin
    do i = 1, nktot
       read ( 99, * ) ii, wnam( i, ispin )
    end do
  end do
  close( unit=99 )
  !
  ! detect jdftx
  inquire( file='dft', exist=is_jdftx )
  if( is_jdftx ) then
    open( unit=99, file='dft', form=f9, status='old' )
    read(99,*) DFT
    if( DFT .eq. 'jdf' ) then
      is_jdftx = .true.
    else
      is_jdftx = .false. 
    endif
  endif
  
  write ( 6, * ) 'looping over k to load...'
  ii = 0
  nq = 0
  if( metal ) then
    open(unit=20,file='ibeg.h')
    allocate( ibeg_array( nktot, nspin ) )
  endif
!!  do nq = 1, nktot
!!    ik3 = mod( (nq - 1 ), zn( 3 ) )
!!    ik2 = mod( floor( dble(nq - 1 )/ dble(zn( 3 ) )), zn( 2 ) )
!!    ik1 = floor( dble( nq - 1 ) / dble( zn(3 ) * zn( 2 ) ) )
!!    write(6,*) nq, ik3, ik2, ik1
!!  enddo
!
!$OMP PARALLEL DO COLLAPSE( 2 ) &
!$OMP& SCHEDULE( STATIC  ) &
!$OMP& PRIVATE(ik1, ik2, ik3, qraw, nq, i, ng, g, zzr, zzi, ibeg, w, nrm, su, j, betot, qsqd, kvc, &
!$OMP&         ck, ii, ibd, itau, ip, ilm, l, m, iproj, fh, coeff, flip_g, ispin ) &
!$OMP& SHARED(zn, qbase, dbeta, wnam, nbtot, metal, ww, efermi, ivh, ivl, bmet, rm1, is_jdftx, &
!$OMP&        ntau, tau, lmin, lmax, nproj, npmax, nqproj, dqproj, fttab, prefs, edge, sc, nspin, &
!$OMP&        eshift, nbd, nlm, lml, lmm, pcoefr, pcoefi, e0, temperature, conduct, bvec, ibeg_array )   
!  do ik1 = 0, zn( 1 ) - 1
!     do ik2 = 0, zn( 2 ) - 1
!        do ik3 = 0, zn( 3 ) - 1
  do ispin = 1, nspin
   do nq = 1, nktot
    ik3 = mod( (nq - 1 ), zn( 3 ) ) 
    ik2 = mod( floor( dble(nq - 1 )/ dble(zn( 3 ) )), zn( 2 ) )
    ik1 = floor( dble( nq - 1 ) / dble( zn(3 ) * zn( 2 ) ) )
!           open( unit=99, file='kdpprog', form=f9, status=u7 )
!           rewind 99
!           write ( 99, '(2x,3i5,5x,3i5)' ) ik1, ik2, ik3, zn( : )
!           close( unit=99 )              
           !
           qraw( 1 ) = ( qbase( 1 ) + dble( ik1 ) ) / dble( zn( 1 ) ) + dbeta( 1 )
           qraw( 2 ) = ( qbase( 2 ) + dble( ik2 ) ) / dble( zn( 2 ) ) + dbeta( 2 )
           qraw( 3 ) = ( qbase( 3 ) + dble( ik3 ) ) / dble( zn( 3 ) ) + dbeta( 3 )
           !
!           nq = nq + 1
!!           nq = 1 + ik3 + zn(3)*ik2 + zn(3)*zn(2)*ik1
           !
           fh = 1000 !+ OMP_GET_THREAD_NUM()
!$         fh = fh + OMP_GET_THREAD_NUM()
           if( is_jdftx ) then
!JTV need to switch to more compliant, possibly access='stream'??
            stop
!             open( unit=fh, file=wnam( nq ), form='binary', status='old' )
           else
             open( unit=fh, file=wnam( nq, ispin ), form='unformatted', status='old' )
           endif
           rewind fh
           read ( fh ) ng
           allocate( g( ng, 3 ), zzr( ng, nbtot ), zzi( ng, nbtot ), w(nbd) )

           if( is_jdftx ) then
             allocate( flip_g( 3, ng ) )
             read( fh ) flip_g
             g = transpose( flip_g )
             deallocate( flip_g )
           else
             read ( fh ) g  
           endif

           read ( fh ) zzr  
           read ( fh ) zzi  
           close( unit=fh )
           !
           if ( conduct ) then
             if ( metal ) then
               if( temperature .gt. 0.000001 ) then
!                 write(6,*) temperature, ivh - ivl + 2
                 ibeg = ivh -ivl + 2
               else
                ! ww is stored with first the valence bands ( ivh -ivl + 1 of them), then the conduction
                do i = nbtot, 2 + ivh - ivl, -1 ! Changed 7 Aug to 2
!!                do i = nbtot, icl, -1 ! Changed 7 Aug to 2
                   if ( ww( i, nq, ispin ) .gt. efermi ) ibeg = i
                end do      
               endif
!               write(20,*) nq, ibeg
                ibeg_array( nq, ispin ) = ibeg
             else
                ibeg = 1 + ( 1 + ivh - ivl )
            end if
           else
             if ( metal ) stop 'metal rixs not implemented'
             ibeg = ivl
             if (ivh .lt. nbd) stop 'not enough bands for valence...'
           endif
           if ( nbtot .lt. ibeg + nbd - 1 ) stop 'not enough bands...'
           do i = 1, nbd
!              write(6,*) nq, ibeg, nbtot, nbd
              w( i ) = ww( ibeg + i - 1, nq, ispin )
              zzr( :, i ) = zzr( :, ibeg + i - 1 )
              zzi( :, i ) = zzi( :, ibeg + i - 1 )
              nrm = sum( zzr( :, i ) ** 2 + zzi( :, i ) ** 2 ) 
              su = 0      
              do j = 1, ng
                 betot( : ) = qraw( : ) + dble( g( j, : ) )
                 qsqd = dot_product( betot, matmul( bmet, betot ) )
                 su = su + qsqd * ( zzr( j, i ) ** 2 + zzi( j, i ) ** 2 ) 
              end do    
           end do
           allocate( kvc( 3, ng ) )
           do i = 1, ng
              kvc( :, i ) = g( i, : ) 
           end do
           allocate( ck( ng, 1 : nbd ) )
           ck( :, 1 : nbd ) = zzr( :, 1 : nbd ) + rm1 * zzi( :, 1 : nbd )
           allocate( coeff( -lmax : lmax, 1 : nbd, npmax, lmin : lmax, ntau ) )
           call nbsecoeffs( ng, kvc, bmet, bvec, ck, 1, nbd, qraw, ntau, tau, lmin, lmax, nproj, npmax, &
                nqproj, dqproj, fttab, coeff, prefs, temperature, efermi, w )
           !
           ! on LHS, w is in Hartree
           ! on RHS, sc is unitless, edge & cs are in eV, w is in Rydberg
           w( : ) = ( edge + sc * ( w( : ) * 13.6057d0 + eshift - edge ) ) / 27.2114d0
           !
           ! tabulate coeffs for OBFs and projector per Bloch state
           ii = ( nq - 1 ) * nbd + ( ispin - 1 ) * nktot * nbd
!           write(6, *) nq, ii
           do ibd = 1, nbd
              ii = ii + 1
           !   if (nq .eq. 1) then
                e0( ii ) = w( ibd )
           !   endif
!             eraw( ii ) = wraw( ibd )
              do itau = 1, ntau
                 ip = 0
                 do ilm = 1, nlm
                    l = lml( ilm )
                    m = lmm( ilm )
                    do iproj = 1, nproj( l )
                       ip = ip + 1
                       pcoefr( ip, ii, itau ) = coeff( m, ibd, iproj, l, itau )
                       pcoefi( ip, ii, itau ) = -rm1 * coeff( m, ibd, iproj, l, itau )
                    end do
                 end do
              end do
           end do
           !
           deallocate( g, zzr, zzi, ck, kvc, w, coeff )
           !  
!        end do
!     end do
   end do
  end do
!$OMP END PARALLEL DO 
  if( metal ) then
    do ispin = 1, nspin
      do nq = 1, nktot
        write(20,*) nq, ibeg_array( nq, ispin )
      enddo
    enddo
    close(20)
  endif
  write ( 6, * ) 'band states are solved'
  ! 
  ! output results generic to all cases
  if (conduct) then
    infoname = 'wvfcninfo'
  else
    infoname = 'wvfvainfo'
  endif
  open( unit=99, file=infoname, form=u11, status=u7 )
  rewind 99
  write ( 99 ) nbd, zn(1)*zn(2)*zn(3), nspin
  write ( 99 ) e0( : )
  close( unit=99 )
  if (conduct) then
    infoname = 'wvfcninfo2'
  else
    infoname = 'wvfvainfo2'
  endif
  open( unit=99, file=infoname, form=f9, status=u7 )
  rewind 99
  write ( 99, * ) nbd, zn(1)*zn(2)*zn(3) 
  write ( 99, * ) e0( : )
  close( unit=99 )

  !
  ! For easy accounting spin was glommed on to ntot, but we need to pull it out now
  ntot = ntot / nspin
  ! output information site specific, i.e. Bloch fcn overlap with projectors
  do itau = 1, ntau
     write ( 6, * ) ' ... ', itau, ntau, fntau( itau )
     open( unit=99, file=fntau( itau ), form=u11, status=u7 )
     rewind 99
     write ( 99 ) nptot, ntot, nspin
     write ( 99 ) tau( :, itau )
     write ( 99 ) pcoefr( :, :, itau )
     write ( 99 ) pcoefi( :, :, itau )
     close( unit=99 )
  end do
  !
  ! output raw band energies
! open( unit=99, file='rawband', form='unformatted', status='unknown' )
! rewind 99
! write ( 99 ) eraw
! close( unit=99 ) 
  !  
  return
end subroutine cainkset
