! program gbuilder
!
! Original subroutie by E. Shirley
! Modified by J. Vinson  feb. 09
!
! Purpose: makes ximat and !!!! not rbfile.bin
!
! Inputs: brange.ipt, avecsinbohr.ipt, bvecs, specpnt, pquadrature
!   enkfile, listwfile, wavefunctions
! Inputs from stdin: nkpts, bands
!
!!!!
!subroutine gbuilder( nv, ng, avec, bvec, bmet, uu, &
!     & kvc, tr, ti, vr, vi, kr, ki, cfr, cfi )
program builder
  implicit none
  !     
  integer nv, ng
  double precision avec( 3, 3 ), bvec( 3, 3 ), bmet( 3, 3 )
  !
  character * 7, parameter :: u7 = 'unknown'
  character * 9, parameter :: f9 = 'formatted'
  character * 11, parameter :: u11 = 'unformatted'
  !
  integer, parameter :: stdin = 5
  integer, parameter :: stdout = 6
  !
  integer :: i, j, ibl, ibh, nk1, nk2, nk3, ik1, ik2, ik3, it, nt, npt, ii, overlap
  !
  real( kind = kind( 1.0d0 ) ) :: celvol, pi, fr, fi, phse, su, pref, su2
  real( kind = kind( 1.0d0 ) ) :: denr, deni, iden2, mu, x, s, qcart( 3 ), qin( 3 ), qout( 3 )
  real( kind = kind( 1.0d0 ) ) :: gr, gi, su0, rmax
  real( kind = kind( 1.0d0 ) ) :: vlev, vhev, clev, chev, muev, sev, mindif, maxdif, vec( 3 )
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: zr( :, : ), zi( :, : ), w( : ), w2( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vipt( : ), drel( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: bre( :, : ), bim( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: t( : ), wgt( : ), newwgt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: wfr( : ), wfi( : ), wfp( :, : ), dmat( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: cosqr( : ), sinqr( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ure( :, : ), uim( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: gre( :, :, : ), gim( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vind( : ), nind( : ), xirow( : )
  !
  real( kind = kind( 1.0d0 ) ) :: term1, term2, rscal, dsqdmod1, dsqdmod2, dsqdmod3, exppref, nav
  real( kind = kind( 1.0d0 ) ) :: efermi
  integer :: dumint, gtot,bandtot,brange(4), nang, nr, indx
  integer, allocatable :: gvec( : , : )
  character*12 :: wfname
  character*2 :: element
  real( kind = kind( 1.0d0 ) ) :: qinb(3), kshift(3)
  !
  ! file units by name
  integer, parameter :: enkfile = 20
  integer, parameter :: listwfile = 21
  !
  !
  pi = 4.0d0 * atan( 1.0d0 )
  ! brange
  open(unit=99,file='brange.ipt',form='formatted',status='old')
  read(99,*) brange(:)
  close( 99 )
  bandtot = brange(2)-brange(1) + brange(4) - brange(3) + 2
  write(6,*) brange(:)
  !
  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  read(99,*) avec(:,:)
  close(99)
  write(6,*) 'avecs loaded'
  !
  open(unit=99,file='efermiinrydberg.ipt',form='formatted',status='old')
  read(99,*) efermi
  close(99)
  !
  open(unit=99,file='bvecs',form='formatted',status='old')
  read(99,*) bvec(:,:)
  close(99)
  write(6,*) 'bvecs loaded'
  call getomega( avec, celvol )
  write(6,*) 'celvol =', celvol
  !
  open(unit=99,file='kmesh.ipt',form=f9,status='old')
  read(99,*) nk1, nk2, nk3
  close(99)
  !
  open(unit=99,file='k0.ipt',form=f9,status='old')
  read(99,*) kshift(:)
  close(99)
  !
  open(unit=99,file='clipbands',form=f9,status='old') 
  read(99,*) ibl,ibh
  close(99)
  !
!!!!
  overlap = brange(2) - brange(3) + 1
  if (overlap .lt. 0) stop "!!!"
  !! ibh = ibh + overlap
!!!!
  !
  open( unit=99, file='rbfile.bin', form=u11, status=u7 )
  rewind 99
  read ( 99 ) npt
  allocate(posn(3,npt), wpt(npt), drel(npt))
  read ( 99 ) posn
  read ( 99 ) wpt
  read ( 99 ) drel 
  close( unit =99 )

  nav = 0
  ! need to get the points on the sphere and their weights
  !
  open( unit=99, file='specpnt', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) nang
  close(99)

  allocate(wfp(npt,npt),dmat(npt,npt))
  allocate(vipt(npt),vind(npt),nind(npt))
  call mkvipt( npt, drel, vipt )
  !
  open( unit=99, file='Pquadrature', form=f9, status=u7 )
  rewind 99
  read ( 99, * ) nt
  allocate( t( nt ), wgt( nt ), newwgt( nt ) )
  su = 0
  do it = 1, nt
     read ( 99, * ) t( it ), wgt( it )
     t( it ) = ( 1 + t( it ) ) / 2
     su = su + wgt( it )
  end do
  wgt = wgt / su
  do it = 1, nt
     newwgt( it ) = wgt( it ) / ( 1.d0 - t( it ) ) ** 2
  end do
  close( unit=99 )
  !
  ! What follows if outlined in Ultramicroscopy 106 (2006) 986 by E.L. Shirley
  ! to summarize;
  ! To simplify the expression of \Chi_0 we will calculate the Green's function integrating 
  !  over a dummy variable which is scaled by the gemoetric mean of the smallest and largest
  !  values of abs( \epsilon_i - \mu )
  write(6,*) "calling vcbder"
  call vcbder2( ibl, ibh, vlev, vhev, clev, chev )
  write(6,*) "done with vcbder"
  muev = ( clev + vhev ) / 2.0d0
  mindif = min( muev - vhev, clev - muev )
  maxdif = max( muev - vlev, chev - muev )
  sev = sqrt( mindif * maxdif )
  write( 6, '(4(1x,1e15.8))' ) vlev, vhev, clev, chev
  write ( 6, '(4(1x,1e15.8))' ) muev, mindif, maxdif, sev
  mu = muev / 13.6057d0    ! s & mu are in Ry!
  s = sev / 13.6057d0      ! s & mu are in Ry!
  !
  !
  allocate( ure( npt, ibl : ibh + overlap ), uim( npt, ibl : ibh + overlap ) )
  allocate( gre( npt, npt, 2 * nt ), gim( npt, npt, 2 * nt ) )
  gre = 0; gim = 0; dmat = 0
  allocate( wfr( npt ), wfi( npt ), cosqr( npt ), sinqr( npt ) )
  pref = 1 / ( nk1 * nk2 * nk3 * celvol )
  ! 

  ! open up the enkfile, read through while in kpt loop
  open(unit=enkfile,file='enkfile',form='formatted',status='old')
  rewind(enkfile)
  ! w is the energy list
  allocate(w(bandtot), w2(brange(4)-brange(1)+1))
  !
  ! open up the listwfile to get the wf names
  open(unit=listwfile,file='listwfile',form='formatted',status='old')
  rewind(listwfile)
  ! need to make sure kpt naming scheme and enkfile reflect this kpoint order
  ! otherwise might need to interchange loop order, won't change anything else
  do ik1=1,nk1
     qin(1) = ( kshift(1) + dble( ik1 - 1 )) / dble( nk1 )
     do ik2=1,nk2
        qin(2) = ( kshift(2) + dble( ik2 - 1 )) / dble( nk2 )
        do ik3=1,nk3
           qin(3) = ( kshift(3) + dble( ik3 - 1 )) / dble( nk3 )
           qcart(:) = 0.0
           do i=1,3
              ! changed 16 april
              ! changed july 7
              qcart(:) = qcart(:) + bvec(:,i)*qin(i)
           enddo
           write(6,*) qin(:), qcart(:)
           !
           ! get energies for this kpoint, read and store all bands
           read(enkfile,*) w(:)
           ! why even go to w2?
           w2(1:brange(2))=w(1:brange(2))
           w2(brange(2)+1:brange(4)) = w(brange(2)+overlap+1:brange(4)+overlap)
           ! get the wfks for this kpt
           read(listwfile,*) dumint, wfname
           ! open wf file and read in wavefunction info
           open(unit=22,file=wfname,form='unformatted',status='old')
           rewind(22)
           read(22) gtot
           allocate(gvec(gtot,3),zr(gtot,bandtot),zi(gtot,bandtot))
           read(22) gvec
           read(22) zr
           read(22) zi
           close(22)

!           open( unit=97 )
!           do i = ibl, ibh
!              do j = ibl, ibh
!                 su = 0
!                 do it = 1, nt
!                    x = 1.0d0 / ( 1.0d0 - t( it ) )
!                    deni = s * t( it ) * x
!                    denr = mu - w2( i )
!                    iden2 = 1.0d0 / ( denr ** 2 + deni ** 2 )
!                    fr = iden2 * denr; fi = - iden2 * deni
!                    denr = mu - w2( j )
!                    iden2 = 1.0d0 / ( denr ** 2 + deni ** 2 )
!                    gr = iden2 * denr; gi = - iden2 * deni
!                    su = su + wgt( it ) * x ** 2 * ( fr * gr - fi * gi )
!                 end do
!                 su = su * s / pi
!                 su0 = 0
!                 ! if ( ( ( w( i ) - mu ) * ( w( j ) - mu ) ) .lt. 0 ) su0 = - 1.d0 / abs( w( i ) - w( j ) )
!                 write ( 97, '(2i5,5(1x,1e15.8))' ) i, j, w2( i ), w2( j ), mu, su, su0
!              end do
!           end do
!           close( unit=97 )
           !
           ! this will probably fail for overlapping brange, but should be fine for the test case
           call realu2( gtot, npt, ibl, ibh, brange(2), overlap, zr( 1, ibl ), zi( 1, ibl ), gvec, bvec, posn, ure, uim )
           do i = 1, npt
              phse = 0
              do j = 1, 3
                 phse = phse + qcart( j ) * posn( j, i )
              end do
              cosqr( i ) = dcos( phse )
              sinqr( i ) = dsin( phse )
           end do
           do j = ibl, brange( 2 )
              !open( unit=99, file='kdpprog', form='formatted', status='unknown' )
              !rewind 99
              !write ( 99, '(8i5)' ) ik1, ik2, ik3, 2 * nk1, 2 * nk2, 2 * nk3, j, ibh
              !close( unit=99 )
              if( w( j ) .le. efermi ) then
                do i = 1, npt
                   wfr( i ) = ure( i, j ) * cosqr( i ) - uim( i, j ) * sinqr( i )
                   wfi( i ) = ure( i, j ) * sinqr( i ) + uim( i, j ) * cosqr( i )
                end do
                do i = 1, npt ! for each psi k, conjg. at -k.  hence, we are real!
                   do ii = 1, npt
                      wfp( ii, i ) = pref * ( wfr( i ) * wfr( ii ) + wfi( i ) * wfi( ii ) )
                   end do
                end do
                dmat( :, : ) = dmat( :, : ) + 2.0d0 * wfp( :, : )
                do it = 1, nt
                   x = 1.0d0 / ( 1.0d0 - t( it ) )
                   deni = s * t( it ) * x
                   denr = mu - w( j )
                   iden2 = 1.0d0 / ( denr ** 2 + deni ** 2 )
                   fr = iden2 * denr; fi = - iden2 * deni
                   gre( :, :, it ) = gre( :, :, it ) + fr * wfp( :, : )
                  gim( :, :, it ) = gim( :, :, it ) + fi * wfp( :, : )
                end do
              endif
           end do
           !
           do j = brange( 2 ) + 1, ibh + overlap
              !open( unit=99, file='kdpprog', form='formatted', status='unknown' )
              !rewind 99
              !write ( 99, '(8i5)' ) ik1, ik2, ik3, 2 * nk1, 2 * nk2, 2 * nk3, j, ibh
              !close( unit=99 )
              if( w( j ) .gt. efermi ) then
                do i = 1, npt
                   wfr( i ) = ure( i, j ) * cosqr( i ) - uim( i, j ) * sinqr( i )
                   wfi( i ) = ure( i, j ) * sinqr( i ) + uim( i, j ) * cosqr( i )
                end do
                do i = 1, npt ! for each psi k, conjg. at -k.  hence, we are real!
                   do ii = 1, npt
                      wfp( ii, i ) = pref * ( wfr( i ) * wfr( ii ) + wfi( i ) * wfi( ii ) )
                   end do
                end do
                dmat( :, : ) = dmat( :, : ) + 2.0d0 * wfp( :, : )
                do it = 1, nt
                   x = 1.0d0 / ( 1.0d0 - t( it ) )
                   deni = s * t( it ) * x
                   denr = mu - w( j )
                   iden2 = 1.0d0 / ( denr ** 2 + deni ** 2 )
                   fr = iden2 * denr; fi = - iden2 * deni
                   gre( :, :, it ) = gre( :, :, it ) + fr * wfp( :, : )
                   gim( :, :, it ) = gim( :, :, it ) + fi * wfp( :, : )
                end do
              endif
           end do
           !
           deallocate(gvec,zr,zi)
           !
        end do ! end k runs
     end do ! end k runs
  end do ! end k runs
  !
  exppref = pi / ( 2.0d0 ** ( 5.0d0 / 3.0d0 ) )
  !
  write(6,*)"starting ximat"
  allocate( xirow( npt ) )
  open( unit=99, file='ximat', form='unformatted', status='unknown' )
  rewind 99
  vind = 0
  do i = 1, npt
     su2 = 0
     do j = 1, npt
        su = 0
        do it = 1, nt
           su = su + ( gre( i, j, it ) ** 2 - gim( i, j, it ) ** 2 ) * newwgt( it )
        end do
        su = 4 * su * s / pi     ! 2 for spin, 2 for Ry
        xirow( j ) = su
        su2 = su2 + vipt( j ) * su * wpt( j )
     end do
     write ( 99 ) xirow

     do j = 1, npt
        vind( j ) = vind( j ) + wpt( i ) * su2 / max( drel( j ), drel( i ) )
     end do
     nind( i ) = su2
  end do
  close( unit=99 )
  !
  write(6,*)"starting nopt"
  open( unit=99, file='nopt', form='formatted', status='unknown' )
  rewind 99
  do i = 1, npt
     write ( 99, '(4(1x,1e15.8))' ) drel( i ), vipt( i ), nind( i ), vind( i )
  end do
  close( unit=99 )
  !
  !  return
  !
end program builder
