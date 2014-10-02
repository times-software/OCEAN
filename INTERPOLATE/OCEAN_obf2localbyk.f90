  subroutine OCEAN_obf2localbyk
! ----------------------------------------------------------------------

  USE constants, ONLY : pi
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE cell_base, ONLY : tpiba, tpiba2, bg
  USE gvect
  USE klist, ONLY : nkstot, xk, nks, ngk
  USE wvfct
  USE io_files, ONLY: nd_nmbr, prefix, nwordwfc, iunwfc, tmp_dir, diropn
  USE wavefunctions_module, ONLY: evc
  use mp, only : mp_sum, mp_bcast, mp_barrier
  use mp_global, only : intra_pool_comm
  use shirley_ham_input, only : band_subset

  use OCEAN_timer

  implicit none

  integer,external :: freeunit
  complex(dp), parameter :: zero = ( 0.d0, 0.d0 )

  integer :: i, j, itau, ik1, ik2, ik3, l, ibnd, ip, m, iproj, ilm, k, ii, jj, kk
  integer :: iuntmp, nwordo2l, iunout, iuntxt
  character(len=6) :: nd_nmbr_tmp
  logical :: exst, have_kshift

  real(dp) :: qbase( 3 ), dqproj, bvec( 3, 3 ), bmet( 3, 3 ), qraw( 3 ), prefs(0 : 1000 ), avec( 3, 3 ), kshift( 3 )

  integer :: atno, nc, lc, ntau, lmin, lmax, npmax, nqproj, nptot, zn( 3 ), indx, nprj, numlm, ig
  integer :: nshift, ishift
  real(dp), allocatable :: tau( :, : ), fttab( :, :, : ), myg( :, : )
  integer, allocatable :: nproj( : ), lml( : ), lmm( : ), ibnd_indx( : ), sub_mill( :, : )
  complex(dp), allocatable :: o2l( :, :, : ), coeff( :, :, :, :, : )

  character(25) :: element, fnroot, add04, add10
  character(6), allocatable :: fntau( : )


  call OCEAN_t_reset

  WRITE( stdout, '(/5x,"Calling obf2localbyk .... ",/)')
  write(stdout,*) ' npwx  = ', npwx
  write(stdout,*) ' npw   = ', npw
  write(stdout,*) ' nbnd  = ', nbnd
  write(stdout,*) ' nbndx = ', nbndx

  ! sort out which band subset we will work with
  if( band_subset(1) > band_subset(2) ) then
    i=band_subset(2)
    band_subset(2) = band_subset(1)
    band_subset(1) = i
  endif
  if( band_subset(1)>=nbnd .or. band_subset(1)<=0 ) band_subset(1) = 1
  if( band_subset(2)>=nbnd .or. band_subset(2)<=0 ) band_subset(2) = nbnd

  if( band_subset(2)-band_subset(1)+1 < nbnd ) then
    write(stdout,*) ' Requested band subset:', band_subset(1), &
                    ' ... ', band_subset(2)
    write(stdout,*) ' Reducing total number of bands from ', nbnd, &
                    ' to ', band_subset(2)-band_subset(1)+1
  endif

  nbnd = band_subset(2)-band_subset(1)+1
  allocate( ibnd_indx(nbnd) )
  do i=1,nbnd
    ibnd_indx(i) = band_subset(1)+i-1
  enddo

  !
  call summary
  !
  ! ======================================================================
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  ! ======================================================================
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)
!  !call sum_band
  write(stdout,*) ' ecutwfc = ', ecutwfc
  write(stdout,*) ' tpiba2 = ', tpiba2
  write(stdout,*) ' nks, nktot = ', nks, nkstot
  write(stdout,*) ' xk = ', xk(1:3,1:nkstot)
  write(stdout,*) '     npw = ', ngk(1:nks)
  write(stdout,*) '    npwx = ', npwx




  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  write(stdout,*)
  write(stdout,*) ' construct interaction0:'
  write(stdout,*)

  ! ======================================================================
  ! Set up the Ylm prefactors; opening files so do it only on ionode
  ! ======================================================================
  if( ionode ) then
    lmax = 5
    call make_prefs( prefs, lmax )
  endif
  call mp_bcast( prefs, ionode_id )

  ! ======================================================================
  ! This is a bunch of house keeping crap that need to be integrated in to 
  ! the input parser at some point.
  !  ======================================================================

  if( ionode ) then

    iuntxt = freeunit()
    open(unit=iuntxt,file='o2l_test.txt', form='formatted' )
    rewind(iuntxt)

    iuntmp = freeunit()
!    open( unit=iuntmp, file='nspn.ipt', form='formatted', status='old' )
!    read( iuntmp, * ) nspn
!    close( iuntmp )
    open( unit=iuntmp, file='avecsinbohr.ipt', form='formatted', status='unknown' )
    rewind( iuntmp )
    read( iuntmp, * ) avec( :, : )
    close( iuntmp )
    !
    do i = 1, 3
       j = 1 + mod( i, 3 )
       k = 1 + mod( j, 3 )
       do ii = 1, 3
          jj = 1 + mod( ii, 3 )
          kk = 1 + mod( jj, 3 )
          bvec( ii, i ) = avec( jj, j ) * avec( kk, k ) - avec( kk, j ) * avec( jj, k )
       end do
       bvec( :, i ) = bvec( :, i ) * 2.0d0 * pi / dot_product( bvec( :, i ), avec( :, i ) )
    end do

    open( unit=iuntmp, file='scaledkzero.ipt', form='formatted', status='old' )
    rewind iuntmp
    read ( iuntmp, * ) qbase( : )
    close( unit=iuntmp )

    open( unit=iuntmp, file='kmesh.ipt', form='formatted', status='old' )
    rewind iuntmp
    read ( iuntmp, * ) zn( : )
    close( unit=iuntmp )

    kshift = 0.0d0
    inquire(file='qinunitsofbvectors.ipt',exist=have_kshift)
    if( have_kshift) then
      open(unit=iuntmp,file='qinunitsofbvectors.ipt',form='formatted',status='old')
      read(iuntmp,*) kshift( : )
      close(iuntmp)
      if( sum(abs(kshift( : ) ) ) .lt. 1.d-14 ) have_kshift = .false.
    endif
    if( have_kshift ) then
      write(stdout,*) 'Have k shift'
    endif

    open(iuntmp,file='ZNL', form='formatted', status='old' )
    read ( iuntmp, * ) atno, nc, lc
    close( iuntmp )
    write ( add04, '(1a1,1i3.3)' ) 'z', atno
    write ( add10, '(1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'z', atno, 'n', nc, 'l', lc
!JTV conduction/valence
    fnroot = 'cksc.'

    open(iuntmp,file='cks.in', form='formatted', status='old' )
    read(iuntmp,*) ntau
    allocate( tau( 3, ntau ) )
    allocate( fntau( ntau ) )
    do itau=1,ntau
      read(iuntmp,*) element, indx, fnroot
      call snatch( element, indx, tau( 1, itau ) )
      write ( fntau( itau ), '(1a2,1i4.4)' ) element, indx
    enddo
    close(iuntmp)
  !======================================================================
  ! Need to load up and initialize all of the projectors
  ! These will end up repeated over all processes
  ! ======================================================================
    call nbseprjprep( lmin, lmax, npmax, dqproj, nqproj, add04 )
    allocate( nproj( lmin : lmax ) )
    call nbseprjnproj( lmin, lmax, nproj, add04 )
    nprj = 0
    do l = lmin, lmax
       nprj = nprj + nproj( l ) * ( 2 * l + 1 )
    end do
    allocate( fttab( nqproj, npmax, lmin : lmax ) )
    call nbseftread( nqproj, npmax, nproj, lmin, lmax, fttab, add04 )
    numlm = nint( ( 1 + lmin + lmax ) * dble( 1 + lmax - lmin ) )
    allocate( lml( numlm ), lmm( numlm ) )
    call setlm( numlm, lml, lmm, lmin, lmax, nproj, nptot )
  endif

  write(stdout,*) ' Read in little data files'
  call mp_barrier
  write(stdout,*) ' Preparing to share little data files'

!  mp_bcast( nspn, ionode_id )
  call mp_bcast( qbase, ionode_id )
  call mp_bcast( zn, ionode_id )
  call mp_bcast( have_kshift, ionode_id )
  call mp_bcast( kshift, ionode_id )
  call mp_bcast( ntau, ionode_id )
  call mp_bcast( lmin, ionode_id )
  call mp_bcast( lmax, ionode_id )
  call mp_bcast( npmax, ionode_id )
  call mp_bcast( dqproj, ionode_id )
  call mp_bcast( nqproj, ionode_id )
  call mp_bcast( nprj, ionode_id )
  call mp_bcast( numlm, ionode_id )
  call mp_bcast( nptot, ionode_id )
  if( .not. ionode ) then
    allocate( tau( 3, ntau ) )
    allocate( nproj( lmin : lmax ) )
    allocate( fttab( nqproj, npmax, lmin : lmax ) )
    allocate( lml( numlm ), lmm( numlm ) )
  endif
  call mp_bcast( tau, ionode_id )
  call mp_bcast( nproj, ionode_id )
  call mp_bcast( fttab, ionode_id )
  call mp_bcast( lml, ionode_id )
  call mp_bcast( lmm, ionode_id )
  call mp_barrier
  write(stdout,*) ' Done sharing little data files'
  call mp_bcast( bvec, ionode_id )


  ! ======================================================================
  ! Now we are ready to allocate the output and prep the output file
  ! ======================================================================
  if( ionode ) then
    nwordo2l = 2 * nbnd * nptot * ntau
    iunout = freeunit()
    open( unit=iunout, form='formatted')
    write(iunout,* ) mill( : , : )
    close( iunout )
    iunout = freeunit()

    ! since only the ionode is going to write to this file
    ! there is no need for the nd_nmbr suffix
    nd_nmbr_tmp = nd_nmbr
    nd_nmbr=''    
    call diropn( iunout, 'o2l', nwordo2l, exst )
    ! restore in case needed again
    nd_nmbr = nd_nmbr_tmp
    write(stdout,*) ' will be saved to file: ', trim(tmp_dir)//trim(prefix)//'.o2l'
  endif

  allocate( o2l( nbnd, nptot, ntau ) )


! ! for simplicity right now 
  allocate( myg( 3, npw ) )
  allocate( sub_mill( 3, npw ) )
  sub_mill = 0
  do ig = 1, npw
    myg( :, ig ) = tpiba * g( :, igk( ig ) )
    sub_mill( :, ig ) = matmul( g( :, igk( ig ) ), bg( :, : ) )
  enddo
  if( ionode ) then
    open(unit=99,form='formatted')
      do ig = 1, npw
        write( 99,*) sub_mill( :, : ) 
      enddo
    close( 99 )
  endif


!JTV not sure about these guys
  !bvec( :, : ) = bg( :, : )
  do i = 1, 3
    do j = 1, 3
      bmet( i, j ) =dot_product( bvec( :, i ), bvec( :, j ) )
    end do
  end do

  call OCEAN_t_printtime( "Stupid prep", stdout )

  nshift = 1
  if( have_kshift ) nshift = 2
  allocate( coeff( -lmax : lmax, 1 : nbnd, npmax, lmin : lmax, ntau ) )
  ! loop over all kpts
  i = 0
  do ishift = 1, nshift
    do ik1 = 0, zn( 1 ) - 1
      do ik2 = 0, zn( 2 ) - 1
         do ik3 = 0, zn( 3 ) - 1
          o2l = zero
          coeff = zero
          i = i + 1
          write(stdout,*) i
          
          qraw( 1 ) = ( qbase( 1 ) + dble( ik1 ) ) / dble( zn( 1 ) )
          qraw( 2 ) = ( qbase( 2 ) + dble( ik2 ) ) / dble( zn( 2 ) )
          qraw( 3 ) = ( qbase( 3 ) + dble( ik3 ) ) / dble( zn( 3 ) )
           
          call OCEAN_t_reset
          call nbsecoeffs( npw, mill, myg, bmet, bvec, evc, 1, nbnd, qraw, ntau, tau, lmin, &
                           lmax, nproj, npmax, nqproj, dqproj, fttab, o2l, prefs, lml, lmm, numlm, nptot )
          call OCEAN_t_printtime( "nbsecoeffs", stdout )

  !JTV I'm assuming here that the basis functions will be distributed by G
  !      such that every proc has every B_i but only B_i( local G )
          call OCEAN_t_reset
!          call mp_sum( coeff )!, intra_pool_comm )
          call mp_sum( o2l )
          call OCEAN_t_printtime( "mp_sum", stdout )

          call OCEAN_t_reset
!          do itau = 1, ntau
!            do ibnd=1,nbnd
!            ip = 0
!              do ilm = 1, numlm
!                l = lml( ilm )
!                m = lmm( ilm )
!                do iproj = 1, nproj( l )
!                  ip = ip + 1
!                  ! interleavig here for later
!                  o2l( ibnd, ip, itau ) = coeff( m, ibnd, iproj, l, itau )
!                end do
!              end do
!             end do
!          enddo
          if( ionode ) call davcio( o2l, nwordo2l, iunout, i, +1 )
                  if( ionode ) write(iuntxt,*)  o2l(  1, 1, 1 )
          call OCEAN_t_printtime( "Reorder & write", stdout )

        enddo
      enddo
    enddo
    qbase( : ) = qbase( : ) + kshift( : )
  enddo

  if( ionode ) then 
    close(iunout )
    iunout = freeunit()
    ! Unformatted so we don't round tau at all
    call seqopn(iuntmp, 'o2li', 'unformatted', exst)
    if( exst ) rewind(iuntmp)
    write(iunout) nptot, ntau 
    write(iunout) tau( :, : )
    write(iunout) fntau( : )
    close(iunout)
  endif

  deallocate( o2l )


  deallocate( tau )
  deallocate( nproj )
  deallocate( fttab )
  deallocate( lml, lmm )

  return


  contains

!
! ng = number of Gvecs for the writing of each wavevector
! kvc = inverse of the g list
! bmet 
! bvec
! ck = complex wavefunction
! ibl = start band
! ibh = stop band
  subroutine nbsecoeffs( ng, kvc, gvecs, bmet, bvec, ck, ibl, ibh, q, ntau, tau,  &
                         lmin, lmax, nproj, npmax, nqproj, dqproj, fttab, o2l, prefs, & 
                         lml, lmm, numlm, nptot )
!    use constants, ONLY : tpi
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, ibl, ibh, ntau, numlm, nptot
    integer, intent( in ) :: lmin, lmax, npmax, nqproj
    integer, intent( in ) :: kvc( 3, ng ), nproj( lmin : lmax ), lml( numlm ), lmm( numlm )
    real(dp), intent( in ) :: dqproj, bvec( 3, 3 ), bmet( 3, 3 ), q( 3 ), prefs( 0 : 1000 )
    real(dp), intent( in ) :: tau( 3, ntau ), gvecs( 3, ng )
    real(dp), intent( in ) :: fttab( nqproj, npmax, lmin : lmax )
    complex(dp), intent( in ) :: ck( ng, ibl : ibh )
!    complex(dp), intent( out ) :: coeff( -lmax:lmax, ibl:ibh, npmax, lmin:lmax, ntau )
    complex(dp), intent( out ) :: o2l( nbnd, nptot, ntau )
    
    !
    integer itau, l, ig, ip
    real(dp) :: qphase, phase
    complex(dp), allocatable :: tauphs( : , : )
    complex(dp), allocatable :: ylmfac( : , : )
    real(dp), allocatable :: sfq( : , : , : )
    complex(dp) :: rm1
    real(dp) :: pi 
    !
    pi = 4.0d0 * atan( 1.0d0 )
    rm1 = -1.d0
    rm1 = sqrt( rm1 )
    allocate( tauphs( ng, ntau ) )
    allocate( sfq( ng, npmax, lmin : lmax ) )
!    allocate( ylmfac( ng * ( 2 * lmax + 1 ), lmin : lmax ) )
    allocate( ylmfac( ng, numlm ) )
    !
    do itau = 1, ntau
!       call getphase( ng, kvc, q, tau( 1, itau ), tauphs( 1, itau ) )
      qphase = 2.0d0 * pi * sum( tau( :, itau ) * q( : ) )
      do ig = 1, ng
        phase = qphase + 2.0d0 * pi * sum( tau( :, itau ) * kvc( :, ig ) )
        tauphs( ig, itau ) = cos( phase ) + rm1 * sin( phase ) !cmplx( cos(phase), sin(phase) )
      enddo
    end do
    call OCEAN_t_printtime( "Phase", stdout )
    call OCEAN_t_reset

    do l = lmin, lmax
       call seanitup( ng, q, kvc, bmet, l, nproj( l ), sfq( 1, 1, l ), &
                      npmax, nqproj, dqproj, fttab( 1, 1, l ) )
!       call getylmfac( ng, kvc, q, bvec, l, ylmfac( 1, l ), prefs )
    end do
    call newgetylmfac( ng, kvc, q, bvec, numlm, lml, lmm, ylmfac, prefs )
    call OCEAN_t_printtime( "Seanitup", stdout )
    call OCEAN_t_reset

    if( .false. ) then
    do itau = 1, ntau
       ip = 0
       do l = lmin, lmax
          call getcoeff( l, ng, 1 + ibh - ibl, ck, tauphs( 1, itau ), ylmfac,  &
!                         coeff( -lmax, ibl, 1, l, itau ), lmax, npmax, sfq( 1, 1, l ), &
                         o2l( 1, 1, itau ), lmax, npmax, sfq( 1, 1, l ), &
                         nproj( l ), lml, lmm, numlm, nptot, ip )
       end do
    end do
    else
      call fullgetcoeff( ng, 1 + ibh - ibl, ntau, lmin, lmax, npmax, nproj, ck, tauphs, ylmfac, &
                           o2l, sfq, lml, lmm, numlm, nptot )
    endif
    call OCEAN_t_printtime( "Get Coeff", stdout )
    !
    deallocate( tauphs, ylmfac, sfq )
    !
    return

  end subroutine nbsecoeffs


  subroutine seanitup( ng, q, gvec, bmet, l, nproj, &
       &               seanfq, npmax, nq, dq, fttab )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, l, nproj, nq, npmax
    integer, intent( in ) :: gvec( 3, ng )
    !
    real(dp), intent( in ) :: dq, q( 3 ), bmet( 3, 3 )
    real(dp), intent( in ) :: fttab( nq, nproj )
    real(dp), intent( out ) :: seanfq( ng, nproj )
    !
    integer :: i, j, i1, i2, ii
    real(dp) :: gc( 3 ), qbar, qii, x
    real(dp) :: wm1, wz0, wp1, wp2
    !
    real(dp), parameter :: pm1 = - 1.d0 / 6.d0
    real(dp), parameter :: pz0 = 1.d0 / 2.d0
    real(dp), parameter :: pp1 = - 1.d0 / 2.d0
    real(dp), parameter :: pp2 = 1.d0 / 6.d0


    do j = 1, ng
       gc( : ) = q( : ) + dble(gvec( :, j ))
       qbar = 0.d0
       do i1 = 1, 3
          do i2 = 1, 3
             qbar = qbar + gc( i1 ) * gc( i2 ) * bmet( i1, i2 )
          end do
       end do
       qbar = sqrt( qbar )
       ii = 1 + qbar / dq
       if ( ii .lt. 2 ) ii = 2
       if ( ii .gt. nq - 2 ) ii = nq - 2
       qii = dq * ( ii - 1 )
       x = ( qbar - qii ) / dq
       wm1 = pm1 * x * ( x - 1.d0 ) * ( x - 2.d0 )
       wz0 = pz0 * ( x + 1.d0 ) * ( x - 1.d0 ) * ( x - 2.d0 ) 
       wp1 = pp1 * ( x + 1.d0 ) * x * ( x - 2.d0 )
       wp2 = pp2 * ( x + 1.d0 ) * x * ( x - 1.d0 )
       do i = 1, nproj
          seanfq( j, i ) = fttab( ii - 1, i ) * wm1 + &
               &           fttab( ii, i ) * wz0 + &
               &           fttab( ii + 1, i ) * wp1 + &
               &           fttab( ii + 2, i ) * wp2
       end do
    end do
    !
    return
  end subroutine seanitup

  subroutine getylmfac( ng, gvec, q, bvec, l, ylmfac, prefs )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, l
    integer, intent( in ) :: gvec( 3, ng )
    real(dp), intent( in ) :: q( 3 ), bvec( 3, 3 ), prefs( 0 : 1000 )
    complex(dp), intent( out ) :: ylmfac( -l : l, ng )
    !
    integer :: jj, ig, m
    real(dp) :: x( 3 ), pi
    complex(dp) :: pref, ylm, rm1
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    pref = 4.0d0 * pi * rm1 ** l
    !
    do ig = 1, ng
       x = 0
       do jj = 1, 3
          x( : ) = x( : ) + bvec( :, jj ) * ( q( jj ) + gvec( jj, ig ) )
       end do
       do m = -l, l
          call getylm( l, m, x( 1 ), x( 2 ), x( 3 ), ylm, prefs )
          ylmfac( m, ig ) = pref * conjg( ylm )
       end do
    end do
    !
    return
  end subroutine getylmfac

  subroutine newgetylmfac(ng, gvec, q, bvec, numlm, lml, lmm, ylmfac, prefs )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, numlm
    integer, intent( in ) :: gvec( 3, ng ), lml( numlm ), lmm( numlm )
    real(dp), intent( in ) :: q( 3 ), bvec( 3, 3 ), prefs( 0 : 1000 ) 
    complex(dp), intent( out ) :: ylmfac( ng, numlm )
    !
    integer :: jj, ig, m, l
    real(dp) :: x( 3 ), pi
    complex(dp) :: pref, ylm, rm1
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    pref = 4.0d0 * pi * rm1 ** l
    !
    do ilm = 1, numlm
      l = lml( ilm )
      m = lmm( ilm )
      do ig = 1, ng
        x = 0
        do jj = 1, 3
          x( : ) = x( : ) + bvec( :, jj ) * ( q( jj ) + gvec( jj, ig ) )
        enddo
        call getylm( l, m, x( 1 ), x( 2 ), x( 3 ), ylm, prefs )
        ylmfac( ig, ilm ) = pref * conjg( ylm )
      enddo
    enddo


  end subroutine newgetylmfac

  subroutine getcoeff( l, ng, nbd, ck, tauphs, ylmfac, o2l, lmax, npmax, seanfq, nproj, &
                       lml, lmm, numlm, nptot, ip )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: l, ng, nbd, nproj, lmax, npmax, numlm, nptot
    integer, intent( in ) :: lml(numlm), lmm(numlm)
    real(dp), intent( in ) :: seanfq( ng, nproj )
    complex(dp), intent( in ) :: ck( ng, nbd )
    complex(dp), intent( in ) :: tauphs( ng )!, ylmfac( -l : l, ng )
    complex(dp), intent( in ) :: ylmfac( ng, numlm )
!    complex(dp), intent( out ) :: coeff( -lmax : lmax, nbd, npmax )
    complex(dp), intent( inout ) :: o2l( nbd, nptot )
    integer, intent(inout) :: ip
    !
    integer :: m, ig, ibd, iproj
    complex(dp) :: su
    !
    do ilm = 1, numlm
      if( lml( ilm ) .eq. l ) then
        goto 111
      endif
    enddo
111 continue
    ilm = ilm - 1

      do m = -l, l
        ilm = ilm + 1
        do iproj = 1, nproj
          ip = ip + 1
            do ibd = 1, nbd
             su = 0
             do ig = 1, ng
                su = su + ylmfac( ig, ilm ) * tauphs( ig ) * ck( ig, ibd ) * seanfq( ig, iproj )
             end do
!             coeff( m, ibd, iproj ) = su
             o2l( ibd, ip ) = su
            end do
          end do
       end do
    !
    return
  end subroutine getcoeff

  subroutine fullgetcoeff( ng, nbd, ntau, lmin, lmax, npmax, nproj, ck, tauphs, ylmfac, & 
                           o2l, seanfq, lml, lmm, numlm, nptot )
    use kinds, only : DP
    implicit none
    integer, intent( in ) :: ng, nbd, lmin, lmax, nproj( lmin:lmax ), npmax, lml(numlm), lmm(numlm), numlm, nptot, ntau
    real(dp), intent( in ) :: seanfq( ng, npmax, lmin : lmax )
    complex(dp), intent( in ) :: ck( ng, nbd )
    complex(DP), intent( in ) :: tauphs( ng, ntau ), ylmfac( ng, numlm )
    complex(dp), intent( out ) :: o2l( nbd, nptot, ntau )

    integer :: itau, ibd, ip, ilm, l, iproj, ig, ig_big, ig_by64, ib_big, ib_by64, ig_remain


!   L2 should be larger than tile_ng * tile_nb * sizeof(complex_dp)
    integer, parameter :: tile_ng = 128
    integer, parameter :: tile_nb = 64
    integer, parameter :: tile_ngm1 = tile_ng - 1
    integer, parameter :: tile_nbm1 = tile_nb - 1

    complex(dp) :: su, su2( tile_ng ), su3( tile_ng )



! !$OMP PARALLEL DEFAULT(NONE) &
! !$OMP& PRIVATE( ig_big, itau, ip, ilm, l, iproj, su, ig, ig_by64, ib_big, ib_by64 ) &
! !$OMP& SHARED( ylmfac, tauphs, ck, seanfq, o2l ) &
! !$OMP& FIRSTPRIVATE( ntau, numlm, lml, nproj, nbd, ng )

    
    ig_by64 = (ng/tile_ng) * tile_ng
    ib_by64 = (nbd/tile_nb) * tile_nb

    if( ig_by64 .eq. 0 ) goto 1110

    do ib_big = 1, ib_by64, tile_nb
    do ig_big = 1, ig_by64, tile_ng

! !$OMP DO SCHEDULE(STATIC)
    do itau = 1, ntau
      ip = 0
      do ilm = 1, numlm
        l = lml( ilm )
!          m = lmm( ilm )
        su2( 1 : tile_ng ) = ylmfac( ig_big : ig_big + tile_ngm1, ilm ) * tauphs( ig_big : ig_big + tile_ngm1, itau )
        do iproj = 1, nproj( l )
          ip = ip + 1
          su3( 1 : tile_ng ) = su2( 1 : tile_ng ) * seanfq( ig_big : ig_big + tile_ngm1, iproj, l )
          do ibd = ib_big, ib_big+tile_nbm1
!          do ibd = 1, nbd
            su = 0.0_dp
!            do ig = 1, ng
!            do ig = ig_big, ig_big+31
!              su = su &
!                 + ylmfac( ig, ilm ) * tauphs( ig, itau ) * ck( ig, ibd ) * seanfq( ig, iproj, l )
!            enddo
            su = sum( su3( 1 : tile_ng ) * ck( ig_big : ig_big + tile_ngm1, ibd ) )
            o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) + su
          enddo
        enddo
      enddo
    enddo
! !$OMP ENDDO NOWAIT
    enddo
    enddo

    if( ig_by64 .eq. ng ) then
      goto 1111
    endif

1110 continue

    ig_remain = ng - ig_by64

! ig remainder, still subset of ib
    do ib_big = 1, ib_by64, tile_nb
! !$OMP DO SCHEDULE(STATIC)
    do itau = 1, ntau
        ip = 0
        do ilm = 1, numlm
          l = lml( ilm )
          su2( 1 : ig_remain ) = ylmfac( ig_by64+1:ng, ilm ) * tauphs( ig_by64+1:ng, itau )
!          m = lmm( ilm )
          do iproj = 1, nproj( l )
            ip = ip + 1
            su3( 1 : ig_remain ) = su2( 1 : ig_remain ) * seanfq( ig_by64+1:ng, iproj, l )
!            do ibd = 1, nbd
          do ibd = ib_big, ib_big+tile_nbm1
!            su = 0.0_dp
!!            do ig = 1, ng
!            do ig = ig_by64 + 1, ng
!!              o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) &
!              su = su &
!          + ylmfac( ig, ilm ) * tauphs( ig, itau ) * ck( ig, ibd ) * seanfq( ig, iproj, l )
!            enddo
            su = sum( su3( 1 : ig_remain ) * ck( ig_by64+1:ng, ibd ) )
            o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) + su
          enddo
        enddo
      enddo
    enddo
! !$OMP ENDDO

    enddo

1111 continue
    if( ib_by64 .eq. nbd ) then
      goto 1112
    endif

! remainder of ib, all of ig

! !$OMP DO SCHEDULE(STATIC)
    do itau = 1, ntau
      ip = 0 
      do ilm = 1, numlm
        l = lml( ilm )
!          m = lmm( ilm )
        do iproj = 1, nproj( l )
          ip = ip + 1
          do ibd = ib_by64 + 1, nbd
            su = 0.0_dp
            do ig = 1, ng
!              o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) &
              su = su &
          + ylmfac( ig, ilm ) * tauphs( ig, itau ) * ck( ig, ibd ) * seanfq( ig, iproj, l )
            enddo
            o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) + su
          enddo
        enddo
      enddo
    enddo
! !$OMP ENDDO

1112 continue

! !$OMP END PARALLEL

!
!
  end subroutine fullgetcoeff


  subroutine make_prefs( prefs, lmax )
    use kinds, only : dp
    implicit none
    integer, intent( in ) :: lmax
    real(dp), intent( out ) :: prefs(0 : 1000 )

    integer, external :: freeunit

    integer :: nsphpt, isphpt, iuntmp
    real(dp) :: sphsu
    real(dp), allocatable, dimension( : ) :: xsph, ysph, zsph, wsph

    iuntmp = freeunit()
   
    open( unit=iuntmp, file='sphpts', form='formatted', status='old' )
    rewind iuntmp
    read ( iuntmp, * ) nsphpt
    allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
    do isphpt = 1, nsphpt
       read ( iuntmp, * ) xsph( isphpt ), ysph( isphpt ), zsph( isphpt ), wsph( isphpt )
    end do
    close( unit=iuntmp )
    sphsu = sum( wsph( : ) )
    wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / sphsu )
    call getprefs( prefs, lmax, nsphpt, wsph, xsph, ysph, zsph )    
  end subroutine make_prefs

  end subroutine OCEAN_obf2localbyk
