! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the GPL 2 License. See the file `License' in the current subdirectory.
!

!  con_coeff is coeffs on the shifted grid for all bands 1 - Nbands
!  val_coeff is coeffs on the un-shifted grid for valence only
module OCEAN_obf2loc

  use kinds, only : dp
  use io_global,  ONLY : stdout, ionode, ionode_id
  use mp, only : mp_sum, mp_bcast, mp_barrier
  
  implicit none

  private
  save

  real(dp) :: prefs( 0 : 1000 )
  real(dp) :: bvec( 3, 3 ), bmet( 3, 3 ), qbase( 3 ), kshift( 3 ), dqproj
  real(dp), allocatable :: tau(:,:), fttab( :, :, : )

  complex(dp), allocatable :: o2l( :, :, : ), con_coeff(:,:,:,:,:), val_coeff(:,:,:,:,:)

!  integer, parameter :: lmax = 5
  integer :: nbuse_xes, lmin, lmax, npmax, nqproj, nptot, ntau, nband, kpts( 3 ), nkpt, nspin, nshift, numlm
  integer, allocatable :: nproj( : ), lml( : ), lmm( : )
  
  logical :: have_kshift

  character(len=6), allocatable :: fntau(:)


  public :: OCEAN_obf2loc_init, OCEAN_obf2loc_alloc, OCEAN_obf2loc_coeffs, OCEAN_obf2loc_sumO2L
  public :: OCEAN_obf2loc_write, OCEAN_obf2loc_reorder

  contains

  subroutine OCEAN_obf2loc_init( nband_, kpts_, nspin_, nshift_, ierr )

    use OCEAN_timer
!    use mp, only : mp_bcast, mp_barrier
!    USE io_global,  ONLY : stdout, ionode, ionode_id

    implicit none

    integer, intent( in ) :: nband_, kpts_(3), nspin_, nshift_
    integer, intent( inout ) :: ierr
!    logical, intent( in ) :: ionode

    integer,external :: freeunit

    integer :: i, j, itau, l, k, ii, jj, kk
    integer :: iuntmp

    real(dp) :: avec( 3, 3 )

    integer :: atno, nc, lc, indx, nprj

    character(25) :: element, fnroot, add04, add10

    real(dp) :: pi

    pi = 4.0d0 * atan( 1.0d0 )



    nband = nband_
    kpts(:) = kpts_(:)
    nkpt = product( kpts(:) )
    nspin = nspin_
    nshift = nshift_

    ! ======================================================================
    ! Set up the Ylm prefactors; opening files so do it only on ionode
    ! ======================================================================
    if( ionode ) then
      lmax = 5
      call Amake_prefs

    ! ======================================================================
    ! This is a bunch of house keeping crap that need to be integrated in to 
    ! the input parser at some point.
    !  ======================================================================

      iuntmp = freeunit()

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


!      open(unit=iuntmp,file='nbuse.ipt',form='formatted',status='old')
!      read(iuntmp,*) nbuse
!      close(iuntmp)

      open(unit=iuntmp,file='nbuse_xes.ipt',form='formatted',status='old')
      read(iuntmp,*) nbuse_xes
      close(iuntmp)

    endif

    write(stdout,*) ' Read in little data files'
!    call MPI_BARRIER( intra_pool_comm, ierr )
    call mp_barrier
    write(stdout,*) ' Preparing to share little data files'

    call mp_bcast( prefs, ionode_id )
!    call mp_bcast( nbuse, ionode_id )
    call mp_bcast( nbuse_xes, ionode_id )


  !  mp_bcast( nspn, ionode_id )
    call mp_bcast( qbase, ionode_id )
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

    call mp_bcast( bvec, ionode_id )

!     call MPI_BARRIER( intra_pool_comm, ierr )
    call mp_barrier
    write(stdout,*) ' Done sharing little data files'


    do i = 1, 3
      do j = 1, 3
        bmet( i, j ) =dot_product( bvec( :, i ), bvec( :, j ) )
      end do
    end do

  end subroutine

  



  subroutine OCEAN_obf2loc_alloc( ierr )

    implicit none
    
    integer, intent( inout ) :: ierr
!    logical, intent( in ) :: ionode

    allocate( o2l( nband, nptot, ntau ) )
    if( ionode ) then
      allocate( con_coeff( nptot, nband, nkpt, nspin, ntau ), &
                val_coeff( nptot, nbuse_xes, nkpt, nspin, ntau ) )
    else
      allocate( con_coeff(1,1,1,1,1), val_coeff(1,1,1,1,1) )
    endif

  end subroutine OCEAN_obf2loc_alloc

  subroutine OCEAN_obf2loc_dealloc( ierr )
    implicit none
    integer, intent( inout ) :: ierr

    deallocate( o2l, con_coeff, val_coeff )
  end subroutine OCEAN_obf2loc_dealloc



  subroutine OCEAN_obf2loc_sumO2L
!    use mpi, only : MPI_IN_PLACE, MPI_SUM, MPI_DOUBLE_COMPLEX, MPI_BARRIER, MPI_COMM_WORLD
!    use mp_global, only : intra_pool_comm, root_pool, me_pool
!    use mp, only : mp_sum
    implicit none
    integer :: nelements, ierr

    write(stdout,*) nband, nptot, ntau
!    call MPI_BARRIER( intra_pool_comm, ierr )
    nelements = nband * nptot * ntau

!    if( me_pool == root_pool ) then
!      call MPI_REDUCE( MPI_IN_PLACE, o2l, nelements, MPI_DOUBLE_COMPLEX, MPI_SUM, & 
!                       root_pool, intra_pool_comm )
!    else
!      call MPI_REDUCE( o2l, o2l, nelements, MPI_DOUBLE_COMPLEX, MPI_SUM, & 
!                       root_pool, intra_pool_comm )
!    endif
!    call MPI_ALLREDUCE( MPI_IN_PLACE, o2l, nelements, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr )
    call mp_sum( o2l )

  end subroutine OCEAN_obf2loc_sumO2L


  subroutine OCEAN_obf2loc_reorder( ispin, ishift, ikpt, ibeg )
    implicit none
!    logical, intent(in) :: ionode
    integer, intent(in) :: ispin, ishift, ikpt
    integer, intent(in) :: ibeg

    integer :: itau, j, k

    if( ionode ) then

      write(6,*) nkpt, nspin

      if( ishift .eq. 1 ) then
        do itau = 1, ntau 
          val_coeff( :, 1:ibeg-1, ikpt, ispin, itau ) = transpose( o2l( 1:ibeg-1, :, itau ) )
        enddo
      endif

      if( ishift .eq. nshift ) then
        do itau = 1, ntau
          do j = 1, nband
            do k = 1, nptot
!JTV switch to transpose call
              con_coeff( k, j, ikpt, ispin, itau ) = o2l( j, k, itau )
            enddo
          enddo
        enddo
      endif

    endif

  end subroutine OCEAN_obf2loc_reorder


  subroutine OCEAN_obf2loc_write( ierr )
    implicit none
!    logical, intent(in) :: ionode
!    integer, intent(in) :: stdout
    integer, intent(inout) :: ierr
    integer, external :: freeunit

    integer :: iuntmp, itau
    character(len=11) :: filout
    real(kind=dp), allocatable :: out_coeffs(:,:,:,:)


    iuntmp = freeunit()

    if( ionode ) then
      write(stdout,*) nptot, nband*nkpt, nspin, ntau
      allocate( out_coeffs( nptot, nband, nkpt, nspin ) )
      do itau = 1, ntau
        write( filout, '(1a5,1a6)' ) 'cksc.', fntau( itau )
        open(unit=iuntmp,file=filout,form='unformatted',status='unknown')!,buffered='yes',blocksize=1048576,buffercount=128)
        rewind(iuntmp)
        write(iuntmp) nptot, nband*nkpt, nspin
        write(iuntmp) tau( :, itau )

        ! Getting strange size-dependent crashes when attempting to combine type conversion and write
        !   There shouldn't be any time costs to two-stepping it
        out_coeffs( :,:,:,:) = real(con_coeff(:,:,:,:,itau))
        write(iuntmp) out_coeffs
        out_coeffs( :,:,:,:) = aimag(con_coeff(:,:,:,:,itau))
        write(iuntmp) out_coeffs

        close(iuntmp) 
      enddo

      deallocate( out_coeffs )
      allocate( out_coeffs( nptot, nbuse_xes, nkpt, nspin ) )

      do itau = 1, ntau
        write( filout, '(1a5,1a6)' ) 'cksv.', fntau( itau )
        open(unit=iuntmp,file=filout,form='unformatted',status='unknown')!,buffered='yes',blocksize=1048576,buffercount=128)
        rewind(iuntmp)
        write(iuntmp) nptot, nbuse_xes*nkpt, nspin
        write(iuntmp) tau( :, itau )

        ! Getting strange size-dependent crashes when attempting to combine type conversion and write
        !   There shouldn't be any time costs to two-stepping it
        out_coeffs( :,:,:,:) = real(val_coeff(:,:,:,:,itau))
        write(iuntmp) out_coeffs
        out_coeffs( :,:,:,:) = aimag(val_coeff(:,:,:,:,itau))
        write(iuntmp) out_coeffs

        close(iuntmp)
      enddo

      deallocate( out_coeffs )

    endif

  end subroutine OCEAN_obf2loc_write
!
! ng = number of Gvecs for the writing of each wavevector
! kvc = inverse of the g list
! bmet 
! bvec
! ck = complex wavefunction
! ibl = start band
! ibh = stop band
  subroutine OCEAN_obf2loc_coeffs( ng, kvc, ck, iq, ispin, ishift, ibeg, loud )
!    use constants, ONLY : tpi
    use OCEAN_timer
!    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, iq, ispin
    integer, intent( in ) :: kvc( 3, ng ), ibeg(nkpt,nspin)
    complex(dp), intent( in ) :: ck( ng, nband )
    integer, intent( inout ) :: ishift
    logical, intent( in ) :: loud
    
    complex(dp), parameter :: zero = 0.0_dp
    !
    integer itau, l, ig, ik1, ik2, ik3, i, j, k
    real(dp) :: qphase, phase
    real(dp) :: q( 3 )
    complex(dp), allocatable :: tauphs( : , : )
    complex(dp), allocatable :: ylmfac( : , : )
    real(dp), allocatable :: sfq( : , : , : )
    complex(dp) :: rm1
    real(dp) :: pi 
    !
    
    ! determine q(3)
!    ik1 = iq / (nkpt(2) * nkpt(3))
!    ik2 = 
!    ik3 = mod( iq, nkpt(3) )
    i = 0
    do ik1 = 0, kpts(1)-1
      do ik2 = 0, kpts(2)-1
        do ik3 = 0, kpts(3)-1
          i = i + 1
          if( i .eq. iq ) then
            q( 1 ) = ( qbase( 1 ) + dble( ik1 ) ) / dble( kpts( 1 ) )
            q( 2 ) = ( qbase( 2 ) + dble( ik2 ) ) / dble( kpts( 2 ) )
            q( 3 ) = ( qbase( 3 ) + dble( ik3 ) ) / dble( kpts( 3 ) )
            goto 11
          endif
        enddo
      enddo
    enddo

    

11 continue
    

    do ishift = 1, nshift

    if( loud ) write(stdout,*) q(:)
    o2l = zero

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
        phase = qphase + 2.0d0 * pi * sum( tau( :, itau ) * real( kvc( :, ig ), dp ) )
        tauphs( ig, itau ) = cos( phase ) + rm1 * sin( phase ) !cmplx( cos(phase), sin(phase) )
      enddo
    end do
    if( loud ) call OCEAN_t_printtime( "Phase", stdout )
    call OCEAN_t_reset

    do l = lmin, lmax
       call Aseanitup( ng, q, kvc, l, sfq( 1, 1, l ) )
    end do
    call Anewgetylmfac( ng, kvc, q, ylmfac )
    if( loud ) call OCEAN_t_printtime( "Seanitup", stdout )
    call OCEAN_t_reset

    call Afullgetcoeff( ng, ck, tauphs, ylmfac, sfq )

    if( loud ) call OCEAN_t_printtime( "Get Coeff", stdout )
    !
    deallocate( tauphs, ylmfac, sfq )
    !

!    call OCEAN_t_reset
!    call OCEAN_obf2loc_sumO2L()
!    call OCEAN_t_printtime( "Sum O2L", stdout )
!
!    call OCEAN_t_reset
!    call OCEAN_obf2loc_reorder( ispin, ishift, iq, ibeg(iq,ispin) )
!    call OCEAN_t_printtime( "Reorder O2L", stdout )
      call mp_sum( o2l )
        if( ionode ) then
          if( ishift .eq. 1 ) then
            do itau = 1, ntau
              val_coeff( :, 1:ibeg(i,ispin)-1, i, ispin, itau ) = transpose( o2l( 1:ibeg(i,ispin)-1, :, itau ) )
            enddo
          endif

          if( ishift .eq. nshift ) then
            do itau = 1, ntau
              do j = 1, nband
                do k = 1, nptot
                  con_coeff( k, j, i, ispin, itau ) = o2l( j, k, itau )
!                con_coeff( 1:nptot, 1:nbuse, i, ispin, itau ) = &
!                           transpose( o2l( ibeg(i,ispin):ibeg(i,ispin)+nbuse-1, 1:nptot, itau ) )
                enddo
              enddo
            enddo
          endif

        endif

      q(:) = q(:) + kshift(:)
    enddo

    return

  end subroutine OCEAN_obf2loc_coeffs


  subroutine Aseanitup( ng, q, gvec, l, seanfq )

!    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, l
    integer, intent( in ) :: gvec( 3, ng )
    !
    real(dp), intent( in ) :: q( 3 )
    real(dp), intent( inout ) :: seanfq( ng, nproj(l) )
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
       gc( : ) = q( : ) + real(gvec( :, j ), dp)
       qbar = 0.d0
       do i1 = 1, 3
          do i2 = 1, 3
             qbar = qbar + gc( i1 ) * gc( i2 ) * bmet( i1, i2 )
          end do
       end do
       qbar = sqrt( qbar )
       ii = 1 + qbar / dqproj
       if ( ii .lt. 2 ) ii = 2
       if ( ii .gt. nqproj - 2 ) ii = nqproj - 2
       qii = dqproj * ( ii - 1 )
       x = ( qbar - qii ) / dqproj
       wm1 = pm1 * x * ( x - 1.d0 ) * ( x - 2.d0 )
       wz0 = pz0 * ( x + 1.d0 ) * ( x - 1.d0 ) * ( x - 2.d0 ) 
       wp1 = pp1 * ( x + 1.d0 ) * x * ( x - 2.d0 )
       wp2 = pp2 * ( x + 1.d0 ) * x * ( x - 1.d0 )
       do i = 1, nproj(l)
          seanfq( j, i ) = fttab( ii - 1, i, l ) * wm1 + &
               &           fttab( ii, i, l ) * wz0 + &
               &           fttab( ii + 1, i, l ) * wp1 + &
               &           fttab( ii + 2, i, l ) * wp2
       end do
    end do
    !
    return
  end subroutine Aseanitup

  subroutine getylmfac( ng, gvec, q, l, ylmfac )
!    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, l
    integer, intent( in ) :: gvec( 3, ng )
    real(dp), intent( in ) :: q( 3 )
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

  subroutine Anewgetylmfac(ng, gvec, q, ylmfac )
!    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng
    integer, intent( in ) :: gvec( 3, ng )
    real(dp), intent( in ) :: q( 3 )
    complex(dp), intent( out ) :: ylmfac( ng, numlm )
    !
    integer :: jj, ig, m, l, ilm
    real(dp) :: x( 3 ), pi
    complex(dp) :: pref, ylm, rm1
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    !
    do ilm = 1, numlm
      l = lml( ilm )
      m = lmm( ilm )
      pref = 4.0d0 * pi * rm1 ** l
      do ig = 1, ng
        x = 0
        do jj = 1, 3
          x( : ) = x( : ) + bvec( :, jj ) * ( q( jj ) + real(gvec( jj, ig ), dp ) )
        enddo
        call getylm( l, m, x( 1 ), x( 2 ), x( 3 ), ylm, prefs )
        ylmfac( ig, ilm ) = pref * conjg( ylm )
      enddo
    enddo


  end subroutine Anewgetylmfac


  subroutine Afullgetcoeff( ng, ck, tauphs, ylmfac, seanfq  )
!    use kinds, only : DP
    implicit none
    integer, intent( in ) :: ng
    real(dp), intent( in ) :: seanfq( ng, npmax, lmin : lmax )
    complex(dp), intent( in ) :: ck( ng, nband )
    complex(DP), intent( in ) :: tauphs( ng, ntau ), ylmfac( ng, numlm )

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
    ib_by64 = (nband/tile_nb) * tile_nb

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
!          do ibd = 1, nband
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
!            do ibd = 1, nband
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
    if( ib_by64 .eq. nband ) then
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
          do ibd = ib_by64 + 1, nband
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
  end subroutine Afullgetcoeff


  subroutine Amake_prefs
!    use kinds, only : dp
    implicit none

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
  end subroutine Amake_prefs


end module OCEAN_obf2loc
