subroutine OCEAN_bofr_multi( )
! the idea here is to project the basis functions onti our wonky OCEAN chi real-space grid
!    The grid is defined by the file rbfile.bin
  use constants, only : tpi
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE cell_base
  USE gvect
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  USE io_files, ONLY: prefix, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY: evc
  USE mp, ONLY : mp_sum, mp_bcast, mp_root_sum, mp_barrier
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset
  use OCEAN_timer

  implicit none

  real(dp),parameter ::  eps=1.d-10
  
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)


  integer :: ibnd  
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  
  integer :: i, ipt, ig
  integer,allocatable :: igk_l2g(:,:)


  integer :: iuntmp, iunbofr
  integer :: npt, ntau, ntau_offset, itau
  real(dp),allocatable :: posn( :, :, : ), drel( :, : ), wpt( :, : )

  complex(dp),allocatable :: expiGr( :, : ), bofr( :, : ), bofr_out( : )
  real(dp) :: gcart( 3 ), bvec(3,3), phase

  integer,external :: freeunit

  write(stdout,*) 'BOFR_MULTI'
  call OCEAN_t_reset


  WRITE( stdout, '(/5x,"Calling realspacebasis .... ",/)')
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

  if( ionode ) then
    write(stdout,*) ' generating plane-wave matrix elements in this basis'
    write(stdout,*) "   < B_i | exp(-iG.r) | B_j > exp(iG'.r')"
  endif


  call flush_unit( stdout )

  call summary
  !
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  call n_plane_waves (ecutwfc, tpiba2, nkstot, xk, g, ngm, npwx, ngk)

  write(stdout,*) ' ecutwfc = ', ecutwfc
  write(stdout,*) ' tpiba2 = ', tpiba2
  write(stdout,*) ' nks, nktot = ', nks, nkstot
  write(stdout,*) ' xk = ', xk(1:3,1:nkstot)
  write(stdout,*) '     npw = ', ngk(1:nks)
  write(stdout,*) '    npwx = ', npwx


!  ! gamma-point only
!  qvec = 0.d0
!  current_k = 1
!  if( lsda ) current_spin = isk(1)

  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )

  ! I no longer have access to igk_l2g - need to regenerate
  if( allocated(igk_l2g) ) deallocate( igk_l2g )
  allocate( igk_l2g( npwx, nks ) )
  igk_l2g = 0

  call gk_l2gmap (ngm, ig_l2g(1), npw, igk, igk_l2g(1,1))
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! report norms
  allocate( norm(nbnd) )
  do ibnd=1,nbnd
    norm(ibnd) = dot_product( evc(:,ibnd_indx(ibnd)), evc(:,ibnd_indx(ibnd)) )
  enddo
  call mp_sum( norm, intra_pool_comm )
  do ibnd=1,nbnd
    if( abs(norm(ibnd)-1.d0) > eps ) then
      write(stdout,'(a,i6,a,f14.10)') ' band ', ibnd, ' norm = ', norm(ibnd)
      call errore('hamq','wave function norm is not 1',1)
    endif
  enddo
  deallocate( norm )


  write(stdout,*) ' Loading the real-space points'
  ! Read in real-space grid info
  if( ionode ) then
    iuntmp = freeunit()
    open(iuntmp,file='bofr_tau_control',form='formatted',status='old')
    read(iuntmp,*) ntau, ntau_offset
    close(iuntmp)
    open(iuntmp,file='rbfile.bin',form='unformatted',status='old')
    rewind(iuntmp)
    read(iuntmp) npt
    allocate( posn(3,npt,ntau), wpt(npt,ntau), drel(npt,ntau) )
    do itau = 1, ntau
      write(stdout,*) itau
      read(iuntmp) posn(:,:,itau)
      read(iuntmp) wpt(:,itau)
      read(iuntmp) drel(:,itau)
    enddo
    close(iuntmp)
  endif
  ! share grid info
  write(stdout,*) ' Sharing the real-space points'
  call mp_bcast( ntau, ionode_id )
  call mp_bcast( npt, ionode_id )
  if( .not. ionode ) allocate( posn(3,npt,ntau), wpt(npt,ntau), drel(npt,ntau) )
  call mp_bcast( posn, ionode_id )
  call mp_bcast( wpt, ionode_id )
  call mp_bcast( drel, ionode_id )


  ! Prep output file
  if( ionode ) then
    iunbofr = freeunit()
    open(iunbofr,file=trim(prefix)//'.bofr',form='unformatted',buffered='yes')
    write(stdout,*) ' will be saved to file: ', trim(prefix)//'.bofr' 

    open(unit=99,file='bvecs',form='formatted',status='old')
    read(99,*) bvec(:,:)
    close(99)

  endif
  call mp_bcast( bvec, ionode_id )



!  ! find total number of G vectors over processes
!  igwx = maxval( igk_l2g(1:npw,1) )
!  call mp_max( igwx )
!  allocate( psid(nrxxs), wtmp(igwx), jtmp(npw) )


  !!!!!!
  allocate( expiGr( npw, npt ), bofr( npt, nbnd ), bofr_out( npt ) )


  call mp_barrier
  call OCEAN_t_printtime( "Init", stdout )

  do itau = 1, ntau
    write(stdout,*) 'Site: ', itau, ntau

    
  !  expiGr = 0.d0
  !  bofr = 0.d0
    ! precalculate exp( G . r ) for all of the local gvecs
  !JTV at the moment I don't know what units these guys are in
    do ipt=1,npt
      do ig=1,npw
  !JTV indexing of g is suspect
        gcart( : ) = 0.d0
        do i = 1, 3
          gcart( : ) = gcart( : ) + bvec( :, i ) * mill( i, ig )
        enddo
        phase = dot_product( gcart( : ), posn( :, ipt, itau ) )
  !JTV need to overload sincos from C
        expiGr( ig, ipt ) = cmplx( cos(phase), sin(phase) )
  !      expiGr( ig, ipt ) = exp( iota * cmplx( dot_product( gcart( : ), posn( :, ipt ) ) ) )
  !      expiGr( ig, ipt ) = exp( iota * tpiba * dot_product( g( :, igk( ig ) ), posn( :, ipt ) ) )
      enddo
    enddo
    write(stdout,*) mill( 1:3, 1:2 )
    write(stdout,*) nbnd, npw, npt
  !  write(stdout,*) expiGr( 1:3, 1 )

    call OCEAN_t_reset

    if( .true. ) then
    do ibnd=1,nbnd
      bofr( :,1 ) = 0.0d0
      do ipt = 1, npt
        bofr( ipt,1 ) = sum( evc( :, ibnd_indx(ibnd) ) * expiGr( :, ipt ) )
  !      bofr( ipt ) = dot_product( evc( :, ibnd_indx(ibnd) ), expiGr( :, ipt ) )
      enddo
  !JTV should only end up on ionode
      call mp_sum( bofr(:,1), intra_pool_comm )
  !    call mp_root_sum( bofr, bofr_out, ionode_id )
      if( ionode ) write(iunbofr) bofr(:,1)
!      if( mod( ibnd, 20 ) .eq. 0 ) &
!      write(stdout,*) ibnd, nbnd
    enddo

    else
    
    call ZGEMM( 'T', 'N', npt, nbnd, npw, one, evc, npw, expiGr, npw, zero, bofr, npt )
    do ibnd = 1, nbnd
      call mp_sum( bofr(:,ibnd), intra_pool_comm)
      if( ionode ) write(iunbofr) bofr(:,ibnd)
    enddo


    endif
    call OCEAN_t_printtime( "matmul", stdout )
  enddo

!  write(stdout,*) '======================================'
!  write(stdout,*) bofr( : )

  if( ionode ) close(iunbofr)

  deallocate( expiGr, bofr, bofr_out )
  deallocate( posn, wpt, drel )
  deallocate( igk_l2g, ibnd_indx )


  return

  end subroutine OCEAN_bofr_multi
