  program OCEAN_sphere_basis

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Nov 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
  use parallel_include
  use hamq_shirley
  use scalapack_module, only : CTXT_, N_, NB_, DLEN_, scalapack_blocksize
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, root
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use mpio
  use kpt_module
  use corerepair_module
  use shirley_input_module 
  use diag_module
!  USE wavefunctions_module, ONLY: evc
  USE io_files, ONLY: nwordwfc, iunwfc, prefix_io=>prefix, &
                                  tmpdir_io=>tmp_dir, nd_nmbr, diropn
  USE control_flags,        ONLY : lscf
  USE basis,                ONLY : starting_pot, starting_wfc
  USE gvect
  USE wvfct, only : npwx, npw
  use constants, only : pi
  use hamq_pool, only : nproc_per_pool, npool, &
                        mypool, rootpool, mypoolid, mypoolroot, &
                        cross_pool_comm, intra_pool_comm, &
                        desc_cyclic, desc_striped, context_cyclic, &
                        local_striped_dim, cyclic_localindex
!                        context_global
!  use hamq_pool, only : nbasis_global

  implicit none

  integer,external :: freeunit

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)
  complex(dp),parameter :: iota=(0.d0,1.d0)
  real(dp), parameter :: eps = 1.0d-14
  real(dp), parameter :: real_one = 1_dp
  real(dp), parameter :: real_zero = 0_dp

  character(255) :: info_file, string

  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr, errorcode, fmode, finfo, resultlen


  character(len=6) :: nd_nmbr_tmp


  integer :: ik, jk, lrindx, lcindx, rsrc, csrc
  integer :: i
  integer :: nbasis_subset

  integer :: nk, nbnd, ispin, itau, icol
  real(dp) :: nelec_, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin
  real(dp) :: fermi_energy
  complex(dp),allocatable :: ztmp(:,:)

  integer :: iuntmp, ntot,nptot, ibd, iunrbf, bi, bj
  real(dp),allocatable :: tau(:,:)

  complex(dp),allocatable :: bofr( :, : ), eikr( : ), single_bofr(:,:)
  complex(dp),allocatable :: phased_bofr(:,:), uofrandb(:,:), &
                             local_uofrandb(:,:), uofrandb2( :, : )
  real(dp),allocatable :: posn( :, : ), wpt( : ), drel( : ), wpt_in(:), &
         xirow(:), vind(:),vipt(:),phase(:)
  real(dp),allocatable :: wgt(:)
  real(dp) :: pref, spinfac, su, omega, avec(3,3), bvec(3,3), qin(3),qcart(3)
  real(dp) :: mindif

  complex(DP), allocatable :: S_l(:,:), B_l(:,:), work(:)
  real(DP), allocatable :: eigU(:), rwork(:)
  real(DP) :: pctrc, trace, trace_tol
  complex(DP), external :: PZLATRA
  integer :: big_nbasis, bnr_l, bnc_l, nblock,  lwork, lrwork, liwork, nbasis_trunc
  integer, allocatable :: irev(:), iwork(:)


  integer :: ntau, ntau_offset, fh_rad
  integer :: nb, nb2, local_npt, local_npt2
  integer :: nprow, npcol, myrow, mycol
  integer, dimension( DLEN_ ) :: desc_bofr_in, desc_bofr, &
                                 desc_local_uofrandb, desc_S, desc_wpt, desc_wpt_in

  integer :: npt, iunbofr, nbnd_small

  integer, external :: numroc
  integer :: iuninf
  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u


! Read in the input file, shirley_input sets up what 
!  is needed for the interpolation
  call shirley_input
    write(6,*) 'Here'

! Still need to read in the actual OBFs. We will be 
!   calculating overlaps and matrix elements and stuff
  prefix_io = prefix
  tmpdir_io = outdir
  lscf = .false.
  starting_pot = 'file'
  starting_wfc = 'file'
  call dump_system( nelec_, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
  write(stdout,*) ' tpiba = ', tpiba
  write(stdout,*) ' nspin = ', nspin
  write(stdout,*) ' lda_plus_u = ', lda_plus_u

  if( nspin .ne. 1 ) then
    write(stdout,*) 'nspin: ', nspin
    write(stdout,*) 'Spin ne 1 is not yet implemented. Trying to quit ... '
    goto 111
  endif
    write(6,*) 'Here'



  if( band_subset(1) < 1 ) band_subset(1)=1
!  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis_global
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

!  call diagx_init( band_subset(1), band_subset(2) )
  call diag_init

! moved earlier
!  call dump_system( nelec_, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
!  write(stdout,*) ' nspin = ', nspin
!  write(stdout,*) ' lda_plus_u = ', lda_plus_u

  ! band structure is given in units of tpiba
  if( trim(kpt%param%grid_type) == 'bandstructure' ) then
    kpt%list%kvec = kpt%list%kvec*tpiba
    kpt%bandstructure%kpathlen = kpt%bandstructure%kpathlen*tpiba
    kpt%bandstructure%kpathlensp = kpt%bandstructure%kpathlensp*tpiba
    kpt%bandstructure%xksp = kpt%bandstructure%xksp*tpiba
  endif
 
  ! store total number of k-points
  nk=kpt%list%nk

  if( mpime==root ) then
    info_file=trim(outfile)//'.info'
    nbnd=nbasis_subset  ! important difference

    fermi_energy = efermi

    write(stdout,*) ' Fermi energy = ', fermi_energy

    iuninf=freeunit()
    open(iuninf,file=info_file,form='formatted')
    write(iuninf,nml=info)
    write(iuninf,*) kpt%param%cartesian
    write(iuninf,*) trim(kpt%param%grid_type)
    write(iuninf,*) kpt%list%wk
    write(iuninf,*) kpt%list%kvec
    if( trim(kpt%param%grid_type) == 'tetrahedra' ) then
      write(iuninf,*) kpt%tetra%ntetra
      write(iuninf,*) kpt%tetra%tetra
    else if( trim(kpt%param%grid_type) == 'bandstructure' ) then
      write(iuninf,*) kpt%bandstructure%nksp
      do i=1,kpt%bandstructure%nksp
        write(iuninf,*) kpt%bandstructure%kpathlensp(i), &
                        trim(kpt%bandstructure%labelsp(i))
      enddo
      write(iuninf,*) kpt%bandstructure%kpathlen
    endif
    close(iuninf)

    write(stdout,*) info_file
    
  endif



  if( ionode ) then
    iuntmp=freeunit()
    
!    open(unit=iuntmp,file='bofr_tau_control',form='formatted',status='old')
!    read(iuntmp,*) ntau, ntau_offset
!    close(iuntmp)
    ntau = 1
    ntau_offset = 0
    
    iunrbf=97 !freeunit()
    open(unit=iunrbf,file='rbfile.bin',form='unformatted',status='old')
    read(iunrbf) npt
    write(stdout,*) npt, ntau
    
    allocate( posn(3,npt), wpt_in(npt), drel(npt) )
    allocate(vipt(npt))

!    do itau = 1, ntau
!      read(iuntmp) posn(:,:,itau)
!      read(iuntmp) wpt(:,itau)
!      read(iuntmp) drel(:,itau)
!      call mkvipt( npt, drel(:,itau), vipt(:,itau) )
!    enddo
!    close(iuntmp)



    open( unit=iuntmp, file='avecsinbohr.ipt',form='formatted',status='old')
    read(iuntmp,*) avec(:,:)
    close(iuntmp)
    call getomega( avec,omega)

    open( unit=iuntmp, file='bvecs',form='formatted',status='old')
    read(iuntmp,*) bvec(:,:)
    close(iuntmp)

    open( unit=iuntmp,file='efermiinrydberg.ipt',form='formatted',status='old')
    read( iuntmp, * ) fermi_energy
    close( iuntmp )

  
    write(stdout,*) npt, ntau
      
  endif
  ! share grid info
  call mp_bcast( npt, ionode_id )

  call mp_bcast( ntau, ionode_id )


!  call mp_bcast( nt, ionode_id )
!  if( .not. ionode ) allocate( t( nt ) )
!  call mp_bcast( t, ionode_id )
!  write(stdout,*) ' nt: ', nt
  call mp_bcast( bvec, ionode_id )
  call mp_bcast( omega, ionode_id )

  ! Prep bofr
  call BLACS_GRIDINFO( context_cyclic, nprow, npcol, myrow, mycol )
  write(stdout,*) nprow, npcol
  nb = 900 / nprow
  nb = min( nb, 32 )
  if( nb .lt. 1 ) nb = 1
  local_npt = numroc( npt, nb, myrow, 0, nprow )
  

  call mp_barrier
  write(stdout,*) npt, nbasis, nb, local_npt, desc_cyclic(NB_)

  allocate( wpt( local_npt ) )

  call descinit( desc_bofr, npt, nbasis, nb, desc_cyclic(NB_), 0, 0, context_cyclic, &
                 local_npt, ierr )
  if( ierr .ne. 0 ) then
    write(1000+mpime,*) 'desc_bofr'
    write(1000+mpime,*) ierr
    write(1000+mpime,*) desc_bofr(:)
    write(1000+mpime,*) npt, nbasis, nb, desc_cyclic(NB_), context_cyclic, local_npt
    flush(1000+mpime)
    stop
  endif

  call mp_barrier
  write(stdout,*) 'desc_bofr complete'
  call descinit( desc_bofr_in, npt, nbasis, npt, nbasis, 0, 0, context_cyclic, npt, ierr )
  if( ierr .ne. 0 ) stop
  call mp_barrier
  write(stdout,*) 'desc_bofr_in complete'

  call descinit( desc_local_uofrandb, npt, nbasis_subset, npt, nbasis_subset, 0, 0, context_cyclic, npt, ierr )



  ! open bofr file 
  call mp_barrier
  write(stdout,*) 'Read in bofr'
  if( ionode ) then
    iunbofr = 98 !freeunit()
    open(iunbofr,file=trim(prefix)//'.bofr',form='unformatted')
    write(stdout,*) ' will read from file: ', trim(prefix)//'.bofr'
    write(stdout,*) ' npt, nbasis: ', npt, nbasis
  endif
  if( mypoolid .eq. 0 ) then
    allocate( single_bofr( npt, nbasis ), phased_bofr( npt, nbasis ) )
    allocate( local_uofrandb( npt, nbasis_subset ) )
  else
    allocate( single_bofr( 1, 1 ), phased_bofr( 1, 1 ) )
    allocate( local_uofrandb( 1, 1 ) )
  endif
    

  if( mypoolid == mypoolroot .and. .false.) then
    call MPI_INFO_CREATE( finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    call MPI_FILE_OPEN( cross_pool_comm, 'radial.dat', fmode, MPI_INFO_NULL, fh_rad, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    !JTV At this point it would be good to create a custom MPI_DATATYPE
    !  so that we can get optimized file writing
    call MPI_FILE_SET_VIEW( fh_rad, offset, MPI_DOUBLE_COMPLEX, &
                          MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

  endif


!  if( ( .not. ionode ) .and. ( mypoolid .eq. 0 ) ) allocate( posn(3,npt), wpt_in( npt ) )
  if( .not. ionode ) then
    if( mypoolid .eq. 0 ) then
      allocate( posn(3,npt), wpt_in( npt ) )
    else
      allocate( posn(1,1), wpt_in(npt) )
    endif
  endif
  

  nb2 = 900 /npcol
  nb2 = min( nb2, 32 )
  if( nb2 .lt. 1 ) nb2 = 1
  local_npt2 =  numroc( npt, nb2, mycol, 0, npcol )


  call mp_bcast( omega, ionode_id )
  pref = 1.d0 / ( kpt%list%nk * omega )
  write(stdout,*) ' pref: ', pref, kpt%list%nk, omega


  if( nspin .eq. 2 ) then
    spinfac = 1.d0
  else
    spinfac = 2.d0 
  endif

  allocate( phase( npt ), eikr( npt ) )

! =====================================================================
  if( npool /= 1 ) stop
  big_nbasis = nbasis_subset * kpt%list%nk

  call scalapack_blocksize( nblock, 32, big_nbasis, big_nbasis, nprow, npcol )
  ! padding constants
  bnr_l = numroc( big_nbasis, nblock, myrow, 0, nprow )
  bnc_l = numroc( big_nbasis, nblock, mycol, 0, npcol )

  call descinit( desc_S, big_nbasis, big_nbasis, nblock, nblock, 0, 0, context_cyclic, max(1,bnr_l), ierr )

  allocate( S_l( bnr_l, bnc_l ), B_l( bnr_l, bnc_l ), eigU( big_nbasis ) )
  lwork=-1
  lrwork = -1
  liwork = -1
  allocate( work(1), rwork(1), iwork(1) )
  call pzheevd( 'V', 'U', big_nbasis, S_l, 1, 1, desc_S, eigU, B_l, 1, 1, desc_S, &
                work, lwork, rwork, lrwork, iwork, liwork, ierr )
  lwork = work(1)
  lrwork = rwork(1) 
  liwork = iwork(1) 
  deallocate( work, rwork, iwork )
  allocate(  work(lwork), rwork(lrwork), iwork(liwork), stat=ierr )
  call errore('scalapack_diag','unable to allocate workspace',abs(ierr))

  
  trace_tol = 1.0d-10

  do itau = 1, ntau
    write(stdout,*) 'Now treating core :', itau, ' of ', ntau
    ! read in posn
    if( ionode ) then
      read(iunrbf) posn
      read(iunrbf) wpt_in
      read(iunrbf) drel
    endif
    if( mypoolid .eq. 0 ) then
      call mp_bcast( posn, ionode_id, cross_pool_comm )
!      call mp_bcast( wpt_in, ionode_id, cross_pool_comm )
    endif
    if( ionode ) call mkvipt( npt, drel, vipt )

    ! read in bofr
    if( ionode ) then
      do i = 1, nbasis
        read(iunbofr) single_bofr(:,i)
      enddo
    endif

    if( mypoolid .eq. 0 ) then
      call mp_bcast( single_bofr, ionode_id, cross_pool_comm )
    endif


    ! get wpt on all procs
    call mp_bcast( wpt_in, ionode_id )
    call descinit( desc_wpt, npt, 1, nb, 1, 0, 0, context_cyclic, local_npt, ierr )

    ! get wpt into the same layout as the other matrices, spec. uofrandb
    do i = 1, npt
      call INFOG2L( i, 1, desc_wpt, nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc )
      if( myrow .eq. rsrc )  wpt( lrindx ) = wpt_in( i )
    enddo
!    call descinit( desc_wpt_in, npt, 1, npt, 1, 0, 0, context_cyclic, npt, ierr )
!    do icol = 0, npcol - 1
!      call descinit( desc_wpt, npt, 1, nb, 1, 0, icol, context_cyclic, local_npt, ierr )
!      call PZGEMR2D( npt, 1, wpt_in, 1, 1, desc_wpt_in, wpt, 1, 1, desc_wpt, context_cyclic )
!    enddo

    
! ======================================================================
  do ispin=1,nspin

    do ik=1,kpt%list%nk
! ======================================================================

  !    if( mod(ik-1,npool)/=mypool ) cycle

      write(stdout,'(a,3f12.5,i6,a,i6,a,i3)') ' k-point ', &
        kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                    ' on node ', mpime

      ! build the Hamiltonian for this q-point
      call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin )

      if( kinetic_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_kink( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else if( local_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_vlock( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else if( nonlocal_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_vnlk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else if( smatrix_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_sk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else

  !    call diagx_ham
        call diag_ham

      endif

      write(stdout,*) ik, eigval(band_subset(1))*rytoev, eigval(band_subset(2))*rytoev
  !    write(stdout,*) kpt%list%kvec( 1 : 3, ik )
  !    write(stdout,*) posn( 1 : 3, 1 )
      qin( : ) =  kpt%list%kvec( 1 : 3, ik ) 
      qcart(:) = 0.d0
      if( kpt%param%cartesian ) then
        qcart(:) = qin(:)*tpiba
      else
        do i=1,3
          qcart(:) = qcart(:) + bvec(:,i)*qin(i)
        enddo
      endif
      write(stdout,*) qin(:)
      write(stdout,*) qcart(:)
      write(stdout,*) kpt%list%kvec( 1 : 3, ik ) * tpiba


      allocate( bofr( local_npt, desc_cyclic(N_) ) )
      allocate( uofrandb( local_npt, desc_cyclic(N_) ) )
      allocate( uofrandb2( local_npt, desc_cyclic(N_) ) )



! This is probably sub-optimal, but requires the least thinking, I think.
! Delay bofr distribution
! Add in phase for the k-point we are at
      if( mypoolid .eq. 0 ) then
        call DGEMV( 'T', 3, npt, real_one, posn(1,1), 3, qcart, 1, real_zero, phase, 1 )
        eikr( : ) = exp( iota * phase( : ) )
        do ibd = 1, nbasis
          phased_bofr( :, ibd ) = single_bofr( :, ibd ) * eikr( : )
        enddo
      endif
! Now distribute bofr, or don't
!      write(stdout,*) 'Distribute bofr'
      call PZGEMR2D( npt, nbasis, phased_bofr, 1, 1, desc_bofr_in, bofr, 1, 1, desc_bofr, context_cyclic )
!      call mp_barrier
!      write(stdout,*) 'bofr distributed'

      ! can probably change to nbasis_subset
      call PZGEMM( 'N', 'N', npt, nbasis, nbasis, &
                   one, bofr, 1, 1, desc_bofr, eigvec, 1, 1, desc_cyclic, &
                   zero, uofrandb, 1, 1, desc_bofr )

      do ibd = 1, nbasis_subset 
        su = REAL( SUM( CONJG( uofrandb( :, ibd ) ) * wpt( : ) * uofrandb( :, ibd ) ) )
        su = sqrt( 1.0_DP / su )
        uofrandb( :, ibd ) = uofrandb( :, ibd ) * su * wpt( : )
      enddo

!===========!
      ! Now we have data for kpt i !
      ! Build up kpt j !

      do jk = ik,kpt%list%nk
        ! build the Hamiltonian for this q-point
        call diag_build_hamk( kpt%list%kvec(1:3,jk), kpt%param%cartesian, ispin )

        if( kinetic_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_kink( kpt%list%kvec(1:3,jk), kpt%param%cartesian, ztmp )
          forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
          deallocate( ztmp )
        else if( local_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_vlock( kpt%list%kvec(1:3,jk), kpt%param%cartesian, ispin, ztmp )
          forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
          deallocate( ztmp )
        else if( nonlocal_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_vnlk( kpt%list%kvec(1:3,jk), kpt%param%cartesian, ispin, ztmp )
          forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
          deallocate( ztmp )
        else if( smatrix_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_sk( kpt%list%kvec(1:3,jk), kpt%param%cartesian, ztmp )
          forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
          deallocate( ztmp )
        else

    !    call diagx_ham
          call diag_ham

        endif

        qin( : ) =  kpt%list%kvec( 1 : 3, jk )
        qcart(:) = 0.d0
        if( kpt%param%cartesian ) then
          qcart(:) = qin(:)*tpiba
        else
          do i=1,3
            qcart(:) = qcart(:) + bvec(:,i)*qin(i)
          enddo
        endif
!        write(stdout,*) qin(:)
!        write(stdout,*) qcart(:)
!        write(stdout,*) kpt%list%kvec( 1 : 3, jk ) * tpiba



  ! This is probably sub-optimal, but requires the least thinking, I think.
  ! Delay bofr distribution
  ! Add in phase for the k-point we are at
        if( mypoolid .eq. 0 ) then
          call DGEMV( 'T', 3, npt, real_one, posn(1,1), 3, qcart, 1, real_zero, phase, 1 )
          eikr( : ) = exp( iota * phase( : ) )
          do ibd = 1, nbasis
            phased_bofr( :, ibd ) = single_bofr( :, ibd ) * eikr( : )
          enddo
        endif
  ! Now distribute bofr, or don't
  !      write(stdout,*) 'Distribute bofr'
        call PZGEMR2D( npt, nbasis, phased_bofr, 1, 1, desc_bofr_in, bofr, 1, 1, desc_bofr, context_cyclic )
  !      call mp_barrier
  !      write(stdout,*) 'bofr distributed'

        call PZGEMM( 'N', 'N', npt, nbasis, nbasis, &
                     one, bofr, 1, 1, desc_bofr, eigvec, 1, 1, desc_cyclic, &
                     zero, uofrandb2, 1, 1, desc_bofr )

        do ibd = 1, nbasis_subset
          su = REAL( SUM( CONJG( uofrandb2( :, ibd ) ) * wpt( : ) * uofrandb2( :, ibd ) ) )
          su = sqrt( 1.0_DP / su )
          uofrandb2( :, ibd ) = uofrandb2( :, ibd ) * su 
        enddo


        ! now overlap the two

!        do ibd = 1, nbasis_subset
!          do jbd = 1, nbasis_subset 
!            S( (ik-1)*nbasis_subset + ibd, (jk-1)*nbasis_subset + jbd ) = &
!              DOT_PRODUCT( uofrandb( :, ibd ), uofrandb2( :, jbd ) )
!          enddo
!        enddo
        ! or PZGEMM
        call PZGEMM( 'C', 'N', nbasis_subset, nbasis_subset, npt, one, &
                     uofrandb, 1, 1, desc_bofr, uofrandb2, 1, 1, desc_bofr, &
                     zero, S_l, 1+(ik-1)*nbasis_subset, 1+(jk-1)*nbasis_subset, desc_S )
        

      enddo




!      call mp_barrier
!      write(stdout,*) 'uofrandb constructed'
!      call mp_barrier

!      ! Condense uofrandb down to io node per pool
!
!      call PZGEMR2D( npt, nbasis_subset, uofrandb, 1, band_subset(1), desc_bofr, &
!                     local_uofrandb, 1, 1, desc_local_uofrandb, context_cyclic )
!
!      if( mypoolid .eq. 0 ) then
!        offset = ((ispin-1)*kpt%list%nk + ik-1) * nbasis_subset * npt
!        call MPI_FILE_WRITE_AT( fh_rad, offset, local_uofrandb, nbasis_subset * npt, &
!                                MPI_DOUBLE_COMPLEX, status, ierr )        
!      endif

      deallocate( uofrandb, bofr, uofrandb2 )

    enddo ! ik

    write(6,*) 'Done with loop over k-points'
    write(stdout,'(a,i20)') ' construct_overlap: true trace         (nbnd*nkpt) = ', big_nbasis

    trace = real( pzlatra( big_nbasis, S_l, 1, 1, desc_S ) )
    write(*,*) ' scalapack_diag: trace (matrix estimate) = ', trace 

    ! Diagonalize
    !  call scalapack_diag( big_nbasis, S_l, eigU, B_l )
    call pzheevd( 'V', 'U', big_nbasis, S_l, 1, 1, desc_S, eigU, B_l, 1, 1, desc_S, &
                work, lwork, rwork, lrwork, iwork, liwork, ierr )

    ! check trace
    trace = sum(eigU)
    write(stdout,'(a,f15.10)') ' scalapack_diag: trace (eigval estimate) = ', trace

    

    ! Truncate
    ! redefine the eigenvalues
    eigU = - eigU**2
    write(stdout,*) ' done diagonalize overlap'

    ! report and truncate
    ! note that diagonalizing S instead of U means that the eigenvalues
    ! are in increasing order, so we need to reverse them
    allocate( irev(big_nbasis) )
    do i=1,big_nbasis
      irev(i)=big_nbasis-i+1
    enddo
    trace = sum(eigU)
    write(stdout,*)
    write(stdout,*) ' -trace(U=-S^2) = ', -trace
    pctrc = 0.d0
    write(stdout,'(2x,a10,3a19)') 'i', 'eigU(i)', 'coverage', '%cover'
    do i=1,big_nbasis
      pctrc=pctrc+eigU(irev(i))
      write(stdout,'(2x,i10,2e19.10,3f19.10)') i, eigU(irev(i)), -pctrc, pctrc/trace*100.d0
    enddo

    ! truncation
    nbasis_trunc = big_nbasis
    if( trace_tol > 0 ) then
      pctrc=0.d0
      do i=1,big_nbasis
        pctrc=pctrc+eigU(irev(i))
        if( abs(pctrc/trace) > 1.d0-abs(trace_tol) ) exit
      enddo
      nbasis_trunc = min( i, big_nbasis )
      write(stdout,*) 'truncating to ', nbasis_trunc, ' basis functions'
      write(stdout,*)
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now construct B_i(r) on this grid and dump to file for posterity

    b2nc_l = numroc( nbasis_trunc, nblock, mycol, 0, npcol )
    allocate( basis_r( local_npt, b2nc_l ) )
    call descinit( desc_B, npt, nbasis_trunk, nb, nblock, 0, 0, &
                   context_cyclic, local_npt, ierr )

    basis_r = 0.0_DP

    do ik=1,kpt%list%nk
! ======================================================================

  !    if( mod(ik-1,npool)/=mypool ) cycle

      write(stdout,'(a,3f12.5,i6,a,i6,a,i3)') ' k-point ', &
        kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                    ' on node ', mpime

      ! build the Hamiltonian for this q-point
      call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin )

      if( kinetic_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_kink( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else if( local_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_vlock( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else if( nonlocal_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_vnlk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else if( smatrix_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_sk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
        forall( i=1:nbasis ) eigval(i)=real(ztmp(i,i))
        deallocate( ztmp )
      else

  !    call diagx_ham
        call diag_ham

      endif
      qin( : ) =  kpt%list%kvec( 1 : 3, ik )
      qcart(:) = 0.d0
      if( kpt%param%cartesian ) then
        qcart(:) = qin(:)*tpiba
      else
        do i=1,3
          qcart(:) = qcart(:) + bvec(:,i)*qin(i)
        enddo
      endif

      allocate( bofr( local_npt, desc_cyclic(N_) ) )
      allocate( uofrandb( local_npt, desc_cyclic(N_) ) )


! This is probably sub-optimal, but requires the least thinking, I think.
! Delay bofr distribution
! Add in phase for the k-point we are at
      if( mypoolid .eq. 0 ) then
        call DGEMV( 'T', 3, npt, real_one, posn(1,1), 3, qcart, 1, real_zero, phase, 1 )
        eikr( : ) = exp( iota * phase( : ) )
        do ibd = 1, nbasis
          phased_bofr( :, ibd ) = single_bofr( :, ibd ) * eikr( : )
        enddo
      endif
! Now distribute bofr, or don't
!      write(stdout,*) 'Distribute bofr'
      call PZGEMR2D( npt, nbasis, phased_bofr, 1, 1, desc_bofr_in, bofr, 1, 1, desc_bofr, context_cyclic )
!      call mp_barrier
!      write(stdout,*) 'bofr distributed'

      ! can probably change to nbasis_subset
      call PZGEMM( 'N', 'N', npt, nbasis, nbasis, &
                   one, bofr, 1, 1, desc_bofr, eigvec, 1, 1, desc_cyclic, &
                   zero, uofrandb, 1, 1, desc_bofr )

      do ibd = 1, nbasis_subset
        su = REAL( SUM( CONJG( uofrandb( :, ibd ) ) * wpt( : ) * uofrandb( :, ibd ) ) )
        su = sqrt( 1.0_DP / su )
        uofrandb( :, ibd ) = CONJG(uofrandb( :, ibd )) * su !WPT?
      enddo

!     Eigenvalues are in ascending order, want only the top nbasis_trunc
      bi = (ik-1)*nbasis_subset+1
      bj = 1+big_nbasis-nbasis_trunc
      call PZGEMM( 'N', 'N', npt, nbasis_trunc, nbasis_subset, &
                   one, uofrandb, 1, 1, desc_bofr, B_l, bi, bj, desc_S, &
                   one, basis_r, 1, 1, desc_B )

    enddo ! ik              

    ! Need to dump to file here

    
    !JTV




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now construct W(i,j) using these basis functions
    if( mypoolid .eq. mypoolroot ) then
      ! Open up W( r )
    endif
    ! share with the pool
    call mp_bcast( W_r, mypoolroot, intra_pool_comm )
    !
    do ipt = 1, local_npt
      
    enddo
    allocate( Wij( nbasis_trunc, nbasis_trunc ) )
    
    


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now create B_i^{nk} for use in BSE



  enddo ! ispin

  enddo ! itau

  call mp_barrier
  write(stdout,*) 'Loop over k-points and bands complete.'

!  if( mypoolid .eq. mypoolroot ) then
!    call MPI_FILE_CLOSE( fh_rad, ierr )
!  endif
 
  call mp_end
!  stop
111 continue
  
end program 
