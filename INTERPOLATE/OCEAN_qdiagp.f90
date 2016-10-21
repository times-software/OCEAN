! Copyright (C) 2015,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the GPL 2 License. See the file `License' in the current subdirectory.
!
  program ocean_qdiag  

  ! Originally called shirley_qdiag

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Nov 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
!  use parallel_include !JTV
!  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : mpime, root  !nproc
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use mpio
!  use kpt_module !JTV
!  use corerepair_module ! JTV
  use shirley_input_module 
  use diag_module
!  USE wavefunctions_module, ONLY: evc
!  USE io_files, ONLY: diropn !nwordwfc, iunwfc, prefix_io=>prefix, &
                             !     tmpdir_io=>tmp_dir, nd_nmbr, diropn, wfc_dir
!  USE control_flags,        ONLY : lscf
!  USE basis,                ONLY : starting_pot, starting_wfc
  use scalapack_module, only : DLEN_!,
  use hamq_pool, only : npool, &  !nproc_per_pool,rootpool, intra_pool_comm, desc_striped, context_global,
                        mypool, mypoolid, mypoolroot, &
                        cross_pool_comm, &
                        desc_cyclic, context_cyclic, &
                        local_striped_dim, cyclic_localindex, &
                        local_cyclic_dims
  use OCEAN_timer

  implicit none

  integer,external :: freeunit

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)
  complex(dp),parameter :: iota=(0.d0,1.d0)
  real(dp), parameter :: eps = 1.0d-14
  real(dp), parameter :: real_one = 1_dp
  real(dp), parameter :: real_zero = 0_dp


  character(255) :: info_file
  logical :: ex, mpiio_workaround
  logical :: legacy_zero_lumo, noshiftlumo

  integer(kind=MPI_OFFSET_KIND) :: offset, locsize_x, basis_vector_size_x
  integer :: fheigval, fheig
!  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  character(255) :: string

  real(dp) :: dumf

!  real(dp) :: kvec(3)
!  real(dp) :: qpathlen, cqvec(3), dqvec(3)
!  character(255) :: ic
!  character(6), allocatable :: fntau(:)
  character (len=5) :: band_style = 'defau'

  integer :: ik
  integer :: i
  integer :: nbasis_subset

  integer :: nk, nbnd, ispin
  real(dp) :: nelec_, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin,  nxpts, xmesh( 3 )
  real(dp) :: fermi_energy
!  real(dp),allocatable :: kr( :, : )
!  integer,allocatable :: gflip( :, : )
  complex(dp),allocatable :: ztmp(:,:)

  integer :: iuntmp, ntot
!  complex(dp),allocatable :: o2l(:,:,:,:,:), coeff( :, :, :, :, :, : ), coeff_small( :, :, :, :, : )
  real(dp),allocatable :: e0(:,:,:,:)
  logical :: exst

  complex(dp),allocatable :: tmels(:,:,: )

  integer,allocatable :: start_band(:), new_start_band(:,:,:)
  real(dp),allocatable :: lumo(:), homo(:)
  integer :: brange( 4 ), nelectron, max_val
  real(dp) :: lumo_point, homo_point,dft_energy_range, kshift( 3 ), kplusq(3), lumo_shift
  logical :: have_kshift

  integer :: iuninf, fhu2, errorcode, fmode, resultlen, nshift, fhtmels, fh_val_energies, fh_con_energies
  integer :: fh_val_eigvecs, fh_con_eigvecs

  integer :: nband, val_band
!  real(dp) :: qin(3), qcart(3), bvec(3,3)
  complex(dp), allocatable :: eigvec_single(:,:)

!  integer :: u1_MB, u1_M, u1_NB, u1_N
!  integer :: nprow, npcol, myrow, mycol


  integer :: nr_eigvec, nc_eigvec, pool_tot_k, pool_ik, eigvec_info

  integer :: pool_rank
!  integer(8) :: long_nbasis, long_nx
!  integer(8) :: long_2g = 2147483648

  integer, dimension( DLEN_ ) :: desc_tmels,  desc_eigvec_single

!  integer, external :: OMP_GET_NUM_THREADS, NUMROC
!  complex(dp), external :: ZDOTC
!  real(dp), external :: DZNRM2


  integer :: dims(2), ndims, array_of_gsizes(2), array_of_distribs(2), array_of_dargs(2), file_type, nelements, mpistatus, &
             locsize, basis_vector_type
  integer, allocatable :: energy_request(:), pool_root_map(:)
  complex(dp), allocatable :: store_eigvec(:,:,:,:), out_eigvec(:,:,:)
!  real(dp), allocatable :: store_energy(:,:)
 
  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u


! Read in the input file, shirley_input sets up what 
!  is needed for the interpolation
  call shirley_input

#ifdef FALSE
! for right now we don't though, everything is taken care of earlier
  tmpdir_io = outdir
  prefix_io = prefix
#ifdef __NIST
  wfc_dir = wfcdir
#endif
  if( .true. ) then
  lscf = .false.
  starting_pot = 'file'
  starting_wfc = 'file'
  call dump_system( nelec_, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
  write(stdout,*) ' nspin = ', nspin
  write(stdout,*) ' lda_plus_u = ', lda_plus_u

  endif
#endif


  if( band_subset(1) < 1 ) band_subset(1)=1
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

  call diag_init

  call dump_system( nelec_, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
  write(stdout,*) ' nspin = ', nspin
  write(stdout,*) ' lda_plus_u = ', lda_plus_u



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
    iuntmp = freeunit()

! OCEAN energy offset. See tutorial on alignment for more complete thoughts !
! The legacy scheme sets the LUMO to be 0    
    legacy_zero_lumo = .true.

! If we have core_offset then we (almost certainly) want to not touch the LUMO energy
! core_offset either contains .false. or a number
    inquire(file='core_offset',exist=ex)
    if( ex ) then
      open(iuntmp,file='core_offset',form='formatted',status='old')
      read(iuntmp,*,err=1001) dumf
      legacy_zero_lumo = .false.
1001 continue
      close(iuntmp) 
      write(stdout,*) 'LUMO shift from core_offset:', legacy_zero_lumo
    endif

! But if you want to override then allow noshift_lumo to be set
    inquire(file='noshift_lumo',exist=ex)
    if( ex ) then
      open(unit=iuntmp,file='noshift_lumo',form='formatted',status='old')
      read(iuntmp,*) noshiftlumo
      close(iuntmp)
      if( noshiftlumo ) then
        legacy_zero_lumo = .false.
      else
        legacy_zero_lumo = .true.
      endif
      write(stdout,*) 'LUMO shift from no_lumoshiftt:', legacy_zero_lumo
    endif

    inquire(file='band_style.inp',exist=ex)
    if( ex ) then
      open(unit=iuntmp,file='band_style.inp',form='formatted',status='old')
      read(iuntmp,*) band_style
      close(iuntmp)
    else
      band_style = 'defau'
    endif
      

    kshift = 0.0d0
    inquire(file='qinunitsofbvectors.ipt',exist=have_kshift)
    if( have_kshift) then
      open(unit=iuntmp,file='qinunitsofbvectors.ipt',form='formatted',status='old')
      read(iuntmp,*) kshift( : )
      close(iuntmp)
      if( sum(abs(kshift( : ) ) ) .lt. 1.d-14 ) have_kshift = .false.
    endif
  ! Do we store k + q ? or just k
    nshift = 1
    if( have_kshift ) nshift = 2


    iuntmp = freeunit()
    open(unit=iuntmp,file='xmesh.ipt',form='formatted',status='old')
    read(iuntmp,*) xmesh(:)
    close(iuntmp)
    nxpts = product( xmesh )


!    open( unit=iuntmp, file='bvecs',form='formatted',status='old')
!    read(iuntmp,*) bvec(:,:)
!    close(iuntmp)


    open(unit=iuntmp,file='efermiinrydberg.ipt',form='formatted',status='old')
    read(iuntmp,*) fermi_energy
    close(iuntmp)

    inquire(file='dft_energy_range.ipt',exist=exst)
    if( exst ) then
      open(unit=iuntmp,file='dft_energy_range.ipt',form='formatted',status='old')
      read(iuntmp,*) dft_energy_range 
      close(iuntmp)
      dft_energy_range = dft_energy_range / rytoev
    else
      dft_energy_range = -1
    endif

    open(unit=iuntmp,file='nelectron',form='formatted',status='old')
    read(iuntmp,*) nelectron
    close(iuntmp)

  endif

  
  call mp_bcast( legacy_zero_lumo, ionode_id )
  call mp_bcast( kshift, ionode_id )
  call mp_bcast( have_kshift, ionode_id )
  call mp_bcast( nshift, ionode_id )

  call mp_bcast( fermi_energy, ionode_id )
  call mp_bcast( dft_energy_range, ionode_id)
  
  call mp_bcast( nelectron, ionode_id )


  ntot = nbasis_subset * kpt%list%nk

  nband = band_subset(2) - band_subset(1) + 1

  call descinit( desc_eigvec_single, nbasis, nbasis, nbasis, nbasis, 0, 0, context_cyclic, nbasis, ierr )

  if( mypoolid .eq. mypoolroot ) then
!!???    allocate( u2( nxpts, band_subset(1) : band_subset(2), kpt%list%nk, nshift ) )
!    allocate( u2( nxpts, band_subset(1) : band_subset(2), 1, nshift ) )
    allocate( eigvec_single( nbasis, nbasis ) )
  else  ! For error checking nonsense
!!!!    allocate( u2( 1, band_subset(1), kpt%list%nk, nshift ) )
!     ! JTV memleak 29 July
!!    allocate( u2( 1, band_subset(1), 1, nshift ) )
!    allocate( u2( nxpts, band_subset(1) : band_subset(2), 1, nshift ) )
    allocate( eigvec_single( 1, 1:band_subset(1) ) )
  endif
!  u2 = 0.d0
!  call descinit( desc_u2, nxpts, nband, nxpts, nband, 0, 0, context_cyclic, nxpts, ierr )


  if( ionode ) write(6,*) 'nshift =', nshift


  allocate(e0( band_subset(1):band_subset(2), kpt%list%nk, nspin, nshift) )
  e0 = 0.d0

  allocate( lumo(kpt%list%nk), start_band(kpt%list%nk), homo(kpt%list%nk) )
  lumo = 0.d0
  start_band = 0

  ! Set up different options for valence/conduction split
  select case (band_style)
    case( 'metal' )
      max_val = nelectron / 2 + 10
    case( 'bands' )
      max_val = nelectron / 2
    case default
      max_val = nelectron / 2 + 10
  end select
  max_val = min( max_val, band_subset(2) )
!  if( have_kshift ) then
!    allocate( tmels( max_val, band_subset( 1 ) : band_subset( 2 ), kpt%list%nk ) )
!    tmels = 0.0d0
!!    allocate( tmp_eigvec( desc_cyclic(M_), desc_cyclic(N_) ) )
!    call descinit( desc_tmels,  max_val, band_subset( 2 )-band_subset( 1 )+1, &
!                                max_val, band_subset( 2 )-band_subset( 1 )+1, &
!                   0, 0, context_cyclic, max_val, ierr )

! dec 3
!    tmp_eig_n = numroc( max_val, desc_cyclic(NB_), mycol, 0, npcol )
!    tmp_eig_m = numroc( nbasis, desc_cyclic(MB_), myrow, 0, nprow )
!    call descinit( desc_tmpeig, nbasis, max_val, desc_cyclic(MB_), desc_cyclic(NB_), 0, 0, context_cyclic, tmp_eig_m, ierr )
!    if( ierr .ne. 0 ) then
!      write(stdout,*) desc_tmpeig(:)
!      stop
!    endif
!
!
!    allocate( tmp_eigvec(  tmp_eig_m, tmp_eig_n ) )

!    allocate( tmp_eigvec( nbasis, max_val ) )
!  endif



  ! ========= pre-fetch all of the energies? =========== !
  ! Determine number of k-points each processor's pool will do
  ! Get an array of the ids for each pool head
  ! Energy_request is sequential for ionode and by kpt/spin for everyone else
  allocate( energy_request( kpt%list%nk * nspin ), pool_root_map( 0:npool-1 ) )

  pool_root_map(:) = 0
  if( mypoolid .eq. mypoolroot ) then
    ! mpime is wrong!
!    pool_root_map(mypool) = mpime
    call MPI_COMM_SIZE( cross_pool_comm, i, ierr )
    write(stdout,*) 'Cross pool size:', i
    call MPI_COMM_RANK( cross_pool_comm, pool_rank, ierr )
    write(stdout,*) 'Cross pool id:', pool_rank, ionode_id
    pool_root_map(mypool) = pool_rank
    call MPI_ALLREDUCE( MPI_IN_PLACE, pool_root_map, npool, MPI_INTEGER, MPI_SUM, cross_pool_comm, ierr )
  endif

  do i = 0, npool - 1 
    write(stdout,*) i, pool_root_map( i )
  enddo

!  if( have_kshift ) then
!    nband = max_val
!  else
    nband = band_subset(2) - band_subset(1) + 1
!  endif

  pool_ik = 0

!  i_energy_request = 0
  do ispin=1,nspin
    do ik=1,kpt%list%nk
!      if( ionode .and. (mod((ik-1)+(ispin-1)*kpt%list%nk,npool) .ne. mypool ) ) then
!        i_energy_request = i_energy_request + 1
!        call MPI_IRECV( e0( 1, ik, ispin, 1 ), nband, MPI_DOUBLE_PRECISION,  &
!                        pool_root_map(mod((ik-1)+(ispin-1)*kpt%list%nk,npool)), ik+(ispin-1)*kpt%list%nk, &
!                        cross_pool_comm, energy_request( i_energy_request ), ierr )
        write(stdout,*) pool_root_map(mod((ik-1)+(ispin-1)*kpt%list%nk,npool)), ik+(ispin-1)*kpt%list%nk
!      endif
      if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle
      pool_ik = pool_ik + 1
    enddo
  enddo
  pool_tot_k = max( pool_ik, 1 )
!  n_energy_request = i_energy_request


  ! Allocate storage for eigenvectors
  ! Create descriptors  
  ! For laziness start by copying *all* of eigvector which is way larger than we need 
  call local_cyclic_dims( nr_eigvec, nc_eigvec )
  allocate( store_eigvec( nr_eigvec, nc_eigvec, pool_tot_k, nshift ) ) !, &
!            store_energy( nband, pool_tot_k ) )
  



  pool_ik = 0
  do ispin=1,nspin
    do ik=1,kpt%list%nk
! ======================================================================

      ! (ik-1)+(ispin-1)*kpt%list%nk
      if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle


      pool_ik = pool_ik + 1

  
      write(stdout,'(a,3f12.5,i6,a,i6,a,i3)') ' k-point ', &
        kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                    ' on node ', mpime

      ! build the Hamiltonian for this q-point
      call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin )

      if( kinetic_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_kink( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
        forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
        deallocate( ztmp )
      else if( local_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_vlock( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
        forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
        deallocate( ztmp )
      else if( nonlocal_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_vnlk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
        forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
        deallocate( ztmp )
      else if( smatrix_only ) then
        allocate( ztmp(nbasis,nbasis) )
        call diag_build_sk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
        forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
        deallocate( ztmp )
      else
        call diag_ham
      endif


      write(stdout,*) ik,  eigval(band_subset(1))*rytoev, eigval(band_subset(2))*rytoev
      e0( :, ik, ispin, 1 ) = 0.5d0 * eigval(band_subset(1):band_subset(2))


      !!!!!!! store away eigenvectors
      store_eigvec( :, :, pool_ik, 1 ) = eigvec( :, : )

      if( have_kshift ) then
         kplusq(:) = kpt%list%kvec(1:3,ik) + kshift(:)
        ! build the Hamiltonian for this q-point
        call diag_build_hamk( kplusq, kpt%param%cartesian, ispin )

        if( kinetic_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_kink( kplusq, kpt%param%cartesian, ztmp )
          forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
          deallocate( ztmp )
        else if( local_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_vlock( kplusq, kpt%param%cartesian, ispin, ztmp )
          forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
          deallocate( ztmp )
        else if( nonlocal_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_vnlk( kplusq, kpt%param%cartesian, ispin, ztmp )
          forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
          deallocate( ztmp )
        else if( smatrix_only ) then
          allocate( ztmp(nbasis,nbasis) )
          call diag_build_sk( kplusq, kpt%param%cartesian, ztmp )
          forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
          deallocate( ztmp )
        else
          call diag_ham
        endif


        write(stdout,*) ik,  eigval(band_subset(1))*rytoev, eigval(band_subset(2))*rytoev
        e0( :, ik, ispin, 2 ) = 0.5d0 * eigval(band_subset(1):band_subset(2))

        store_eigvec( :, :, pool_ik, 2 ) = eigvec( :, : )

      endif

    enddo ! ik
  enddo ! ispin
  ! done with interpolation here

  call mp_barrier
  write(stdout,*) 'Done with diag'



  ! Clean up energy communications and find fermi, lumo, homo 

  ! Find start bands and save out energyfile
  
#ifdef FALSE
  if( ionode ) then
    do ispin = 1, nspin
      do ik = 1, kpt%list%nk
        if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)==mypool ) cycle
        call MPI_RECV( e0( :, ik, ispin, 1 ), nband, MPI_DOUBLE_PRECISION, &
!                       pool_root_map(mod((ik-1)+(ispin-1)*kpt%list%nk,npool)), ik+(ispin-1)*kpt%list%nk, &
                       MPI_ANY_SOURCE, ik+(ispin-1)*kpt%list%nk, &
                       cross_pool_comm, MPI_STATUS_IGNORE, ierr )
      enddo
    enddo
  elseif( mypoolid .eq. mypoolroot ) then
    do ispin = 1, nspin
      do ik = 1, kpt%list%nk
        if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle
        call MPI_SEND( e0( :, ik, ispin, 1 ), nband, MPI_DOUBLE_PRECISION, &
                       ionode_id, ik+(ispin-1)*kpt%list%nk, &
                       cross_pool_comm, MPI_STATUS_IGNORE, ierr )
      enddo
    enddo
  endif
#else
  ! Possibly having problems with the non-blocking
  if( mypoolid .eq. mypoolroot ) then
    call mp_sum( e0, cross_pool_comm )
  endif
#endif

  if( ionode ) then
    write(stdout,*) 'Sharing energies'

!    call MPI_WAITALL( n_energy_request, energy_request, MPI_STATUSES_IGNORE, ierr )

    call fix_fermi( nband, kpt%list%nk, nspin, nshift, max_val, nelectron, 0, &
                    e0, homo_point, lumo_point, fermi_energy )

    
    allocate( new_start_band( kpt%list%nk, nspin, nshift ) )
    call find_startband( nband, kpt%list%nk, nspin, nshift, nelectron, fermi_energy, &
                         band_style, e0, new_start_band )


    if( legacy_zero_lumo ) then
      lumo_shift = lumo_point * 0.5_DP
    else
      lumo_shift = 0.0_DP
    endif

    
    call dump_energies( band_subset, nband, kpt%list%nk, nspin, nshift, e0, lumo_shift, new_start_band, brange, ierr )

    write(stdout,*) 'Sharing energies'
    write(stdout,*) brange(1), brange(2)
    write(stdout,*) brange(3), brange(4)
  endif

  call mp_bcast( brange, ionode_id )
  if( .not. ionode ) then
    allocate( new_start_band( kpt%list%nk, nspin, nshift ) )
  endif
  call mp_bcast( new_start_band, ionode_id )

!  if( mypoolid .eq. 0 ) then
!    if( .not. ionode ) then
!      do ispin=1,nspin
!        do ik=1,kpt%list%nk
!          if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle
!          call MPI_WAIT( energy_request( ik+(ispin-1)*kpt%list%nk ), MPI_STATUS_IGNORE, ierr )
!        enddo
!      enddo
!    endif
!    call MPI_BCAST( fermi_energy, 1, MPI_DOUBLE_PRECISION, ionode_id, cross_pool_comm, ierr )
!  endif
!  call MPI_BCAST( fermi_energy, 1, MPI_DOUBLE_PRECISION, mypoolroot, intra_pool_comm, ierr )

  deallocate( energy_request )
  write(stdout,*) 'Done sharing energies'

!! Everything needs to be C ordering
!  ndims = 3
!  dims = ( npool, nprow, npcol )
!  array_of_gsizes = (/ kpt%list%nk, nbasis, nbasis /)
!  array_of_distribs = (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC MPI_DISTRIBUTE_CYCLIC /)
!  array_of_dargs = (/ 1, desc_cyclic[MB_], desc_cyclic[NB_] /)
!  periods = (/.false.,.false.,.false./)
!  reorder=.false.
 ! call MPI_CART_CREATE( mpi_comm_world, ndims, dims, periods, reorder, file_comm, ierr )
!  CALL MPI_TYPE_CREATE_DARRAY( nproc, mpime, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
!                               dims, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, file_type )


!  nband = band_subset(2) - band_subset(1) + 1
  nband = brange( 2 ) - brange( 1 ) + 1 + brange( 4 ) - brange( 3 ) + 1
  if( mypoolid .eq. mypoolroot ) then
    allocate( out_eigvec( nbasis, nband, pool_tot_k ) )
  else
    allocate( out_eigvec( 1, 1, pool_tot_k ) )
  endif

  pool_ik = 0
  do ispin=1,nspin
    do ik=1,kpt%list%nk
! ======================================================================

      ! (ik-1)+(ispin-1)*kpt%list%nk
      if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle

      pool_ik = pool_ik + 1

      if( have_kshift ) then

        nband = brange( 2 ) - brange( 1 ) + 1
        call PZGEMR2D( nbasis, nband, store_eigvec( 1, 1, pool_ik, 1 ), 1, brange( 1 ), desc_cyclic, &
                                      out_eigvec( 1, 1, pool_ik ), 1, 1, desc_eigvec_single, &
                       context_cyclic, ierr )

        ! val_band gives out offset of out_eigvec
        val_band = nband + 1
        nband = brange( 4 ) - brange( 3 ) + 1
        call PZGEMR2D( nbasis, nband, store_eigvec( 1, 1, pool_ik, 2 ), 1, new_start_band( ik, ispin, 2 ), desc_cyclic, &
                                      out_eigvec( 1, 1, pool_ik ), 1, val_band, desc_eigvec_single, &
                       context_cyclic, ierr )
        nband = brange( 2 ) - brange( 1 ) + 1 + brange( 4 ) - brange( 3 ) + 1

      else

        call PZGEMR2D( nbasis, nband, store_eigvec( 1, 1, pool_ik, 1 ), 1, 1, desc_cyclic, &
                                      out_eigvec( 1, 1, pool_ik ), 1, 1, desc_eigvec_single, &
                       context_cyclic, ierr )
      endif
    enddo
  enddo

  call mp_barrier
  write(stdout,*) 'Writing out eigvecs'

!  if( nbasis > 32000 .or. nbasis * nband > 2**27/npool )  then

  ! If each k-point takes up less than 256MB maybe try and use good mpi mapping
  if( int( nbasis, MPI_OFFSET_KIND ) * int( nband, MPI_OFFSET_KIND ) < 16777216_MPI_OFFSET_KIND ) then
    mpiio_workaround = .false.
  else

    mpiio_workaround = .true.
    write(stdout,*) 'Using mpi/io workaround for eigvecs.dat'
    write(stdout,*) '   nbasis =', nbasis
    write(stdout,*) '   ', nbasis * nband, 2**27/npool
  endif
!  else
!    mpiio_workaround = .false.
!  endif

  mpiio_workaround = .true.
  if( mypoolid .eq. mypoolroot ) then 
    if( .not. mpiio_workaround ) then

    ndims = 2
    dims = (/ 1, npool /)
    array_of_gsizes = (/ nbasis * nband, nspin * kpt%list%nk /)
    array_of_distribs = (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC /)
    array_of_dargs = (/ nbasis * nband, 1 /)
    !
    !CALL MPI_TYPE_CREATE_DARRAY( npool, mypool, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
    CALL MPI_TYPE_CREATE_DARRAY( npool, pool_rank, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
                                 dims, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, file_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Create DARRAY",string,errorcode)
    endif

    CALL MPI_TYPE_COMMIT( file_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Commit filetype",string,errorcode)
    endif

    ! clear out existing
!    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_DELETE_ON_CLOSE)
!    fmode = MPI_MODE_CREATE
!    call MPI_FILE_OPEN( cross_pool_comm, 'eigvecs.dat', fmode, MPI_INFO_NULL, fheig, ierr )
!    if( ierr/=0 ) then
!      errorcode=ierr
!      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
!      write(6,*) string
!      call errore(string,errorcode)
!    endif
!    call MPI_FILE_CLOSE( fheig, ierr )
    if( ionode .and. .false. ) then
      inquire(file='eigvec.dat',exist=ex)
      if( ex ) then
        fheigval=freeunit()
        open(fheigval,file='eigvec.dat',form='unformatted')
        close(fheigval,status='delete')
      endif
    endif
    call MPI_BARRIER( cross_pool_comm, ierr )



    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    fmode = IOR(fmode, MPI_MODE_UNIQUE_OPEN)

    call MPI_INFO_CREATE( eigvec_info, ierr )
!    call MPI_INFO_SET( eigvec_info, 'key', 'value', ierr )
    call MPI_INFO_SET( eigvec_info, 'access_style', 'write_once', ierr )
!    call MPI_INFO_SET( eigvec_info, 'cb_nodes', '4', ierr )
!    call MPI_INFO_SET( eigvec_info, 'striping_factor', '4', ierr )

    call MPI_FILE_OPEN( cross_pool_comm, 'eigvecs.dat', fmode, eigvec_info, fheig, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Open eigvecs.dat",string,errorcode)
    endif
    offset=0
    !JTV At this point it would be good to create a custom MPI_DATATYPE
    !  so that we can get optimized file writing
    call MPI_FILE_SET_VIEW( fheig, offset, MPI_DOUBLE_COMPLEX, &
                            file_type, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Set view fheig",string,errorcode)
    endif

    call MPI_TYPE_SIZE( file_type, locsize, ierr )
    nelements = locsize / 16
!    if( nelements .ne.  nbasis * nband * pool_tot_k ) then
!    endif
    call MPI_BARRIER( cross_pool_comm, ierr )
    call OCEAN_t_reset
    call MPI_FILE_WRITE_ALL( fheig, out_eigvec, nelements, MPI_DOUBLE_COMPLEX, mpistatus, ierr )
    call OCEAN_t_printtime( "File out", stdout)
    write(6,*) nbasis * nband * nspin * kpt%list%nk / 64
    
    call MPI_FILE_CLOSE( fheig, ierr )

    call MPI_TYPE_FREE( file_type, ierr )

  elseif( .true. ) then

    ! This shouldn't break until either nbasis * 16 > 2GB or nband > 2GB

    CALL MPI_TYPE_CONTIGUOUS( nbasis, MPI_DOUBLE_COMPLEX, file_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Create DARRAY",string,errorcode)
    endif

    CALL MPI_TYPE_COMMIT( file_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Commit filetype",string,errorcode)
    endif


    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
!    fmode = IOR(fmode, MPI_MODE_UNIQUE_OPEN)

!    call MPI_INFO_CREATE( eigvec_info, ierr )
!    call MPI_INFO_SET( eigvec_info, 'access_style', 'write_once', ierr )
  
!    call MPI_FILE_OPEN( cross_pool_comm, 'eigvecs.dat', fmode, eigvec_info, fheig, ierr )
    call MPI_FILE_OPEN( cross_pool_comm, 'eigvecs.dat', fmode, MPI_INFO_NULL, fheig, ierr )
    if( ierr/=0 ) then
      errorcode=ierr 
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Open eigvecs.dat",string,errorcode)
    endif
    offset=0
    !JTV At this point it would be good to create a custom MPI_DATATYPE
    !  so that we can get optimized file writing
!    call MPI_FILE_SET_VIEW( fheig, offset, MPI_DOUBLE_COMPLEX, &
!                            file_type, 'native', MPI_INFO_NULL, ierr )
    call MPI_FILE_SET_VIEW( fheig, offset, file_type, &
                            file_type, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Set view fheig",string,errorcode)
    endif

    call MPI_TYPE_SIZE( file_type, locsize, ierr )
    nelements = locsize / 16
    if( nelements .ne.  nbasis ) then
      write(stdout,*) 'Something is wrong'
    endif
    nelements = nband
    call MPI_BARRIER( cross_pool_comm, ierr )
    call OCEAN_t_reset

    pool_ik = 0
    if( .false. ) then
      do ik = 0, nspin * kpt%list%nk - 1, npool
        offset = ik + mypool
        pool_ik = pool_ik + 1

        if( pool_ik .gt. pool_tot_k ) then 
          nelements = 0                ! write nothing
          offset = nspin * kpt%list%nk ! do it at the end of the file
          pool_ik = pool_tot_k         ! avoid bounds fails
        endif

        offset = offset * int(nband,MPI_OFFSET_KIND)

        call MPI_FILE_WRITE_AT_ALL( fheig, offset, out_eigvec(1,1,pool_ik), nelements, file_type, mpistatus, ierr )
      enddo
    else

      do ispin = 1, nspin
      do ik = 1, kpt%list%nk

        if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle

        offset = int( ik-1 + (ispin-1)*kpt%list%nk, MPI_OFFSET_KIND ) * int( nband, MPI_OFFSET_KIND )
        pool_ik = pool_ik + 1
    
        call MPI_FILE_WRITE_AT( fheig, offset, out_eigvec(1,1,pool_ik), nelements, file_type, mpistatus, ierr )

      enddo
      enddo
    endif

    call OCEAN_t_printtime( "File out", stdout)
    write(stdout,*) int(nbasis,MPI_OFFSET_KIND) * int(nband,MPI_OFFSET_KIND) * (nspin * kpt%list%nk)

    call MPI_FILE_CLOSE( fheig, ierr )

  else  ! Trying a completely different method here

    write(stdout,*) 'Using new approach for eigvecs'
  
    CALL MPI_TYPE_CONTIGUOUS( nbasis, MPI_DOUBLE_COMPLEX, basis_vector_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Create CONTIGUOUS",string,errorcode)
    endif

    CALL MPI_TYPE_COMMIT( basis_vector_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Commit basis_vector_type",string,errorcode)
    endif


    ndims = 2
    dims = (/ 1, npool /)
    array_of_gsizes = (/ nband, nspin * kpt%list%nk /)
    array_of_distribs = (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC /)
    array_of_dargs = (/ nband, 1 /)
    !
    !CALL MPI_TYPE_CREATE_DARRAY( npool, mypool, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
    CALL MPI_TYPE_CREATE_DARRAY( npool, pool_rank, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
                                 dims, MPI_ORDER_FORTRAN, basis_vector_type, file_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Create DARRAY",string,errorcode)
    endif

    CALL MPI_TYPE_COMMIT( file_type, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Commit filetype",string,errorcode)
    endif

    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    fmode = IOR(fmode, MPI_MODE_UNIQUE_OPEN)

    call MPI_INFO_CREATE( eigvec_info, ierr )
    call MPI_INFO_SET( eigvec_info, 'access_style', 'write_once', ierr )

    call MPI_FILE_OPEN( cross_pool_comm, 'eigvecs.dat', fmode, eigvec_info, fheig, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Open eigvecs.dat",string,errorcode)
    endif
    offset=0
    !JTV At this point it would be good to create a custom MPI_DATATYPE
    !  so that we can get optimized file writing
    call MPI_FILE_SET_VIEW( fheig, offset, basis_vector_type, &
                            file_type, 'native', MPI_INFO_NULL, ierr )

    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore("Set view fheig",string,errorcode)
    endif

!    call MPI_TYPE_SIZE( file_type, locsize, ierr )
!    nelements = locsize / 16
!    call MPI_TYPE_SIZE( basis_vector_type, basis_vector_size, ierr )

    call MPI_TYPE_SIZE_X( file_type, locsize_x, ierr )
    call MPI_TYPE_SIZE_X( basis_vector_type, basis_vector_size_x, ierr )

    offset = locsize_x / basis_vector_size_x
    nelements = offset

!    nelements = locsize / basis_vector_size
    write(stdout,*) nelements, locsize_x, basis_vector_size_x
!    if( nelements .ne.  nbasis * nband * pool_tot_k ) then
!    endif
    call MPI_BARRIER( cross_pool_comm, ierr )
    call OCEAN_t_reset
    call MPI_FILE_WRITE_ALL( fheig, out_eigvec, nelements, basis_vector_type, mpistatus, ierr )
    call OCEAN_t_printtime( "File out", stdout)
    call MPI_GET_COUNT( mpistatus, basis_vector_type, nelements, ierr )
    write(stdout,*) nelements

    call MPI_FILE_CLOSE( fheig, ierr )



  endif
  endif


  write(stdout,*) 'qdiag.info'
  if( ionode ) then
    iuntmp = freeunit()
    open(iuntmp,file='qdiag.info',form='formatted',status='unknown')
    write(iuntmp,*) nbasis, nband, kpt%list%nk, nspin, nshift
    close(iuntmp)
  endif

  write(stdout,*) 'obf_control'
  if( ionode ) then
    open(unit=iuntmp,file='obf_control',form='formatted',status='unknown')
    rewind(iuntmp)
    write(iuntmp,*) nbasis
    write(iuntmp,*) max_val
    write(iuntmp,*) nbasis_subset
  endif



! Need to calculate and write tmels
  if( have_kshift ) then

    write(stdout,*) 'Calculating tmels'

    val_band = brange(2) - brange(1) + 1
    nband = brange(4) - brange(3) + 1

    write(stdout,*) val_band, nband

    if( mypoolid == mypoolroot ) then
      allocate( tmels( val_band, nband, pool_tot_k ) )
    else
      allocate( tmels( 1, 1, pool_tot_k ) )
    endif
!    tmels = 0.0d0
    call descinit( desc_tmels, val_band, nband, val_band, nband, 0, 0, context_cyclic, val_band, ierr )

    pool_ik = 0
    do ispin=1,nspin
      do ik = 1, kpt%list%nk
        if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle

        pool_ik = pool_ik + 1

        call PZGEMM( 'C', 'N', val_band, nband, nbasis, one, &
                     store_eigvec( 1, 1, pool_ik, 1 ), 1, brange(1), desc_cyclic, &
                     store_eigvec( 1, 1, pool_ik, 2 ), 1, brange(3), desc_cyclic, &
                     zero, tmels( 1, 1, pool_ik ), 1, 1, desc_tmels )

      enddo
    enddo

    if( mypoolid == mypoolroot ) then
      write(stdout,*) "Writing out tmels"

      if( .false. ) then 
        ndims = 2
        dims = (/ 1, npool /)
        array_of_gsizes = (/ val_band * nband, nspin * kpt%list%nk /)
        array_of_distribs = (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC /)
        array_of_dargs = (/ val_band * nband, 1 /)
        !
        CALL MPI_TYPE_CREATE_DARRAY( npool, pool_rank, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
                                     dims, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, file_type, ierr )
        if( ierr/=0 ) then
          errorcode=ierr
          call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
          call errore("Create DARRAY",string,errorcode)
        endif

        CALL MPI_TYPE_COMMIT( file_type, ierr )
        if( ierr/=0 ) then
          errorcode=ierr
          call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
          call errore("Commit filetype",string,errorcode)
        endif

        fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        fmode = IOR(fmode, MPI_MODE_UNIQUE_OPEN)

        call MPI_FILE_OPEN( cross_pool_comm, 'ptmels.dat', fmode, MPI_INFO_NULL, fhtmels, ierr )
        if( ierr/=0 ) then
          errorcode=ierr
          call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
          call errore("Open eigvecs.dat",string,errorcode)
        endif
        offset=0
        call MPI_FILE_SET_VIEW( fhtmels, offset, MPI_DOUBLE_COMPLEX, &
                                file_type, 'native', MPI_INFO_NULL, ierr )
        if( ierr/=0 ) then
          errorcode=ierr
          call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
          call errore("Set view fhtmels",string,errorcode)
        endif

        call MPI_TYPE_SIZE( file_type, locsize, ierr )
        nelements = locsize / 16

        call OCEAN_t_reset
        call MPI_FILE_WRITE_ALL( fhtmels, tmels, nelements, MPI_DOUBLE_COMPLEX, mpistatus, ierr )
        call OCEAN_t_printtime( "File out", stdout)
        write(6,*) nbasis * nband * nspin * kpt%list%nk / 64

        call MPI_FILE_CLOSE( fhtmels, ierr )

        call MPI_TYPE_FREE( file_type, ierr )

      else


        fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
        fmode = IOR(fmode, MPI_MODE_UNIQUE_OPEN)
        call MPI_FILE_OPEN( cross_pool_comm, 'ptmels.dat', fmode, MPI_INFO_NULL, fhtmels, ierr )
        if( ierr/=0 ) then
          errorcode=ierr
          call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
          call errore("Open ptmels.dat",string,errorcode)
        endif
        offset=0
        !JTV At this point it would be good to create a custom MPI_DATATYPE
        !  so that we can get optimized file writing
        call MPI_FILE_SET_VIEW( fhtmels, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, &
                                'native', MPI_INFO_NULL, ierr )
        if( ierr/=0 ) then
          errorcode=ierr
          call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
          call errore("Set view fheig",string,errorcode)
        endif

        nelements = val_band * nband

        pool_ik = 0
        do ispin = 1, nspin
          do ik = 1, kpt%list%nk

            if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle

            offset = int( ik-1 + (ispin-1)*kpt%list%nk, MPI_OFFSET_KIND ) * int( nband, MPI_OFFSET_KIND ) & 
                   * int( val_band, MPI_OFFSET_KIND )
            pool_ik = pool_ik + 1

            call MPI_FILE_WRITE_AT( fhtmels, offset, tmels(1,1,pool_ik), nelements, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )

          enddo
        enddo
        call MPI_FILE_CLOSE( fhtmels, ierr )
      endif
    endif
    
    if( ionode ) then

      open(unit=iuntmp,file='tmels.info',form='formatted',status='unknown')
      rewind(iuntmp)
      write(iuntmp,*) val_band, brange(3), brange(4), kpt%list%nk
      close(iuntmp)

      open(unit=iuntmp,file='tmels_partial.txt',form='formatted',status='unknown')
      rewind(iuntmp)
      write(iuntmp,*) tmels(:,:,1)
      close(iuntmp)
    endif

    
  endif



  if( mypoolid==mypoolroot ) then
    call MPI_FILE_CLOSE( fhu2, ierr )
    call MPI_FILE_CLOSE( fhtmels, ierr )
    call MPI_FILE_CLOSE( fh_val_energies, ierr )
    call MPI_FILE_CLOSE( fh_con_energies, ierr )

    call MPI_FILE_CLOSE( fh_val_eigvecs, ierr )
    call MPI_FILE_CLOSE( fh_con_eigvecs, ierr )
  endif

111 continue 
  write(stdout,*) ' end ocean_qdiag'
  call mp_end
  stop
  
  end program ocean_qdiag
