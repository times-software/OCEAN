  program OCEAN_builder

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Nov 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
  use parallel_include
  use hamq_shirley
  use scalapack_module, only : CTXT_, N_, NB_, DLEN_
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, root
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use mpio
  use kpt_module
  use corerepair_module
  use shirley_input_module 
  use diag_module
!  USE wavefunctions_module, ONLY: evc
  USE io_files, ONLY: nwordwfc, prefix_io=>prefix, &
                                  tmpdir_io=>tmp_dir, nd_nmbr, diropn
  USE control_flags,        ONLY : lscf
  USE basis,                ONLY : starting_pot, starting_wfc
  USE gvect
!  USE wvfct, only : npwx, npw
  use constants, only : pi
  use hamq_pool, only : nproc_per_pool, npool, &
                        mypool, rootpool, mypoolid, mypoolroot, &
                        cross_pool_comm, intra_pool_comm, &
                        desc_cyclic, context_cyclic, &
                        local_striped_dim, cyclic_localindex
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

!  integer(kind=MPI_OFFSET_KIND) :: offset
!  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  character(20 ) :: ximat_name, nopt_name

!  character(len=6) :: nd_nmbr_tmp


  integer :: ik
  integer :: i,j
  integer :: nbasis_subset

  integer :: nk, nbnd, ispin, itau
  real(dp) :: nelec_, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin
  real(dp) :: fermi_energy
  complex(dp),allocatable :: ztmp(:,:)

  integer :: iuntmp, nptot, ibd, n_se, iunrbf
  real(dp),allocatable :: tau(:,:), se_list(:)
  logical :: se_exist

  complex(dp),allocatable :: gre(:,:,:), bofr( :, : ), eikr( : ), gre_small(:,:,:), single_bofr(:,:)
  complex(dp),allocatable :: gre_local(:,:,:), phased_bofr(:,:), uofrandb(:,:)
  real(dp),allocatable :: posn( :, : ), wpt( : ), drel( : ), t(:), xirow(:), nind(:), vind(:),vipt(:),phase(:)
  real(dp),allocatable :: wgt(:), newwgt(:), full_xi(:,:), xi_local(:,:)
  real(dp) :: pref, spinfac, denr, deni, sigma, su, su2, omega, avec(3,3), bvec(3,3), qin(3),qcart(3)
  real(dp) :: vlev, vhev, clev, chev, mindif, maxdif, absdiff, newdiff, ktmp( 3 ), eshift, shifted_eig
  complex(dp) :: scalar 

  integer :: ntau, ntau_offset
  integer :: nb, nb2, local_npt, local_npt2
  integer :: nprow, npcol, myrow, mycol
  integer, dimension( DLEN_ ) :: desc_bofr_in, desc_bofr, desc_gre, desc_gre_local

  integer :: npt, it, nt, iunbofr, nbnd_small
!  integer :: num_threads, npt_start, npt_left, npt_chunk
  integer, external :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, numroc
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
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

  call diag_init


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
    
    open(unit=iuntmp,file='bofr_tau_control',form='formatted',status='old')
    read(iuntmp,*) ntau, ntau_offset
    close(iuntmp)
    
    iunrbf=97 !freeunit()
    open(unit=iunrbf,file='rbfile.bin',form='unformatted',status='old')
    read(iunrbf) npt
    write(6,*) npt, ntau
    
    allocate( posn(3,npt), wpt(npt), drel(npt) )
    allocate(vipt(npt))

    open( unit=iuntmp, file='avecsinbohr.ipt',form='formatted',status='old')
    read(iuntmp,*) avec(:,:)
    close(iuntmp)
    call getomega( avec,omega)

    open( unit=iuntmp, file='bvecs',form='formatted',status='old')
    read(iuntmp,*) bvec(:,:)
    close(iuntmp)

    open( unit=iuntmp, file='Pquadrature', form='formatted', status='old' )
    rewind iuntmp
    read ( iuntmp, * ) nt
    allocate( t( nt ), wgt( nt ), newwgt( nt ) )
    su = 0
    do it = 1, nt
       read ( iuntmp, * ) t( it ), wgt( it )
       t( it ) = ( 1 + t( it ) ) / 2
       su = su + wgt( it )
    end do
    wgt = wgt / su
    do it = 1, nt
       newwgt( it ) = wgt( it ) / ( 1.d0 - t( it ) ) ** 2
    end do
    close( unit=iuntmp )

    open( unit=iuntmp,file='efermiinrydberg.ipt',form='formatted',status='old')
    read( iuntmp, * ) fermi_energy
    close( iuntmp )

    n_se = 0
    inquire(file='screeningenergies.ipt',exist=se_exist)
    if( se_exist ) then
      open( unit=iuntmp,file='screeningenergies.ipt',form='formatted',status='old')
      read( iuntmp, * ) n_se
      allocate( se_list( n_se ) )
      read( iuntmp, * ) se_list( : )
      se_list( : ) = se_list( : ) /rytoev
      write(stdout,*) se_list(:)
      close( iuntmp )
    endif

    eshift = 0.0d0
    inquire(file='scissor',exist=se_exist)
    if( se_exist ) then
      ! eshift -> eshift / 2 because we are shifting both occ and unocc away from Efermi
      open(unit=iuntmp,file='scissor',form='formatted',status='old')
      read( iuntmp, * ) eshift
      eshift = eshift * 0.5 / rytoev
      close(iuntmp)
      write(stdout,*) eshift
    endif
      
  endif



  ! share grid info
  call mp_bcast( npt, ionode_id )

  call mp_bcast( ntau, ionode_id )

  call mp_bcast( nt, ionode_id )
  if( .not. ionode ) allocate( t( nt ), newwgt( nt ) )
  call mp_bcast( t, ionode_id )
  call mp_bcast( newwgt, ionode_id )
  write(stdout,*) ' nt: ', nt
  call mp_bcast( bvec, ionode_id )
  call mp_bcast( omega, ionode_id )

  call mp_bcast( n_se, ionode_id )
  if( n_se .gt. 0 ) then
    if( .not. ionode ) allocate( se_list( n_se ) )
    call mp_bcast( se_list, ionode_id )
  endif
  call mp_bcast( eshift, ionode_id )

  ! Prep bofr
  call BLACS_GRIDINFO( context_cyclic, nprow, npcol, myrow, mycol )
  write(stdout,*) nprow, npcol
  nb = 900 / nprow
  nb = min( nb, 32 )
  if( nb .lt. 1 ) nb = 1
  local_npt = numroc( npt, nb, myrow, 0, nprow )
  

  call mp_barrier
  write(stdout,*) npt, nbasis, nb, local_npt, desc_cyclic(NB_)
  call descinit( desc_bofr, npt, nbasis, nb, desc_cyclic(NB_), 0, 0, context_cyclic, local_npt, ierr )
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
  else
    allocate( single_bofr( 1, 1 ), phased_bofr( 1, 1 ) )
  endif
    



!!!!!!!!!!!!!!!!!!!!!!!!!

  if( ( .not. ionode ) .and. ( mypoolid .eq. 0 ) ) allocate( posn(3,npt) )


  

  nb2 = 900 /npcol
  nb2 = min( nb2, 32 )
  if( nb2 .lt. 1 ) nb2 = 1
  local_npt2 =  numroc( npt, nb2, mycol, 0, npcol )

  call descinit( desc_gre, npt, npt, nb, nb2, 0, 0, context_cyclic, local_npt, ierr )
  if( ierr .ne. 0 ) stop
  allocate(gre(local_npt,local_npt2, nt), gre_small(local_npt,local_npt2, nt) )
  gre = 0.d0
  gre_small = 0.d0

  call mp_bcast( omega, ionode_id )
  pref = 1.d0 / ( kpt%list%nk * omega )
  write(stdout,*) ' pref: ', pref, kpt%list%nk, omega


  if( nspin .eq. 2 ) then
    spinfac = 1.d0
  else
    spinfac = 2.d0 
  endif

  allocate( phase( npt ), eikr( npt ) )

! ======================================================================
! Need to set up geometric mean thingy for energy integral
! For now do this at Gamma
  call mp_bcast( fermi_energy, ionode_id )
  ktmp( : ) = 0.d0
  call diag_build_hamk( ktmp, .true., 1 )
  call diag_ham
  vhev = eigval( 1 )
  vlev = eigval( 1 )
  chev = eigval( band_subset(2) )
  clev = eigval( band_subset(2) )
  do ibd = band_subset(1),band_subset(2)
    if( eigval( ibd ) .le. fermi_energy ) then
      if( eigval( ibd ) .gt. vhev ) vhev = eigval( ibd )
    else
      if( eigval( ibd ) .lt. clev ) clev = eigval( ibd )
    endif
  enddo
  write(stdout,*) "done with vcbder"
  write(stdout, '(4(1x,1e15.8))' ) vlev*rytoev, vhev*rytoev, clev*rytoev, chev*rytoev
  mindif = min( fermi_energy - vhev, clev - fermi_energy )
  maxdif = max( fermi_energy - vlev, chev - fermi_energy )
  sigma = sqrt( mindif * maxdif )
  write (stdout,'(4(1x,1e15.8))') fermi_energy*rytoev, mindif*rytoev, maxdif*rytoev, sigma*rytoev

  nbnd_small = floor( dble(nbasis_subset) * 0.75 )
  write(stdout,*) 'Nbands        Nbands small'
  write(stdout,*) nbasis_subset , nbnd_small


! =====================================================================

  do itau = 1, ntau
    gre = 0.d0
    gre_small = 0.d0
    write(stdout,*) 'Now treating core :', itau, ' of ', ntau
    ! read in posn
    if( ionode ) then
      read(iunrbf) posn
      read(iunrbf) wpt
      read(iunrbf) drel
    endif
    if( mypoolid .eq. 0 ) then
      call mp_bcast( posn, ionode_id, cross_pool_comm )
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
    

    call OCEAN_t_reset 

! ======================================================================
    do ispin=1,nspin

    do ik=1,kpt%list%nk
! ======================================================================

      if( mod(ik-1,npool)/=mypool ) cycle

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
        call diag_ham
      endif

      write(stdout,'(a10,x,a10,x,a10,x,a10)') 'K-point', 'Min (eV)', 'Max (eV)', 'Small (eV)'
      write(stdout,'(i10,x,f10.3,x,f10.3,x,f10.3)') ik, eigval(band_subset(1))*rytoev, &
                              eigval(band_subset(2))*rytoev, eigval(nbnd_small)*rytoev

      qin( : ) =  kpt%list%kvec( 1 : 3, ik ) 
      qcart(:) = 0.d0
      if( kpt%param%cartesian ) then
        qcart(:) = qin(:)*tpiba
      else
        do i=1,3
          qcart(:) = qcart(:) + bvec(:,i)*qin(i)
        enddo
      endif
  !    write(stdout,*) qin(:)
  !    write(stdout,*) qcart(:)
  !    write(stdout,*) kpt%list%kvec( 1 : 3, ik ) * tpiba


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
      call PZGEMR2D( npt, nbasis, phased_bofr, 1, 1, desc_bofr_in, &
                     bofr, 1, 1, desc_bofr, context_cyclic )
      

      call PZGEMM( 'N', 'N', npt, nbasis, nbasis, &
                   one, bofr, 1, 1, desc_bofr, eigvec, 1, 1, desc_cyclic, &
                   zero, uofrandb, 1, 1, desc_bofr )


      do ibd = band_subset(1),band_subset(2)

        ! We shift conduction up AND valence down to get away from eFermi
        if( eigval( ibd ) .gt. fermi_energy ) then
          shifted_eig = eigval( ibd ) + eshift
          absdiff = shifted_eig - fermi_energy
          denr = -sqrt( absdiff**2 + 1.0d-12 )
        else
          shifted_eig = eigval( ibd ) - eshift
          absdiff = shifted_eig - fermi_energy
          denr = sqrt( absdiff**2 + 1.0d-12 )
        endif


      ! I think we may be bandwidth limited
        do it = 1, nt
          deni = sigma * t( it ) / ( 1.0_dp - t( it ) )

  !          absdiff = abs( fermi_energy - shifted_eig )
  !          newdiff = sqrt( absdiff**2 + 1.d0*10**(-12) )
  !          denr = sign( newdiff, fermi_energy - shifted_eig )
          scalar = pref / cmplx( denr, deni )


          call PZGERC( npt, npt, scalar, uofrandb, 1, ibd, desc_bofr, 1,  &
                                         uofrandb, 1, ibd, desc_bofr, 1,  &
                                         gre(1,1,it), 1, 1, desc_gre )

          if( ibd .lt.  nbnd_small + band_subset(1)-1 ) then 
            call PZGERC( npt, npt, scalar, uofrandb, 1, ibd, desc_bofr, 1, &
                                           uofrandb, 1, ibd, desc_bofr, 1, &
                                           gre_small(1,1,it), 1, 1, desc_gre )
          endif
        enddo

      enddo ! ibd


      deallocate( uofrandb, bofr )

    enddo ! ik

    enddo ! ispin

    call mp_barrier
    call OCEAN_t_printtime( 'K-point loop', stdout )
    write(stdout,*) 'Loop over k-points and bands complete.'
    write(stdout,*) 'Gather chi to ionode'
 
    call descinit( desc_gre_local, npt, npt, npt, npt, 0, 0,  context_cyclic, npt, ierr )

    allocate( full_xi( local_npt, local_npt2 ) )
    allocate( vind( npt ), nind( npt ), xirow( npt ) )
    if( mypoolid .eq. 0 ) then
      allocate( xi_local( npt, npt ) )
    else
      allocate( xi_local( 1, 1 ) )
    endif
    full_xi = 0.0_DP

    ! now each pool has gre( npt_local, npt_local, nt )
    ! better feature is non-blocking MPI calls
    ! then in the loop below check for the finished call before each loop?
    ! Almost certqainly not worth it, but, would be good to call non-blocking to share
    !  gre_small
    call OCEAN_t_reset
    call mp_sum( gre, cross_pool_comm )
    call OCEAN_t_printtime( 'GRE share', stdout )

    ! later have each pool work on subset of i ?
    if( mypool .eq. 0 ) then
      do it = 1, nt
        su = newwgt( it ) * 4.0_DP * sigma / pi ! 2 for spin 2 for Ry
        do j = 1, local_npt2
          do i = 1, local_npt
            full_xi( i, j ) = full_xi( i, j ) &
                            + su * (real(gre( i, j, it )) ** 2 - aimag( gre( i, j, it ) )**2 )
          enddo
        enddo
      enddo


      call PDGEMR2D( npt, npt, full_xi, 1, 1, desc_gre, & 
                     xi_local, 1, 1, desc_gre_local, desc_gre_local(CTXT_) )

      if( mypoolid .eq. 0 ) then
        write(ximat_name,'(a5,i4.4)') 'ximat', itau+ntau_offset
        open(unit=99,file=ximat_name,form='unformatted', status='unknown' )
        rewind(99)

        vind = 0.0_DP
        do i = 1, npt
          su2 = 0
          do j = 1, npt
            xirow( j ) = xi_local( i, j )
            su2 = su2 + vipt( j ) * wpt( j )
          enddo
          write(99) xirow

          do j = 1, npt
             vind( j ) = vind( j ) + wpt( i ) * su2 / max( drel( j ), drel( i ) )
          enddo
          nind( i ) = su2
        enddo
        close( unit=99 )
        !
        write(6,*)"starting nopt"
        write(nopt_name,'(a4,i4.4)') 'nopt', itau+ntau_offset
        open( unit=99, file=nopt_name, form='formatted', status='unknown' )
        rewind 99
        do i = 1, npt
           write ( 99, '(4(1x,1e15.8))' ) drel( i ), vipt( i ), nind( i ), vind( i )
        end do
        close( unit=99 )
      endif
    endif
    !
    ! this is slow
    call mp_sum( gre_small, cross_pool_comm )
    full_xi = 0.0_DP
    ! later have each pool work on subset of i ?
    if( mypool .eq. 0 ) then
      do it = 1, nt
        su = newwgt( it ) * 4.0_DP * sigma / pi ! 2 for spin 2 for Ry
        do j = 1, local_npt2
          do i = 1, local_npt
            full_xi( i, j ) = full_xi( i, j ) &
                  + su * (real(gre_small( i, j, it )) ** 2 - aimag( gre_small( i, j, it ) )**2 )
          enddo
        enddo
      enddo

      call PDGEMR2D( npt, npt, full_xi, 1, 1, desc_gre, &
                     xi_local, 1, 1, desc_gre_local, desc_gre_local(CTXT_) )

      if( mypoolid .eq. 0 ) then
        write(ximat_name,'(a11,i4.4)') 'ximat_small', itau+ntau_offset
        open(unit=99,file=ximat_name,form='unformatted', status='unknown' )
        rewind(99)

        vind = 0.0_DP
        do i = 1, npt
          su2 = 0
          do j = 1, npt
            xirow( j ) = xi_local( i, j )
            su2 = su2 + vipt( j ) * wpt( j )
          enddo
          write(99) xirow

          do j = 1, npt
             vind( j ) = vind( j ) + wpt( i ) * su2 / max( drel( j ), drel( i ) )
          enddo
          nind( i ) = su2
        enddo
        close( unit=99 )
        !
        write(6,*)"starting nopt"
        write(nopt_name,'(a11,i4.4)') 'nopt_small', itau+ntau_offset
        open( unit=99, file=nopt_name, form='formatted', status='unknown' )
        rewind 99
        do i = 1, npt
           write ( 99, '(4(1x,1e15.8))' ) drel( i ), vipt( i ), nind( i ), vind( i )
        end do
        close( unit=99 )
      endif
    endif

    deallocate( full_xi, vind, nind, xirow, xi_local )

  enddo ! itau

  111 continue 
  write(stdout,*) ' end OCEAN_builder'
  call mp_end
!  stop
  
end program OCEAN_builder
