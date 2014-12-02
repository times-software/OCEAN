  program shirley_qdiag

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Nov 2007

  ! now parallelized
#include "f_defs.h"
  use kinds, only : dp
  use parallel_include
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_global, only : nproc, mpime, root
  use mp, only : mp_bcast, mp_end, mp_barrier, mp_sum
  use mpio
  use kpt_module
  use corerepair_module
  use shirley_input_module 
  use diag_module
  USE wavefunctions_module, ONLY: evc
  USE io_files, ONLY: nwordwfc, iunwfc, prefix_io=>prefix, &
                                  tmpdir_io=>tmp_dir, nd_nmbr, diropn, wfc_dir
  USE control_flags,        ONLY : lscf
  USE basis,                ONLY : starting_pot, starting_wfc
  USE gvect
  USE wvfct, only : npwx, npw
  use scalapack_module, only : DLEN_, M_, N_, MB_, NB_
  use hamq_pool, only : nproc_per_pool, npool, &
                        mypool, rootpool, mypoolid, mypoolroot, &
                        cross_pool_comm, intra_pool_comm, &
                        desc_cyclic, desc_striped, context_cyclic, &
                        local_striped_dim, cyclic_localindex, &
                        context_global, local_cyclic_dims
 use OCEAN_timer
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


  character(255) :: eigval_file, eigvec_file, info_file, filout
  logical :: ex
  logical :: legacy_zero_lumo

  integer(kind=MPI_OFFSET_KIND) :: offset, off1, off2, off3, off4
  integer :: fheigval, fheigvec, fhenk, fhpsi, fheig
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  character(255) :: fmtstr, string
  character(20 ) :: filnam

  real(dp) :: lambda
  character(len=6) :: nd_nmbr_tmp


  real(dp) :: kvec(3)
  real(dp) :: qpathlen, cqvec(3), dqvec(3)
  character(255) :: ic
  character(6), allocatable :: fntau(:)
  character (len=5) :: band_style = 'defau'

  integer :: ik
  integer :: i,j,k, ig
  integer :: nbasis_subset

  integer :: nk, nbnd, ispin, itau, ip, ibnd
  real(dp) :: nelec_, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin, ix, nxpts, xmesh( 3 )
  real(dp) :: fermi_energy
  real(dp),allocatable :: kr( :, : )
  integer,allocatable :: gflip( :, : )
  complex(dp),allocatable :: ztmp(:,:)

  integer :: nwordo2l, iuntmp, ntot,nptot, ntau, ibd, ibp, iunbofr, tmp_eig_m, tmp_eig_n 
  complex(dp),allocatable :: o2l(:,:,:,:,:), coeff( :, :, :, :, :, : ), coeff_small( :, :, :, :, : )
  real(dp),allocatable :: tau(:,:),e0(:,:,:,:), e0_small(:,:,:)
  logical :: exst

  complex(dp),allocatable :: u1(:,:), u2(:,:,:,:), tmp_eigvec(:,:), tmels(:,:,: ), u1_single(:,:)
  complex(dp) :: x, w

  integer,allocatable :: start_band(:),stop_band(:), new_start_band(:,:,:)
  real(dp),allocatable :: lumo(:), homo(:)
  integer :: nbuse, brange( 4 ), nelectron, nbuse_xes, max_val, ishift
  real(dp) :: lumo_point, homo_point,dft_energy_range, cs, vs, newgap, sef, egw, elda, eshift, kshift( 3 ), kplusq(3), lumo_shift
  logical :: have_gwipt, have_kshift, radial_ufunc

  integer :: iuninf, fhu2, finfo, errorcode, fmode, resultlen, nshift, fhtmels, fh_val_energies, fh_con_energies, fho2l
  integer :: fh_val_eigvecs, fh_con_eigvecs, ibeg_unit

  integer :: nr, nr2, ir, nang, nang2, iang, lmax, nthreads, nband, npt
  real(dp) :: qin(3), qcart(3), deltaL, rmax, rmin, bvec(3,3), norm
  real(dp), allocatable ::  rad_grid(:), angular_grid(:,:), rpts(:,:,:), posn(:,:), wpt(:), bofr(:,:), phased_bofr(:,:),phase(:)
  complex(dp), allocatable :: BofRad(:,:,:), u_radial(:,:,:,:), eikr(:), eigvec_single(:,:)

  integer :: u1_MB, u1_M, u1_NB, u1_N
  integer :: nprow, npcol, myrow, mycol


  integer :: nr_eigvec, nc_eigvec, pool_tot_k, pool_ik, eigvec_info

  integer :: block_start, block_stop, block_width, block, niter, iter
  integer(8) :: long_nbasis, long_nx
  integer(8) :: long_2g = 2147483648

  integer, dimension( DLEN_ ) :: desc_coeff, desc_o2l, desc_u1, desc_u2, desc_tmels, single_desc_u1, desc_eigvec_single, desc_tmpeig

  integer, external :: OMP_GET_NUM_THREADS, NUMROC
  complex(dp), external :: ZDOTC
  real(dp), external :: DZNRM2


  integer :: dims(2), ndims, array_of_gsizes(2), array_of_distribs(2), array_of_dargs(2), file_type, nelements, mpistatus, &
             locsize, i_energy_request, n_energy_request
  integer, allocatable :: energy_request(:), pool_root_map(:)
  complex(dp), allocatable :: store_eigvec(:,:,:), out_eigvec(:,:,:)
 
  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u


! Read in the input file, shirley_input sets up what 
!  is needed for the interpolation
  call shirley_input


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


! At this point we can run the OBF2cks



  if( band_subset(1) < 1 ) band_subset(1)=1
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

!  call diagx_init( band_subset(1), band_subset(2) )
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

  ! MPI-IO
  if( mypoolid==mypoolroot ) then
    eigval_file=trim(outfile)//'.eigval'
    inquire(file=trim(eigval_file),exist=ex)
    if( ex ) then
      ! delete pre-existing files
 !     fheigval=freeunit()
 !     open(fheigval,file=trim(eigval_file),form='unformatted')
!      close(fheigval,status='delete')
    endif

    call mp_file_open_dp( eigval_file, fheigval, rootpool, cross_pool_comm )

    if( eigvec_output ) then
      eigvec_file=trim(outfile)//'.eigvec'
      call mp_file_open_dp( eigvec_file, fheigvec, rootpool, cross_pool_comm )
    endif


    call MPI_INFO_CREATE( finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    call MPI_FILE_OPEN( cross_pool_comm, 'u2par.dat', fmode, MPI_INFO_NULL, fhu2, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    !JTV At this point it would be good to create a custom MPI_DATATYPE
    !  so that we can get optimized file writing
    call MPI_FILE_SET_VIEW( fhu2, offset, MPI_DOUBLE_COMPLEX, &
                          MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

    call MPI_FILE_OPEN( cross_pool_comm, 'val_energies.dat', fmode, finfo, fh_val_energies, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    call MPI_FILE_SET_VIEW( fh_val_energies, offset, MPI_DOUBLE_PRECISION, &
                          MPI_DOUBLE_PRECISION, 'native', finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

    call MPI_FILE_OPEN( cross_pool_comm, 'con_energies.dat', fmode, finfo, fh_con_energies, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    call MPI_FILE_SET_VIEW( fh_con_energies, offset, MPI_DOUBLE_PRECISION, &
                          MPI_DOUBLE_PRECISION, 'native', finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

    call MPI_FILE_OPEN( cross_pool_comm, 'val_eigvecs.dat', fmode, finfo, fh_val_eigvecs, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    call MPI_FILE_SET_VIEW( fh_val_eigvecs, offset, MPI_DOUBLE_COMPLEX, &
                          MPI_DOUBLE_COMPLEX, 'native', finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

    call MPI_FILE_OPEN( cross_pool_comm, 'con_eigvecs.dat', fmode, finfo, fh_con_eigvecs, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    call MPI_FILE_SET_VIEW( fh_con_eigvecs, offset, MPI_DOUBLE_COMPLEX, &
                          MPI_DOUBLE_COMPLEX, 'native', finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif




    call MPI_INFO_CREATE( finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif

    call MPI_FILE_OPEN( cross_pool_comm, 'ptmels.dat', fmode, finfo, fhtmels, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset = 0
    call MPI_FILE_SET_VIEW( fhtmels, offset, MPI_DOUBLE_COMPLEX, &
                          MPI_DOUBLE_COMPLEX, 'native', finfo, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif


!    ! Need to add in UW-BSE type files
!    inquire(file='enkfile',exist=ex)
!    if(ex) then
!      ! delete pre-existing file
!      fhenk=freeunit()
!      open(fhenk,file='enkfile')
!      close(fhenk,status='delete')
!    endif

    ! Open enkfile as native, formatted file    
!!    fhenk=freeunit()
!!    open(unit=fhenk,file='enkfile',form='formatted',status='unknown')
!!    rewind(fhenk)

    
  endif


  if( ionode ) then
    iuntmp = freeunit()

    inquire(file='noshift_lumo',exist=ex)
    if( ex ) then
      open(unit=iuntmp,file='noshift_lumo',form='formatted',status='old')
      read(iuntmp,*) legacy_zero_lumo
      close(iuntmp)
    else
      legacy_zero_lumo = .true.
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


    call seqopn(iuntmp, 'o2li', 'unformatted', exst)
    if( .not. exst ) call errore('ocean', 'file does not exist', 'o2li')
    read(iuntmp) nptot, ntau
    allocate( tau( 3, ntau ) )
    allocate( fntau( ntau ) )
    read(iuntmp) tau( :, : )
    read(iuntmp) fntau( : )
    close(iuntmp)
!test
!    nwordo2l = 2 * nbasis * nptot * ntau
    write(stdout,*) nbasis, nptot, ntau
!test    allocate( o2l( nbasis, nptot, ntau, kpt%list%nk, nshift ) )
!    allocate( o2l( nbasis, nptot, ntau, 1, 1 ) )


    nang = 0
    nr = 0

    iuntmp = freeunit()
    open(unit=iuntmp,file='xmesh.ipt',form='formatted',status='old')
    read(iuntmp,*) xmesh(:)
    close(iuntmp)
    nxpts = product( xmesh )


    open( unit=iuntmp, file='bvecs',form='formatted',status='old')
    read(iuntmp,*) bvec(:,:)
    close(iuntmp)


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
    
    allocate( u1_single( nxpts, nbasis ) )
    iuntmp = freeunit()
    open(unit=iuntmp,file='u1.dat',form='unformatted', status='old')
    rewind(iuntmp)
    do i = 1, nbasis
      read(iuntmp) u1_single( :, i )
    enddo
    close(iuntmp)
!    u1( :, : ) = conjg( u1( :,: ) )


    open(unit=iuntmp,file='nelectron',form='formatted',status='old')
    read(iuntmp,*) nelectron
    close(iuntmp)

  endif

  call mp_bcast( nang, ionode_id )
  call mp_bcast( nr, ionode_id )
  radial_ufunc = .false. 
  
  call mp_bcast( legacy_zero_lumo, ionode_id )
  call mp_bcast( kshift, ionode_id )
  call mp_bcast( have_kshift, ionode_id )
  call mp_bcast( nshift, ionode_id )

  call mp_bcast( fermi_energy, ionode_id )
  call mp_bcast( dft_energy_range, ionode_id)

  call mp_bcast( nptot, ionode_id )
  call mp_bcast( ntau, ionode_id )
!test
!  if( .not. ionode ) allocate( o2l( nbasis, nptot, ntau, kpt%list%nk, nshift ) )
!  call mp_bcast( o2l, ionode_id )
  call mp_barrier
  write(stdout,*) ' shared o2l array'
  

!test
  nwordo2l = 2 * nbasis * nptot * ntau
  if( mypoolid == mypoolroot ) then
    allocate( o2l(  nbasis, nptot, ntau, 1, 1 ) )
  else
    allocate( o2l( 1, 1, ntau, 1, 1 ) )
  endif
  fho2l = freeunit()
  ! since only the ionode is going to write to this file
  ! there is no need for the nd_nmbr suffix
  nd_nmbr_tmp = nd_nmbr
  nd_nmbr=''
  call diropn( fho2l, 'o2l', nwordo2l, exst )
  ! restore in case needed again
  nd_nmbr = nd_nmbr_tmp
  write(stdout,*) ' reading from file: ', trim(tmpdir_io)//trim(prefix)//'.o2l'
  if( .not. exst ) write(stdout,*) exst
!  do ishift = 1, nshift
!    do ik = 1,kpt%list%nk
!      call davcio( o2l( 1, 1, 1, ik, ishift ), nwordo2l, fho2l, ik + (ishift-1)*kpt%list%nk, -1 )
!    enddo
!  enddo
!  close( iuntmp )


  call descinit( desc_o2l, nbasis, nptot, nbasis, nptot, 0, 0, context_cyclic, nbasis, ierr )

  call mp_bcast(nxpts, ionode_id )

  if( mypoolid .eq. mypoolroot ) then
    if( .not. ionode ) allocate( u1_single( nxpts, nbasis ) )
    call mp_bcast( u1_single, ionode_id, cross_pool_comm )
  else
    allocate( u1_single( 1, 1 ) )
  endif

  call BLACS_GRIDINFO( context_cyclic, nprow, npcol, myrow, mycol )
  u1_MB = min( 32, nxpts / nprow )
  if( u1_MB .lt. 1 ) u1_MB = 1
  u1_M = numroc( nxpts, u1_MB, myrow, 0, nprow )


!  u1_NB = min( 32, nbasis / npcol )
!  if( u1_NB .lt. 1 ) u1_NB = 1
!  u1_N = numroc( nbasis, u1_NB, mycol, 0, npcol )
  u1_NB = desc_cyclic( NB_ )
  u1_N = desc_cyclic( N_ )

  if( u1_N .lt. 1 ) u1_N = 1


  write(stdout,*) u1_M, u1_N, u1_MB, u1_NB, nxpts, nbasis
 

  call descinit( single_desc_u1, nxpts, nbasis, nxpts, nbasis, 0, 0, context_cyclic, nxpts, ierr )
  call descinit( desc_u1, nxpts, nbasis, u1_MB, u1_NB, 0, 0, context_cyclic, nxpts, ierr )

  allocate( u1( u1_M, u1_N ) )


! !!!!!!!
! Currently there is a bug in SCALAPACK
! This matrix cannot be > 2GB

! No idea why, but I'm crashing WAY short of this limit
  long_nx = nxpts
  long_nbasis = nbasis
  long_nx = long_nx * long_nbasis

  long_2g = long_2g / 64

  if( long_nx .ge. long_2g ) then
    niter = 1 + ( long_nx/long_2g )
    block = nbasis / niter
    

    write(stdout,*) 'Breaking share into sections', niter
    do iter = 1, niter
      write(stdout,*) iter
      block_start = 1 + (iter-1) * block
      block_stop = min( block_start + block - 1, nbasis )
      block_width = block_stop - block_start + 1
      call mkl_free_buffers()
      write(stdout,*) nxpts, block_width, nxpts * block_width * 16

      call PZGEMR2D( nxpts, block_width, u1_single, 1, block_start, single_desc_u1, u1, 1, &
                     block_start, desc_u1, context_cyclic, ierr )
    enddo
  
  else 
  
    call mkl_free_buffers()
    call PZGEMR2D( nxpts, nbasis, u1_single, 1, 1, single_desc_u1, u1, 1, 1, desc_u1, context_cyclic, ierr )

  endif

  deallocate( u1_single )


  call mp_bcast( nelectron, ionode_id )


  ntot = nbasis_subset * kpt%list%nk

  nband = band_subset(2) - band_subset(1) + 1
  if( mypoolid .eq. mypoolroot ) then
    allocate( coeff( nptot, band_subset(1):band_subset(2), kpt%list%nk, nspin, ntau, nshift ) )
  else ! for error checking nonsense
    allocate( coeff( 1, band_subset(1), kpt%list%nk, nspin, ntau, nshift ) )
  endif
  coeff = 0.d0
  call descinit( desc_coeff, nptot, nband, nptot, nband, 0, 0, context_cyclic, nptot, ierr )
  if( ierr .ne. 0 ) then
    stop
  endif

  call descinit( desc_eigvec_single, nbasis, nbasis, nbasis, nbasis, 0, 0, context_cyclic, nbasis, ierr )
  if( mypoolid .eq. mypoolroot ) then
!???    allocate( u2( nxpts, band_subset(1) : band_subset(2), kpt%list%nk, nshift ) )
    allocate( u2( nxpts, band_subset(1) : band_subset(2), 1, nshift ) )
    allocate( eigvec_single( nbasis, nbasis ) )
  else  ! For error checking nonsense
!!!!    allocate( u2( 1, band_subset(1), kpt%list%nk, nshift ) )
     ! JTV memleak 29 July
!    allocate( u2( 1, band_subset(1), 1, nshift ) )
    allocate( u2( nxpts, band_subset(1) : band_subset(2), 1, nshift ) )
    allocate( eigvec_single( 1, 1:band_subset(1) ) )
  endif
  u2 = 0.d0
  call descinit( desc_u2, nxpts, nband, nxpts, nband, 0, 0, context_cyclic, nxpts, ierr )

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
  if( have_kshift ) then
    allocate( tmels( max_val, band_subset( 1 ) : band_subset( 2 ), kpt%list%nk ) )
    tmels = 0.0d0
!    allocate( tmp_eigvec( desc_cyclic(M_), desc_cyclic(N_) ) )
    call descinit( desc_tmels,  max_val, band_subset( 2 )-band_subset( 1 )+1, &
                                max_val, band_subset( 2 )-band_subset( 1 )+1, &
                   0, 0, context_cyclic, max_val, ierr )


    tmp_eig_n = numroc( max_val, desc_cyclic(NB_), mycol, 0, npcol )
    tmp_eig_m = numroc( nbasis, desc_cyclic(MB_), myrow, 0, nprow )
    call descinit( desc_tmpeig, nbasis, max_val, desc_cyclic(MB_), desc_cyclic(NB_), 0, 0, context_cyclic, tmp_eig_m, ierr )
    if( ierr .ne. 0 ) then
      write(stdout,*) desc_tmpeig(:)
      stop
    endif


    allocate( tmp_eigvec(  tmp_eig_m, tmp_eig_n ) )

!    allocate( tmp_eigvec( nbasis, max_val ) )
  endif



  ! ========= pre-fetch all of the energies? =========== !
  ! Determine number of k-points each processor's pool will do
  ! Get an array of the ids for each pool head
  ! Energy_request is sequential for ionode and by kpt/spin for everyone else
  allocate( energy_request( kpt%list%nk * nspin ), pool_root_map( 0:npool-1 ) )

  pool_root_map(:) = 0
  if( mypoolid .eq. mypoolroot ) then
    pool_root_map(mypool) = mpime
    call MPI_ALLREDUCE( MPI_IN_PLACE, pool_root_map, npool, MPI_INTEGER, MPI_SUM, cross_pool_comm, ierr )
  endif

!  if( have_kshift ) then
!    nband = max_val
!  else
    nband = band_subset(2) - band_subset(1) + 1
!  endif

  pool_ik = 0
  i_energy_request = 0
  do ispin=1,nspin
    do ik=1,kpt%list%nk
      if( ionode .and. (mod((ik-1)+(ispin-1)*kpt%list%nk,npool) .ne. mypool ) ) then
        i_energy_request = i_energy_request + 1
        call MPI_IRECV( e0( 1, ik, ispin, 1 ), nband, MPI_DOUBLE_PRECISION,  &
                        pool_root_map(mod((ik-1)+(ispin-1)*kpt%list%nk,npool)), ik+(ispin-1)*kpt%list%nk, &
                        cross_pool_comm, energy_request( i_energy_request ), ierr )
        write(stdout,*) pool_root_map(mod((ik-1)+(ispin-1)*kpt%list%nk,npool)), ik+(ispin-1)*kpt%list%nk
      endif
      if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle
      pool_ik = pool_ik + 1
    enddo
  enddo
  pool_tot_k = max( pool_ik, 1 )
  n_energy_request = i_energy_request


  ! Allocate storage for eigenvectors
  ! Create descriptors  
  ! For laziness start by copying *all* of eigvector which is way larger than we need 
  call local_cyclic_dims( nr_eigvec, nc_eigvec )
  allocate( store_eigvec( nr_eigvec, nc_eigvec, pool_tot_k ) )
  



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

  !    call diagx_ham
        call diag_ham

      endif


      write(stdout,*) ik,  eigval(band_subset(1))*rytoev, eigval(band_subset(2))*rytoev
      e0( :, ik, ispin, 1 ) = 0.5d0 * eigval(band_subset(1):band_subset(2))


      ! Send energies to ionode
      if( .not. ionode .and. mypoolid .eq. 0 ) then
        if( have_kshift ) then 
          nband = max_val
        else
          nband = band_subset(2) - band_subset(1) + 1
        endif
   
        write(stdout,*) 
        call MPI_ISEND( e0( :, ik, ispin, 1 ), nband, MPI_DOUBLE_PRECISION, ionode_id, ik+(ispin-1)*kpt%list%nk, &
                        cross_pool_comm, energy_request( ik+(ispin-1)*kpt%list%nk ), ierr )
      endif



      !!!!!!! store away eigenvectors
      store_eigvec( :, :, pool_ik ) = eigvec( :, : )

    enddo ! ik
  enddo ! ispin
  ! done with interpolation here

!#ifdef __NIST
!  call diag_free
!#endif



  ! Clean up energy communications and find fermi, lumo, homo 

  ! Find start bands and save out energyfile

  if( ionode ) then
    write(stdout,*) 'Sharing energies'

    call MPI_WAITALL( n_energy_request, energy_request, MPI_STATUSES_IGNORE, ierr )

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

    
    call dump_energies( band_subset, nband, kpt%list%nk, nspin, nshift, e0, lumo_shift, new_start_band, ierr )

    write(stdout,*) 'Sharing energies'


  endif

  if( mypoolid .eq. 0 ) then
    if( .not. ionode ) then
      do ispin=1,nspin
        do ik=1,kpt%list%nk
          if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle
!          write(stdout,*)  pool_root_map(mod((ik-1)+(ispin-1)*kpt%list%nk,npool)), ik+(ispin-1)*kpt%list%nk, mpime 
!          flush(stdout)
          call MPI_WAIT( energy_request( ik+(ispin-1)*kpt%list%nk ), MPI_STATUS_IGNORE, ierr )
        enddo
      enddo
    endif
    call MPI_BCAST( fermi_energy, 1, MPI_DOUBLE_PRECISION, ionode_id, cross_pool_comm, ierr )
  endif
  call MPI_BCAST( fermi_energy, 1, MPI_DOUBLE_PRECISION, mypoolroot, intra_pool_comm, ierr )

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


  nband = band_subset(2) - band_subset(1) + 1
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


      call PZGEMR2D( nbasis, nband, store_eigvec( 1, 1, pool_ik ), 1, 1, desc_cyclic, &
                                      out_eigvec( 1, 1, pool_ik ), 1, 1, desc_eigvec_single, &
                     context_cyclic, ierr )
    enddo
  enddo

  if( mypoolid .eq. mypoolroot ) then 
  
    ndims = 2
    dims = (/ 1, npool /)
    array_of_gsizes = (/ nbasis * nband, nspin * kpt%list%nk /)
    array_of_distribs = (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC /)
    array_of_dargs = (/ nbasis * nband, 1 /)
    CALL MPI_TYPE_CREATE_DARRAY( npool, mypool, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, &
                                 dims, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, file_type, ierr )
    CALL MPI_TYPE_COMMIT( file_type, ierr )

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
      call errore(string,errorcode)
    endif
    offset=0
    !JTV At this point it would be good to create a custom MPI_DATATYPE
    !  so that we can get optimized file writing
    call MPI_FILE_SET_VIEW( fheig, offset, MPI_DOUBLE_COMPLEX, &
                            file_type, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
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


  endif

  if( ionode ) then
    iuntmp = freeunit()
    open(iuntmp,file='qdiag.info',form='formatted',status='unknown')
    write(iuntmp,*) nbasis, nband, kpt%list%nk, nspin, nshift
    close(iuntmp)
  endif



#ifdef FALSE

  if( mypoolid .eq. mypoolroot ) then

    call MPI_FILE_OPEN( cross_pool_comm, 'eigvecs.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fheig, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    call MPI_FILE_SET_VIEW( fheig, offset, MPI_DOUBLE_COMPLEX, &
                            MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif


  endif
  
  pool_ik = 0
  do ispin=1,nspin
    do ik=1,kpt%list%nk
! ======================================================================

      ! (ik-1)+(ispin-1)*kpt%list%nk
      if( mod((ik-1)+(ispin-1)*kpt%list%nk,npool)/=mypool ) cycle

      pool_ik = pool_ik + 1

      if( mypoolid .eq. mypoolroot ) then
        offset =  ((ispin-1)*kpt%list%nk + ik-1) * nband * nbasis
        call mpi_file_read_at( fheig, offset, out_eigvec( 1, 1, 1 ), &
                                nband * nbasis, &
                                MPI_DOUBLE_COMPLEX, status, ierr )
      endif
      call PZGEMR2D( nbasis, nband, out_eigvec( 1, 1, 1 ), 1, 1, desc_eigvec_single, &
                     eigvec, 1, 1, desc_cyclic, context_cyclic, ierr )






! Calculate u2, the U functions on the xmesh grid
! u1 and eigvec are distributed
! u2 is localized
    nband = band_subset(2) - band_subset(1) + 1
    call PZGEMM( 'N', 'N', nxpts, nband, nbasis, & 
                 one, u1, 1, 1, desc_u1, &
                 eigvec, 1, 1, desc_cyclic, &
                 zero, u2(1, band_subset(1), 1, 1 ), 1, 1, desc_u2 )
!?????                 zero, u2(1, band_subset(1), ik, 1 ), 1, 1, desc_u2 )

!    call ZGEMM( 'N', 'N', nxpts, nband, nbasis, one, u1, nxpts, eigvec(1,band_subset(1)), nbasis, &
!                zero, u2( 1, band_subset(1), ik, 1 ), nxpts )


!JTV 
! For now to get it working, u2 exists only on the primary node in each pool

! Get eigvec all onto root node for writing
    call PZGEMR2D( nbasis, nbasis, eigvec, 1, 1, desc_cyclic, eigvec_single, 1, 1, & 
                   desc_eigvec_single, context_cyclic, ierr )


    if( mypoolid == mypoolroot ) then
      do ibd = band_subset(1), band_subset( 2 )
!???? u2( ik -> 1 )
        norm =  DZNRM2( nxpts, u2( 1, ibd, 1, 1 ), 1 )
        norm = 1_dp / norm
        call ZDSCAL( nxpts, norm, u2( 1, ibd, 1, 1 ), 1 )

        do ibp = band_subset(1), ibd-1
           w = -ZDOTC( nxpts,  u2( 1, ibd, 1, 1 ), 1, u2( 1, ibp, 1, 1 ), 1 )
           CALL ZAXPY( nxpts, w, u2(1, ibp, 1, 1 ), 1, u2( 1, ibd, 1, 1 ), 1 )
        enddo
        norm =  DZNRM2( nxpts, u2( 1, ibd, 1, 1 ), 1 )
        norm = 1_dp / norm
        call ZDSCAL( nxpts, norm, u2( 1, ibd, 1, 1 ), 1 )
      enddo


      ! write val energies
      offset = ((ispin-1)*kpt%list%nk + ik-1) * max_val
      call mpi_file_write_at( fh_val_energies, offset, eigval(1:max_val), max_val, &
                              MPI_DOUBLE_PRECISION, status, ierr )

      ! write val eigvecs 
      offset = ((ispin-1)*kpt%list%nk + ik-1) * nbasis * max_val
! disabling for now
!      call mpi_file_write_at( fh_val_eigvecs, offset, eigvec_single, nbasis*max_val, &
!                              MPI_DOUBLE_COMPLEX, status, ierr )

    endif

! Calculate the local projectors
    
    do itau=1,ntau
      call PZGEMM( 'T', 'N', nptot, nband, nbasis, &
                   one, o2l(1,1,itau,1,1), 1, 1, desc_o2l, &
                   eigvec(1,band_subset(1)), 1, 1, desc_cyclic, &
                   zero, coeff(1,band_subset(1),ik,ispin,itau,1), 1, 1, desc_coeff )
!      CALL ZGEMM( 'T', 'N', nptot, nband, nbasis, one, o2l(1,1,itau,1,1), nbasis, &
!                    eigvec(1,band_subset(1)), nbasis, zero, coeff(1,band_subset(1),ik,ispin,itau,1), nptot )
    enddo

    if( have_kshift ) then
      if( mypoolid == mypoolroot ) then
        call davcio( o2l, nwordo2l, fho2l, ik + kpt%list%nk, -1 )
      endif

!      tmp_eigvec( :, 1 : max_val ) = conjg( eigvec( :, band_subset( 1 ) : band_subset( 1 ) + max_val - 1) )
!      tmp_eigvec( :, 1 : max_val ) = eigvec( :, band_subset( 1 ) : band_subset( 1 ) + max_val - 1)

!JTV
! Can probably be done with PZGEMR2D
!      tmp_eigvec( :, : ) = eigvec( :, : )
       call PZGEMR2D( nbasis, max_val, eigvec, 1, 1, desc_cyclic, &
                      tmp_eigvec, 1, 1, desc_tmpeig, context_cyclic, ierr )



      ! build the Hamiltonian for this q-point !!! k + q
      kplusq(:) = kpt%list%kvec(1:3,ik) + kshift(:)
 
      call diag_build_hamk( kplusq, kpt%param%cartesian, ispin )

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
      
  !    call diagx_ham
        call diag_ham

      endif
      
      write(stdout,*) ik, eigval(band_subset(1):band_subset(2))*rytoev
      e0( :, ik, ispin, 2 ) = 0.5d0 * eigval(band_subset(1):band_subset(2))
  ! Find the start_band
      do ibd = band_subset(1), band_subset( 2 )
        if( eigval( ibd ) .gt. fermi_energy ) then
          start_band( ik ) = ibd
          lumo( ik ) = eigval( ibd )
          goto 11
        endif
      enddo
11    continue
  !    write(stdout,*) ik, start_band(ik)
    

    
!      u2( :, band_subset(1) : band_subset(2), ik, 2 ) = matmul( u1, eigvec( :, band_subset(1) : band_subset(2) ) )
      call PZGEMM( 'N', 'N', nxpts, nband, nbasis, &
                 one, u1, 1, 1, desc_u1, &
                 eigvec, 1, 1, desc_cyclic, &
                 zero, u2(1, band_subset(1), 1, 2 ), 1, 1, desc_u2 )
!!!????                 zero, u2(1, band_subset(1), ik, 2 ), 1, 1, desc_u2 )

!!???? ik -> 1 
      if( mypoolid == mypoolroot ) then
        do ibd = band_subset(1), band_subset( 2 )
          norm =  DZNRM2( nxpts, u2( 1, ibd, 1, 2 ), 1 )
          norm = 1_dp / norm
          call ZDSCAL( nxpts, norm, u2( 1, ibd, 1, 2 ), 1 )

          do ibp = band_subset(1), ibd-1
            w = -ZDOTC( nxpts,  u2( 1, ibd, 1, 2 ), 1, u2( 1, ibp, 1, 2 ), 1 )
            CALL ZAXPY( nxpts, w, u2(1, ibp, 1, 2 ), 1, u2( 1, ibd, 1, 2 ), 1 )
          enddo
          norm =  DZNRM2( nxpts, u2( 1, ibd, 1, 2 ), 1 )
          norm = 1_dp / norm
          call ZDSCAL( nxpts, norm, u2( 1, ibd, 1, 2 ), 1 )
        enddo
      endif

      call BLACS_BARRIER( context_cyclic, 'A')

      do itau=1,ntau
        call PZGEMM( 'T', 'N', nptot, nband, nbasis, &
                   one, o2l(1,1,itau,1,1), 1, 1, desc_o2l, &
                   eigvec(1,band_subset(1)), 1, 1, desc_cyclic, &
                   zero, coeff(1,band_subset(1),ik,ispin,itau,2), 1, 1, desc_coeff )
      enddo

      call BLACS_BARRIER( context_cyclic, 'A')
!      do itau=1,ntau
!        do ibnd=band_subset(1),band_subset(2)
!          do ip=1,nptot
!            coeff( ip, ibnd, ik, ispin, itau, 2 ) = & 
!              dot_product(conjg( o2l( 1: nbasis, ip, itau, 1, 1 )), eigvec( 1 : nbasis, ibnd ) ) 
!          enddo
!        enddo
!      enddo

    
!      call PZGEMM( 'C', 'N', max_val, nband, nbasis, &
!                   one, tmp_eigvec, 1, 1, desc_cyclic, &
!                   eigvec, 1, 1, desc_cyclic, &
!                   zero, tmels, 1, 1, desc_tmels )

!      do ibd = band_subset(1), band_subset( 2 )
!        do i = 1, max_val
!          tmels( i, ibd, ik ) = dot_product( eigvec( :, ibd ), tmp_eigvec( :, i ) ) 
!        enddo
!      enddo

      call PZGEMR2D( nbasis, nbasis, eigvec, 1, 1, desc_cyclic, eigvec_single, 1, 1, &
                     desc_eigvec_single, context_cyclic, ierr )
      
      call BLACS_BARRIER( context_cyclic, 'A')

    endif 
    ! write con energies (if shift has new eigenvalues, otherwise it doesn't )
    offset = ((ispin-1)*kpt%list%nk + ik-1) * nbasis_subset


    if( mypoolid==mypoolroot ) then
      call mpi_file_write_at( fh_con_energies, offset, eigval(band_subset(1):band_subset(2)), nbasis_subset,  &
                              MPI_DOUBLE_PRECISION, status, ierr )

      ! write val eigvecs 
      offset = ((ispin-1)*kpt%list%nk + ik-1) * nbasis * nbasis_subset
! disabling for now
!      call mpi_file_write_at( fh_con_eigvecs, offset, eigvec_single(1,band_subset(1)), &
!                              nbasis*nbasis_subset, MPI_DOUBLE_COMPLEX, status, ierr )

!      u2( :, :, ik, : ) = conjg( u2( :,:,ik,:) )
!      tmels( :, :, ik ) = conjg( tmels( :,:,ik) )
      if( nshift .eq. 1 ) then
        offset =  ((ispin-1)*kpt%list%nk + ik-1) * nbasis_subset * nxpts
!????? ik -> 1
        call mpi_file_write_at( fhu2, offset, u2(1,1,1,1), & 
                                nbasis_subset * nxpts, &
!                                (band_subset(2) - band_subset(1) + 1 ) * nxpts, &
                                MPI_DOUBLE_COMPLEX, status, ierr )
      else
        off1 = ((ispin-1)*kpt%list%nk + ik-1)
        off2 = ( nbasis_subset + max_val )
        off3 = nxpts
        offset = off1 * off2 * off3
!        offset =  ((ispin-1)*kpt%list%nk + ik-1)*( nbasis_subset + max_val ) * nxpts
        
!????? ik -> 1
        call mpi_file_write_at( fhu2, offset, u2(1,1,1,1), &
                                max_val * nxpts, &
                                MPI_DOUBLE_COMPLEX, status, ierr )
        offset = offset + max_val * nxpts
!????? ik -> 1
        call mpi_file_write_at( fhu2, offset, u2(1,1,1,2), &
                                nbasis_subset * nxpts, &
                                MPI_DOUBLE_COMPLEX, status, ierr )
      endif  
    

      if( have_kshift ) then
        offset =  ((ispin-1)*kpt%list%nk + ik-1)*(band_subset(2) - band_subset(1) + 1 ) * max_val
        call mpi_file_write_at( fhtmels, offset, tmels(1,1,ik), &
                                (band_subset(2) - band_subset(1) + 1 ) * max_val, &
                                MPI_DOUBLE_COMPLEX, status, ierr )
      endif



    endif
      

  enddo ! ik

  enddo ! ispin

  if( mypoolid .eq. mypoolroot ) then
    
    call MPI_FILE_CLOSE( fheig, ierr )
  endif

!#ifdef FALSE



  deallocate( u1, o2l )

  if( mypoolid == mypoolroot ) then
    call mp_sum( homo, cross_pool_comm )
    homo_point = minval( homo )
    call mp_sum( lumo, cross_pool_comm )
    lumo_point = minval( lumo )
    call mp_sum( start_band, cross_pool_comm )

    call mp_sum( e0, cross_pool_comm )

    ibeg_unit = freeunit()
    open(unit=ibeg_unit,file='ibeg.h',form='formatted',status='unknown')
    rewind(ibeg_unit)

!    select case (band_style)
!    if( band_style .eq. 'band' ) then
!      case( 'band' )
!        start_band(:) = nelectron / 2 + 1
        do ispin = 1, nspin
          do ik = 1, kpt%list%nk
            start_band( ik ) = nelectron / 2 + 1
            write(ibeg_unit,*) ik, start_band( ik )
          enddo
        enddo

!      case default
    else
    
    ! Sort energies to determine true Fermi, homo, lumo
!    fermi_energy = fermi_energy / 2.0_DP



!    call fix_fermi( nbasis_subset, kpt%list%nk, nspin, nshift, max_val, nelectron, 0, &
!                    e0, homo_point, lumo_point, fermi_energy )
!    
!    do ispin = 1, nspin
!      do ik = 1, kpt%list%nk
!        do ibd = band_subset(1), band_subset( 2 )
!        if( 2.0_DP * e0( ibd, ik, ispin, nshift ) .gt. fermi_energy ) then
!          start_band( ik ) = ibd
!          lumo( ik ) = 2.0_DP * e0( ibd, ik, ispin, nshift )
!          goto 21
!        endif
!        enddo
!21    continue
!        write(ibeg_unit,*) ik, start_band( ik ), -1
!      enddo
!    enddo
!
!    endif
!    end select
    close(ibeg_unit)
    !


!  endif
!  call mp_bcast( homo, ionode_id )
!  call mp_bcast( lumo, ionode_id )
!  call mp_bcast( start_band, ionode_id )
!  call mp_bcast( e0, ionode_id )
    write(stdout,*) start_band(1:8)

  allocate(stop_band(kpt%list%nk))
!  goto 13
  if( dft_energy_range .le. 0.0 ) then
    nbuse = band_subset( 2 )
    do ik = 1,kpt%list%nk
      nbuse = min( nbuse, band_subset( 2 ) - start_band( ik ) + 1 )
      write(6,*) ik, nbuse
    enddo
  else
    nbuse = 0
    do ispin=1,nspin
       do ik = 1,kpt%list%nk
          do ibd = band_subset( 2 ) - 1, start_band( ik ), -1
             if( 2.0d0* e0( ibd, ik, ispin, 1 ) .le. dft_energy_range + lumo_point ) then
                stop_band( ik ) = ibd + 1
                goto 12
             endif
          enddo
12        continue
!      if( stop_band( ik ) .eq. band_subset( 2 ) ) call errore('OCEAN', 'DFT_ENERGY_RANGE too large' )
          if( stop_band( ik ) .eq. band_subset( 2 ) ) then 
             write(stdout,*) 'DFT_ENERGY_RANGE too large: ', dft_energy_range*rytoev
             dft_energy_range = 2.0d0* e0(  band_subset( 2 ) - 1, ik, ispin, 1 ) - lumo_point
             write(stdout,*) 'New energy range: ', dft_energy_range*rytoev
          endif
          nbuse = max( nbuse, stop_band( ik ) - start_band( ik ) + 1 )
       enddo
    enddo
    do ispin=1,nspin
       do ik = 1,kpt%list%nk
          if( nbuse + start_band(ik) .gt. band_subset( 2 ) ) then
             write(stdout,*) 'Reducing range to fit', ik, start_band(ik)
             nbuse = band_subset( 2 ) - start_band(ik)
             !        write(stdout,*) 'New energy range: ', dft_energy_range*rytoev
          endif
          dft_energy_range = min(dft_energy_range, 2.0d0 * e0( nbuse + start_band( ik ), ik, ispin, 1 ) -lumo_point )
       enddo
    enddo
    write(stdout,*) 'Final energy range: ', dft_energy_range * rytoev
  endif
  write(stdout,*) 'nbuse: ', nbuse
13 continue

  ! close MPI-IO files
!  if( mypoolid==mypoolroot ) then
    ! close binary dump
    call mpi_file_close( fheigval, ierr )
    if( ierr/=0 ) &
      call errore('shirley_qdiagp','problem closing eigval file',abs(ierr))

    if( eigvec_output ) then
      ! close eigvec file
      call mpi_file_close( fheigvec, ierr )
      if( ierr/=0 ) &
        call errore('shirley_qdiagp','problem closing eigvec file',abs(ierr))
    endif
!!    close(fhenk)
  endif

  call mp_bcast( nbuse, ionode_id )
  call mp_bcast( start_band, ionode_id )

! Downsize coeff
  
  nbuse_xes = maxval( start_band ) - 1
!  call mp_sum( coeff, intra_pool_comm )

  if( mypoolid == mypoolroot ) then
    allocate( coeff_small( nptot, nbuse, kpt%list%nk, nspin, ntau ) )
    if( have_kshift ) then
      do ik=1,kpt%list%nk
        coeff_small( :, 1:nbuse, ik, :, : ) = coeff( :, start_band(ik):start_band(ik)+nbuse-1, ik, :, :, 2 )
      enddo
    else
      do ik=1,kpt%list%nk
        coeff_small( :, 1:nbuse, ik, :, : ) = coeff( :, start_band(ik):start_band(ik)+nbuse-1, ik, :, :, 1 )
    enddo
    endif
    
    call mp_sum( coeff_small, cross_pool_comm )
  endif

  if( ionode ) then
    iuntmp = freeunit()
    do itau = 1, ntau
      write( filout, '(1a5,1a6)' ) 'cksc.', fntau( itau )
  !    open(unit=iuntmp,file='test_cks.bin',form='unformatted',status='unknown')
      open(unit=iuntmp,file=filout,form='unformatted',status='unknown')
      rewind(iuntmp)
      write(iuntmp) nptot, nbuse*kpt%list%nk, nspin
  !     write(iuntmp) nptot, ntot, nspin
      write(iuntmp) tau( :, itau )
      write(iuntmp) real(coeff_small(:,:,:,:,itau))
      write(iuntmp) aimag(coeff_small(:,:,:,:,itau))
      close(iuntmp)
    enddo
  endif

  if( mypoolid == mypoolroot ) then
    deallocate( coeff_small )
  
    allocate( coeff_small( nptot, nbuse_xes, kpt%list%nk, nspin, ntau ) )
    coeff_small = 0.0d0
    do ik=1,kpt%list%nk
      coeff_small( :, 1:start_band(ik)-1, ik, :, : ) = coeff( :, 1:start_band(ik)-1, ik, :, :, 1 )
    enddo

    call mp_sum( coeff_small, cross_pool_comm )
  endif

  if( ionode ) then
    iuntmp = freeunit()
    do itau = 1, ntau
      write( filout, '(1a5,1a6)' ) 'cksv.', fntau( itau )
!    open(unit=iuntmp,file='test_cks_val.bin',form='unformatted',status='unknown')
      open(unit=iuntmp,file=filout,form='unformatted',status='unknown')
      rewind(iuntmp)
      write(iuntmp) nptot, nbuse_xes*kpt%list%nk, nspin
  !     write(iuntmp) nptot, ntot, nspin
      write(iuntmp) tau( :, itau )
      write(iuntmp) real(coeff_small(:,:,:,:,itau))
      write(iuntmp) aimag(coeff_small(:,:,:,:,itau))
      close(iuntmp)
    enddo
  endif

  if( mypoolid == mypoolroot)   deallocate( coeff_small )
  
  deallocate( coeff )

  if( ionode ) then 
    brange( 1 ) = band_subset( 1 )
    brange( 2 ) = maxval( start_band( : ) ) - 1
    brange( 3 ) = minval( start_band( : ) )
!JTV
!    brange( 4 ) = brange( 3 ) + nbuse - 1
    brange( 4 ) = nbasis_subset
    iuntmp = freeunit()
    open(unit=iuntmp,file='brange.ipt',form='formatted',status='unknown')
    rewind(iuntmp)
    write(iuntmp,*) brange( 1 ), brange( 2 )
    write(iuntmp,*) brange( 3 ), brange( 4 )
    close( iuntmp )
  endif

!?????  Would have to go back and fix u2
  if( mypoolid == mypoolroot .and. .false. ) then
    call mp_sum( u2, cross_pool_comm )
!  endif
!    do ik = 1, kpt%list%nk
!      call mp_put( ionode_id, mod(ik-1,npool), mypool, u2(:,:,ik,1), u2(:,:,ik,1),ik,cross_pool_comm )
!    enddo
  if( ionode ) then
    brange( 1 ) = band_subset( 1 )
    brange( 2 ) = maxval( start_band( : ) ) - 1
    brange( 3 ) = minval( start_band( : ) )
!JTV
    brange( 4 ) = brange( 3 ) + nbuse - 1
    iuntmp = freeunit()
    open(unit=iuntmp,file='brange.ipt',form='formatted',status='unknown')
    rewind(iuntmp)
    write(iuntmp,*) brange( 1 ), brange( 2 )
    write(iuntmp,*) brange( 3 ), brange( 4 )
    close( iuntmp )
    
    open(unit=iuntmp,file='u2.dat',form='unformatted',status='unknown')
    rewind(iuntmp)
    do ik=1,kpt%list%nk
!      do ibd = 1, band_subset(2)
      do ibd = brange( 1 ), brange( 2 )
        do ix = 1, nxpts
          write ( iuntmp ) ibd, ik, ix, u2( ix, ibd, ik, 1 )
        end do
      end do
      do ibd = brange( 3 ), brange( 4 )
        do ix = 1, nxpts
          write ( iuntmp ) ibd, ik, ix, u2( ix, ibd, ik, nshift )
        end do
      end do
    end do
    close(iuntmp)
  endif
  endif


    

! Downsize e0
  ! lumo to Ha
  lumo_point = 0.5d0 * lumo_point
  homo_point = 0.5d0 * homo_point

  if( legacy_zero_lumo ) then
    lumo_shift = lumo_point
  else
    lumo_shift = 0.0_DP
  endif

  allocate( e0_small( nbuse, kpt%list%nk, nspin ) )
  do ispin=1,nspin
     if( have_kshift ) then
        do ik=1,kpt%list%nk
           e0_small( 1 : nbuse, ik, ispin ) = e0( start_band(ik) : start_band(ik) + nbuse - 1, ik, ispin, 2 ) &
                - lumo_shift
        enddo
     else
        do ik=1,kpt%list%nk
           e0_small( 1 : nbuse, ik, ispin ) = e0( start_band(ik) : start_band(ik) + nbuse - 1, ik, ispin, 1 ) &
                - lumo_shift
        enddo
     endif
  enddo


!  call mp_sum( e0_small )
!  call mp_sum( e0 )
  if( ionode ) then 
    write(stdout,*) 'LUMO: ', lumo_point
    open(unit=iuntmp,file='wvfcninfo',form='unformatted',status='unknown')
    rewind(iuntmp)
!    write(iuntmp) nbasis_subset,kpt%list%nk,nspin
    write(iuntmp) nbuse,kpt%list%nk,nspin
    write(iuntmp) e0_small(:,:,:)
    close(iuntmp)
  endif

  deallocate( e0_small )
  if( ionode ) then
    open(unit=iuntmp,file='nbuse.ipt',form='formatted',status='unknown')
    write(iuntmp,*) nbuse
    close(iuntmp)
  endif

  allocate( e0_small( nbuse_xes, kpt%list%nk, nspin ) )
  e0_small = 0.0d0
  do ispin=1,nspin
     do ik=1,kpt%list%nk
        e0_small( 1 : start_band(ik) - 1, ik, ispin ) = e0( 1 : start_band(ik) - 1, ik, ispin, 1 ) - lumo_shift
     enddo
  enddo
  if( ionode ) then
    open(unit=iuntmp,file='wvfvainfo',form='unformatted',status='unknown')
    rewind(iuntmp)
!    write(iuntmp) nbasis_subset,kpt%list%nk,nspin
    write(iuntmp) nbuse_xes,kpt%list%nk,nspin
    write(iuntmp) e0_small(:,:,:)
    close(iuntmp)
  endif

  deallocate( e0_small )
  if( ionode ) then
    open(unit=iuntmp,file='nbuse_xes.ipt',form='formatted',status='unknown')
    write(iuntmp,*) nbuse_xes
    close(iuntmp)
  endif



  ! Legacy AI2NBSE energy file
  ! indicies are reversed for no reason
  if( ionode .and. have_kshift) then
    e0( :,:,:,: ) = e0( :, :,:,: ) - homo_point
    sef = efermi - homo_point
    inquire( file='gwipt', exist=have_gwipt)
    if( have_gwipt ) then
      open(unit=iuntmp,file='gwipt',form='formatted',status='old')
      read(iuntmp,*) egw, elda, vs, cs
      eshift = egw - elda
      newgap = lumo_point - homo_point + eshift
    else
      eshift = 0.0d0
      newgap = lumo_point - homo_point
      vs = 0.0d0
      cs = 0.0d0
    endif
    open(unit=iuntmp,file='ebdat',form='formatted',status='unknown')
    rewind(iuntmp)
    do ispin=1, nspin
       do ibd = brange(1),brange(2)
          do ik=1,kpt%list%nk
             if( e0( ibd, ik, ispin, 1 ) .lt. sef ) then
                write(iuntmp,*) (e0( ibd, ik, ispin, 1 ) * (1.d0 + vs ) ),  e0( ibd, ik, ispin, 1 ) + homo_point
             else
                write(iuntmp,*) (newgap + (1.d0 + cs) * (e0( ibd, ik, ispin, 1 ) - newgap ) ),  e0( ibd, ik, ispin, 1 ) + homo_point
             endif
          enddo
       enddo
       do ibd = brange(3),brange(4)
          do ik=1,kpt%list%nk
             if( e0( ibd, ik, ispin, 2 ) .lt. sef ) then
                write(iuntmp,*) (e0( ibd, ik, ispin, 2 ) * (1.d0 + vs ) ),  e0( ibd, ik, ispin, 2 ) + homo_point
             else
                write(iuntmp,*) (newgap + (1.d0 + cs) * (e0( ibd, ik, ispin, 2 ) - newgap ) ),  e0( ibd, ik, ispin, 2 ) + homo_point
             endif
          enddo
       enddo
    enddo ! ispin
    close(iuntmp)
  endif

  deallocate( e0 )



  if( have_kshift .and. ( mypoolid .eq. mypoolroot ) ) then
    call mp_sum( tmels, cross_pool_comm )
    if( ionode ) then
      open(unit=iuntmp,file='pdadat',form='formatted',status='unknown')
      rewind(iuntmp)
      do ik = 1, kpt%list%nk
          do i = brange( 1 ), brange( 2 )
            do j = brange( 3 ), brange( 4 )
              write( iuntmp, '(2(1x,1e22.15))' ) real(tmels( i, j, ik )), aimag(tmels( i, j, ik ))
            enddo
          enddo
        do j = 1, 3
          do ibd = brange( 1 ), brange( 2 )
            do i = brange( 3 ), brange( 4 )
              write( iuntmp, '(2(1x,1e22.15))' ) 0.0, 0.0
            enddo
          enddo
        enddo
      enddo
      close(iuntmp)
      open(unit=iuntmp,file='tmels.info',form='formatted',status='unknown')
      rewind(iuntmp)
      write(iuntmp,*) max_val, band_subset( 1 ),  band_subset( 2 ), kpt%list%nk
      close(iuntmp)
    endif
    deallocate(tmels)
  endif

  if( ionode ) then
    open(unit=iuntmp,file='obf_control',form='formatted',status='unknown')
    rewind(iuntmp)
    write(iuntmp,*) nbasis
    write(iuntmp,*) max_val
    write(iuntmp,*) nbasis_subset
  endif

#endif

  if( mypoolid==mypoolroot ) then
    call MPI_FILE_CLOSE( fhu2, ierr )
    call MPI_FILE_CLOSE( fhtmels, ierr )
    call MPI_FILE_CLOSE( fh_val_energies, ierr )
    call MPI_FILE_CLOSE( fh_con_energies, ierr )

    call MPI_FILE_CLOSE( fh_val_eigvecs, ierr )
    call MPI_FILE_CLOSE( fh_con_eigvecs, ierr )
  endif

111 continue 
  write(stdout,*) ' end shirley_qdiag'
  call mp_end
  stop
  
  end program shirley_qdiag
