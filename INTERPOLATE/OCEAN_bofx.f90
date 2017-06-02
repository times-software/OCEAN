! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the GPL 2 License. See the file `License' in the current subdirectory.
!
subroutine OCEAN_bofx( )
! This is a stopgap, going to probably try and evaluate the matrix elements of W
!   in the optimal basis. For now, just using the ai2nbse stylings
!
! Here we project the Basis functions onto a grid specified by xmesh.ipt 
! This is done through fft, currently using the old ai2nbse/ocean coding
!
! Now with option to have normal ordering of x,y,z elements, default is legacy ordering
!
! Output: u1.dat
!
! Aux input: xmesh.ipt
!
  use constants, only : tpi
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE cell_base
  USE gvect
  USE klist, ONLY: xk, nks, nkstot, ngk
  USE wvfct
  USE io_files, ONLY: prefix, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY: evc
  USE mp, ONLY : mp_sum, mp_bcast,  mp_max, mp_min
  USE mp_global, ONLY : me_pool, nproc_pool, root_pool, mpime, intra_pool_comm, & 
                        mypoolid => my_pool_id, mypoolroot => root_pool
  USE buffers, ONLY : get_buffer
  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset
!  use hamq_pool, only : mypool, mypoolid, mypoolroot, cross_pool_comm, intra_pool_comm
  
  use mpi
  use OCEAN_bofx_mod
  use OCEAN_obf2loc
  use OCEAN_timer

  implicit none

  real(dp),parameter ::  eps=1.d-10
  
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)

  integer(kind=MPI_OFFSET_KIND), parameter :: oneGBinComplex = 67108864

  character(255) :: string
  logical :: loud

  integer :: ibnd  
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  
  integer :: i, ig, dumi
  integer,allocatable :: igk_l2g(:,:), ibeg(:,:)

  integer :: xmesh( 3 )
  logical :: bt_logical, invert_xmesh
  character(len=3) :: dum_c

  integer :: iuntmp, iobegin, tick, tock, ishift

  integer :: min_mill(3), uni_min_mill(3), max_mill(3), uni_max_mill(3), igl, igh, nfft(3)

  integer :: nshift, nspin, nkpt, ispin, ikpt, kpts(3)
  integer :: nbasis_, nband, ierr, errorcode, fmode, resultlen, band_block, band_length

  complex(dp), allocatable :: eigvec(:,:), unk(:,:)

  integer :: ibd, nelement, writeEveryNKpt
  integer :: fheig, fhu2
  integer :: u2_type
  integer(kind=MPI_OFFSET_KIND) :: offset, long_s, long_k, long_x, long_b

  integer,external :: freeunit

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
#if defined( __NIST ) && defined( __BUFFER )
  CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
#else
  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )
#endif

!  CALL davcio( evc, 2*nwordwfc, iunwfc, 1, - 1 )

  ! report norms
  allocate( norm(nbnd) )
  do ibnd=1,nbnd
    norm(ibnd) = dot_product( evc(:,ibnd_indx(ibnd)), evc(:,ibnd_indx(ibnd)) )
  enddo
  call mpi_barrier( intra_pool_comm, ierr )
  call mp_sum( norm, intra_pool_comm )
  do ibnd=1,nbnd
    if( abs(norm(ibnd)-1.d0) > eps ) then
      write(stdout,'(a,i6,a,f14.10)') ' band ', ibnd, ' norm = ', norm(ibnd)
      call errore('hamq','wave function norm is not 1',1)
    endif
  enddo
  deallocate( norm )

  if( ionode ) then
    iuntmp = freeunit()
    open(iuntmp,file='xmesh.ipt',form='formatted',status='old')
    read(iuntmp,*) xmesh( : )
    close(iuntmp)

    inquire(file="bloch_type",exist=bt_logical)
    if( bt_logical ) then
      open(iuntmp,file="bloch_type",form='formatted',status='old')
      read(iuntmp,*) dum_c
      read(iuntmp,*) invert_xmesh
      close(iuntmp)
    else
      invert_xmesh = .true.
    endif

    open(iuntmp,file='qdiag.info',form='formatted',status='old')
    read(iuntmp,*)  nbasis_, nband, nkpt, nspin, nshift
    close(iuntmp)

    open(iuntmp,file='kmesh.ipt',form='formatted',status='old')
    read(iuntmp,*) kpts(:)
    close(iuntmp)

    iobegin = freeunit()
    open(iobegin,file='ibeg.h',form='formatted',status='old')
    allocate(ibeg(nkpt,nspin) )
    do ispin = 1, nspin
      do ikpt = 1, nkpt
        read(iobegin,*) dumi, ibeg( ikpt, ispin )
      enddo
    enddo
    close(iobegin)

  endif
  call mp_bcast( xmesh, ionode_id )
  call mp_bcast( invert_xmesh, ionode_id )
  call mp_bcast( nbasis_, ionode_id )
  call mp_bcast( nband,  ionode_id )
  call mp_bcast( nkpt,   ionode_id )
  call mp_bcast( nspin,  ionode_id )
  call mp_bcast( nshift, ionode_id )
  call mp_bcast( kpts, ionode_id )

  if( .not. ionode ) allocate(ibeg(nkpt,nspin) )
  call mp_bcast( ibeg, ionode_id )

  write(stdout,*) nbnd, band_subset(1), band_subset(2)




  !!!!!!!!!  Eigvecs !!!!!!!!!!!
  CALL MPI_FILE_OPEN( intra_pool_comm, 'eigvecs.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fheig, ierr )
  if( ierr/=0 ) then
    errorcode=ierr
    call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
    call errore(string,errorcode)
  endif
  offset=0

  call MPI_FILE_SET_VIEW( fheig, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, 'native', & 
                          MPI_INFO_NULL, ierr )
  if( ierr/=0 ) then
    errorcode=ierr
    call MPI_ERROR_STRING( errorcode, string, resultlen, ierr ) 
    call errore(string,errorcode)
  endif
  !!!!!!!!!  Eigvecs !!!!!!!!!!!


  !!!!!!!!!  u2par.dat !!!!!!!!!
!  if( mypoolid .eq. mypoolroot ) then
  if( .false. ) then
    call MPI_TYPE_CONTIGUOUS( product(xmesh), MPI_DOUBLE_COMPLEX, u2_type, ierr )
    call MPI_TYPE_COMMIT( u2_type, ierr )

    fmode = IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY)
    call MPI_FILE_OPEN( intra_pool_comm, 'u2par.dat', fmode, MPI_INFO_NULL, fhu2, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
    offset=0
    call MPI_FILE_SET_VIEW( fhu2, offset, u2_type, &
                            u2_type, 'native', MPI_INFO_NULL, ierr )
!    call MPI_FILE_SET_VIEW( fhu2, offset, MPI_DOUBLE_COMPLEX, &
!                          MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr )
    if( ierr/=0 ) then
      errorcode=ierr
      call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
      call errore(string,errorcode)
    endif
  else


  call OCEAN_bofx_open_file( product(xmesh), ierr )
  if( ierr .ne. 0 ) then
    write(stdout,*) 'Failed calling OCEAN_bofx_open_file', ierr
  endif

  endif


#define OBF
#ifdef OBF
!!!! obf?
  call OCEAN_obf2loc_init( nband, kpts, nspin, nshift, ierr )
  if( ierr .ne. 0 ) then
    write(stdout,*) 'Failed to run OCEAN_obf2loc_init'
  endif

  call OCEAN_obf2loc_alloc( ierr )
  if( ierr .ne. 0 ) then
    write(stdout,*) 'Failed to run OCEAN_obf2loc_alloc'
  endif

!!!!
#endif


  allocate( unk( npw, nband ) )

  long_b = nband
  long_x = product(xmesh)


  ! Block nband for large systems
  if( int(nband,MPI_OFFSET_KIND) * int( nbasis_, MPI_OFFSET_KIND ) > oneGBinComplex ) then
    if( nbasis_ > oneGBinComplex ) then
      band_block = min( 1, nband )
    elseif( nbasis_ > 65536 ) then
      band_block = min( 128, nband )
    elseif( nbasis_ > 16382 ) then
      band_block = min( 512, nband )
    else
      band_block = min( 2048, nband ) 
    endif
  else
    band_block = nband
  endif
  !
  write(stdout,*) 'Blocking eigvec:', band_block, nband


  ! Need to reign in amount of stdout
  if( nkpt .gt. 100 ) then
    writeEveryNKpt = max( 10, nkpt / 100 )
  else
    writeEveryNKpt = 1
  endif
  
  
  

  loud = .true.

! Set up the fft grid
  call par_gen_setup( xmesh, npw, mill, loud )


  ishift = 1

  do ispin = 1 , nspin
    do ikpt = 1, nkpt

      allocate( eigvec( nbasis_, band_block ) )


      offset = ( ispin - 1 ) * nkpt + (ikpt - 1 )
      offset = offset * int(nbasis_, kind( offset ) )
      offset = offset * int(nband, kind( offset ) )

      call OCEAN_t_reset

      do ibd = 1, nband, band_block
        nelement = min( band_block, nband - ibd + 1 )

        call MPI_FILE_READ_AT_ALL( fheig, offset, eigvec, nelement*nbasis_, MPI_DOUBLE_COMPLEX, &
                                   MPI_STATUS_IGNORE, ierr )

      
        ! still zero because we are doing a complete chunk in bands
        call ZGEMM( 'N', 'N', npw, nelement, nbasis_, one, evc, npw, eigvec, nbasis_, zero, &
                     unk(1,ibd), npw )
        offset = offset + int( nelement, MPI_OFFSET_KIND ) * int( nbasis_, MPI_OFFSET_KIND )

      enddo
      if( loud ) call OCEAN_t_printtime( "Read eigvec", stdout )

      deallocate( eigvec )

      offset = ( ispin - 1 ) * nkpt + (ikpt - 1 )
!      offset = offset * product( xmesh )
!      offset = offset * nband
      offset = offset * long_b
!      offset = offset * long_x
      call par_gentoreal( xmesh, nband, unk, npw, mill, offset, invert_xmesh, loud, ibeg(ikpt,ispin) )  
      loud = .false.
    
      if( mod( ikpt, writeEveryNKpt ) .eq. 0 ) then
        write(stdout,*) ikpt, offset
      endif
  

#ifdef OBF
!!!!!!
      call OCEAN_obf2loc_coeffs( npw, mill, unk, ikpt, ispin, ishift, ibeg, loud )

!!!!!!
#endif

    enddo
  enddo

#ifdef OBF
!!!!
  call OCEAN_obf2loc_write( ierr )
  if( ierr .ne. 0 ) then
    write(stdout,*) 'Failed obf2loc_write'
  endif
  ! catch errors from obf2loc
  call MPI_BARRIER( intra_pool_comm, ierr )
  if( ierr .ne. 0 ) then
    write(stdout,*) 'Failed barrier after OCEAN_obf2loc_write'
  else
    write(stdout,*) 'Finished OCEAN_obf2loc_write'
  endif

!!!!
#endif
  

  call MPI_FILE_CLOSE( fheig, ierr )


  call par_gen_shutdown( ierr )
  call OCEAN_bofx_close_file( ierr )

  deallocate( unk )


  call MPI_BARRIER( intra_pool_comm, ierr )

  return

end subroutine OCEAN_bofx
