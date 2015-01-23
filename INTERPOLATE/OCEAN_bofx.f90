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
  USE mp_global, ONLY : me_pool, nproc_pool, root_pool, mpime
  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset
  use hamq_pool, only : mypool, mypoolid, mypoolroot, cross_pool_comm, intra_pool_comm
  use mpi
  use OCEAN_bofx_mod
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

  integer :: iuntmp, iobegin

  integer :: min_mill(3), uni_min_mill(3), max_mill(3), uni_max_mill(3), igl, igh, nfft(3)

  integer :: nshift, nspin, nkpt, ispin, ikpt
  integer :: nbasis_, nband, ierr, errorcode, fmode, resultlen, band_block, band_length

  complex(dp), allocatable :: eigvec(:,:), unk(:,:)

  integer :: ibd, nelement
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
    read(iuntmp,*),  nbasis_, nband, nkpt, nspin, nshift
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

  if( .not. ionode ) allocate(ibeg(nkpt,nspin) )
  call mp_bcast( ibeg, ionode_id )

  write(stdout,*) nbnd, band_subset(1), band_subset(2)
!  allocate( sub_mill( 3, npw ) )
!  sub_mill = 0
!  do ig = 1, npw
!    sub_mill( :, ig ) = matmul( g( :, igk( ig ) ), bg( :, : ) )
!  enddo


  !!!!!!!!!
  ! Find maximum extent of mill
  do i = 1, 3
    igl = mill( i, 1 )
    igh = mill( i, 1 )
    do ig = 2, npw
      igl = min( igl, mill( i, ig ) )
      igh = max( igh, mill( i, ig ) )
    end do

    min_mill( i ) = igl !minval( mill( i, : ) )
    max_mill( i ) = igh !maxval( mill( i, : ) )

    uni_min_mill(i) = min_mill(i)
    uni_max_mill(i) = max_mill(i)

    call mp_min( uni_min_mill(i) )
    call mp_max( uni_max_mill(i) )

  enddo

  write(stdout,*) 'loc min:', min_mill(:)
  write(stdout,*) 'loc max:', max_mill(:)

  write(stdout,*) 'uni min:', uni_min_mill(:)
  write(stdout,*) 'uni max:', uni_max_mill(:)

  do i = 0, nproc_pool-1
    if( i .eq. mpime ) then
      write(1000+mpime,*) mpime, npw
    endif
  enddo


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
!  else
!  endif



  allocate( unk( npw, nband ) )

  long_b = nband
  long_x = product(xmesh)


  ! Block nband for large systems
  if( int(nband,MPI_OFFSET_KIND) * int( nbasis_, MPI_OFFSET_KIND ) > oneGBinComplex ) then
    if( nbasis_ > oneGBinComplex ) then
      band_block = min( 1, nband )
    elseif( nbasis_ > 65536 ) then
      band_block = min( 256, nband )
    elseif( nbasis_ > 16382 ) then
      band_block = min( 1024, nband )
    else
      band_block = min( 4096, nband ) 
    endif
  else
    band_block = nband
  endif
  !
  write(stdout,*) 'Blocking eigvec:', band_block, nband
  
  

  loud = .true.

! Set up the fft grid
  call par_gen_setup( xmesh, npw, mill, loud )


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
      call OCEAN_t_printtime( "Read eigvec", stdout )
  
      deallocate( eigvec )

      offset = ( ispin - 1 ) * nkpt + (ikpt - 1 )
!      offset = offset * product( xmesh )
!      offset = offset * nband
      offset = offset * long_b
!      offset = offset * long_x
      call par_gentoreal( xmesh, nband, unk, npw, mill, fhu2, offset, invert_xmesh, loud, ibeg(ikpt,ispin), u2_type )  
    !  loud = .false.
      write(stdout,*) ikpt, offset
  

    enddo
  enddo


  call MPI_FILE_CLOSE( fheig, ierr )
!  if(  mypoolid .eq. mypoolroot ) then
  call MPI_FILE_CLOSE( fhu2, ierr )
!  endif

  deallocate( unk )


  return

  contains


  subroutine gentoreal( nx, nfcn, fcn, ng, gvec, iu, offset, invert_xmesh, loud, ibeg)
    use kinds, only : dp
    USE io_global,  ONLY : stdout, ionode
    USE mp, ONLY : mp_sum, mp_max, mp_min
    implicit none
    !
    integer, intent( in ) :: nx( 3 ), nfcn, ng, iu
    integer, intent( in ) :: gvec( 3, ng )
    complex(dp), intent( in ) :: fcn( ng, nfcn )
    logical, intent( in ) :: invert_xmesh
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: offset
    logical, intent( in ) :: loud
    integer, intent( in ) :: ibeg
!    logical, parameter :: loud = .true.
    !
    integer :: ix, iy, iz, i1, i2, i3, nfft( 3 ), idwrk, igl, igh, nmin, ii
    integer :: fac( 3 ), j, ig, i, locap( 3 ), hicap( 3 ), toreal, torecp, nftot
    real(dp) :: normreal, normrecp
    complex(dp) :: rm1, w
    character * 80 :: fstr
    !
    integer, parameter :: nfac = 3
    integer, allocatable :: ilist( :, : )
    real(dp), allocatable :: zr( :, :, : ), zi( :, :, : ), wrk( : )
    complex(dp), allocatable :: cres( :, :, :, : )
    real(dp), external :: DZNRM2
    complex(dp), external :: ZDOTC
    integer, external :: optim
    integer :: start_band
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    fac( 1 ) = 2; fac( 2 ) = 3; fac( 3 ) = 5
    hicap( 1 ) = 20; hicap( 2 ) = 3; hicap( 3 ) = 1
    !
    allocate( ilist( max( nx( 1 ), nx( 2 ), nx( 3 ) ), 3 ) )
    do j = 1, 3
       igl = gvec( j, 1 )
       igh = gvec( j, 1 )
       do ig = 2, ng
          igl = min( igl, gvec( j, ig ) )
          igh = max( igh, gvec( j, ig ) )
       end do
       call mp_max( igh )
       call mp_min( igl )
       nmin = 1 + igh - igl
       call facpowfind( nx( j ), nfac, fac, locap ) 
       nfft( j ) = optim( nmin, nfac, fac, locap, hicap )
       fstr = '(1i1,1a1,5x,1a4,1i4,5x,1a4,1i4,5x,1a3,1i4,5x,1a5,1i4)'
       if ( loud ) write ( stdout, fstr ) j, ':', 'igl=', igl, 'igh=', igh, 'nx=', nx( j ), 'nfft=', nfft( j )
       if ( igl * igh .ge. 0 ) stop 'zero not included!'
       do i = 1, nx( j )
          ilist( i, j ) = 1 + ( i - 1 ) * nfft( j ) / nx( j )
       end do
       if ( loud ) write ( stdout, '(15i5)' ) ilist( 1 : nx( j ), j )
    end do
    ! 
    i = max( nfft( 1 ), nfft( 2 ), nfft( 3 ) )
    idwrk = 2 * i * ( i + 1 )
    allocate( zr( nfft( 1 ), nfft( 2 ), nfft( 3 ) ) )
    allocate( zi( nfft( 1 ), nfft( 2 ), nfft( 3 ) ) )
    allocate( wrk( idwrk ) )
    nftot = nfft( 1 ) * nfft( 2 ) * nfft( 3 ) 
    !
    if( invert_xmesh ) then
      ! the indices only of the output are reversed here.
      allocate( cres( nx( 3 ), nx( 2 ), nx( 1 ), nfcn ) )
    else
      allocate( cres( nx( 1 ), nx( 2 ), nx( 3 ), nfcn ) )
    endif
    call chkfftreal( toreal, normreal, loud )
    call chkfftrecp( torecp, normrecp, loud )
    !
    do i = 1, nfcn
       zr( :, :, : ) = 0.0d0; zi( :, :, : ) = 0.0d0
       fstr = '(3(1a1,2i8,2x),1a5,2i8)'
       do ig = 1, ng
          i1 = 1 + gvec( 1, ig )
          if ( i1 .le. 0 ) i1 = i1 + nfft( 1 )
          i2 = 1 + gvec( 2, ig )
          if ( i2 .le. 0 ) i2 = i2 + nfft( 2 )
          i3 = 1 + gvec( 3, ig )
          if ( i3 .le. 0 ) i3 = i3 + nfft( 3 )
          zr( i1, i2, i3 ) = real( fcn( ig, i ) )
          zi( i1, i2, i3 ) = aimag( fcn( ig, i ) )
          if ( loud .and. ( ig .le. 10 ) .and. ( i .le. 3 ) ) then
             if ( loud ) write ( stdout, fstr ) 'x', gvec( 1, ig ), i1,  &
                               'y', gvec( 2, ig ), i2, 'z', gvec( 3, ig ), i3, 'ig, i', ig, i
          end if
       end do
!       call mp_root_sum( zr, zr, ionode_id, intra_pool_comm )
!       call mp_root_sum( zi, zi, ionode_id, intra_pool_comm )
       call mp_sum( zr, intra_pool_comm )
       call mp_sum( zi, intra_pool_comm )
       if( mypoolid .eq. mypoolroot ) then
         call cfft( zr, zi, nfft( 1 ), nfft( 1 ), nfft( 2 ), nfft( 3 ), toreal, wrk, idwrk )
         zr = zr / dble( nftot ) ** normreal
         zi = zi / dble( nftot ) ** normreal
         ii = 0
         fstr = '(2(1a9,3i5,5x),1a9,2(1x,1e15.8))'
         do iz = 1, nx( 3 )
            i3 = ilist( iz, 3 )
            do iy = 1, nx( 2 )
               i2 = ilist( iy, 2 )
               do ix = 1, nx( 1 )
                  i1 = ilist( ix, 1 )
                  ii = ii + 1
                  if ( loud .and. ( ii .le. 10 ) .and. ( i .le. 3 ) ) then
                     write ( stdout, fstr ) 'mesh ind.', i1, i2, i3,  &
                                            'cell ind.', ix, iy, iz,  &
                                            'value = ', zr( i1, i2, i3 ), zi( i1, i2, i3 )
                  end if
                  if( invert_xmesh ) then
                    cres( iz, iy, ix, i ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
                  else
                    cres( ix, iy, iz, i ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
                  endif
               end do
            end do
         end do 
!         write ( iu ) cres
!         if( mypoolid .eq. mypoolroot ) then
!           call MPI_FILE_WRITE_AT( iu, offset, cres, product(nx), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
!           offset = offset + product(nx)
!         endif
       endif
    end do

    if( mypoolid == mypoolroot ) then

!      write(stdout,*) ibeg
      ! This now mimics orthog better. 
      do i = 1, nfcn
        normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
        normreal = 1.0_dp / normreal
        call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )

        if( i .lt. ibeg ) then
          start_band = 1
        else
          start_band = ibeg
        endif

        do j = start_band, i-1
          w = -ZDOTC( product(nx), cres( 1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
! Added to better mimic orthog.f90
          normreal = DZNRM2( product(nx), cres( 1, 1, 1, j ), 1 )
          w = w / normreal
!\\
          CALL ZAXPY( product(nx), w, cres(1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
        enddo
        normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
        normreal = 1.0_dp / normreal
        call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )
      enddo


      call MPI_FILE_WRITE_AT( iu, offset, cres, product(nx)*nfcn, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
    endif


    deallocate( zr, zi, wrk, cres, ilist )
    !
    return
  end subroutine gentoreal
  end subroutine OCEAN_bofx
