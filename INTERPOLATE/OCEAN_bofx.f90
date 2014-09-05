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
  USE mp, ONLY : mp_sum, mp_bcast
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  use hamq_shirley
  use shirley_ham_input, only : debug, band_subset

  implicit none

  real(dp),parameter ::  eps=1.d-10
  
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)


  integer :: ibnd  
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  
  integer :: i, ipt, ig
  integer,allocatable :: igk_l2g(:,:), sub_mill(:,:)

  integer :: xmesh( 3 )
  logical :: bt_logical, invert_xmesh
  character(len=3) :: dum_c

  integer :: iuntmp, iunout

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

    iunout = freeunit()
    open(iunout,file='u1.dat', form='unformatted', status='unknown' )
    rewind iunout
  endif
  call mp_bcast( xmesh, ionode_id )
  call mp_bcast( invert_xmesh, ionode_id )

  write(stdout,*) nbnd, band_subset(1), band_subset(2)
  allocate( sub_mill( 3, npw ) )
  sub_mill = 0
  do ig = 1, npw
    sub_mill( :, ig ) = matmul( g( :, igk( ig ) ), bg( :, : ) )
  enddo

  call gentoreal( xmesh, nbnd, evc, npw, mill, iunout, invert_xmesh )

  if( ionode ) close( iunout )
  deallocate( sub_mill )


  return

  contains


  subroutine gentoreal( nx, nfcn, fcn, ng, gvec, iu, invert_xmesh )
    use kinds, only : dp
    USE io_global,  ONLY : stdout, ionode
    USE mp, ONLY : mp_sum, mp_max, mp_min
    implicit none
    !
    integer, intent( in ) :: nx( 3 ), nfcn, ng, iu
    integer, intent( in ) :: gvec( 3, ng )
    complex(dp), intent( in ) :: fcn( ng, nfcn )
    logical, intent( in ) :: invert_xmesh
    logical, parameter :: loud = .true.
    !
    integer :: ix, iy, iz, i1, i2, i3, nfft( 3 ), idwrk, igl, igh, nmin, ii
    integer :: fac( 3 ), j, ig, i, locap( 3 ), hicap( 3 ), toreal, torecp, nftot
    real(dp) :: normreal, normrecp
    complex(dp) :: rm1
    character * 80 :: fstr
    !
    integer, parameter :: nfac = 3
    integer, allocatable :: ilist( :, : )
    real(dp), allocatable :: zr( :, :, : ), zi( :, :, : ), wrk( : )
    complex(dp), allocatable :: cres( :, :, : )
    integer, external :: optim
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
      allocate( cres( nx( 3 ), nx( 2 ), nx( 1 ) ) )
    else
      allocate( cres( nx( 1 ), nx( 2 ), nx( 3 ) ) )
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
       if( ionode ) then
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
                    cres( iz, iy, ix ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
                  else
                    cres( ix, iy, iz ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
                  endif
               end do
            end do
         end do 
         write ( iu ) cres
       endif
    end do
    deallocate( zr, zi, wrk, cres, ilist )
    !
    return
  end subroutine gentoreal
  end subroutine OCEAN_bofx