! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 06-2019
!
!
module ocean_tmels
  use AI_kinds, only : DP

  implicit none
  private

  logical, save :: fancyFileView = .false.
  logical, save :: is_init = .false.
  integer, save :: tmelsFH

  public :: ocean_tmels_open, ocean_tmels_close, ocean_tmels_calc, ocean_tmels_legacy

  contains


  subroutine ocean_tmels_calc( ikpt, ispin, nValBands, nConBands, valGvecs, conGvecs, &
                               ValUofG, conUofG, nkpts, totValBand, totConBand, startVal, startCon, ierr )
    integer, intent( in ) :: ikpt, ispin, nValBands, nConBands, valGvecs(:,:), conGvecs(:,:), &
                             nkpts, totValBand, totConBand, startVal, startCon
    complex(DP), intent( in ) :: ValUofG(:,:), conUofG(:,:)
    integer, intent(inout) :: ierr

    if( .not. is_init ) then
      ierr = 10123
      return
    endif

    if( .true. ) then
      call tmelsCalc( ikpt, ispin, nValBands, nConBanDs, valGvecs, conGvecs, &
                       ValUofG, conUofG, nkpts, totValBand, totConBand, startVal, startCon, ierr )

    else
!      call tmelsCalcUmklapp
    endif

  end subroutine ocean_tmels_calc


  subroutine tmelsCalc( ikpt, ispin, nValBands, nConBands, valGvecs, conGvecs, &
                        ValUofG, conUofG, nkpts, totValBand, totConBand, startVal, startCon, ierr )
    use ocean_mpi, only : MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, MPI_OFFSET_KIND
    use prep_system, only : psys, params, physical_system, system_parameters
    integer, intent( in ) :: ikpt, ispin, nValBands, nConBands, valGvecs(:,:), conGvecs(:,:), &
                             nkpts, totValBand, totConBand, startVal, startCon
    complex(DP), intent( in ) :: ValUofG(:,:), conUofG(:,:)
    integer, intent(inout) :: ierr

    complex(DP), allocatable :: tmels(:,:)
    integer, allocatable :: gVecMap(:)
    integer :: i, j, k, nGV, nGC
    integer(MPI_OFFSET_KIND) :: offset

    if( ikpt .lt. 1 ) then
      i = 0
      ! zero length write, but do need a valid buffer
      if( fancyFileView ) then
        call MPI_FILE_WRITE_ALL( tmelsFH, valUofG, i, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
      endif
      return
    endif

    allocate( tmels( nValBands, nConBands ) )
    ! for openmp move to first touch
    tmels = 0.0_DP
    nGV = size( valGvecs, 2 )
    nGC = size( conGvecs, 2 )

    allocate( gVecMap( nGV ) )
    do i = 1, nGV
      do j = i, nGC
        if( valGvecs( 1, i ) .eq. conGvecs( 1, j ) &
            .and. valGvecs( 2, i ) .eq. conGvecs( 2, j ) & 
            .and. valGvecs( 3, i ) .eq. conGvecs( 3, j ) ) then
          gVecMap( i ) = j
          goto 11
        endif
      enddo
      do j = 1, min( nGC, i - 1 )
        if( valGvecs( 1, i ) .eq. conGvecs( 1, j ) .and. valGvecs( 2, i ) .eq. conGvecs( 2, j ) &
            .and. valGvecs( 3, i ) .eq. conGvecs( 3, j ) ) then
          gVecMap( i ) = j
          goto 11
        endif
      enddo
      gVecMap( i ) = 0
11      continue
    enddo
    
    do i = 1, nConBands
      do j = 1, nValBands
        do k = 1, nGV
          if( gVecMap( k ) .ne. 0 ) then
            tmels( j, i ) = tmels( j, i ) + conUofG( gVecMap( k ), i ) * conjg( valUofG( k, j ) )
!            tmels( j, i ) = tmels( j, i ) + conjg( conUofG( gVecMap( k ), i ) ) * valUofG( k, j )
          endif
        enddo
      enddo
    enddo
  

    deallocate( gVecMap )
    i = nValBands * nConBands
    if( fancyFileView ) then
      call MPI_FILE_WRITE_ALL( tmelsFH, tmels, i, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
    else
      offset = ( ( ispin - 1 ) * nkpts + ikpt - 1 ) * totValBand * totConBand
      offset = offset + ( startCon - 1 ) * totValBand
!      offset = offset * 16
      write(6,*) ikpt, offset
      call MPI_FILE_WRITE_AT( tmelsFH, offset, tmels, i, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )

#if 0
      do i = 1, nConBands
        do j = 1, nValBands
          write(999,*) real(tmels(j,i)), aimag(tmels(j,i))
        enddo
      enddo
#endif

    endif
  
    deallocate( tmels ) 

  end subroutine tmelsCalc

    
  ! nValBands = number of valence bands in calc
  ! nConBands = number of conduction bands in calc
  ! myConBands = number that this processor is use
  ! myConBandStart => The block of bands this processor is treating might start somewhere other than the first
  ! myKStart -> Like the previous might not start on the first k-point
  ! kStride = this is the number of k-point reading pools
  !
  ! If there are 4 i/o pools, I might read the 1, 5, 9, k-points
  subroutine ocean_tmels_open( nValBands, nConBands, myConBands, myConBandStart, myNK, myKStart, kStride, ierr )
    use ocean_mpi, only : myid, comm, &
                          MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &
                          MPI_OFFSET_KIND, MPI_DOUBLE_COMPLEX
    use prep_system, only : system_parameters, params

    integer, intent( in ) :: nValBands, nConBands, myConBands, myConBandStart, myNK, myKStart, kStride
    integer, intent( inout ) :: ierr

#ifdef MPI_F08
    type( MPI_DATATYPE ):: fileType
#else
    integer :: fileType
#endif
    integer(MPI_OFFSET_KIND) :: offset
    integer :: fflags, sizeofcomplex, vcBandBlock, vectorStride
    complex(DP) :: dumz


    if( myid .eq. 0 ) then
      open(unit=99, file='tmels.info', form='formatted', status='unknown' )
      write(99,*) params%brange(2)-params%brange(1)+1, params%brange(3), params%brange(4), params%nkpts, params%nspin
      close(99)
    endif
#ifdef MPI
    fflags = IOR( MPI_MODE_WRONLY, MPI_MODE_CREATE )
    fflags = IOR( fflags, MPI_MODE_UNIQUE_OPEN )
    call MPI_FILE_OPEN( comm, "ptmels.dat", fflags, MPI_INFO_NULL, tmelsFH, ierr )
    if( ierr .ne. 0 ) return

    vcBandBlock = nValBands * myConBands
    vectorStride = nValBands * nConBands * kStride
    write(1000+myid, * ) myNK, vcBandBlock, vectorStride
    flush(1000+myid)
    call MPI_BARRIER( comm, ierr )
    call MPI_TYPE_VECTOR( myNK, vcBandBlock, vectorStride, MPI_DOUBLE_COMPLEX, fileType, ierr )
    if( ierr .ne. 0 ) return
    call MPI_TYPE_COMMIT( fileType, ierr )
    if( ierr .ne. 0 ) return


    offset = ( myKStart - 1 ) * int( nValBands, MPI_OFFSET_KIND ) * int( nConBands, MPI_OFFSET_KIND )
    offset = offset + ( myConBandStart - 1 ) * nValBands 
!    call MPI_SIZEOF( dumz, sizeofcomplex, ierr )
!    if( ierr .ne. 0 ) return
    sizeofcomplex = 16
    offset = offset * int(sizeofcomplex, MPI_OFFSET_KIND )
    write(1000+myid, * ) 'offset', offset
    flush(1000+myid)

    if( fancyFileView ) then
      call MPI_FILE_SET_VIEW( tmelsFH, offset, MPI_DOUBLE_COMPLEX, fileType, "native", MPI_INFO_NULL, ierr )
    else
      offset = 0
      call MPI_FILE_SET_VIEW( tmelsFH, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, "native", MPI_INFO_NULL, ierr )
    endif
    if( ierr .ne. 0 ) return

    call MPI_TYPE_FREE( fileType, ierr )

#else
    open( file='ptmels', form='unformatted', status='unknown', newunit=fh, access='stream' )
#endif
    

    is_init = .true.
  end subroutine ocean_tmels_open

  subroutine ocean_tmels_close( ierr )
    integer, intent( inout ) :: ierr
#ifdef MPI
    call MPI_FILE_CLOSE( tmelsFH, ierr )
#else
    close(fh)
#endif
    is_init = .false.
  end subroutine


!# define TESTMPI
  subroutine ocean_tmels_legacy( nv, nc, nk, nspin, ierr )
    use ocean_mpi
    integer, intent( in ) :: nv, nc, nk, nspin
    integer, intent( inout ) :: ierr

    complex(dp), allocatable :: ttt( :, : )
    integer :: ispn, ik, fh, elements
    integer(8) :: offset

#ifdef TESTMPI
    call MPI_FILE_OPEN( comm, 'ptmels.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    offset = 0
    call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, &
                            'native', MPI_INFO_NULL, ierr)
    if( ierr .ne. MPI_SUCCESS ) return
#else
    if( myid .eq. root ) open(unit=98, form='unformatted', access='stream', file='ptmels.dat')
#endif
    allocate( ttt( nv, nc ) )

    if( myid .eq. root ) then
    open( unit=99, file='tmels', form='formatted' )
!JTV if we want to allow spin-flipping transitions we'll need to do it here
    do ispn = 1, nspin
      do ik = 1, nk
        ! Read in the tmels
        offset = nv * nc  * ( ik - 1 )
        elements = nv * nc
#ifdef TESTMPI
        call MPI_FILE_READ_AT( FH, offset, ttt, elements, &
                                   MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
        if( ierr .ne. MPI_SUCCESS) return
#else
        read(98) ttt
#endif

        write(99, '(2E23.15)' ) ttt(:,:)

      enddo
    enddo

    close(99)
    endif

#ifdef TESTMPI
    call MPI_FILE_CLOSE( fh, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#else
    if( myid .eq. root )close(98)
#endif

  end subroutine ocean_tmels_legacy

end module ocean_tmels
