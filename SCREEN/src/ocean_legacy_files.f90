! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module ocean_legacy_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif


  implicit none
  private
  save


  logical :: is_init = .false.
  integer :: nfiles
  integer, allocatable :: file_indx( : )
  character( len=12 ), allocatable :: file_names( : )

  integer :: brange(4)
  integer :: kpts(3)
  integer :: nspin 


#ifdef MPI_F08
  type( MPI_COMM ) :: inter_comm
  type( MPI_COMM ) :: pool_comm
#else
  integer :: inter_comm
  integer :: pool_comm
#endif

  integer :: inter_myid
  integer :: inter_nproc
  integer, parameter :: inter_root = 0

  integer :: npool
  integer :: mypool
  integer :: pool_myid
  integer :: pool_nproc
  integer, parameter :: pool_root = 0
  
  integer :: pool_nbands

  public :: olf_read_init, olf_read_at_kpt, olf_clean, olf_read_energies, olf_get_ngvecs_at_kpt, &
            olf_read_energies_single
  public :: olf_kpts_and_spins, olf_return_my_bands, olf_is_my_kpt
  public :: olf_nprocPerPool, olf_getPoolIndex, olf_getBandsForPoolID, olf_returnGlobalID

  contains

  pure function olf_nprocPerPool() result( nproc )
    integer :: nproc
    
    nproc = pool_nproc
  end function olf_nprocPerPool

  pure function olf_getPoolIndex( ispin, ikpt ) result( poolIndex )
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    integer :: kptCounter
    !
    kptCounter = ikpt + ( ispin - 1 ) * product(kpts(:))
    poolIndex = mod( kptCounter, npool )
  end function olf_getPoolIndex

  pure function olf_getBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    bands_remain = brange(4)-brange(3)+brange(2)-brange(1)+2

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function olf_getBandsForPoolID

  pure function olf_returnGlobalID( poolIndex, poolID ) result( globalID )
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID

    globalID = poolIndex * pool_nproc + poolID
  end function olf_returnGlobalID

  subroutine olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    integer, intent( in ) :: ikpt, ispin
    logical, intent( out ) :: is_kpt
    integer, intent( inout ) :: ierr
    !
    if( .not. is_init ) then
      ierr = 2
      return
    endif
    if( ispin .lt. 1 .or. ispin .gt. nspin ) then
      ierr = 3 
      return
    endif
    if( ikpt .lt. 1 .or. ikpt .gt. product(kpts(:)) ) then
      ierr = 4
      return
    endif

!    if( mod( ikpt + (ispin-1)*product(kpts(:)), npool ) .eq. mypool ) then
    if( olf_getPoolIndex( ispin, ikpt ) .eq. mypool ) then
      is_kpt = .true.
    else
      is_kpt = .false.
    endif

  end subroutine olf_is_my_kpt

  subroutine olf_return_my_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr
    
    if( .not. is_init ) ierr = 1
    nbands = pool_nbands
  end subroutine olf_return_my_bands

  integer function olf_kpts_and_spins()
    olf_kpts_and_spins = product(kpts(:)) * nspin / npool
    return
  end function olf_kpts_and_spins

  subroutine olf_read_energies_single( myid, root, comm, energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_BCAST, MPI_DOUBLE_PRECISION
#endif
    integer, intent( in ) :: myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    real( DP ), intent( out ) :: energies( :, :, : )
    integer, intent( inout ) :: ierr
    !
    integer :: ispn, ikpt, bstop1, bstart2, bstop2

    if( is_init .eqv. .false. ) then
      write(6,*) 'olf_read_energies_single called but was not initialized'
      ierr = 4
      return
    endif 

    if( size( energies, 1 ) .lt. ( brange(4)-brange(1)+1 ) ) then
      write(6,*) 'Error! Band dimension of energies too small in olf_read_energies_single'
      ierr = 1
      return
    endif
    if( size( energies, 2 ) .lt. product(kpts(:)) ) then
      write(6,*) 'Error! K-point dimension of energies too small in olf_read_energies_single'
      ierr = 2
      return
    endif
    if( size( energies, 3 ) .lt. nspin ) then
      write(6,*) 'Error! Spin dimension of energies too small in olf_read_energies_single'
      ierr = 3
      return
    endif


    if( myid .eq. root ) then
      bstop1  = brange(2) - brange(1) + 1
      bstart2 = brange(3) - brange(1) + 1
      bstop2  = brange(4) - brange(1) + 1

      open( unit=99, file='enkfile', form='formatted', status='old' )
      do ispn = 1, nspin
        do ikpt = 1, product(kpts(:))
          read( 99, * ) energies( 1 : bstop1, ikpt, ispn )
          read( 99, * ) energies( bstart2 : bstop2, ikpt, ispn )
        enddo
      enddo
      close( 99 )
    endif

#ifdef MPI
    call MPI_BCAST( energies, size(energies), MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine olf_read_energies_single



  subroutine olf_read_energies( myid, root, comm, nbv, nbc, nkpts, nspns, val_energies, con_energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_BCAST, MPI_DOUBLE_PRECISION
#endif
    integer, intent( in ) :: myid, root, nbv, nbc, nkpts, nspns
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    real( DP ), intent( out ) :: val_energies( nbv, nkpts, nspns ), con_energies( nbv, nkpts, nspns )
    integer, intent( inout ) :: ierr
    !
    integer :: ispn, ikpt

    if( myid .eq. root ) then
      open( unit=99, file='enkfile', form='formatted', status='old' )
      do ispn = 1, nspns
        do ikpt = 1, nkpts
          read( 99, * ) val_energies( :, ikpt, ispn )
          read( 99, * ) con_energies( :, ikpt, ispn )
        enddo
      enddo
      close( 99 )
    endif

#ifdef MPI
    call MPI_BCAST( val_energies, nbv*nkpts*nspns, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( con_energies, nbc*nkpts*nspns, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine olf_read_energies



  subroutine olf_clean( ierr )
    integer, intent( inout ) :: ierr
    !
    nfiles = 0
    if( allocated( file_indx ) ) deallocate( file_indx )
    if( allocated( file_names ) ) deallocate( file_names )

#ifdef MPI
    call MPI_COMM_FREE( pool_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_FREE( inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine olf_clean

  ! Read the universal little files
  subroutine olf_read_init( comm, ierr )
    use ocean_mpi, only :  MPI_BCAST, MPI_INTEGER, MPI_CHARACTER
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr
    !
    integer :: i
    !
    ! Set the comms for file handling
    call MPI_COMM_DUP( comm, inter_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( inter_comm, inter_nproc, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( inter_comm, inter_myid, ierr )
    if( ierr .ne. 0 ) return


    if( inter_myid .eq. inter_root ) then
      open( unit=99, file='masterwfile', form='formatted', status='old' )
      read( 99, * ) nfiles
      close( 99 )
      allocate( file_names( nfiles ), file_indx( nfiles ) )
      open( unit=99, file='listwfile', form='formatted', status='old' )
      do i = 1, nfiles
        read( 99, * ) file_indx( i ), file_names( i )
      enddo
      close( 99 )

      open(unit=99,file='brange.ipt', form='formatted', status='old' )
      read(99,*) brange(:)
      close(99)
      open(unit=99,file='kmesh.ipt', form='formatted', status='old' )
      read(99,*) kpts(:)
      close(99)
      open(unit=99,file='nspin', form='formatted', status='old' )
      read(99,*) nspin
      close(99)
    endif

#ifdef MPI
    call MPI_BCAST( nfiles, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    if( inter_myid .ne. inter_root ) then
      allocate( file_names( nfiles ), file_indx( nfiles ) )
    endif

    call MPI_BCAST( file_indx, nfiles,  MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( file_names, 12 * nfiles, MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( brange, 4, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( kpts, 3, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( nspin, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    call set_pools( ierr )
    if( ierr .ne. 0 ) return

    call writeDiagnostics( )

!    write(6,*) 'olf_read_init was successful'
    is_init = .true.
  
  end subroutine  olf_read_init 

  subroutine writeDiagnostics( )
    if( inter_myid .eq. inter_root ) then
      write( 6, '(A)' ) "    #############################"
      write( 6, '(A,I8)' ) " Npools:      ", npool
      write( 6, '(A,I8)' ) " Nprocs/pool: ", pool_nproc
    endif

    write(1000+inter_myid, '(A)' ) "    #############################"
    write(1000+inter_myid, '(A,I8)' ) " Npools:      ", npool
    write(1000+inter_myid, '(A,I8)' ) " Nprocs/pool: ", pool_nproc
    write(1000+inter_myid, '(A,I8)' ) " My pool:     ", mypool
    write(1000+inter_myid, '(A,I8)' ) " My pool id:  ", pool_myid
    write(1000+inter_myid, '(A,I8)' ) " My bands:    ", pool_nbands


    write(1000+inter_myid, '(A)' ) "    #############################"
    flush(1000+inter_myid)
  end subroutine 

  subroutine set_pools( ierr )
    integer, intent( inout ) :: ierr
    !
    integer :: i

    if( nfiles .ge. inter_nproc ) then
      mypool = inter_myid
      npool = inter_nproc

    else
      do i = 2, inter_nproc
        if( mod( inter_nproc, i ) .eq. 0 ) then
          if( inter_myid .eq. 0 ) write(6,*) i, inter_nproc
          if( nfiles .ge. (inter_nproc/i) ) then
            npool = inter_nproc/i
            mypool = inter_myid/ i 
#ifdef DEBUG
            write(6,*) '*', inter_myid, npool, mypool
#endif
            goto 11
          endif
        endif
      enddo
11    continue
    endif
    
    call MPI_COMM_SPLIT( inter_comm, mypool, inter_myid, pool_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_COMM_RANK( pool_comm, pool_myid, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( pool_comm, pool_nproc, ierr )
    if( ierr .ne. 0 ) return


    pool_nbands = olf_getBandsForPoolID( pool_myid )
!    nbands_left = brange(4)-brange(3)+brange(2)-brange(1)+2
!    do i = 0, pool_nproc-1
!      nbands = nbands_left / ( pool_nproc - i )
!      if( i .eq. pool_myid ) pool_nbands = nbands
!      nbands_left = nbands_left - nbands
!    enddo

  end subroutine set_pools


  subroutine olf_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
    use OCEAN_mpi
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr
    !
    integer :: i, itarg
    logical :: is_kpt

    call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    if( ierr .ne. 0 ) return
    if( is_kpt .eqv. .false. ) then 
      gvecs = 0
      return
    endif

    if( pool_myid .eq. pool_root ) then
      ! determine which file we want (usually we have 1 kpt per file)
      itarg = 0
      if( file_indx( nfiles ) .gt. ikpt ) then
        do i = min( ikpt, nfiles-1), 1, -1
          if( file_indx( i ) .le. ikpt .and. file_indx( i + 1 ) .gt. ikpt ) then
            itarg = i
            exit
          endif
        enddo
      else
        itarg = nfiles
      endif

      if( itarg .eq. 0 ) then
        write( 6, * ) "Couldn't figure out which file to open", ikpt
        ierr = 12
        goto 111
      endif
      open( unit=99,file=file_names(itarg), form='unformatted', status='old' )
      read( 99 ) gvecs
      close( 99 )
    endif

!    write(6,*) 'gvecs', pool_root, pool_comm
111 continue
    
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, i )
    if( ierr .ne. 0 ) return
    if( i .ne. 0 ) then
      ierr = i
      return
    endif

    call MPI_BCAST( gvecs, 1, MPI_INTEGER, pool_root, pool_comm, ierr )
    if( ierr .ne. 0 ) return
#endif
!    write(6,*) 'gvecs', pool_root, pool_comm

  end subroutine olf_get_ngvecs_at_kpt

  subroutine olf_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_IBCAST, MPI_INTEGER, MPI_BCAST, MPI_IRSEND, MPI_IRECV, &
                          MPI_DOUBLE_PRECISION, MPI_STATUSES_IGNORE, MPI_CANCEL, myid, MPI_REQUEST_NULL
!    use OCEAN_mpi
#endif
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
    real( DP ), allocatable, dimension( :, : ) :: re_wvfn, im_wvfn
    integer, allocatable, dimension( :, : ) :: trans_gvecs
    integer :: test_gvec, itarg, nbands_to_send, nr, ierr_, nbands, id, start_band
#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: requests( : )
#else
    integer, allocatable :: requests( : )
#endif

    if( olf_getPoolIndex( ispin, ikpt ) .ne. mypool ) return
    !
    if( olf_getBandsForPoolID( pool_myid ) .ne. my_bands ) then
      ierr = 1
      return
    endif

    nbands = brange(4)-brange(3)+brange(2)-brange(1)+2

    if( pool_myid .eq. pool_root ) then
      itarg = wvfn_file_indx( ikpt )
      if( itarg .lt. 1 ) then
        ierr = -1
        goto 10
      endif

      open( unit=99, file=file_names( itarg ), form='unformatted', status='old' )
      read( 99 ) test_gvec
      if( test_gvec .ne. ngvecs ) then
        ierr = -2
      endif

      ! Error synch. Also ensures that the other procs have posted their recvs
10    continue
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif

      allocate( trans_gvecs( ngvecs, 3 ) )
      read( 99 ) trans_gvecs
      gvecs = transpose( trans_gvecs )
      deallocate( trans_gvecs )

      nr = 2 * ( pool_nproc - 1 ) + 1 
!      nr = 2 * pool_nproc
      allocate( requests( 0:nr ) )
      requests(:) = MPI_REQUEST_NULL
#ifdef MPI
      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, requests( nr ), ierr )
#endif
      allocate( re_wvfn( ngvecs, nbands ), im_wvfn( ngvecs, nbands ) )

      write(1000+myid,*) '***Reading k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands


      read( 99 ) re_wvfn
      start_band = 1

#ifdef MPI
      ! loop over each proc in this pool to send wavefunctions
      do id = 0, pool_nproc - 1
        nbands_to_send = olf_getBandsForPoolID( id )

        ! don't send if I am me
        if( id .ne. pool_myid ) then
          write(1000+myid,'(A,3(1X,I8))') '   Sending ...', id, start_band, nbands_to_send
          call MPI_IRSEND( re_wvfn( 1, start_band ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         id, 1, pool_comm, requests( id ), ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
        else
          write(1000+myid,'(A,3(1X,I8))') "   Don't Send: ", start_band, nbands_to_send, my_bands
        endif

        start_band = start_band + nbands_to_send
      enddo
#endif

      read( 99 ) im_wvfn
      start_band = 1

#ifdef MPI
     ! loop over each proc in this pool to send imag wavefunctions
      do id = 0, pool_nproc - 1
        nbands_to_send = olf_getBandsForPoolID( id )

        ! don't send if I am me
        if( id .ne. pool_myid ) then
          call MPI_IRSEND( im_wvfn( 1, start_band ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         id, 2, pool_comm, requests( id + pool_nproc - 1), ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
        endif

        start_band = start_band + nbands_to_send
      enddo
#endif

      close( 99 )

    else
      nr = 3
      allocate( requests( nr ) )
      requests(:) = MPI_REQUEST_NULL
      write(1000+myid,*) '***Receiving k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands


      allocate( re_wvfn( ngvecs, my_bands ), im_wvfn( ngvecs, my_bands ) )
#ifdef MPI
      call MPI_IRECV( re_wvfn, ngvecs*my_bands, MPI_DOUBLE_PRECISION, pool_root, 1, pool_comm, & 
                      requests( 1 ), ierr )
      call MPI_IRECV( im_wvfn, ngvecs*my_bands, MPI_DOUBLE_PRECISION, pool_root, 2, pool_comm, &
                      requests( 2 ), ierr )

      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) then
        call MPI_CANCEL( requests( 1 ) , ierr )
        call MPI_CANCEL( requests( 2 ) , ierr )
        ierr = 5
        return
      endif

      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, requests( 3 ), ierr )
#endif


    endif


    call MPI_WAITALL( nr, requests, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    wfns( :, : ) = cmplx( re_wvfn( :, 1:my_bands ), im_wvfn( :, 1:my_bands ), DP )
    deallocate( re_wvfn, im_wvfn, requests )

    write(1000+myid,*) '***Finishing k-point: ', ikpt, ispin
    call MPI_BARRIER( pool_comm, ierr )

  end subroutine olf_read_at_kpt


  integer function wvfn_file_indx( ikpt )
    integer, intent( in ) :: ikpt
    !
    integer :: i
    !
    wvfn_file_indx = 0
    if( file_indx( nfiles ) .gt. ikpt ) then
      do i = min( ikpt, nfiles-1), 1, -1
        if( file_indx( i ) .le. ikpt .and. file_indx( i + 1 ) .gt. ikpt ) then
          wvfn_file_indx = i
          return
        endif
      enddo
    else
      wvfn_file_indx = nfiles
    endif

    return
  end function

end module ocean_legacy_files
