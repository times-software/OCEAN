! Copyright (C) 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module ocean_qe54_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif


  implicit none
  private
  save


  logical :: is_init = .false.
  character( len=128 ) :: prefix

  integer :: bands(2)
  integer :: kpts(3)
  integer :: nspin 
  integer :: nfiles


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

  public :: qe54_read_init, qe54_read_at_kpt, qe54_clean, qe54_read_energies, qe54_get_ngvecs_at_kpt, &
            qe54_read_energies_single
  public :: qe54_kpts_and_spins, qe54_return_my_bands, qe54_is_my_kpt
  public :: qe54_nprocPerPool, qe54_getPoolIndex, qe54_getBandsForPoolID, qe54_returnGlobalID

  contains

  pure function qe54_nprocPerPool() result( nproc )
    integer :: nproc
    
    nproc = pool_nproc
  end function qe54_nprocPerPool

  pure function qe54_getPoolIndex( ispin, ikpt ) result( poolIndex )
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    integer :: kptCounter
    !
    kptCounter = ikpt + ( ispin - 1 ) * product(kpts(:))
    poolIndex = mod( kptCounter, npool )
  end function qe54_getPoolIndex

  pure function qe54_getBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    bands_remain = bands(2)-bands(1)+1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function qe54_getBandsForPoolID

  pure function qe54_returnGlobalID( poolIndex, poolID ) result( globalID )
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID

    globalID = poolIndex * pool_nproc + poolID
  end function qe54_returnGlobalID

  subroutine qe54_is_my_kpt( ikpt, ispin, is_kpt, ierr )
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
    if( qe54_getPoolIndex( ispin, ikpt ) .eq. mypool ) then
      is_kpt = .true.
    else
      is_kpt = .false.
    endif

  end subroutine qe54_is_my_kpt

  subroutine qe54_return_my_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr
    
    if( .not. is_init ) ierr = 1
    nbands = pool_nbands
  end subroutine qe54_return_my_bands

  integer function qe54_kpts_and_spins()
    qe54_kpts_and_spins = product(kpts(:)) * nspin / npool
    return
  end function qe54_kpts_and_spins

  subroutine qe54_read_energies_single( myid, root, comm, energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_DOUBLE_PRECISION
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
    integer :: ispn, ikpt, i, bstop2 !bstop1, bstart2, bstop2
    real(DP) :: dumf
    character(len=128) :: lineBurn

    if( is_init .eqv. .false. ) then
      write(6,*) 'qe54_read_energies_single called but was not initialized'
      ierr = 4
      return
    endif 

    if( size( energies, 1 ) .lt. ( bands(2)-bands(1)+1 ) ) then
      write(6,*) 'Error! Band dimension of energies too small in qe54_read_energies_single'
      ierr = 1
      return
    endif
    if( size( energies, 2 ) .lt. product(kpts(:)) ) then
      write(6,*) 'Error! K-point dimension of energies too small in qe54_read_energies_single'
      ierr = 2
      return
    endif
    if( size( energies, 3 ) .lt. nspin ) then
      write(6,*) 'Error! Spin dimension of energies too small in qe54_read_energies_single'
      ierr = 3
      return
    endif

    if( myid .eq. root ) then
      bstop2  = bands(2) - bands(1) + 1

      do ispn = 1, nspin
        do ikpt = 1, product(kpts(:))

          open( unit=99, file=qe54_eigFile(ikpt,ispn), form='formatted', status='old' )
          do i = 1, 9
            read(99,*) !lineBurn
!            write(6,*) lineBurn
          enddo
!          do i = 1, bands( 1 ) - 1
!            read(99,*) dumf
!            write(6,*) dumf
!          enddo
          do i = 1, bstop2
            read(99,*) energies( i, ikpt, ispn )
          enddo

          close( 99 )
        enddo
      enddo
      energies(:,:,:) = energies(:,:,:) * 2.0_DP
    endif

#ifdef MPI
    call MPI_BCAST( energies, size(energies), MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine qe54_read_energies_single

  pure function qe54_eigFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin
    character(len=512) :: fileName

    if( nspin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'eigenval.xml' 
    elseif( ispin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'eigenval1.xml'
    else
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'eigenval2.xml'
    endif

  end function qe54_eigFile

  pure function qe54_gkvFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin
    character(len=512) :: fileName

    write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'gkvectors.dat'

  end function qe54_gkvFile

  pure function qe54_evcFile( ikpt, ispin ) result( fileName )
    integer, intent ( in ) :: ikpt, ispin
    character(len=512) :: fileName

    if( nspin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'evc.dat'
    elseif( ispin .eq. 1 ) then
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'evc1.dat'
    else
      write( fileName, '(a,a1,i5.5,a1,a)' ) trim( prefix ), 'K', ikpt, '/', 'evc2.dat'
    endif

  end function qe54_evcFile

  subroutine qe54_read_energies( myid, root, comm, nbv, nbc, nkpts, nspns, val_energies, con_energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_DOUBLE_PRECISION
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

  end subroutine qe54_read_energies



  subroutine qe54_clean( ierr )
    integer, intent( inout ) :: ierr
    !
    nfiles = 0

#ifdef MPI
    call MPI_COMM_FREE( pool_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_FREE( inter_comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine qe54_clean

  ! Read the universal little files
  subroutine qe54_read_init( comm, ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_CHARACTER
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr
    !
    integer :: i, brange(4)
    logical :: ex
    character(len=128) :: tmp
    !
    ! Set the comms for file handling
    call MPI_COMM_DUP( comm, inter_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( inter_comm, inter_nproc, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( inter_comm, inter_myid, ierr )
    if( ierr .ne. 0 ) return


    if( inter_myid .eq. inter_root ) then
      open( unit=99, file='prefix', form='formatted', status='old' )
      read(99,*)  tmp
      close( 99 )
      write( prefix, '(a,a,a)' ) 'Out/', trim(tmp), '.save/'

      inquire( file='bands.ipt', exist=ex )
      if( ex ) then
        open(unit=99,file='bands.ipt', form='formatted', status='old' )
        read(99,*) bands(:)
        close(99)
      else
        open(unit=99,file='brange.ipt', form='formatted', status='old' )
        read(99,*) brange(:)
        close(99)
        bands(2) = brange(4)
        bands(1) = brange(1)
      endif
      open(unit=99,file='kmesh.ipt', form='formatted', status='old' )
      read(99,*) kpts(:)
      close(99)
      open(unit=99,file='nspin', form='formatted', status='old' )
      read(99,*) nspin
      close(99)

      nfiles = nspin * product(kpts(:) )
    endif

#ifdef MPI
    call MPI_BCAST( nfiles, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( bands, 2, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( kpts, 3, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( nspin, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( prefix, len(prefix), MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

#endif

    call set_pools( ierr )
    if( ierr .ne. 0 ) return

    call writeDiagnostics( )

!    write(6,*) 'qe54_read_init was successful'
    is_init = .true.
  
  end subroutine  qe54_read_init 

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


    pool_nbands = qe54_getBandsForPoolID( pool_myid )
!    nbands_left = brange(4)-brange(3)+brange(2)-brange(1)+2
!    do i = 0, pool_nproc-1
!      nbands = nbands_left / ( pool_nproc - i )
!      if( i .eq. pool_myid ) pool_nbands = nbands
!      nbands_left = nbands_left - nbands
!    enddo

  end subroutine set_pools


  subroutine qe54_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
    use OCEAN_mpi
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr
    !
    integer :: i, crap
    logical :: is_kpt

    call qe54_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    if( ierr .ne. 0 ) return
    if( is_kpt .eqv. .false. ) then 
      gvecs = 0
      return
    endif

    if( pool_myid .eq. pool_root ) then
      open( unit=99,file=qe54_gkvFile( ikpt, ispin ), form='unformatted', status='old' )
      do i = 1, 12
        read( 99 )
      enddo
      read(99) crap, gvecs
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

  end subroutine qe54_get_ngvecs_at_kpt

  subroutine qe54_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_STATUSES_IGNORE, myid, MPI_REQUEST_NULL
!    use OCEAN_mpi
#endif
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
    complex( DP ), allocatable, dimension( :, : ) :: cmplx_wvfn
    integer :: test_gvec, itarg, nbands_to_send, nr, ierr_, nbands, id, start_band, crap, i
#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: requests( : )
#else
    integer, allocatable :: requests( : )
#endif

    if( qe54_getPoolIndex( ispin, ikpt ) .ne. mypool ) return
    !
    if( qe54_getBandsForPoolID( pool_myid ) .ne. my_bands ) then
      ierr = 1
      return
    endif

    nbands = bands(2)-bands(1)+1

    if( pool_myid .eq. pool_root ) then
    
      ! get gvecs first
      open( unit=99, file=qe54_gkvFile( ikpt, ispin), form='unformatted', status='old' )
      do i = 1, 12
        read( 99 )
      enddo
      read(99) crap, test_gvec
      if( test_gvec .ne. ngvecs ) then
        ierr = -2
        write(6,*) test_gvec, ngvecs
        return
      endif
      do i = 1, 5
        read(99)
      enddo

      ! Error synch. Also ensures that the other procs have posted their recvs
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif

      do i = 1, 19
        read( 99 )
      enddo

      read( 99 ) crap, gvecs( :, 1:ngvecs )
      close( 99 )

      nr = pool_nproc 
!      nr = 2 * pool_nproc
      allocate( requests( 0:nr ) )
      requests(:) = MPI_REQUEST_NULL
#ifdef MPI
      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, requests( nr ), ierr )
#endif
      allocate( cmplx_wvfn( ngvecs, nbands ) )

      write(1000+myid,*) '***Reading k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands

      open( unit=99, file=qe54_evcFile( ikpt, ispin), form='unformatted', status='old' )
      do i = 1, 12
        read(99)
      enddo

      ! if we are skipping bands, do that here
!      do i = 1, bands(1)-1
!        read( 99 ) crap, cmplx_wvfn( 1, 1 )
!      enddo

      start_band = 1

#ifdef MPI
      ! loop over each proc in this pool to send wavefunctions
      do id = 0, pool_nproc - 1
        nbands_to_send = qe54_getBandsForPoolID( id )

        do i = start_band, nbands_to_send + start_band - 1
          read(99)
          read(99)
          read( 99 ) crap, cmplx_wvfn( :, i )
          read(99)
          read(99)
        enddo

        ! don't send if I am me
        if( id .ne. pool_myid ) then
          write(1000+myid,'(A,3(1X,I8))') '   Sending ...', id, start_band, nbands_to_send
          call MPI_IRSEND( cmplx_wvfn( 1, start_band ), nbands_to_send*ngvecs, MPI_DOUBLE_COMPLEX, &
                         id, 1, pool_comm, requests( id ), ierr )
          ! this might not sync up
          if( ierr .ne. 0 ) return
        else
          write(1000+myid,'(A,3(1X,I8))') "   Don't Send: ", start_band, nbands_to_send, my_bands
          wfns( :, : ) = cmplx_wvfn( :, start_band : nbands_to_send + start_band - 1 )
        endif

        start_band = start_band + nbands_to_send
      enddo
#endif

      close( 99 )

      nr=nr+1
    else
      nr = 2
      allocate( requests( nr ), cmplx_wvfn( 1, 1 ) )
      requests(:) = MPI_REQUEST_NULL
      write(1000+myid,*) '***Receiving k-point: ', ikpt, ispin
      write(1000+myid,*) '   Ngvecs: ', ngvecs
      write(1000+myid,*) '   Nbands: ', nbands


#ifdef MPI
      call MPI_IRECV( wfns, ngvecs*my_bands, MPI_DOUBLE_COMPLEX, pool_root, 1, pool_comm, & 
                      requests( 1 ), ierr )

      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, requests( 2 ), ierr )
      call MPI_BCAST( ierr, 1, MPI_INTEGER, pool_root, pool_comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) then
        call MPI_CANCEL( requests( 1 ) , ierr )
        call MPI_CANCEL( requests( 2 ) , ierr )
        ierr = 5
        return
      endif

#endif


    endif


    call MPI_WAITALL( nr, requests, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    deallocate( cmplx_wvfn, requests )

    write(1000+myid,*) '***Finishing k-point: ', ikpt, ispin
    call MPI_BARRIER( pool_comm, ierr )

  end subroutine qe54_read_at_kpt



end module ocean_qe54_files
