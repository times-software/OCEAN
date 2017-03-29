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


  integer :: nfiles
  integer, allocatable :: file_indx( : )
  character( len=12 ), allocatable :: file_names( : )

  public :: olf_read_init, olf_read_at_kpt, olf_clean, olf_read_energies

  contains

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



  subroutine olf_clean()
    nfiles = 0
    if( allocated( file_indx ) ) deallocate( file_indx )
    if( allocated( file_names ) ) deallocate( file_names )
  end subroutine olf_clean

  ! Read the universal little files
  subroutine olf_read_init( myid, root, comm, ierr )
    use ocean_mpi, only :  MPI_BCAST, MPI_INTEGER, MPI_CHARACTER
    integer, intent( in ) :: myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr
    !
    integer :: i
    !
    if( myid .eq. root ) then
      open( unit=99, file='masterwfile', form='formatted', status='old' )
      read( 99, * ) nfiles
      close( 99 )
      allocate( file_names( nfiles ), file_indx( nfiles ) )
      open( unit=99, file='listwfile', form='formatted', status='old' )
      do i = 1, nfiles
        read( 99, * ) file_indx( i ), file_names( i )
      enddo
      close( 99 )
    endif

#ifdef MPI
    call MPI_BCAST( nfiles, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    if( myid .ne. root ) then
      allocate( file_names( nfiles ), file_indx( nfiles ) )
    endif

    call MPI_BCAST( file_indx, nfiles,  MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( file_names, 12 * nfiles, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif
  
  end subroutine  olf_read_init 


  subroutine olf_get_ngvecs_at_kpt( ikpt, myid, root, comm, gvecs, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_BCAST, MPI_INTEGER
#endif
    integer, intent( in ) :: ikpt, myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr
    !
    integer :: i, itarg
    if( myid .eq. root ) then
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
111 continue
    
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, i )
    if( ierr .ne. 0 ) return
    if( i .ne. 0 ) then
      ierr = i
      return
    endif

    call MPI_BCAST( gvecs, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine olf_get_ngvecs_at_kpt

  subroutine olf_read_at_kpt( ikpt, myid, root, nproc, comm, nbands, my_bands, ngvecs, gvecs, wfns, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_IBCAST, MPI_INTEGER, MPI_BCAST, MPI_IRSEND, MPI_IRECV, &
                          MPI_DOUBLE_PRECISION, MPI_STATUSES_IGNORE, MPI_CANCEL
#endif
    integer, intent( in ) :: ikpt, myid, root, nproc, ngvecs, nbands, my_bands
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
    real( DP ), allocatable, dimension( :, : ) :: re_wvfn, im_wvfn, trans_gvecs
    integer :: test_gvec, itarg, nbands_to_send, iband, i, ir, nr, ierr_
#ifdef MPI_F09
    type( MPI_REQUEST ), allocatable :: requests( : )
#else
    integer, allocatable :: requests( : )
#endif

    !
    if( myid .eq. root ) then
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
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif

      allocate( trans_gvecs( ngvecs, 3 ) )
      read( 99 ) trans_gvecs
      gvecs = transpose( trans_gvecs )
      deallocate( trans_gvecs )

      nr = 2 * ( nproc - 1 ) + 1 
      allocate( requests( nr ) )
      ir = 0
#ifdef MPI
      ir = ir + 1
      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, root, comm, requests( ir ), ierr )
#endif
      allocate( re_wvfn( ngvecs, nbands ), im_wvfn( ngvecs, nbands ) )

      read( 99 ) re_wvfn
#ifdef MPI
      nbands_to_send = nbands / nproc + 1
      iband = 1
      do i = 1, mod( nbands, nproc )
        ir = ir + 1
        call MPI_IRSEND( re_wvfn( 1, iband ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         i, 1, comm, requests( ir ), ierr )
        iband = iband + nbands_to_send
      enddo
      nbands_to_send = nbands / nproc
      do i = mod( nbands, nproc ) + 1, nproc - 1
        ir = ir + 1
        call MPI_IRSEND( re_wvfn( 1, iband ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         i, 1, comm, requests( ir ), ierr )
        iband = iband + nbands_to_send
      enddo
#endif

      read( 99 ) im_wvfn
#ifdef MPI
      nbands_to_send = nbands / nproc + 1
      iband = 1
      do i = 1, mod( nbands, nproc )
        ir = ir + 1
        call MPI_IRSEND( im_wvfn( 1, iband ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         i, 2, comm, requests( ir ), ierr )
        iband = iband + nbands_to_send
      enddo
      nbands_to_send = nbands / nproc
      do i = mod( nbands, nproc ) + 1, nproc - 1
        ir = ir + 1
        call MPI_IRSEND( im_wvfn( 1, iband ), nbands_to_send*ngvecs, MPI_DOUBLE_PRECISION, &
                         i, 2, comm, requests( ir ), ierr )
        iband = iband + nbands_to_send
      enddo
#endif
      close( 99 )

    else
      nr = 3
      allocate( requests( nr ) )
      allocate( re_wvfn( ngvecs, my_bands ), im_wvfn( ngvecs, my_bands ) )
#ifdef MPI
      call MPI_IRECV( re_wvfn, ngvecs*my_bands, MPI_DOUBLE_PRECISION, root, 1, comm, requests( 1 ), ierr )
      call MPI_IRECV( im_wvfn, ngvecs*my_bands, MPI_DOUBLE_PRECISION, root, 2, comm, requests( 2 ), ierr )

      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. 0 ) then
        call MPI_CANCEL( requests( 1 ) , ierr )
        call MPI_CANCEL( requests( 2 ) , ierr )
        ierr = 5
        return
      endif

      call MPI_IBCAST( gvecs, 3*ngvecs, MPI_INTEGER, root, comm, requests( 3 ), ierr )
#endif


    endif

    call MPI_WAITALL( nr, requests, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    wfns( :, : ) = cmplx( re_wvfn( :, 1:my_bands ), im_wvfn( :, 1:my_bands ), DP )
    deallocate( re_wvfn, im_wvfn )


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
