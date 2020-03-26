#if 1
! Copyright (C) 2018 - 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module ocean_abi_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif

  use ocean_mpi, only : MPI_OFFSET_KIND



  implicit none
  private
  save


  logical :: is_init = .false.
  logical :: is_gamma
  logical :: is_shift  ! Is there actually a difference between val and con?
  logical :: is_split  ! Are val and con in different dft runs?
  logical, parameter :: GammaFullStorage = .false.
  character( len=128 ) :: prefix
  character( len=128 ) :: prefixSplit

  integer :: bands(2)
  integer :: brange(4)
  integer :: kpts(3)
  integer :: nkpt
  integer :: nspin

  ! Bands by k-point should be uniform in situations ocean encouters
  ! is indexed by kpoint and spin
  integer, allocatable :: bandsByK( :, : )
  ! Planewaves are only indexed by k-point (same for spin up and down)
  integer, allocatable :: planewavesByK( : )
  ! offsets, which should be of size MPI_OFFSET_KIND
  ! gives file location (including offset for header ) for the start of each
  !  indexed by k-point and spin
  integer( MPI_OFFSET_KIND ), allocatable :: gVecOffsets( :, : )
  integer( MPI_OFFSET_KIND ), allocatable :: eigenOffsets( :, : )
  integer( MPI_OFFSET_KIND ), allocatable :: wvfOffsets( :, : )


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
  integer :: pool_val_nbands
  integer :: pool_con_nbands

  integer :: WFK_FH, WFK_splitFH
  
  public :: abi_read_init
  public :: abi_getAllBandsForPoolID, abi_getValenceBandsForPoolID, abi_getConductionBandsForPoolID

  contains

  pure function abi_getAllBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    nbands = 0
    bands_remain = bands(2)-bands(1)+1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function abi_getAllBandsForPoolID

  pure function abi_getValenceBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    bands_remain = brange(2) - brange(1) + 1
    nbands = 0

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function abi_getValenceBandsForPoolID

  pure function abi_getConductionBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    nbands = 0
    bands_remain = brange(4) - brange(3) + 1

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function abi_getConductionBandsForPoolID

! At the moment we aren't parsing the header of the second WFK file if the run was split valence/conduction for kshifted
  subroutine abi_read_init( comm, ierr )
    use ai_kinds, only : sizeChar
    use ocean_mpi, only : MPI_INTEGER, MPI_CHARACTER, MPI_LOGICAL, MPI_DOUBLE_PRECISION, &
                          MPI_MODE_RDONLY, MPI_INFO_NULL, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                          MPI_OFFSET_KIND
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr

    character(len=6) :: codvsn
    character(len=132) :: title
    character(len=132) :: filnam
    integer :: ierr_, fh, idum(8)
    integer( MPI_OFFSET_KIND ) :: pos, pos2

    logical :: isSplit
!    real(dp) :: d1(20)

    ! Set the comms for file handling
    call MPI_COMM_DUP( comm, inter_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( inter_comm, inter_nproc, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( inter_comm, inter_myid, ierr )
    if( ierr .ne. 0 ) return 

    call load_ocean_inputs( ierr )
    if( ierr .ne. 0 ) return

    if( inter_myid .eq. inter_root ) then
      
      ! Get file name(s)
      isSplit = .false.
      call get_fileName( filnam, isSplit )
      
      open( unit=99, file=filnam, form='unformatted', status='old' )

      call parseHeader( 99, pos, ierr )

      close( 99 )

!      call MPI_FILE_OPEN( MPI_COMM_WORLD, "filnam", MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
!      pos2 = 4
!      call MPI_FILE_READ_AT( fh, pos2, codvsn, 6, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr )
!      pos2 = pos2 + 6 * sizeChar
!      call MPI_FILE_READ_AT( fh, pos2, idum, 2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
!      write(6,*) codvsn, idum(1:2)

!!      call MPI_FILE_READ_AT( fh, pos, idum, 2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
!!      write(6,*)  idum(1:2)
!!      call MPI_FILE_READ_AT( fh, pos, d1, 2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
!!      write(6,*) d1(:)
!!      write(6,*) 'ETOT: ', d1(1)
!!      write(6,*) 'EF  : ', d1(2)
!!      call MPI_FILE_READ_AT( fh, pos, title, 132, MPI_CHARACTER,  MPI_STATUS_IGNORE, ierr )
!!      write(6,*) title
!      write(6,*) pos
!      call MPI_FILE_CLOSE( fh, ierr )
    endif

    ! Error sync
111 continue
    call MPI_BCAST( ierr, 1, MPI_INTEGER, inter_root, inter_comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    call share_init( ierr )
    if( ierr .ne. 0 ) return

    call set_pools( ierr )
    if( ierr .ne. 0 ) return

    call open_wvfn_files( ierr )
    if ( ierr .ne. 0 ) return

    call writeDiagnostics( )



  end subroutine abi_read_init


  ! shares the important info from the header, like planewaves and offsets
  subroutine share_init( ierr )
    use ocean_mpi, only : comm, root, myid, MPI_INTEGER, MPI_OFFSET
    
    integer, intent( inout ) :: ierr

    if( myid .ne. root ) then
      allocate( gVecOffsets( nkpt, nspin ), eigenOffsets( nkpt, nspin ), wvfOffsets( nkpt, nspin ), &
                bandsByK( nkpt, nspin ), planewavesByK( nkpt ) )
    endif
    call MPI_BCAST( bands, 2, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( bandsByK, nkpt*nspin, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( planewavesByK, nkpt, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( gVecOffsets, nkpt*nspin, MPI_OFFSET, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( eigenOffsets, nkpt*nspin, MPI_OFFSET, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( wvfOffsets, nkpt*nspin, MPI_OFFSET, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

  end subroutine share_init

  ! Determines the position of the start of info (after the header ) in the file
  ! Also returns the important info about the file 
  subroutine parseHeader( iun, pos, ierr )
    use ai_kinds, only : sizeInt, sizeChar, sizeDouble, sizeRecord

    integer, intent( in ) :: iun
    integer(MPI_OFFSET_KIND), intent( out ) :: pos
    integer, intent( inout ) :: ierr


!    integer :: bantot, date, intxc, ixc, natom, ngfft(3), nkpt, nspden, nspinor, nsppol, nsym, & 
!               npsp, ntypat, occopt, pertcase, usepaw, qptn, usewvl, nshiftk_orig, nshiftk, mband
!    real(DP) :: ecut, ecutdg, ecutsm, ecut_eff, rprimd(3,3), stmbias, tphysel, tsmear
    integer, allocatable, dimension( : ) :: i1, i2, i3
    integer, allocatable, dimension( : ) :: wavefunctionStorageByK
!    integer, allocatable, dimension( : ) :: wavefunctionStorageByK, bandsByK, planewavesByK
    real(DP), allocatable, dimension( : ) :: d1, d2

    character(len=6) :: codvsn
    character(len=132) :: title
    integer :: headform, fform, maxband, noncollinear, numSyms, numPsps, &
               numAtomTypes, natoms, numkshifts, i, k, j
    integer( MPI_OFFSET_KIND ) :: offset

    pos = 0

    read( iun ) codvsn,headform,fform
    write(6,*) 'Abinit version: ', codvsn, headform,fform

    if( headform < 80 ) then
      write(6,*) 'Unsupported ABINIT version (or file reading earlier)'
      ierr = 214907
      return
    endif

    ! 2x sizeRecord per record
    pos = pos + 2 * sizeRecord + 2 * sizeInt + 6 * sizeChar

    allocate( i1( 18 ), d1( 7 ), i2( 1 ), d2( 12 ), i3( 4 ) )
    read( iun ) i1, d1, d2, i3
!bantot, date, intxc, ixc, natom, ngfft(1:3), &
!& nkpt, nspden, nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase,&
!& usepaw, ecut, ecutdg, ecutsm, ecut_eff, qptn, rprimd, &
!& stmbias, tphysel, tsmear, usewvl, nshiftk_orig, nshiftk, mband

    pos = pos + 22 * sizeInt + 19 * sizeDouble + 2 * sizeRecord

    numKshifts = i3(3)
    maxband = i3( 4 )
    if( maxband .lt. bands(2) ) then
      write(6,*) "Abinit file has too few bands"
      ierr = 2147
      return
    endif
 
    if( nkpt .ne. i1( 9 ) ) then
      write(6,*) "K-points in abinit file don't match other inputs"
      ierr = 2148
      return
    endif
!    nkpt = i1( 9 )

    if( nspin .ne. i1( 12 ) ) then
      write(6,*) "Spin in abinit file don't match other inputs"
      ierr = 2149
      return
    endif
!    nspin = i1( 12 )
    noncollinear = i1( 11 ) 
    numSyms = i1( 13 )
    numPsps = i1( 14 )
    numAtomTypes = i1( 15 )
    natoms = i1( 5 )

    if( noncollinear .ne. 1 ) then
      write( 6, * ) 'Non-collinear spins not supported!!'
      ierr = 125908
      return
    endif
    write(6,*) numPsps
    write(6,*) natoms
    write(6,*) numSyms
    write(6,*) numAtomTypes
    write(6,*) nspin
    write(6,*) maxband


    allocate( wavefunctionStorageByK( nkpt ), bandsByK( nkpt, nspin ), planewavesByK( nkpt ) )
   
    read( iun ) wavefunctionStorageByK(:), bandsByK(:,:), planewavesByK(:)!, i1(1:2)
!    write(6,*) i1(1:2)
!    write(6,*) bandsByK(:)
    ! not read
    ! so_psp -- integers, npsp
    ! symafm -- integer, nsym
    ! symrel -- integer, 9 * nsym
    ! typat -- integer, natom
    ! kptns -- real, 3 * nkpt
    ! occ3d -- integer maxBand * nkpt * nspin
    ! tnons --3 * nsym
    ! znucltypat -- ntypat
    ! wtk -- nkpt

    pos = pos + sizeRecord + 2 * nkpt * sizeInt + nkpt * nspin * sizeInt
    
    pos = pos + numPsps * sizeInt + 10 * numSyms * sizeInt + natoms * sizeInt + 3 * nkpt * sizeDouble &
        + maxband * nkpt * nspin * sizeDouble + 3 * numSyms * sizeDouble &
        + numAtomTypes * sizeDouble + nkpt * sizeDouble

    pos = pos + sizeRecord

    deallocate( d1, d2 )
    allocate( d1( 1 + 3 * natoms ), d2( 2 ) )
    read( iun ) d1, d2
    write(6,*) 'ETOT: ', d2(1)
    write(6,*) 'EF  : ', d2(2)
  
    pos = pos + 2 * sizeRecord + 3 * natoms * sizeDouble + 3 * sizeDouble + numAtomTypes * sizeDouble
!    pos = pos + sizeRecord + sizeDouble + 3 * natoms * sizeDouble

    ! kptopt = int
    ! pawcpxocc = int
    ! nelect = real
    ! charge = real
    ! icoulomb = int
    ! kptrlatt_orig = 9 int
    ! kptrlatt = 9 int
    ! shiftk_orig = 3 * nshiftk real
    ! shiftk = 3 * nshiftk real
    pos = pos + 2 * sizeRecord + 21 * sizeInt + 2 * sizeDouble + 6 * numkshifts * sizeDouble

    read( iun )


!    &   hdr%title(ipsp), hdr%znuclpsp(ipsp), hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
!&   hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp), hdr%md5_pseudos(ipsp)
    ! title 132 char
    ! znuclpsp real
    ! zionpsp real
    ! pspso int
    ! pspdat int
    ! psp cod int
    ! pspxc int
    ! lmn_size int
    ! md5 char 32
    ! 132 from title and 32 from MD5
    do i = 1, numPsps
      pos = pos + 2 * sizeRecord + 164 * sizeChar + 5 * sizeInt + 2 * sizeDouble
      read( iun ) title
      write(6,*) title
    enddo

    allocate( gVecOffsets( nkpt, nspin ), eigenOffsets( nkpt, nspin ), wvfOffsets( nkpt, nspin ) )

    offset = pos
    do i = 1, nspin
      do k = 1, nkpt
        offset = offset + 2 * sizeRecord + 3 * sizeInt
        gVecOffsets( k, i ) = offset
        offset = offset + 2 * sizeRecord + 3 * planewavesByK( k ) * sizeInt
        eigenOffsets( k, i ) = offset
        offset = offset + 2 * sizeRecord + 2 * bandsByK( k, i ) * sizeDouble
        wvfOffsets( k, i ) = offset
        do j = 1, bandsByK( k, i )
          offset = offset + 2 * sizeRecord + 2 * planewavesByK( k ) * sizeDouble * noncollinear
        enddo
!        m_wfk.F90 :3435
      enddo
    enddo

  end subroutine parseHeader

  subroutine load_ocean_inputs( ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
    integer, intent( inout ) :: ierr

    real(dp) :: qinb(3)
    logical :: ex
    real(DP), parameter :: tol = 0.0000001_DP

    if( inter_myid .eq. inter_root ) then
!      open( unit=99, file='prefix', form='formatted', status='old' )
!      read(99,*)  tmp
!      close( 99 )
!      write( prefix, '(a,a,a)' ) 'Out/', trim(tmp), '.save/'
!      write( prefixSplit, '(a,a,a)' ) 'Out/', trim(tmp), '_shift.save/'
      prefix = 'RUN0001'
      prefixSplit = 'RUN0002'

      open(unit=99,file='brange.ipt', form='formatted', status='old' )
      read(99,*) brange(:)
      close(99)

      inquire( file='bands.ipt', exist=ex )
      if( ex ) then
        open(unit=99,file='bands.ipt', form='formatted', status='old' )
        read(99,*) bands(:)
        close(99)
      else
        bands(2) = brange(4)
        bands(1) = brange(1)
      endif
      open(unit=99,file='kmesh.ipt', form='formatted', status='old' )
      read(99,*) kpts(:)
      close(99)
      open(unit=99,file='nspin', form='formatted', status='old' )
      read(99,*) nspin
      close(99)

      inquire( file='qinunitsofbvectors.ipt', exist=ex )
      if( ex ) then
        open( unit=99, file='qinunitsofbvectors.ipt', form='formatted', status='old')
        read( 99 , * ) qinb(:)
        close( 99 )
        if( abs( qinb(1) ) + abs( qinb(2) ) + abs( qinb(3) ) > tol ) is_shift = .true.
      else
        is_shift = .false.
      endif

  
      is_split = .false.
      inquire( file='dft.split', exist=ex )
      if( ex ) then
        open( unit=99, file='dft.split', form='formatted', status='old')
        read( 99, * ) is_split
        close( 99 )
      endif
      if( is_split ) is_shift = .true.
    endif

    call MPI_BCAST( bands, 2, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( brange, 4, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( kpts, 3, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( nspin, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( prefix, len(prefix), MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( prefixSplit, len(prefixSplit), MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( is_shift, 1, MPI_LOGICAL, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( is_split, 1, MPI_LOGICAL, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return


    nkpt = product( kpts(:) )

  end subroutine load_ocean_inputs

  subroutine set_pools( ierr )
    integer, intent( inout ) :: ierr
    !
    integer :: i, nks

    nks = nkpt * nspin

    if( nks .ge. inter_nproc ) then
      mypool = inter_myid
      npool = inter_nproc

    else
      do i = 2, inter_nproc
        if( mod( inter_nproc, i ) .eq. 0 ) then
          if( inter_myid .eq. 0 ) write(6,*) i, inter_nproc
          write(1000+inter_myid,*)  i, inter_nproc, nks
          if( nks .ge. (inter_nproc/i) ) then
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


    pool_nbands = abi_getAllBandsForPoolID( pool_myid )
    pool_con_nbands = abi_getConductionBandsForPoolID( pool_myid )
    pool_val_nbands = abi_getValenceBandsForPoolID( pool_myid )

  end subroutine set_pools

  ! The assumption is that the DFT run has at most two WFK files
  !
  ! Each pool opens its own file handle
  !  this is because within a pool we will want to share G-vectors
  subroutine open_wvfn_files( ierr )
    use ocean_mpi, only : MPI_MODE_RDONLY, MPI_INFO_NULL
    integer, intent( inout ) :: ierr
    !
    character( len=128 ) :: filnam

    call get_fileName( filnam, .false. )

    ! should be able to set some info from what was gathered in the header parsing
    call MPI_FILE_OPEN( pool_comm, filnam, MPI_MODE_RDONLY, MPI_INFO_NULL, WFK_FH, ierr )
    if( ierr .ne. 0 ) return

    if( is_split ) then
      call get_fileName( filnam, .true. )
      ! should be able to set some info from what was gathered in the header parsing
      call MPI_FILE_OPEN( pool_comm, filnam, MPI_MODE_RDONLY, MPI_INFO_NULL, WFK_splitFH, ierr )
      if( ierr .ne. 0 ) return
    endif

  end subroutine open_wvfn_files
    

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
    write(1000+inter_myid, '(A,I8)' ) " My con bands:", pool_con_nbands
    write(1000+inter_myid, '(A,I8)' ) " My val bands:", pool_val_nbands


    write(1000+inter_myid, '(A)' ) "    #############################"
!    flush(1000+inter_myid)
  end subroutine

  subroutine get_fileName( filnam, isSplit )
    character( len=*), intent( out ) :: filnam
    logical, intent( in ) :: issplit

    if( isSplit ) then
      write( filnam, '(A,A)' ) trim( prefixSplit ), '_WFK'
    else
      write( filnam, '(A,A)' ) trim( prefix ), '_WFK'
    endif
  end subroutine get_fileName
    

end module ocean_abi_files
#endif
