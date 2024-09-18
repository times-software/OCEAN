! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2019
!
!
module prep_wvfn
  use AI_kinds, only : DP

  implicit none
  private

  logical, parameter :: perKptU2 = .true.

  public :: prep_wvfn_driver


  logical, parameter :: debug = .false.

  contains


  subroutine prep_wvfn_driver( ierr )
    use prep_system, only : system_parameters, params, psys, prep_system_ikpt2kvec, calcParams, &
                            calculation_parameters, prep_system_umklapp
    use ocean_mpi, only : myid, root, nproc, comm
    use ocean_dft_files, only : odf_is_my_kpt, odf_return_my_bands, odf_nprocPerPool, odf_poolID, &
                                ODF_VALENCE, ODF_CONDUCTION, odf_universal2KptandSpin, &
                                odf_get_ngvecs_at_kpt, odf_read_at_kpt_split, odf_npool, &
                                odf_getBandsForPoolID
    use ocean_cks, only : ocean_cks_build, ocean_cks_nsites, ocean_cks_writeCksHolders, & 
                          ocean_cks_makeCksHolders, ocean_cks_doLegacyCks
    use screen_opf, only : screen_opf_largestLMNproj
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    use ocean_tmels, only : ocean_tmels_open, ocean_tmels_close, ocean_tmels_calc, ocean_tmels_legacy
    integer, intent( inout ) :: ierr

    complex(DP), allocatable, target :: valUofG(:,:), conUofG(:,:), wvfn(:,:,:,:), UofX(:,:,:,:) , UofX2(:,:)
    complex(DP), allocatable :: spareValUofG( :, : )
          !valUofX(:,:,:,:), conUofX(:,:,:,:)
    complex(DP), pointer :: UofGPointer(:,:)!, cksPointer(:,:,:) !, UofXpointer(:,:,:,:)
!    complex(DP), allocatable, target :: cksValArray(:,:,:,:), cksConArray(:,:,:,:)

    real(DP) :: kqVec(3), deltaR, kqVecCart(3), memEstimate

    integer, allocatable, target :: valGvecs(:,:), conGvecs(:,:)
    integer, pointer :: gvecPointer(:,:)

    integer :: ispin, ikpt, nprocPool, nuni, iuni, npool, nsites, nprojSpacer
    integer :: valNgvecs, conNgves, valBands, conBands, ngvecs(2), odf_flag, totConBands
    integer :: nG, nbands, fftGrid(3), allBands, nX, vType, cType, totValBands, myConBandStart, myValBandStart
    integer :: umklapp(3)
    integer :: nband_chunk, nchunk, ichunk, nband_todo, ib, ib2

    integer :: conFH, valFH, fileHandle, poolID, i, testFH, iband, bandChunk, omp_threads, kStride
    logical :: is_kpt, wantCKS, wantU2, addShift, wantLegacy, wantTmels
!$  integer, external :: omp_get_max_threads

    wantCKS = calcParams%makeCKS
    wantU2 = calcParams%makeU2
    wantLegacy = .false.
    wantTmels = calcParams%makeTmels


!    if( wantU2 .and. nproc .eq. 1 ) then
!      call prep_wvfn_openLegacy( testFH, ierr )
!      if( ierr .ne. 0 ) return
!    endif

    call odf_return_my_bands( valBands, ierr, ODF_VALENCE )
    if( ierr .ne. 0 ) return
  
    call odf_return_my_bands( conBands, ierr, ODF_CONDUCTION )
    if( ierr .ne. 0 ) return

    write(6,*) valBands, conBands

    if( wantCKS ) then
      nsites = ocean_cks_nsites()
    else
      nsites = 0
    endif

!    isDualFile =  odf_isDualFile()
!!!prefix, nxpts, myBands, totalBands, FH, arrayType, ierr )
!!    nX = product( param%xmesh(:) )
!!    allBands = params%brange(2) - params%brange(1) + 1


    npool = odf_npool()
    nprocPool = odf_nprocPerPool()
    poolID = odf_poolID()

    ! The more correct thing to do would be to make a fancier derived type, such that each 
    ! pool would also have an offset based on k-point and then the offset when writing would 
    ! be only dependent on iuni not ikpt
    allBands = params%nkpts * params%nspin * ( params%brange(2) - params%brange(1) + 1 )
    nx = product( params%xmesh(:) )
    call prep_wvfn_openU2( 'val', valFH, allBands, nx, poolID, nprocPool, ierr, vtype )
    if( ierr .ne. 0 ) return
    allBands = params%nkpts * params%nspin * ( params%brange(4) - params%brange(3) + 1 )
    call prep_wvfn_openU2( 'con', conFH, allBands, nx, poolID, nprocPool, ierr, ctype )
    if( ierr .ne. 0 ) return

    ! need to set the number of x-points proc should epxect locally
    ! probably subroutines/functions in this module
    ! grabbing pool_size and pool_id from ODF
    nX = prep_wvfn_divideXmesh( product( params%xmesh ), nprocPool, poolID )

    nuni = ceiling( real( params%nspin * params%nkpts, DP ) / real( npool, DP ) )

    if( wantCKS ) then
      call ocean_cks_makeCksHolders( valBands, conBands, nuni, ierr )
    endif

  
    kStride = odf_npool()
    call odf_universal2KptandSpin( 1, ikpt, ispin )
    if( ispin .gt. 1 ) ikpt = ikpt + params%nkpts

    myConBandStart = 1
    do i = 0, odf_poolID() - 1
      myConBandStart = myConBandStart + odf_getBandsForPoolID( i, ODF_CONDUCTION )
    enddo
    totValBands = params%brange(2) - params%brange(1) + 1
    totConBAnds = params%brange(4) - params%brange(3) + 1

    if( wantTmels ) then
      write(myid+1000,*) totValBands, totConBands, conBands, myConBandStart, nuni, ikpt, kStride
!      flush(myid+1000)
      call ocean_tmels_open( totValBands, totConBands, conBands, myConBandStart, nuni, ikpt, kStride, ierr )
    endif



!    if( wantCKS ) then
!      nprojSpacer = screen_opf_largestLMNproj()
!      allocate( cksValArray( nprojSpacer, valBands, nsites, nuni ), &
!                cksConArray( nprojSpacer, conBands, nsites, nuni ) )
!    else
!      nprojSpacer = 0
!      allocate( cksValArray(0,0,0,0), cksConArray(0,0,0,0) )
!    endif

!    ! loop over spin and kpt
!    do ispin = 1, params%nspin
!      do ikpt = 1, params%nkpts
    do iuni = 1, nuni

      call odf_universal2KptandSpin( iuni, ikpt, ispin )
      write(1000+myid,*) iuni, ikpt, ispin

      ! universal2KptandSpin returns 0,0 if on the last trip and un-equal distribution of k-points
      ! We still call all the write_all calls for things like tmels, u2.dat, and cks
      !

      if( ikpt .gt. 0 ) then

        if( myid .eq. root ) write(6,*) ikpt, ispin

        ! This is just a precaution while implementing
        call odf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
        if( ierr .ne. 0 ) return

        if( is_kpt .eqv. .false. ) then
          ierr = 101
          return
        endif
        ! end precaution
        
        ! first call gets us either valence or both
        call odf_get_ngvecs_at_kpt( ikpt, ispin, ngvecs, ierr )
        if( ierr .ne. 0 ) then
          write(1000+myid,*) 'Failed to get ngvecs', ierr
          return
        else
          write(1000+myid,*) 'Fetched ngvecs: ', ngvecs(:)
        endif
        
        ! separate logic for shifted and for split

!        ! For things like legacy *or* when doing plain XAS/XES there will only be a single list of G-vecs
!        ! Don't copy/read them twice and instead track in the logic in what follows
!        if( isDualFile ) then
        ! alternatively, just have ODF copy over the g-vecs where needed
          allocate( valGvecs( 3, ngvecs(1) ), conGvecs( 3, ngvecs(2) ), &
                    valUofG( ngvecs(1), valBands ), conUofG( ngvecs(2), conBands ) )
!        else    
!          allocate( valGvecs( 3, ngvecs(1) ), conGvecs( 0,0 ), & 
!                    valUofG( ngvecs(1), valBands ), conUofG( ngvecs(2), conBands ) )
!        endif


    ! read in file, which means distribute bands across procs
      ! need to ensure bounds are appropriately enforced 
      ! this means for shifted calcs keep only occ/unocc

        call odf_read_at_kpt_split( ikpt, ispin, ngvecs(1), ngvecs(2), valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
        if( ierr .ne. 0 ) then
          write(1000+myid,*) 'Failed to read k-point', ierr
          return
        endif


        ! Check for Umklapp -- if the actual k-point was outside [0,1] (or maybe [-1,1])
        ! the calculated k-point was *not* made at the requested k-point, but was shifted 
        ! to be inside the BZ
        call prep_system_umklapp( ikpt, .false., umklapp )
        if( umklapp(1) .ne. 0 .or. umklapp(2) .ne. 0 .or. umklapp(3) .ne. 0 ) then
          do i = 1, ngvecs(1)
          valGvecs(:,i) = valGvecs(:,i) - umklapp(:)
          enddo
        endif
        call prep_system_umklapp( ikpt, .true., umklapp )
        if( umklapp(1) .ne. 0 .or. umklapp(2) .ne. 0 .or. umklapp(3) .ne. 0 ) then
          do i = 1, ngvecs(2)
          conGvecs(:,i) = conGvecs(:,i) - umklapp(:)
          enddo
        endif

    ! If doing TMELS for valence, start by figuring that out
      ! A) possibly at the same time as read in redistribute slices of G-points?
      ! but then also need to ensure uniformity of G
      ! B) Everybody gets the occupied states, possibly band by band, and then calculates
      ! the tmels for only their subset of unoccupied states  <----

        if( wantTmels ) then
          totValBands = params%brange(2) - params%brange(1) + 1
!          totConBands = params%brange(4) - params%brange(3) + 1
          ! If we have more than one proc per k-point, share the valence bands
          if( odf_nprocPerPool() .gt. 1 ) then
            myValBandStart = 1
            allocate( spareValUofG( ngvecs(1), totValBands ) )
            call shareValence( valUofG, spareValUofG, ierr )
            if( ierr .ne. 0 ) return
            call ocean_tmels_calc( ikpt, ispin, totValBands, conBands, valGvecs, conGvecs, &
                                   spareValUofG, conUofG, params%nkpts, totValBands, totConBands, &
                                   myValBandStart, myConBandStart, ierr )
            deallocate( spareValUofG )
          else
            call ocean_tmels_calc( ikpt, ispin, totValBands, conBands, valGvecs, conGvecs, &
                                   ValUofG, conUofG, params%nkpts, totValBands, totConBands, &
                                   myValBandStart, myConBandStart, ierr )
          endif
        endif

      else  ! ikpt == 0, the process isn't doing anything 
        ! This makes sure we don't assign pointers to un-allocated arrays
        allocate( valGvecs( 0, 0 ), conGvecs( 0, 0 ), valUofG( 0, 0 ), conUofG( 0, 0 ) )
        if( wantTmels ) then
          call ocean_tmels_calc( ikpt, ispin, totValBands, conBands, valGvecs, conGvecs, &
                                 spareValUofG, conUofG, params%nkpts, totValBands, totConBands, &
                                 myValBandStart, myConBandStart, ierr )
        endif
      endif

      ! Loop over conversion

      gvecPointer => valGvecs
      UofGPointer => valUofG
!        UofXPointer => valUofX
      nbands = valBands
      nG = ngvecs(1)
      fileHandle = valFH
      odf_flag = ODF_VALENCE
!      cksPointer => cksValArray(:,:,:,iuni)

      !JTV
      ! need to set the complete band range for valence/conduction
      allBands = params%brange(2) - params%brange(1) + 1

!      fileHeader = 'cksv'
        
      do i = 1, 2

#if 0
        if( ikpt .ne. 0 ) then
          call prep_wvfn_checkFFT( nG, gvecPointer, wantCKS, wantU2, fftGrid, ierr )

          allocate( wvfn( fftGrid(1), fftGrid(2), fftGrid(3), nbands ) )

    ! FFT by band, openMP option?

          call prep_wvfn_doFFT( gvecPointer, UofGPointer, wvfn )

        else
          allocate( wvfn( 0, 0, 0, 0 ) )
        endif
#endif

!        ! deallocate as-read wvfn
!        if( i .eq. 1 ) then
!          deallocate( valGvecs, valUofG )
!        else
!          deallocate( conGvecs, conUofG )
!        endif


        if( wantU2 ) then
          
          if( ikpt .ne. 0 ) then

            call prep_wvfn_checkFFT( nG, gvecPointer, .false., wantU2, fftGrid, ierr )


            nband_chunk = 1
! !$          nband_chunk = omp_get_max_threads

            allocate( wvfn( fftGrid(1), fftGrid(2), fftGrid(3), nband_chunk ) )
            !JTV
!            ! This reversal is currently in the code, but we need to remove at some point
            ! This is reversed in the legacy u2.dat. That reversal is removed here, but
            ! taken into account in the legacy writer
            allocate( UofX( params%xmesh(1), params%xmesh(2), params%xmesh(3), nbands ) )

            nchunk = ( nbands - 1 ) / nband_chunk + 1
            do ichunk = 1, nchunk
              nband_todo = min( nband_chunk, nbands - ( ichunk - 1 )*nband_chunk )
              ib = (ichunk-1)*nband_chunk+1
              call prep_wvfn_doFFT( gvecPointer, UofGPointer(:,ib:), wvfn(:,:,:,1:nband_todo) )


!              ib = 1 + (ichunk-1)*nband_chunk
              ib2 = ib + nband_todo - 1
              call prep_wvfn_u1( wvfn(:,:,:,1:nband_todo), UofX(:,:,:,ib:ib2), ierr )
              if( ierr .ne. 0 ) return
            enddo

            deallocate( wvfn )
            allocate( UofX2( nX, allBands ) )

            ! begin orthogonalization
              ! because we should be nearly orthogonal already, can use non-modified Gram-Schmidt
            call prep_wvfn_u2( UofX, UofX2, odf_flag, ierr )
            if( ierr .ne. 0 ) return
          
          else
            allocate( UofX( 0, 0, 0, 0 ), UofX2( 0, 0 ) )
          endif

          ! save subsampled xmesh
          if( perKptU2 ) then
            call prep_wvfn_writeU2_perKpt( ikpt, ispin, UofX2, fileHandle, ierr )
          else
            call prep_wvfn_writeU2( ikpt, ispin, UofX2, fileHandle, ierr )
          endif
          if( ierr .ne. 0 ) return

!          if( wantU2 .and. nproc .eq. 1 ) then
!            call prep_wvfn_writeLegacy( testFH, ikpt, UofX2, ierr )
!            if( ierr .ne. 0 ) return
!          endif

          deallocate( UofX, UofX2 )

        endif

        ! if doing CKS then compute matrix elements here
        if( wantCKS ) then

          if( ikpt .ne. 0 ) then

            call prep_wvfn_checkFFT( nG, gvecPointer, wantCKS, .false., fftGrid, ierr )

            memEstimate = real( fftGrid(1), DP ) * real( fftGrid(2), DP ) * real( fftGrid(3), DP ) 
            bandChunk = max( floor( 67108864_DP / memEstimate ), 1 )
            omp_threads = 1
!$          omp_threads = omp_get_max_threads()
            
            bandChunk = max( bandChunk, omp_threads )
            bandChunk = max( 1, bandChunk/omp_threads ) * omp_threads 

            bandChunk = min( bandChunk, nbands )
            bandChunk = max( bandChunk, 1 )

            write(1000+myid, '(3I8,E24.16)' ) omp_threads, bandChunk, nbands, memEstimate/67108864_DP
            write(1000+myid,*) nbands, bandchunk

!            allocate( wvfn( fftGrid(1), fftGrid(2), fftGrid(3), nbands ) )
            allocate( wvfn( fftGrid(1), fftGrid(2), fftGrid(3), bandChunk ) )

            do iband = 1, nbands, bandChunk

              if( nbands - iband + 1 .lt. bandChunk ) then
                deallocate( wvfn )
                allocate( wvfn( fftGrid(1), fftGrid(2), fftGrid(3), nbands - iband + 1 ) )
              endif

              call prep_wvfn_doFFT( gvecPointer, UofGPointer(:,iband:), wvfn )

          ! allocate CKS holder which requires queries for sites, maybe loop over all sites? need all the FHs

          ! call CKS
!            call prep_wvfn_cks( 
          ! deltaR is the min real-space spacing for the overlap, should make it an input
              deltaR = 0.02_DP
              addShift = (i .eq. 2 )
              call prep_system_ikpt2kvec( ikpt, addShift, kqVec, kqVecCart ) 
              if( iband .eq. 1 ) then
                write(1000+myid, '(A,I8,6(X,E24.16),X,L2,I2)' ) 'kqVec', ikpt, kqVecCart(:), & 
                                                                kqVec(:), addShift, i
              endif
              call ocean_cks_build( wvfn, kqVecCart, deltaR, psys%avecs, (i.eq.1), iuni, iband-1, ierr )
            enddo

            deallocate( wvfn )

          ! call write CKS

          endif
        endif

        ! For realu2 option for cks conversion
        if( i .eq. 1 ) then
          deallocate( valGvecs, valUofG )
        else
          deallocate( conGvecs, conUofG )
        endif

!        ! deallocate FFT wvfn
!        deallocate( wvfn )

        gvecPointer => conGvecs
        UofGPointer => conUofG
!          UofXPointer => conUofX
        nbands = conBands
        nG = ngvecs(2)
        fileHandle = conFH
        allBands = params%brange(4) - params%brange(3) + 1
        odf_flag = ODF_CONDUCTION
!        cksPointer => cksConArray(:,:,:,iuni)
      
      enddo ! i

      

    enddo ! iuni ( combined spin and kpt )


    call prep_wvfn_closeU2( valFH, ierr, vtype )
    call prep_wvfn_closeU2( conFH, ierr, ctype )

    call MPI_BARRIER( comm, ierr )

    call screen_tk_stop( "main" )
    call screen_tk_start( "legacy" )
    
    if( wantU2 ) then
!!      if( nproc .eq. 1 ) then
!!        call prep_wvfn_closeLegacy( testFH, ierr )
!!        if( ierr .ne. 0 ) return
!!        call prep_wvfn_doLegacyParallel( ierr )
!!      else
        if( wantLegacy ) then
          write(1000+myid, * ) 'Calling Legacy Parallel'
          flush(1000+myid)
          call prep_wvfn_doLegacyParallel( ierr )
          call MPI_BARRIER( comm, ierr )
          write(1000+myid, * ) 'Finished Legacy Parallel'
          flush(1000+myid)
        endif
!!      endif
    endif

    if( wantCKS ) then
  !    if( nproc .eq. 1 ) then
  !      call prep_wvfn_cksWriteLegacy( cksValArray, cksConArray )
  !    endif
      call ocean_cks_writeCksHolders( ierr )
      if( ierr .ne. 0 ) return
      call MPI_BARRIER( comm, ierr )
      if( ierr .ne. 0 ) return
      if( wantLegacy ) then
        write(1000+myid, * ) 'Calling Legacy CKS'
        flush(1000+myid)
        call ocean_cks_doLegacyCks( ierr )
        call MPI_BARRIER( comm, ierr )
        write(1000+myid, * ) 'Finished Legacy CKS'
        flush(1000+myid)
      endif
    endif

    if( wantTmels ) then 
      call ocean_tmels_close( ierr )
      if( wantLegacy ) call ocean_tmels_legacy( totValBands, totConBands, params%nkpts, params%nspin, ierr )
    endif

    call screen_tk_stop( "legacy" )
    call screen_tk_start( "main" )

    write(1000+myid, * ) 'Calling Energies'
!    flush(1000+myid)
    call prep_wvfn_writeEnergies( ierr )
    if( ierr .ne. 0 ) return
    call MPI_BARRIER( comm, ierr )
    write(1000+myid, * ) 'Finished Energies'
!    flush(1000+myid)

  end subroutine prep_wvfn_driver

  subroutine prep_wvfn_writeEnergies( ierr )
    use ocean_mpi, only : myid, root
    use ocean_dft_files, only : odf_read_energies_split
    use prep_system, only : params
    use ocean_constants, only : eV2Hartree, Hartree2eV
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: conEnergies(:,:,:), valEnergies(:,:,:) 

    real(DP) :: eshift, efermi
    integer :: nbc, nbv, nk, ioerr, ib, ik, ispin
    logical :: zero_lumo, have_core_offset, ex, noshiftlumo

    nbc = params%brange(4)-params%brange(3)+1
    nbv = params%brange(2)-params%brange(1)+1
    nk = product( params%kmesh(:) )
    allocate( conEnergies( nbc, nk, params%nspin ), valEnergies( nbv, nk, params%nspin ) )

    ! don't bother sharing
    call odf_read_energies_split( myid, root, valEnergies, conEnergies, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) then

      open(unit=99, file='enkfile', form='formatted', status='unknown' )
      do ispin = 1, params%nspin 
        do ik = 1, nk
          write(99, * ) valEnergies(:,ik, ispin)*2.0_DP
          write(99, * ) conEnergies(:,ik, ispin)*2.0_DP
        enddo
      enddo
      close(99)

      open( unit=99, file='efermiinrydberg.ipt', form='formatted', status='old')
      read( 99, * ) efermi
      close( 99 )
      efermi = efermi * 0.5_DP
      eshift = conEnergies( nbc, 1, 1 )
      do ispin = 1, params%nspin 
        do ik = 1, nk
          do ib = 1, nbc
            if( conEnergies( ib, ik, ispin ) .ge. efermi ) then
              if( eshift .gt. conEnergies( ib, ik, ispin ) ) then
                eshift = conEnergies( ib, ik, ispin )
              endif  
              exit
            endif
          enddo
        enddo
      enddo
!      eshift = -eshift

      open( unit=99, file='cbm.out', form='formatted', status='unknown' )
      write(99,*) eshift*Hartree2eV
      close( 99 )

      zero_lumo = .false.

!      open( unit=99, file='core_offset', form='formatted', status='old')
!      rewind 99
!      ! might be true/false, but might be a number (which means true)
!      read( 99, *, IOSTAT=ioerr ) have_core_offset
!      close( 99 )
!      if( ioerr .ne. 0 ) have_core_offset = .true.
!      ! if we are doing a core-level shift, then we don't zero the lumo
!      if( have_core_offset ) zero_lumo = .false.

      ! Allow an override option
      inquire( file='noshift_lumo',exist=ex)
      if( ex ) then
        open(unit=99,file='noshift_lumo',form='formatted',status='old')
        read( 99, * ) noshiftlumo
        close( 99 )
        if( noshiftlumo ) then
          zero_lumo = .false.
        else
          zero_lumo = .true.
        endif
        write(6,*) 'LUMO shift from no_lumoshiftt:', zero_lumo
      endif

      if( zero_lumo ) then
!        open( unit=99, file='eshift.ipt', form='formatted', status='old' )
!        rewind 99
!        read ( 99, * ) eshift
!        close( 99 )
!        eshift = eshift * eV2Hartree

        valEnergies(:,:,:) = valEnergies(:,:,:) + eshift
        conEnergies(:,:,:) = conEnergies(:,:,:) + eshift
      endif

      open( unit=99, file='eshift.ipt', form='formatted', status='unknown' )
      write(99, * ) eshift * Hartree2eV
      close( 99 )

      open(unit=99, file='wvfvainfo', form='unformatted', status='unknown' )
      rewind( 99 )
      write( 99 ) nbv, nk, params%nspin
      write( 99 ) valEnergies
      close( 99 ) 

      open(unit=99, file='wvfcninfo', form='unformatted', status='unknown' )
      rewind( 99 )
      write( 99 ) nbc, nk, params%nspin
      write( 99 ) conEnergies
      close( 99 )

#if 0
      open(unit=99, file='wvfcninfo2', form='formatted', status='unknown' )
      rewind( 99 )
      write( 99, * ) nbc, nk, params%nspin
      write( 99, * ) conEnergies
      close( 99 )
#endif
    endif

    deallocate( valEnergies, conEnergies )
  end subroutine prep_wvfn_writeEnergies

  subroutine prep_wvfn_openLegacy( fh, ierr )
    integer, intent( out ) :: FH 
    integer, intent( inout ) :: ierr

    open( file='u2.dat', form='unformatted', status='unknown', newunit=fh )
  end subroutine

  subroutine prep_wvfn_writeLegacy( FH, ik, UofX2, ierr )
    integer, intent( in ) :: fh, ik
    complex(DP), intent( in ) :: UofX2(:,:)
    integer, intent( inout ) :: ierr

    integer :: nx, nb, ix, ib
    nx = size( UofX2, 1 )
    nb = size( UofX2, 2 )

    do ib = 1, nb
      do ix = 1, nx
        write( fh ) ib, ik, ix, UofX2( ix, ib )
      enddo
    enddo
  end subroutine prep_wvfn_writeLegacy

  subroutine prep_wvfn_closeLegacy( fh, ierr )
    integer, intent( in ) :: FH  
    integer, intent( inout ) :: ierr

    close( fh )
  end subroutine prep_wvfn_closeLegacy



!  subroutine prep_wvfn_openU2( prefix, nxpts, myBands, totalBands, FH, arrayType, ierr )
  subroutine prep_wvfn_openU2( prefix, FH, nb, nx, poolID, nprocPerPool, ierr, fileType )
    use ocean_mpi, only : myid, comm, &
                          MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, & 
                          MPI_OFFSET_KIND, MPI_DOUBLE_COMPLEX
    
    character(len=*), intent( in ) :: prefix
!    integer, intent( in ) :: nxpts, myBands, totalBands
    integer, intent( out ) :: FH  !, arrayType
    integer, intent( in ) :: nb, nx, poolID, nprocPerPool
    integer, intent( inout ) :: ierr

    character(len=128) :: filnam
    integer( MPI_OFFSET_KIND ) :: offset
    integer :: fflags, sizeofcomplex, i, myx
    complex(DP) :: dumz
#ifdef MPI_F08
    type( MPI_DATATYPE ):: fileType
#else
    integer :: fileType
#endif

    write( filnam, '(A,A)' ) trim( prefix ), '.u2.dat'

#ifdef MPI
    fflags = IOR( MPI_MODE_WRONLY, MPI_MODE_CREATE )
    fflags = IOR( fflags, MPI_MODE_UNIQUE_OPEN )
    call MPI_FILE_OPEN( comm, filnam, fflags, MPI_INFO_NULL, fh, ierr )
    if( ierr .ne. 0 ) return

    

    myx = prep_wvfn_divideXmesh( nx, nprocPerPool, poolID )
    call MPI_TYPE_VECTOR( nb, myx, nx, MPI_DOUBLE_COMPLEX, fileType, ierr )
    if( ierr .ne. 0 ) return
    call MPI_TYPE_COMMIT( fileType, ierr )
    if( ierr .ne. 0 ) return
    write(1000+myid,*) 'FILE TYPE'
    write(1000+myid,*) nb, myx, nx
    
    offset = 0
    do i = 0, poolID - 1
      offset = offset + prep_wvfn_divideXmesh( nx, nprocPerPool, i )
    enddo
#ifdef SIZEOF
    call MPI_SIZEOF( dumz, sizeofcomplex, ierr )
#else
    sizeofcomplex = 16
#endif
    if( ierr .ne. 0 ) return
    offset = offset *  sizeofcomplex
    write(1000+myid, * ) 'offset', offset, sizeofcomplex
      
    ! We have two options for writing this out
    if( perKptU2 ) then
      call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, "native", MPI_INFO_NULL, ierr )
    else
      call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_COMPLEX, fileType, "native", MPI_INFO_NULL, ierr )
    endif
    if( ierr .ne. 0 ) return
#else
    open( file=filnam, form='unformatted', status='unknown', newunit=fh )
#endif

  end subroutine prep_wvfn_openU2

  subroutine prep_wvfn_closeU2( FH, ierr, fileType )
    integer, intent( out ) :: FH 
    integer, intent( inout ) :: ierr, fileType

#ifdef MPI
    call MPI_FILE_CLOSE( fh, ierr )
    call MPI_TYPE_FREE( fileType, ierr )
#else
   close( fh )
#endif
  end subroutine prep_wvfn_closeU2
  

  subroutine prep_wvfn_u1( wvfn, UofX, ierr )
    use ocean_mpi, only : myid
    complex(DP), intent( in ) :: wvfn(:,:,:,:)
    complex(DP), intent( out ) :: UofX(:,:,:,:)
    integer, intent( inout ) :: ierr
    !
    integer :: xIn, yIn, zIn, xOut, yOut, Zout, nbands, nbands_
    integer :: xSkip, ySkip, zSkip, ix, iy, iz, iix, iiy, iiz, iband, nOut
    real(DP) :: res
    real(DP), external :: dznrm2

    xIn = size( wvfn, 1 ) 
    yIn = size( wvfn, 2 ) 
    zIn = size( wvfn, 3 ) 
    nbands = size( wvfn, 4 ) 

    xOut = size( UofX, 1 )
    yOut = size( UofX, 2 )
    zOut = size( UofX, 3 )
    nbands_ = size( UofX, 4 )

    nOut = xOut * yOut * zOut

    if( nbands .ne. nbands_ ) then
      ierr = 4
      return
    endif

    if( mod( xIn, xOut ) .ne. 0 ) then
      ierr = 1
      return
    endif
    if( mod( yIn, yOut ) .ne. 0 ) then
      ierr = 2
      return
    endif
    if( mod( zIn, zOut ) .ne. 0 ) then
      ierr = 3
      return
    endif

  
    xSkip = xIn / xOut
    ySkip = yIn / yOut
    zSkip = zIn / zOut

    do iband = 1, nbands
      !JTV
! NO LONGER      ! This reversal is currently in the code, but we need to remove at some point
      iiz = 1
      do iz = 1, zOut

        iiy = 1
        do iy = 1, yOut

          iix = 1
          do ix = 1, xOut
            
            UofX( ix, iy, iz, iband ) = wvfn( iix, iiy, iiz, iband )
            iix = iix + xSkip
          enddo
          
          iiy = iiy + ySkip
        enddo
    
        iiz = iiz + zSkip
      enddo


      ! normalize
      res = dznrm2( nOut, UofX(:,:,:,iband), 1 )
      res = 1.0_DP / res
      call zdscal( nOut, res, UofX(:,:,:,iband), 1 )
      ! end normalize

!      write(1000+myid, * ) ' !', iband, dznrm2( nOut, UofX(:,:,:,iband), 1 )
    enddo

  end subroutine prep_wvfn_u1

  subroutine prep_wvfn_u2( UofX, UofX2, odf_flag, ierr )
    use ocean_mpi, only : MPI_DOUBLE_COMPLEX, MPI_STATUSES_IGNORE, MPI_SUM, MPI_IN_PLACE, myid
    use ocean_dft_files, only : odf_poolComm, odf_nprocPerPool, odf_getBandsForPoolID, odf_poolID

    complex(DP), intent( in ) :: UofX(:,:,:,:)
    complex(DP), intent( out ) :: UofX2(:,:)
    integer, intent( in ) :: odf_flag
    integer, intent(inout ) :: ierr

    complex(DP), allocatable :: coeff( : )
    integer :: iprocPool, nprocPool, mypool, iband, jband, nband, nX, xdim, ydim, zdim, totdim, ix, iy, iz
#ifdef MPI_F08
    type( MPI_COMM ) :: pool_comm
    type( MPI_REQUEST ) :: req(:,:)
    type( MPI_DATATYPE ):: newType
#else
    integer :: pool_comm, newType
    integer, allocatable :: req(:,:)
#endif

    nprocPool = odf_nprocPerPool()
    allocate( req( 0:nprocPool-1, 2) )
    mypool = odf_poolID()

    pool_comm = odf_poolComm()

    nX = size( UofX2, 1 )

    xdim = size( UofX, 1 )
    ydim = size( UofX, 2 )
    zdim = size( UofX, 3 )
    totdim = xdim * ydim * zdim
    
    ! Loop through every pool and post receives
    iband = 1
    do iprocPool = 0, nprocPool - 1
      
      nband = odf_getBandsForPoolID( iprocPool, odf_flag )
      if( debug ) write(1000+myid,*) 'UofX2 recvs:', iprocPool, nx, nband, totdim
      call MPI_IRECV( UofX2(:,iband:), nX*nband, MPI_DOUBLE_COMPLEX, iprocPool, 1, pool_comm, & 
                      req( iprocPool, 1 ), ierr )
      if( ierr .ne. 0 ) return
    
      iband = iband + nband
    enddo

    call MPI_BARRIER( pool_comm, ierr )
    if( ierr .ne. 0 ) return

    iz = 1
    iy = 1
    ix = 1
    nband = odf_getBandsForPoolID( mypool, odf_flag )

    do iprocPool = 0, nprocPool - 1
      nX = prep_wvfn_divideXmesh( totdim, nprocPool, iprocPool )

      call MPI_TYPE_VECTOR( nband, nX, totdim, MPI_DOUBLE_COMPLEX, newType, ierr )
      if( ierr .ne. 0 ) return
      call MPI_TYPE_COMMIT( newType, ierr )
      if( ierr .ne. 0 ) return

      if( debug ) then
        write(1000+myid,*) 'UofX2 sends:', iprocPool, nx, nband
        write(1000+myid,*) '            ', ix, iy, iz
      endif
      call MPI_ISEND( UofX( ix, iy, iz, 1 ), 1, newType, iprocPool, 1, pool_comm, req( iprocPool, 2 ), ierr )
      if( ierr .ne. 0 ) return

      call MPI_TYPE_FREE( newType, ierr )
      if( ierr .ne. 0 ) return

      ix = ix + nX
      do while( ix .gt. xdim )
        ix = ix - xdim
        iy = iy + 1
      enddo
      do while( iy .gt. ydim )
        iy = iy - ydim  
        iz = iz + 1
      enddo
    enddo

!    flush(1000+myid)
    call MPI_WAITALL( nprocPool * 2, req, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    deallocate( req )
    nband = size( UofX2, 2 )
    allocate( coeff( nband ) )

#if 0
!TEST
    do iband = 1, nband
      coeff(iband) = dot_product( UofX2(:,iband), UofX2(:,iband) )
    enddo
    call MPI_ALLREDUCE( MPI_IN_PLACE, coeff, nband, MPI_DOUBLE_COMPLEX, MPI_SUM, pool_comm, ierr )
    do iband = 1, nband
      write(1000+myid, * ) '!!!', real(coeff(iband),DP), size( UofX2,1)
    enddo
!\TEST
#endif

    ! we've already normalized the un-orthogonalized ones, so skip band 1 
    do iband = 2, nband
      do jband = 1, iband - 1
        coeff( jband ) = dot_product( UofX2(:,jband), UofX2(:,iband) )
      enddo
      call MPI_ALLREDUCE( MPI_IN_PLACE, coeff, iband-1, MPI_DOUBLE_COMPLEX, MPI_SUM, pool_comm, ierr )
      if( ierr .ne. 0 ) return 
 !     write(1000+myid,*) coeff(1:iband-1)

      do jband = 1, iband - 1
        UofX2(:,iband) = UofX2(:,iband) - coeff( jband ) * UofX2(:,jband)
      enddo
      coeff( 1 ) = dot_product( UofX2(:,iband), UofX2(:,iband) )

      call MPI_ALLREDUCE( MPI_IN_PLACE, coeff, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, pool_comm, ierr )
      coeff( 1 ) = 1.0_DP / sqrt( coeff(1) )
      UofX2(:,iband) = UofX2(:,iband) * coeff( 1 )

    enddo

    deallocate( coeff )

  end subroutine prep_wvfn_u2

  subroutine prep_wvfn_writeU2_perKpt( ikpt, ispin, UofX2, fileHandle, ierr )
    use prep_system, only : system_parameters, params
    use ocean_dft_files, only : odf_poolID, odf_nprocPerPool, odf_poolComm
    use ocean_mpi, only : &
#ifdef MPI_F08
                          MPI_Datatype, &
#endif
                          MPI_DOUBLE_COMPLEX, MPI_OFFSET_KIND, MPI_STATUS_IGNORE, myid, MPI_STATUS_SIZE, MPI_INFO_NULL

    integer, intent( in ) :: ikpt, ispin, fileHandle
    complex(DP), intent( in ) :: UofX2(:,:)
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: writeBuffer(:,:)

#ifdef MPI_F08
    type( MPI_DATATYPE ):: newType
#else
    integer :: newType, stat( MPI_STATUS_SIZE)
#endif
    integer(MPI_OFFSET_KIND) :: offset
    integer :: ix, i, myPoolID, nprocPerPool, nx, myx, nb, poolComm, id, ib, xstart, xstop, nbchunk, nbwrite, iib



    myPoolID = odf_poolID()
    nprocPerPool = odf_nprocPerPool()
    poolComm = odf_poolComm()

    nb = size( UofX2, 2 )
    myx = size( UofX2, 1)
    nx = product( params%xmesh(:) )

    if( (nx .gt. 32768 ) .or. ( nb .gt. 32768 ) .or. ( nx*nb .gt. 134217728 ) ) then
      nbchunk = max( 1, floor( 134217728.0_DP / real( nx, DP ) ) )
    else 
      nbchunk = max(nb,1) ! in case nb is zero
    endif

    ! offset for k-point
    offset = max( ikpt - 1 + ( ispin - 1 ) * params%nkpts, 0 )
!    offset = offset * int( size( UofX2, 1 ), MPI_OFFSET_KIND ) * int( size( UofX2, 2 ), MPI_OFFSET_KIND )
    offset = offset * int( nb, MPI_OFFSET_KIND ) * int( nx, MPI_OFFSET_KIND )

    if( myPoolID .eq. 0 ) then
      if( nprocPerPool .eq. 0 ) then
        do ib = 1, nb, nbchunk
          nbwrite = min( nb-ib+1, nbchunk )
          call MPI_FILE_WRITE_AT( fileHandle, offset, UofX2( :, ib : ), nx*nbwrite, &
                                  MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          offset = offset + nx * nbwrite
        enddo
      else

        allocate( writeBuffer( nx, nbchunk ) )
        do ib = 1, nb, nbchunk

          do iib = ib, min( ib + nbchunk -1, nb )

            xstart = 1
            xstop = prep_wvfn_divideXmesh( nx, nprocPerPool, 0 )
            writeBuffer( xstart:xstop, iib-ib+1 ) = UofX2( :, iib )
            do id = 1, nprocPerPool-1
              xstart = xstop+1
              myx = prep_wvfn_divideXmesh( nx, nprocPerPool, id )
              xstop = xstart + myx - 1
              call MPI_RECV( writeBuffer( xstart:xstop, iib-ib+1), myx, MPI_DOUBLE_COMPLEX, id, 1, &
                                          poolComm, MPI_STATUS_IGNORE, ierr )
              if( ierr .ne. 0 ) return
            enddo
          enddo

  !JTV !TODO
  ! Here (and above) we should do this in multiple steps if the wave functions are really big
          if( debug ) write(1000+myid,'(A2,X,3(I12))') 'OF', ikpt, ib, offset
          nbwrite = min( nb-ib+1, nbchunk )
          call MPI_FILE_WRITE_AT( fileHandle, offset, writeBuffer, nx*nbwrite, MPI_DOUBLE_COMPLEX, &
                                  MPI_STATUS_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          offset = offset + nx * nbwrite
        enddo      
        deallocate( writeBuffer )
      endif
    else
      do ib = 1, nb
        call MPI_SEND( UofX2( :, ib ), myx, MPI_DOUBLE_COMPLEX, 0, 1, poolComm, MPI_STATUS_IGNORE, ierr )
        if( ierr .ne. 0 ) return
      enddo
    endif
    

  end subroutine prep_wvfn_writeU2_perKpt

  subroutine prep_wvfn_writeU2( ikpt, ispin, UofX2, fileHandle, ierr )
    use prep_system, only : system_parameters, params
    use ocean_dft_files, only : odf_poolID, odf_nprocPerPool
    use ocean_mpi, only : & 
#ifdef MPI_F08
                          MPI_Datatype, &
#endif
                          MPI_DOUBLE_COMPLEX, MPI_OFFSET_KIND, MPI_STATUS_IGNORE, myid, MPI_STATUS_SIZE, MPI_INFO_NULL

    integer, intent( in ) :: ikpt, ispin, fileHandle
    complex(DP), intent( in ) :: UofX2(:,:)
    integer, intent( inout ) :: ierr

#ifdef MPI_F08
    type( MPI_DATATYPE ):: newType
#else
    integer :: newType, stat( MPI_STATUS_SIZE)
#endif
    integer(MPI_OFFSET_KIND) :: offset
    integer :: ix, i, myPoolID, nprocPerPool, nx, myx, nb

    ! offset for k-point
    offset = max( ikpt - 1 + ( ispin - 1 ) * params%nkpts, 0 )
    offset = offset * int( size( UofX2, 1 ), MPI_OFFSET_KIND ) * int( size( UofX2, 2 ), MPI_OFFSET_KIND )
!    offset = offset * int( product( params%xmesh(:) ), MPI_OFFSET_KIND ) & 
!           * int( size( UofX2, 2 ), MPI_OFFSET_KIND )


    myPoolID = odf_poolID()
    nprocPerPool = odf_nprocPerPool()

    ! offset for x position within k-point
    ix = 0
    do i = 0, myPoolID - 1
      ix = ix + prep_wvfn_divideXmesh( product( params%xmesh(:) ), nprocPerPool, i )
    enddo
!    offset = offset + ix

    nb = size( UofX2, 2 )
    myx = size( UofX2, 1)
    nx = product( params%xmesh(:) )

    call MPI_TYPE_VECTOR( nb, myx, nx, MPI_DOUBLE_COMPLEX, newType, ierr )
!    call MPI_TYPE_VECTOR( size( UofX2, 2 ), size( UofX2, 1), product( params%xmesh(:) ), & 
!                          MPI_DOUBLE_COMPLEX, newType, ierr )
    if( ierr .ne. 0 ) return
    call MPI_TYPE_COMMIT( newType, ierr )
    if( ierr .ne. 0 ) return


    if( ikpt .eq. 0 ) then
      i = 0
    else
      i = myx*nb
    endif

    if( debug ) then
      write(1000+myid,'(A2,X,6(I12))') 'OF', ikpt, offset, size( UofX2, 2 ), size( UofX2, 1), & 
            offset*16, ( offset + size( UofX2, 2 ) * size( UofX2, 1) ) * 16
    endif
!    call MPI_FILE_SET_VIEW( fileHandle, offset, MPI_DOUBLE_COMPLEX, newType, "native", MPI_INFO_NULL, ierr )
!    if( ierr .ne. 0 ) return
!    call MPI_FILE_WRITE_ALL( fileHandle,  UofX2, i, newType, stat, ierr )
!
!!    call MPI_FILE_WRITE_AT_ALL( fileHandle, offset, UofX2, i, newType, MPI_STATUS_IGNORE, ierr )
!    call MPI_FILE_WRITE_AT_ALL( fileHandle, offset, UofX2, i, newType, stat, ierr )
    call MPI_FILE_WRITE_AT_ALL( fileHandle, offset, UofX2, i, MPI_DOUBLE_COMPLEX, stat, ierr )
    if( ierr .ne. 0 ) return

    call MPI_GET_COUNT( stat, MPI_DOUBLE_COMPLEX, i, ierr )
    if( debug ) write(1000+myid,*) i
  
    call MPI_TYPE_FREE( newType, ierr )
    if( ierr .ne. 0 ) return

  end subroutine prep_wvfn_writeU2

  
  pure function prep_wvfn_divideXmesh( nx, poolSize, myid ) result( myx )
    integer, intent( in ) :: nx, poolSize, myid
    integer :: myx

    integer :: remain, i

    myx = 0
    remain = nx

    do i = 0, myid 
      myx = remain / ( poolSize - i )
      remain = remain - myx
    enddo
    
  end function prep_wvfn_divideXmesh

  subroutine prep_wvfn_checkFFT( nG, gvecs, wantCKS, wantU2, fftGrid, ierr )
    use ocean_mpi, only : myid
    use prep_system, only : system_parameters, params
    integer, intent( in ) :: nG, gvecs(3,nG)
    logical, intent( in ) :: wantCKS, wantU2
    integer, intent( out ) :: fftGrid( 3 )
    integer, intent( inout ) :: ierr
    !
    integer :: boundaries( 3, 2 )
    integer :: i, j, k, test, tx

    boundaries(:,1) = minval( gvecs, 2 )
    boundaries(:,2) = maxval( gvecs, 2 )
    do i = 1, 3
      if( (-boundaries(i,2)) .lt. boundaries(i,1) ) boundaries(i,1) = -boundaries(i,2)
      if( (-boundaries(i,1)) .gt. boundaries(i,2) ) boundaries(i,2) = -boundaries(i,1)
    enddo
    fftGrid(:) = boundaries(:,2) - boundaries(:,1) + 1
    write(1000+myid,'(A,3(X,I8))') 'FFT grid', fftGrid(:)

    if( wantCKS ) fftGrid(:) = fftGrid(:) * 2
    write(1000+myid,'(A,3(X,I8))') 'FFT grid', fftGrid(:)

!    if( wantU2 ) then
!      do i = 1, 3
!        do j = 0, params%xmesh( i )
!          if( mod( fftGrid( i ), params%xmesh( i ) ) .eq. 0 ) then
!            exit
!          else
!            fftGrid( i ) = fftGrid( i ) + 1
!          endif
!        enddo
!        if( mod( fftGrid( i ), params%xmesh( i ) ) .ne. 0 ) then
!          ierr = i
!          return
!        endif
!      enddo
!    else
      tx = 1
      do i = 1, 3
        if( wantU2 ) tx = params%xmesh( i )
! 60 is the largest gap in the 5-smooth numbers below 1000
! for 7-smooth this drops to 45, but FFTW3 can do any factor
!  with FFTW3 support, only check a few larger options before moving on
#ifdef __FFTW3 
        do k = 0, 5 + tx
#else
        do k = 0, 60 + tx
#endif
          test = fftGrid(i) + k
          if( mod( test, tx ) .ne. 0 ) cycle 
          if( debug ) write(1000+myid,*) 'TESTING: ', test
          do
            if( mod( test, 2 ) .ne. 0 ) exit
            test = test / 2
            if( debug ) write(1000+myid,*) '   ', 2
          enddo
          do j = 5, 3, -2
            do
              if( mod( test, j ) .ne. 0 ) exit
              test = test / j
              if( debug ) write(1000+myid,*) '   ', j
            enddo
          enddo
#ifdef __FFTW3
          do
            if( mod( test, 7 ) .ne. 0 ) exit
            test = test / 7
            if( debug ) write(1000+myid,*) '   ', 7
          enddo
! by default FFTW3 is best on 7-smooth numbers 
!          if( mod( test, 11 ) .eq. 0 ) then
!            test = test / 11
!            if( debug ) write(1000+myid,*) '   ', 11
!          endif
!          if( mod( test, 13 ) .eq. 0 ) then
!            test = test / 13
!            if( debug ) write(1000+myid,*) '   ', 13
!          endif
#endif
          if( test .eq. 1 ) then
            fftGrid(i) = fftGrid(i) + k
            if( debug ) write(1000+myid,*) 'TESTING:    ', fftGrid(i), tx
            exit
          endif

        enddo
        
        ! If we need U2 then we need an FFT grid that is compatible
        !  if we have't found one we can recover for FFTW since
        !  it allows arbitrary lengths. 
        if( mod( fftGrid(i), tx ) .ne. 0 ) then
#ifdef __FFTW3
          do k = 0, params%xmesh( i )
            if( mod( fftGrid( i ), params%xmesh( i ) ) .eq. 0 ) then
              exit
            else
              fftGrid( i ) = fftGrid( i ) + 1
            endif
          enddo
#else
          write(1000+myid,*) 'Failed to find acceptable FFT grid', fftGrid(i), tx
          ierr = 101
          return
#endif
        endif

      enddo
!    endif
    write(1000+myid,'(A,3(X,I8))') 'FFT grid', fftGrid(:)

  end subroutine prep_wvfn_checkFFT


  subroutine prep_wvfn_doFFT( gvecs, UofG, UofX )
    use ocean_dft_files, only : odf_isFullStorage
    use ocean_mpi, only : myid
    use FFT_wrapper, only : fft_wrapper_init, fft_obj, FFT_wrapper_single, OCEAN_BACKWARD, &
                            FFT_wrapper_delete
!#ifdef __FFTW3
!    use iso_c_binding
!    include 'fftw3.f03'
!#endif
    integer, intent( in ) :: gvecs( :,: )
    complex(DP), intent( in ) :: uofg( :, : )
    complex(DP), intent( inout ) :: uofx(:,:,:,:)
    !
!#ifdef __FFTW3
!    type(C_PTR), allocatable :: bplan( : )  
    type(fft_obj), allocatable :: bplan(:)
    integer :: dims(3), nbands, ngvecs
    integer :: i, j, k, ig, ib
    real(DP), external :: dznrm2

    dims(1) = size( uofx, 1 )
    dims(2) = size( uofx, 2 )
    dims(3) = size( uofx, 3 )
    nbands = size( uofx, 4 )
    ngvecs = size( gvecs, 2 )
    if( dims(1) .lt. 1 ) return
    ! almost certainly need to reverse dims
!    bplan = fftw_plan_many_dft( 3, dims, nbands, uofx, 

    allocate( bplan( nbands ) )
    ! should be moved to plan many or something
    ib = 1
    call FFT_wrapper_init( dims, bplan( ib ), uofx(:,:,:,ib), 1000+myid )
    do ib = 2, nbands
!      bplan(ib) = fftw_plan_dft_3d( dims(3), dims(2), dims(1), uofx(:,:,:,ib), uofx(:,:,:,ib), &
!                  FFTW_BACKWARD, FFTW_PATIENT )
      call FFT_wrapper_init( dims, bplan( ib ), uofx(:,:,:,ib) )
    enddo

#ifdef __FFTW3
!$OMP PARALLEL DO DEFAULT( NONE ) &
!$OMP SHARED( nbands, dims, ngvecs, gvecs, uofx, uofg, bplan ) &
!$OMP PRIVATE( ib, ig, i, j, k )
#endif
    do ib = 1, nbands
      uofx(:,:,:,ib) = 0.0_DP
      do ig = 1, ngvecs
        i = 1 + gvecs(1,ig)
        j = 1 + gvecs(2,ig)
        k = 1 + gvecs(3,ig)
        if( i .le. 0 ) i = i + dims(1)
        if( j .le. 0 ) j = j + dims(2)
        if( k .le. 0 ) k = k + dims(3)

        uofx( i, j, k, ib ) = uofg( ig, ib )
      enddo
      if( .not. odf_isFullStorage() ) then
        do ig = 1, ngvecs
          i = 1 - gvecs(1,ig)
          j = 1 - gvecs(2,ig)
          k = 1 - gvecs(3,ig)
          if( i .le. 0 ) i = i + dims(1)
          if( j .le. 0 ) j = j + dims(2)
          if( k .le. 0 ) k = k + dims(3)

          uofx( i, j, k, ib ) = conjg( uofg( ig, ib ) )
        enddo
      endif

!      call fftw_execute_dft( bplan( ib ), uofx(:,:,:,ib), uofx(:,:,:,ib) )
      call  FFT_wrapper_single( uofx(:,:,:,ib), OCEAN_BACKWARD, bplan(ib), .false. )

!      call fftw_destroy_plan( bplan )

!      write(1000+myid, * ) '  ', ib, dznrm2( product( dims(1:3) ), uofx(:,:,:,ib), 1 )
    enddo
#ifdef __FFTW3
!$OMP END PARALLEL DO
#endif

    do ib = 1, nbands
!      call fftw_destroy_plan( bplan( ib ) )
      call FFT_wrapper_delete( bplan( ib ) )
    enddo
    deallocate( bplan )

!#else
!    ! To keep the compiler happy
!    uofx = 0.0_DP
!#endif
  end subroutine prep_wvfn_doFFT

  subroutine prep_wvfn_doLegacyParallel( ierr )
    use prep_system, only : system_parameters, params
    use ocean_mpi, only : myid, root, comm
     
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: val_wvfn(:,:,:,:), con_wvfn(:,:,:,:)
    complex(DP) :: old
    integer :: ik, ix, ib, vb, cb, nb, nx, dumi(3), ny, nz, iy, iz, iix, is
    logical :: ex

    if( myid .eq. root ) then

!    nx = product( params%xmesh(:) )
    nx = params%xmesh(1)
    ny = params%xmesh(2)
    nz = params%xmesh(3)
    vb = params%brange(2) - params%brange(1) + 1
    cb = params%brange(4) - params%brange(3) + 1
    nb = vb + cb
!    write(6,*) vb, cb, nx
!    write(1000,*) 'FINAL'

    open( unit=99, file='val.u2.dat', form='unformatted', access='stream' )
    open( unit=98, file='con.u2.dat', form='unformatted', access='stream' )

!    inquire( file='oldu2.dat', exist = ex )
!    if( ex ) open( unit=96, file='oldu2.dat',form='unformatted' )
    open( unit=97, file='u2.dat', form='unformatted', status='unknown' )
    rewind( 97 )

    allocate( val_wvfn( nx, ny, nz, vb ), con_wvfn( nx, ny, nz, cb ) )
    do is = 1, params%nspin
    do ik = 1, params%nkpts
      read(99) val_wvfn

      do ib = 1, vb
        iix = 0
        do ix = 1, nx
          do iy = 1, ny
            do iz = 1, nz
              iix = iix + 1
              write(97) ib, ik, iix, val_wvfn( ix, iy, iz, ib )
!          if( ex ) then
!            read( 96 ) dumi(:), old
!            write(1000,'(4(E25.16))') val_wvfn( ix, ib ), old
!            write(2000,'(3I8,2E25.8)') ib, ik, ix, val_wvfn( ix, ib )
!            write(3000,'(3I8,2E25.8)') ib, ik, ix, old
!
!          endif
            enddo
          enddo
        enddo
      enddo

      read(98) con_wvfn
      do ib = 1, cb
        iix = 0
        do ix = 1, nx
          do iy = 1, ny
            do iz = 1, nz
              iix = iix + 1
              write(97) ib+vb, ik, iix, con_wvfn( ix, iy, iz, ib )
!          if( ex ) then
!            read( 96 ) dumi(:), old
!            write(1000,'(4(E25.16))') con_wvfn( ix, ib ), old
!            write(2001,'(3I8,2E25.8)') ib+vb, ik, ix, con_wvfn( ix, ib )
!            write(3001,'(3I8,2E25.8)') ib+vb, ik, ix, old
!          endif
            enddo
          enddo
        enddo
      enddo

!      if( ex ) write(1000,*) '!##!', ik

!!      if( ik .eq. 1 ) then
!        do ib = 1, vb
!          write(6,*) ib, dot_product( val_wvfn( :, ib), val_wvfn( :, ib) )
!        enddo
!        do ib = 1, cb
!          write(6,*) ib+vb, dot_product( con_wvfn( :, ib), con_wvfn( :, ib) )
!        enddo
!!      endif
    enddo
    enddo
    deallocate( val_wvfn, con_wvfn )
    close( 99 )
    close( 98 )
    close( 97 )

    endif
    call MPI_BARRIER( comm, ierr )
  end subroutine prep_wvfn_doLegacyParallel

  

  ! everybody (that is working on this k-point) gets a complete copy (all bands) of the valence
  subroutine shareValence( valUofG, spareValUofG, ierr )
#ifdef MPI_F08
    use ocean_mpi, only : MPI_DOUBLE_COMPLEX, MPI_COMM
#else
    use ocean_mpi, only : MPI_DOUBLE_COMPLEX
#endif
    use ocean_dft_files, only : odf_poolComm, odf_nprocPerPool, odf_poolID, odf_getBandsForPoolID, ODF_VALENCE
    complex(DP), intent( in ) :: valUofG( :, : )
    complex(DP), intent( out ) :: spareValUofG( :, : )
    integer, intent( inout ) :: ierr

    integer :: nproc, i, nbands, ng, iband, myid_
#ifdef MPI_F08
    type( MPI_COMM ) :: myComm
#else
    integer :: myComm
#endif

    myComm = odf_poolComm()
    nproc = odf_nprocPerPool()
    ng = size( valUofG, 1 )
    myid_ = odf_poolID()

    iband = 1
    do i = 0, nproc - 1
      nbands = odf_getBandsForPoolID( i, ODF_VALENCE )

      if( i .eq. myid_ ) then
        spareValUofG( :, iband : iband+nbands-1 ) = valUofG(:,:)
      endif

      call MPI_BCAST( spareValUofG( :, iband : iband+nbands-1 ), ng*nbands, MPI_DOUBLE_COMPLEX, & 
                      i, myComm, ierr )
      if( ierr .ne. 0 ) return

    enddo


  end subroutine shareValence

end module prep_wvfn
