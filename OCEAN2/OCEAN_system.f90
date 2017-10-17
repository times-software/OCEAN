! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_system
  use AI_kinds
!  use mpi
  implicit none


  type, public :: o_system
    real(DP)         :: celvol
    real(DP)         :: avec(3,3)
    real(DP)         :: bvec(3,3)
    real(DP)         :: bmet(3,3)
    real(DP)         :: qinunitsofbvectors(3)
    real(DP)         :: epsilon0
    integer( S_INT ) :: nkpts
    integer( S_INT ) :: nxpts
    integer( S_INT ) :: nalpha
    integer( S_INT ) :: nbeta
    integer( S_INT ) :: num_bands
    integer( S_INT ) :: val_bands
    integer          :: brange(4)
    integer          :: nelectron 
    integer( S_INT ) :: xmesh( 3 )
    integer( S_INT ) :: kmesh( 3 )
    integer( S_INT ) :: ZNL(3)
    integer( S_INT ) :: nspn 
    integer          :: valence_ham_spin
    integer( S_INT ) :: nobf = 0
    integer( S_INT ) :: nruns
    integer          :: nedges 
    integer          :: nXES_photon

    integer          :: tmel_selector
    integer          :: enk_selector
    integer          :: bloch_selector

    logical          :: e0
    logical          :: mult
    logical          :: long_range
    logical          :: obf
    logical          :: conduct
    logical          :: valence =.false.
    logical          :: kshift
    logical          :: have_core = .true.
    logical          :: have_val = .false.
    logical          :: backf = .false.
    logical          :: write_rhs 
    logical          :: complex_bse
    character(len=5) :: calc_type

    type(o_run), pointer :: cur_run => null()

  end type o_system


  type :: o_run
    real(DP) :: tau(3)
    integer( S_INT ) :: nalpha
    integer( S_INT ) :: ZNL(3)
    
    integer :: indx
    integer :: photon
    integer :: num_bands
    integer :: val_bands
    integer :: start_band
    integer :: rixs_energy 
    integer :: rixs_pol
    character(len=2) :: elname
    character(len=2) :: corelevel
    character(len=255) :: basename
    character(len=255) :: filename

    character(len=3) :: calc_type

    logical          :: e0
    logical          :: mult
    logical          :: long_range
    logical          :: obf
    logical          :: conduct
    logical          :: valence = .false.
    logical          :: have_core  
    logical          :: have_val  
    logical          :: lflag
    logical          :: bflag
    logical          :: bande = .true.
    logical          :: aldaf
    logical          :: backf
    logical          :: complex_bse
    
    type(o_run), pointer :: prev_run => null()
    type(o_run), pointer :: next_run => null()

  end type o_run


  

  contains 

  subroutine OCEAN_sys_update( sys, ierr )
    implicit none

    type( o_system ), intent( inout ) :: sys
    integer, intent( inout ) :: ierr

    if( .not. associated( sys%cur_run ) ) then
      ierr = -5
      return
    endif

    if( associated(sys%cur_run%next_run ) ) then
      sys%cur_run => sys%cur_run%next_run
    else
      sys%cur_run => null()
    endif
  end subroutine OCEAN_sys_update

  subroutine OCEAN_sys_init( sys, ierr )
    use OCEAN_mpi!, ONLY : myid, comm, root, nproc
    implicit none
     

    type( o_system ), intent( inout ) :: sys
    integer, intent( inout ) :: ierr

    real( DP ) :: inter
    real( DP ), parameter :: inter_min = 0.000001
    integer :: nruns 
    logical :: file_exist

    logical :: exst

    if( myid .eq. root ) then

      open(unit=99,file='xmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%xmesh(:)
      close(99)
      sys%nxpts = product( sys%xmesh(:) )

      open(unit=99,file='kmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%kmesh(:)
      close(99)
      sys%nkpts = product( sys%kmesh(:) )

      open(unit=99,file='nspin',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%nspn
      close(99)
      sys%nbeta = sys%nspn**2
      sys%valence_ham_spin = sys%nspn

      open(unit=99,file='ZNL',form='formatted',status='old')
      rewind(99) 
      read(99,*) sys%ZNL(:)
      write(6,*) sys%ZNL(:)
      close(99) 
      ! nalpha is ( nspin valence ) * ( nspin core ) * ( 2 * l_core + 1 )
      sys%nalpha = 4 * ( 2 * sys%ZNL(3) + 1 )
      write(6,*) sys%nalpha
      if( sys%ZNL(3) .gt. 0 ) then
        if( sys%valence_ham_spin .eq. 1 ) then
          write(6,*) 'Detected L>0 while spin=1'
          write(6,*) '  Upgrading valence Hamiltonian to spin=2'
          sys%nbeta = 4
          sys%valence_ham_spin = 2
        endif
      endif

      inquire( file='force_val_ham_spin.ipt', exist=exst )
      if( exst ) then
        open( unit=99, file='force_val_ham_spin.ipt', form='formatted',status='old')
        read(99,*) sys%valence_ham_spin
        close(99)
        sys%nbeta = sys%valence_ham_spin * sys%valence_ham_spin
        write(6,*) '  Override of valence Hamiltonian spin requested:', sys%valence_ham_spin
      endif

      open(unit=99,file='nbuse.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%num_bands
      close(99)

      open(unit=99,file='qinunitsofbvectors.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%qinunitsofbvectors(:)
      close(99)

      if( sum( abs( sys%qinunitsofbvectors(:) ) ) .gt. 1.0d-14 ) then
        sys%kshift = .true.
      else
        sys%kshift = .false.
      endif

      open(unit=98,file='brange.ipt',form='formatted',status='old')
      rewind(98)
      read(98,*) sys%brange(1:2)
      read(98,*) sys%brange(3:4)
      close(98)

      sys%val_bands = sys%brange(2) - sys%brange(1) + 1

      call getabb( sys%avec, sys%bvec, sys%bmet )
      call getomega( sys%avec, sys%celvol )     


      sys%mult = .true.
      inquire(file="mult.ipt",exist=file_exist)
      if( file_exist ) then
        open(unit=99,file='mult.ipt',form='formatted',status='old')
        rewind(99)
        read(99,*) sys%mult
        close(99)
      endif
 
      sys%long_range = .true.

      open(unit=99,file='cks.normal',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%conduct
      close(99)
      if( .not. sys%conduct ) then
        ! Need to run mult to get spin-orbit
        sys%long_range = .false.
      endif

      open(unit=99,file='mode',form='formatted',status='old')
      rewind(99)
      read(99,*) inter
      close(99)
      if( inter .lt. inter_min ) then
        ! Need to run mult to get spin-orbit
        sys%long_range = .false.
      endif
      

      sys%e0 = .true.
      sys%obf = .false.
      sys%calc_type = 'NaN'
!      sys%conduct = .true.

      open(unit=99,file='nelectron',form='formatted',status='old')
      read(99,*) sys%nelectron
      close(99)

      inquire(file='nedges', exist=exst )
      if( exst ) then
        open(unit=99,file='nedges',form='formatted',status='old')
        read(99,*) sys%nedges
        close(99)
      else
        sys%nedges = 0
      endif

      
      open(unit=99,file='epsilon',form='formatted', status='old' )
      read(99,*) sys%epsilon0
      close( 99 )

      open(unit=99,file='tmel_selector',form='formatted',status='old')
      read(99,*) sys%tmel_selector
      close(99)

      open(unit=99,file='enk_selector',form='formatted',status='old')
      read(99,*) sys%enk_selector
      close(99)

      open(unit=99,file='bloch_selector',form='formatted',status='old')
      read(99,*) sys%bloch_selector
      close(99)

      inquire(file='nXES_photon.ipt', exist=exst )
      if( exst ) then
        open(unit=99,file='nXES_photon.ipt',form='formatted',status='old')
        read(99,*) sys%nXES_photon
        close(99)
      else
        sys%nXES_photon = -1
      endif


      open(unit=99,file='cnbse.write_rhs',form='formatted',status='old')
      read(99,*) sys%write_rhs
      close(99)


      inquire(file='force_complex_bse.ipt', exist=exst )
      if( exst ) then
        open( unit=99, file='force_complex_bse.ipt', form='formatted',status='old')
        read( 99, * ) sys%complex_bse
        close( 99 )
      else
        sys%complex_bse = .false.
      endif
      
    endif
#ifdef MPI
! Could create an mpi_datatype, but probably not worth it

    if( nproc .gt. 1 ) then

    call MPI_BCAST( sys%celvol, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%avec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bvec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bmet, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%qinunitsofbvectors, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%epsilon0, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( sys%nkpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nxpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nalpha, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%num_bands, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%val_bands, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%xmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%kmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%ZNL, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nspn, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%valence_ham_spin, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nbeta, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%brange, 4, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nelectron, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nedges, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%tmel_selector, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%enk_selector, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bloch_selector, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nXES_photon, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111



    call MPI_BCAST( sys%e0, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%mult, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%long_range, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%obf, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%conduct, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%kshift, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%write_rhs, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%complex_bse, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( sys%calc_type, 5, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111


    endif


#endif
111 continue
!    allocate( sys%cur_run, STAT=ierr )
!    if( ierr .ne. 0 ) return
!
!    sys%cur_run%nalpha = sys%nalpha
!    sys%cur_run%ZNL(:) = sys%ZNL(:)
!    sys%cur_run%indx = 1
!    sys%cur_run%photon = 1
!    sys%cur_run%elname = 'F_'
!    sys%cur_run%corelevel = '1s'
!    sys%cur_run%basename = 'xas'
!!    write(sys%cur_run%filename,'(A,A,A,I2.2,A,A,A,I2.2)') sys%cur_run%basename, &
!!          '_', sys%cur_run%elname, '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, &
!!          '_', sys%cur_run%photon
!    sys%cur_run%tau(:) = 0.0_DP
    call OCEAN_runlist_init( sys, nruns, ierr )
    sys%nruns = nruns

  end subroutine OCEAN_sys_init

  subroutine OCEAN_runlist_init( sys, running_total, ierr )
    use OCEAN_mpi!, ONLY : myid, comm, root
    implicit none

    type( o_system ), intent( inout ) :: sys
    integer, intent( out ) :: running_total
    integer, intent( inout ) :: ierr

    real(DP) :: tau(3)
    integer :: ZNL(3), indx, photon
    character(len=2) :: elname, corelevel, ein
    character(len=5) :: calc_type
    type(o_run), pointer :: temp_prev_run, temp_cur_run

    integer :: ntot, nmatch, iter, i, start_band, num_bands, val_bands, val_flag,  &
               rixs_energy, rixs_pol
    logical :: found, have_val, have_core, lflag, bflag
    real(DP) :: tmp(3)

    ! These are optional so should be given defaults
    rixs_pol = 0
    rixs_energy = 0

    running_total = 0 
    temp_prev_run => sys%cur_run

    if( myid .eq. root ) then
      write(6,*) 'Runlist init'
      open(unit=99,file='runlist',form='formatted',status='old')
      rewind(99)
      read(99,*) running_total
    endif

#ifdef MPI
    call MPI_BCAST( running_total, 1, MPI_INTEGER, root, comm, ierr )
#endif

    do i = 1, running_total

      if( myid .eq. root ) then
        read(99,*)  ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type

        have_core = .false.
        have_val = .false.

        select case( calc_type )
        case( 'XAS' )
          start_band = sys%brange(3)
          num_bands = sys%num_bands
          have_core = .true.
          val_bands = sys%brange(2)-sys%brange(1)+1
        case( 'XES' )
          start_band = sys%brange(1)
          num_bands = sys%num_bands
          have_core = .true.
          val_bands = sys%brange(2)-sys%brange(1)+1
        case( 'VAL' )
          num_bands = sys%brange(4)-sys%brange(3)+1
          val_bands = sys%brange(2)-sys%brange(1)+1
          have_val = .true.
        case( 'RXS' )
          backspace( 99 )
          read(99,*)  ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type, &
                      rixs_energy, rixs_pol
          num_bands = sys%brange(4)-sys%brange(3)+1
          val_bands = sys%brange(2)-sys%brange(1)+1
          have_val = .true.
          
        case default
          start_band = sys%brange(3)
        end select 


        if( have_core ) then
          open(unit=98,file='xyz.wyck',form='formatted',status='old')
          rewind(98)
          read(98,*) ntot
          nmatch = 0
          found = .false.
          do iter = 1, ntot
            read ( 98, * ) ein, tmp( : )
            if ( ein .eq. elname ) then
              nmatch = nmatch + 1
              if ( nmatch .eq. indx ) then
                tau( : ) = tmp( : )
                found = .true.
              end if
            end if
            if ( found ) goto 112
          end do
          if( .not. found ) then
            ierr = -1
            return
          endif
112       continue 
          write ( 6, '(1a15,3f10.5)' ) 'snatched alpha=', tau( : )
          close(98)
        endif

        if( have_val ) then
          open(unit=98,file="lflag",form='formatted',status='old')
          rewind(98)
          read(98,*) val_flag
          close(98)
          if( val_flag .gt. 0 ) then
            lflag = .true.
          else
            lflag = .false.
          endif

          open(unit=98,file="bflag",form='formatted',status='old')
          rewind(98)
          read(98,*) val_flag
          close(98)
          if( val_flag .gt. 0 ) then
            bflag = .true.
          else
            bflag = .false.
          endif

        endif
            


      endif


#ifdef MPI
      call MPI_BCAST( tau, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( ZNL, 3, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( elname, 2, MPI_CHARACTER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( corelevel, 2, MPI_CHARACTER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( indx, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( photon, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( calc_type, 5, MPI_CHARACTER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( start_band, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( num_bands, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( val_bands, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( rixs_energy, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( rixs_pol, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( have_core, 1, MPI_LOGICAL, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( have_val, 1, MPI_LOGICAL, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( lflag, 1, MPI_LOGICAL, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
      call MPI_BCAST( bflag, 1, MPI_LOGICAL, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

      !!!!
      if( have_core ) sys%have_core = .true.
      if( have_val ) sys%have_val = .true.


      allocate( temp_cur_run, STAT=ierr )
      if( ierr .ne. 0 ) return
      if( i .eq. 1 ) then
        sys%cur_run => temp_cur_run
      else
        temp_cur_run%prev_run => temp_prev_run
        temp_prev_run%next_run => temp_cur_run
      endif
      temp_cur_run%tau(:) = tau(:)
      temp_cur_run%ZNL(:) = ZNL(:)
      temp_cur_run%nalpha = 4 * ( 2 * temp_cur_run%ZNL(3) + 1 )
      temp_cur_run%elname = elname
      temp_cur_run%indx = indx
      temp_cur_run%corelevel = corelevel
      temp_cur_run%calc_type = calc_type
      temp_cur_run%photon = photon
      temp_cur_run%start_band = start_band
      temp_cur_run%num_bands = num_bands
      temp_cur_run%val_bands = val_bands

      temp_cur_run%have_core = have_core
      temp_cur_run%have_val = have_val
      temp_cur_run%lflag = lflag
      temp_cur_run%bflag = bflag

      temp_cur_run%rixs_energy = rixs_energy
      temp_cur_run%rixs_pol = rixs_pol

      temp_cur_run%complex_bse = sys%complex_bse
      
      temp_cur_run%basename = 'abs'
      write(temp_cur_run%filename,'(A3,A1,A2,A1,I2.2,A1,A2,A1,I2.2)' ) temp_cur_run%basename, '_', temp_cur_run%elname, &
            '.', temp_cur_run%indx, '_', temp_cur_run%corelevel, '_', temp_cur_run%photon

      temp_prev_run => temp_cur_run
!      running_total = running_total + 1
    enddo


    if( running_total .lt. 1 ) then
      ierr = -1
      if(myid .eq. root ) write(6,*) 'Failed to read in any runs'
    endif
    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) then
      close( 99 )
      write(6,*) 'Number of calcs to complete: ', running_total
    endif
111 continue

    
  end subroutine OCEAN_runlist_init


    
end module OCEAN_system
