! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_energies
  use AI_kinds
  use OCEAN_psi, only : OCEAN_vector

  implicit none
  save
  private


  type( OCEAN_vector ) :: p_energy
  type( OCEAN_vector ) :: allow


  LOGICAL :: have_selfenergy
  LOGICAL :: is_init = .false.
  LOGICAL :: is_loaded = .false.
  LOGICAL :: val_init = .false.
  LOGICAL :: val_loaded = .false.

  public :: OCEAN_energies_allow, OCEAN_energies_val_sfact, OCEAN_energies_val_act, &
            OCEAN_energies_val_load, OCEAN_energies_act, OCEAN_energies_init, OCEAN_energies_load
  
  contains

  subroutine OCEAN_energies_allow( sys, psi, ierr )
    use OCEAN_system, only : O_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_2element_mult
    use OCEAN_mpi, only : myid, root
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: psi
    !
    call OCEAN_psi_2element_mult( psi, allow, ierr, is_real_only=.true. )
  end subroutine

  subroutine OCEAN_energies_val_sfact( sys, psi, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_mult
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: psi
    !
    call OCEAN_psi_mult( psi, allow, .false. )
  end subroutine

  subroutine OCEAN_energies_val_act( sys, psi, hpsi, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_f2m_3element_mult
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: hpsi
    !
    logical :: backwards
    backwards = .false.
!    if( psi%val_myid .eq. 0 ) then
!      call OCEAN_psi_cmult( psi, hpsi, p_energy, .false. )
!    endif
    
    call OCEAN_psi_f2m_3element_mult( hpsi, p_energy, psi, ierr, backwards )

  end subroutine OCEAN_energies_val_act


  subroutine OCEAN_energies_val_load( sys, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_new, OCEAN_psi_zero_full
    use OCEAN_val_energy, only : OCEAN_read_energies
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    !
    if( sys%have_val .eqv. .false. ) return
    !

    if( .not. val_init ) then
!      allocate( p_energy, allow )
      call OCEAN_psi_new( p_energy, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_new( allow, ierr )
      if( ierr .ne. 0 ) return

      val_init = .true.
    endif

    if( .not. val_loaded ) then
      call OCEAN_psi_zero_full( p_energy, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_zero_full( allow, ierr )
      if( ierr .ne. 0 ) return
      
      call OCEAN_read_energies( sys, p_energy, allow, ierr )
      if( ierr .ne. 0 ) return

      val_loaded = .true.
    endif

  end subroutine OCEAN_energies_val_load
    

  subroutine OCEAN_energies_init(  sys, ierr )
    use OCEAN_system
    use OCEAN_psi
!    use OCEAN_mpi, only : myid, root

    implicit none
!    include 'fftw3.f03'
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys

    integer, parameter :: cacheline_by_Z = 1


    if( .not. is_init ) then
      call OCEAN_psi_new( p_energy, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_new( allow, ierr )
      if( ierr .ne. 0 ) return

      is_init = .true.
    endif


  end subroutine OCEAN_energies_init


  subroutine OCEAN_energies_kill( ierr )
    use OCEAN_psi, only : OCEAN_psi_kill
    implicit none
    integer, intent(inout) :: ierr

    call OCEAN_psi_kill( p_energy, ierr )
    call OCEAN_psi_kill( allow, ierr )
    is_init = .false.

  end subroutine OCEAN_energies_kill

  subroutine OCEAN_energies_load( sys, complex_bse, ierr )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_psi, only : OCEAN_psi_zero_full, OCEAN_psi_bcast_full
    use OCEAN_constants, only : eV2Hartree

    implicit none
    type(O_system), intent( in ) :: sys
    logical, intent(inout) :: complex_bse
    integer, intent(inout) :: ierr

    real(DP), allocatable :: energies(:,:,:), imag_se(:,:,:)
    real(DP) :: core_offset, efermi
    logical :: metal

    if( myid .eq. root ) then
      allocate( energies( sys%num_bands, sys%nkpts, sys%nspn ), &
                imag_se(  sys%num_bands, sys%nkpts, sys%nspn ) )
    else
      allocate( energies( 1, 1, 1 ), imag_se( 1, 1, 1 ) )
    endif
    energies = 0.0_DP
    imag_se = 0.0_DP

    call OCEAN_psi_zero_full( p_energy, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_full( allow, ierr )
    if( ierr .ne. 0 ) return


    ! Currently all of this is done only on the root proccess 
    !  the energies are only shared at the end, once they are in p_energy.
    ! However, all of the other variables are bcast within each subroutine
    call read_energies( sys, energies, ierr )

    call gw_corrections( sys, energies, imag_se, complex_bse, ierr )

    call read_coreLevelShift( sys, core_offset, ierr )

    call core_find_fermi( sys, energies, complex_bse, metal, efermi, ierr )

    call core_make_allow( sys, energies, metal, efermi )

    if( myid .eq. root ) then
      call stubby( sys, p_energy, energies, imag_se, core_offset )
    endif

    call OCEAN_psi_bcast_full( root, p_energy, ierr )
    call OCEAN_psi_bcast_full( root, allow, ierr )

    deallocate( energies, imag_se )

  end subroutine

!> @brief Wrapper for making allow array for core-levels
  subroutine core_make_allow( sys, ener, metal, efermi )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root
    !
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: ener(:,:,:)
    real(DP), intent( in ) :: efermi
    logical, intent( in ) :: metal
    !
    real(DP), allocatable :: allowArray(:,:,:)
    integer :: n(3)

    if( myid .ne. root ) return

    n = shape( ener )
    allocate( allowArray( n(1), n(2), n(3) ) )

    call core_allow_ff( sys, allowArray, ener, metal, efermi )

    call stubby( sys, allow, allowArray )


    deallocate( allowArray )

  end subroutine core_make_allow


!> brief Makes the core-level allow array using a Fermi-factor
  subroutine core_allow_ff(  sys, allowArray, ener, metal, efermi, temp )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root
    !
    type(O_system), intent( in ) :: sys
    real(DP), intent( out ) :: allowArray(:,:,:)
    real(DP), intent( in ) :: ener(:,:,:)
    logical, intent( in ) :: metal
    real(DP), intent( in ) :: efermi
    real(DP), intent( in ), optional :: temp
    !
    real(DP) :: itemp, scaledFermi, con, val
    integer :: i, j, k, n(3)

    if( present( temp ) ) then
      itemp = 1.0_dp / temp
      scaledFermi = - efermi / temp
    endif

    select case( sys%cur_run%calc_type )

      case( 'XAS' )
        con = 1.0_DP
        val = 0.0_DP

      case( 'XES' )
        con = 0.0_DP
        val = 1.0_DP

      case default
        con = 1.0_DP
        val = 0.0_DP

    end select
    
    n = shape( ener )

    if( myid .eq. root ) write(6,*) 'Dim:', n(:)

    if( present( temp ) ) then
      do k = 1, n(3)
        do j = 1, n(2)
          do i = 1, n(1)
            allowArray(i,j,k) = 1.0_dp / ( exp( ener(i,j,k) * itemp + scaledFermi ) + 1.0_dp )
          enddo
        enddo
      enddo
    else
      do k = 1, n(3)
        do j = 1, n(2)
          do i = 1, n(1)
            if( ener( i, j, k ) .ge. efermi ) then
              allowArray(i,j,k) = con
            else
!              if( myid .eq. root ) write(6,*) i,j,k, ener(i,j,k)
              allowArray(i,j,k) = val
            endif
          enddo
        enddo
      enddo
    endif
  
    open(unit=99,file='allow.txt',form='formatted', status='unknown')
    do k = 1, n(3)
      do j = 1, n(2)
        do i = 1, n(1)
          write(99,*) i,j,k,allowArray(i,j,k)
        enddo
      enddo
    enddo
    close(99)

  end subroutine core_allow_ff
  

  subroutine read_energies( sys, energies, ierr )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root, comm, MPI_SUCCESS, MPI_INTEGER
    
    type(O_system), intent( in ) :: sys
    real(DP), intent( out ) :: energies(:,:,:)
    integer, intent(inout) :: ierr

    real(DP), allocatable :: tmp_e0(:,:,:)
    real(DP) :: core_offset
    integer :: nbd, nq, nspn, iter, i, j, ierr_
    character(len=9) :: infoname
    character(len=4) :: gw_control
    logical :: file_exists, have_gw

    character(len=18) ::clsFile



    if( sys%conduct ) then
      infoname = 'wvfcninfo'
    else
      infoname = 'wvfvainfo'
    endif

    if( myid .eq. root ) then
      write(6,*) 'Loading energies ...'
      open( unit=99, file=infoname, form='unformatted', status='old' )
      rewind(99)
      if( sys%nspn .eq. 2 ) then
         read( 99 ) nbd, nq, nspn
      else
        read ( 99 ) nbd, nq
        nspn = 1
      endif
      if( nbd .ne. sys%num_bands ) then
        write(6,*) 'Band mismatch!', sys%num_bands, nbd
        ierr = 1
        goto 111
      elseif( nq .ne. sys%nkpts ) then
        write(6,*) 'K-point mismatch', sys%nkpts, nq
        ierr = 1
        goto 111
      elseif( nspn .ne. sys%nspn ) then
        write(6,*) 'Spin mismatch', sys%nspn, nspn
        ierr = 1
        goto 111
      endif

      allocate( tmp_e0( sys%num_bands,  sys%nkpts, sys%nspn ), STAT=ierr )
      if( ierr .ne. 0 ) return
      read(99) tmp_e0
      energies( 1 : sys%num_bands, 1 : sys%nkpts, : ) = tmp_e0( :, :, : )
      deallocate( tmp_e0 )
!      read(99) energies( 1 : sys%num_bands, 1 : sys%nkpts, : )


      close( 99 )
    endif

111 continue
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. MPI_SUCCESS ) then
      ierr = ierr_
      return
    endif
#endif

  end subroutine read_energies


  subroutine gw_corrections( sys, energies, imag_se, complex_bse, ierr )
    use OCEAN_system
    use OCEAN_mpi
    use OCEAN_psi
!    use mpi
    use OCEAN_constants, only : eV2Hartree

    implicit none
    type(O_system), intent( in ) :: sys
    real(DP), intent(inout) :: energies(:,:,:)
    real(DP), intent( out ) :: imag_se(:,:,:)
    logical, intent(inout) :: complex_bse
    integer, intent(inout) :: ierr

    real(DP), allocatable :: tmp_e0(:,:,:)
    real(DP) :: core_offset
    integer :: nbd, nq, nspn, iter, i, j, ierr_
    character(len=9) :: infoname
    character(len=4) :: gw_control
    logical :: file_exists, have_gw

    character(len=18) ::clsFile



    if( myid .eq. root ) then


      inquire(file='gw_control',exist=have_gw)
      if( have_gw ) then
        open(unit=99,file='gw_control',form='formatted',status='old')
        rewind(99)
        read(99,*) gw_control
        close(99)
        select case (gw_control)
        case ('full')
          if( myid .eq. root ) write(6,*) 'full'
          complex_bse = .true.
          call OCEAN_abinit_fullgw( sys, energies, imag_se, ierr, .true. )
        case ('real')
          call OCEAN_abinit_fullgw( sys, energies, imag_se, ierr, .false. )
        case( 'band' )
          call OCEAN_gw_by_band( sys, energies, imag_se, ierr, .false. )
        case( 'ibnd' )
          complex_bse = .true.
          call OCEAN_gw_by_band( sys, energies, imag_se, ierr, .true. )
        case( 'cstr' )
          call OCEAN_gw_stretch( sys, energies, ierr )
        case( 'list' )
          call OCEAN_gw_list( sys, energies, imag_se, ierr, .false. )
        case default
          write(6,*) 'Unrecognized gw_control'
        end select
      endif
    
    endif


#ifdef MPI
111 continue
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. MPI_SUCCESS ) then
      ierr = ierr_
      return
    endif

    call MPI_BCAST( complex_bse, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif


  end subroutine gw_corrections

  subroutine read_coreLevelShift( sys, cls, ierr )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root, comm, MPI_DOUBLE_PRECISION
    use OCEAN_constants, only : ev2Hartree, Hartree2eV

    type(O_system), intent( in ) :: sys
    real(DP), intent( out ) :: cls
    integer, intent( inout ) :: ierr

    character(len=18) ::clsFile
    logical :: file_exists, file_exists2, noshift_lumo

    if( myid .eq. root ) then
      write(clsFile,'(1A5,1A2,1I4.4,1A2,1I2.2,1A1,1I2.2)') 'cls.z', sys%cur_run%elname, sys%cur_run%indx, &
            '_n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)
      inquire(file=clsFile,exist=file_exists)
      if( file_exists ) then
        open(unit=99,file=clsFile,form='formatted',status='old')
        read(99,*) cls
        close(99)
      else
        inquire( file='noshift_lumo', exist=file_exists2)
        if( file_exists2 ) then
          open(unit=99,file='noshift_lumo', form='formatted',status='old')
          read(99,*) noshift_lumo
          close(99)
        else
          noshift_lumo = .false.
        endif
        if( noshift_lumo .eqv. .false. ) then
          open( unit=99, file='eshift.ipt', form='formatted', status='old' )
          rewind 99
          read ( 99, * ) cls
          close( unit=99 )
          cls = cls * eV2Hartree
        else
          cls = 0.0_DP
        endif
      endif

      write(6,*) 'Shifting energies by:', cls * Hartree2eV
    endif

#ifdef MPI
    call MPI_BCAST( cls, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine read_coreLevelShift



  subroutine stubby( sys, p, re_array, im_array, shift )
    use OCEAN_system

    implicit none
    type(O_system), intent( in ) :: sys
    type(ocean_vector), intent(inout) :: p
    real(DP), intent( in ) :: re_array(:,:,:)
    real(DP), intent( in ), optional :: im_array(:,:,:)
    real(DP), intent( in ), optional :: shift

    integer :: ialpha, ikpt, ibd, icms, icml, ivms, val_spin( sys%nalpha )


!    if( sys%nspn .ne. 1 ) ierr = -1 

    ! predefine the valence spins
    if( sys%nspn .eq. 1 ) then
      val_spin( : ) = 1
    else
      ialpha = 0
      do icms = 1, 2
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = 1, 2
            ialpha = ialpha + 1
            val_spin( ialpha ) = ivms
          enddo
        enddo
      enddo
    endif

    if( present( shift ) ) then
      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          do ibd = 1, sys%num_bands
            p%r( ibd, ikpt, ialpha ) = re_array( ibd, ikpt, val_spin(ialpha) ) + shift

          enddo
        enddo
      enddo
    else
      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          do ibd = 1, sys%num_bands
            p%r( ibd, ikpt, ialpha ) = re_array( ibd, ikpt, val_spin(ialpha) ) 
          enddo
        enddo
      enddo
    endif

    if( present( im_array ) ) then
      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          do ibd = 1, sys%num_bands
            p%i( ibd, ikpt, ialpha ) = im_array( ibd, ikpt, val_spin(ialpha) )
          enddo
        enddo
      enddo
    endif

  end subroutine stubby

  subroutine OCEAN_energies_act( sys, psi, hpsi, backwards, ierr )
    use OCEAN_system
    use OCEAN_psi

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: hpsi
    logical, intent( in ) :: backwards

!    hpsi%r(:,:,:) = p_energy%r(:,:,:) * psi%r(:,:,:)
!    hpsi%i(:,:,:) = p_energy%r(:,:,:) * psi%i(:,:,:)

    call OCEAN_psi_f2m_3element_mult( hpsi, p_energy, psi, ierr, backwards )

  end subroutine OCEAN_energies_act



#if 0
 !JTV need to move energies into an ocean_vector then:
 ! 1) this can be done w/ min instead of full
 ! 2) we can write an element-wise y(i) = z(i) * x(i) + y(i) in OCEAN_psi
  subroutine OCEAN_energies_act2( sys, psi, hpsi, backwards, ierr )
    use OCEAN_system
    use OCEAN_psi

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: hpsi
    logical, intent( in ) :: backwards

    integer :: ialpha, ikpt, ibd, icms, icml, ivms, val_spin( sys%nalpha )

!    if( sys%nspn .ne. 1 ) ierr = -1 

    ! predefine the valence spins
    if( sys%nspn .eq. 1 ) then
      val_spin( : ) = 1
    else
      ialpha = 0
      do icms = 1, 2
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = 1, 2
            ialpha = ialpha + 1
            val_spin( ialpha ) = ivms
          enddo
        enddo
      enddo
    endif


    ! if back = true then we are acting the Hermitian conjugate, i.e. i -> -i
    if( backwards ) then

      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          do ibd = 1, sys%num_bands
            hpsi%r( ibd, ikpt, ialpha ) = energies( ibd, ikpt, val_spin(ialpha) ) * psi%r( ibd, ikpt, ialpha ) &
                                 + imag_selfenergy( ibd, ikpt, val_spin(ialpha) ) * psi%i( ibd, ikpt, ialpha )
          enddo
          do ibd = 1, sys%num_bands
            hpsi%i( ibd, ikpt, ialpha ) = energies( ibd, ikpt, val_spin(ialpha) ) * psi%i( ibd, ikpt, ialpha ) &
                                 - imag_selfenergy( ibd, ikpt, val_spin(ialpha) ) * psi%r( ibd, ikpt, ialpha )

          enddo
        enddo 
      enddo

    else ! normal 

      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          do ibd = 1, sys%num_bands
            hpsi%r( ibd, ikpt, ialpha ) = energies( ibd, ikpt, val_spin(ialpha) ) * psi%r( ibd, ikpt, ialpha ) &
                                 - imag_selfenergy( ibd, ikpt, val_spin(ialpha) ) * psi%i( ibd, ikpt, ialpha )
          enddo
          do ibd = 1, sys%num_bands
            hpsi%i( ibd, ikpt, ialpha ) = energies( ibd, ikpt, val_spin(ialpha) ) * psi%i( ibd, ikpt, ialpha ) &
                                 + imag_selfenergy( ibd, ikpt, val_spin(ialpha) ) * psi%r( ibd, ikpt, ialpha )

          enddo
        enddo 
      enddo

    endif

  end subroutine OCEAN_energies_act2




  complex(DP) function OCEAN_energies_single( ib, ik, ia )
    use OCEAN_system
    implicit none
    integer, intent( in ) :: ib, ik, ia

    OCEAN_energies_single = CMPLX( energies( ib, ik, 1 ), imag_selfenergy( ib, ik, 1 ), DP )
  end function
#endif

  subroutine OCEAN_gw_stretch( sys, energies, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    real(DP), intent(inout) :: energies(:,:,:)

    real(DP) :: cstr
    logical :: have_gw

    if( myid .ne. root ) return

    write(6,*) 'Attempting GW stretch!'

    inquire(file='gwcstr', exist=have_gw )

    if( .not. have_gw ) then
      write( 6, * ) 'GW corrections requested (stretch style). File gwcstr not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    open( unit=99, file='gwcstr', form='formatted',status='old')
    rewind(99)
    read(99,*) cstr
    close(99)

    if( abs( cstr ) .lt. 0.00000001_DP ) return

    cstr = cstr + 1.0_DP
    energies( :, :, : ) = energies( :, :, : ) * cstr

  end subroutine OCEAN_gw_stretch

  subroutine OCEAN_gw_list( sys, energies, imag_selfenergy, ierr, keep_imag )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_constants, only : eV2Hartree

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    logical, intent( in ) :: keep_imag
    real(DP), intent( inout ) :: energies(:,:,:)
    real(DP), intent( out ) :: imag_selfenergy(:,:,:)

    real( DP ) :: delta_gw, imag_gw
    integer :: band_max, band_loop, ikpt, iband, ispn
    logical :: have_gw

    if( myid .ne. root ) return
    inquire( file="list_val_gw.txt", exist=have_gw )
    if( have_gw ) then
      !
      open( unit=99, file="list_val_gw.txt", form="formatted", status="old" )
      read( 99, * ) band_max
      !

      band_loop = min( band_max, sys%cur_run%num_bands )
      do iband = 1, band_loop
        do ispn = 1, sys%nspn
          do ikpt = 1, sys%nkpts
            if( keep_imag ) then
              read( 99, * ) delta_gw, imag_gw
              energies( iband, ikpt, ispn ) = energies( iband, ikpt, ispn ) + delta_gw * eV2Hartree
              imag_selfenergy( iband, ikpt, ispn ) = imag_gw
            else
              read( 99, * ) delta_gw
              energies( iband, ikpt, ispn ) = energies( iband, ikpt, ispn ) + delta_gw * eV2Hartree
            endif
          enddo
        enddo
      enddo
      write(6,* ) 'GW: list corrections from list_val_gw.txt', sys%cur_run%num_bands
      close( 99 )
    else
      write(6,* ) 'GW: list requested, but list_val_gw.txt not found!'
    endif

  end subroutine OCEAN_gw_list

  subroutine OCEAN_gw_by_band( sys, energies, imag_selfenergy, ierr, keep_imag )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_constants, only : eV2Hartree

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    logical, intent( in ) :: keep_imag
    real(DP), intent( inout ) :: energies(:,:,:)
    real(DP), intent( out ) :: imag_selfenergy(:,:,:)

    integer :: iter, ispn, kiter
    real( DP ), allocatable :: re_se( : ), im_se( : )
    logical :: have_gw

    if( myid .ne. root ) return
    
    inquire(file='GW_band.in',exist=have_gw)
    if( .not. have_gw ) then
      write( 6, * ) 'GW corrections requested (band style). File GW_band.in not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    allocate( re_se( sys%num_bands ), im_se( sys%num_bands ) )
    open(unit=99,file='GW_band.in',form='formatted',status='old')
    rewind(99)

    if( keep_imag ) then
      do iter = 1, sys%num_bands
        read(99,*) re_se( iter ), im_se( iter )
      enddo
    else
      do iter = 1, sys%num_bands
        read(99,*) re_se( iter )
      enddo
      im_se( : ) = 0.0_DP
    endif
    close( 99 )
    re_se( : ) = re_se( : ) * eV2Hartree !/ 27.21138506_DP
    im_se( : ) = -im_se( : ) * eV2Hartree !/ ( 27.21138506_DP ) 


    do ispn = 1, sys%nspn
      do kiter = 1, sys%nkpts
        energies( :, kiter, ispn ) = energies( :, kiter, ispn ) + re_se( : )
        imag_selfenergy( :, kiter, ispn ) = im_se( : )
      enddo
    enddo

    deallocate( re_se, im_se )


  end subroutine OCEAN_gw_by_band

  subroutine OCEAN_abinit_fullgw( sys, energies, imag_selfenergy, ierr, keep_imag )
! Only run on root
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_constants, only : eV2Hartree, Hartree2eV

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    logical, intent( in ) :: keep_imag
    real(DP), intent(inout) :: energies(:,:,:)
    real(DP), intent(out )  :: imag_selfenergy(:,:,:)

    integer :: iter, biter, kiter, ix, iy, iz, r_kiter, ibd, ispn
    integer :: num_band, gw_nkpt, gw_nspn, start_band, stop_band, skip_band
    integer, allocatable :: kpt_map(:), start_b(:), stop_b(:)
    real(DP), allocatable :: re_se(:), im_se(:)
    real(DP) :: re_min, im_min, re_max, im_max, kpt( 3 ), e0
    logical :: have_kpt_map, have_gw

    if( myid .ne. root ) return

    inquire(file='kpt_map',exist=have_kpt_map)
    if( .not. have_kpt_map ) then
      write( 6, * ) 'GW corrections requested. File kpt_map not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    inquire(file='GWx_GW',exist=have_gw)
    if( .not. have_gw ) then 
      write( 6, * ) 'GW corrections requested. File GWx_GW not found.'
      write( 6, * ) 'No corrections will be done'
      return 
    endif

    allocate(kpt_map( sys%nkpts ) )
    open(unit=99,file='kpt_map',form='formatted',status='old')
    do iter = 1, sys%nkpts
      read(99,*) kpt_map( iter )
    enddo

    im_min = 0.0_DP
    re_min = 0.0_DP
    re_max = 0.0_DP
    im_max = 0.0_DP
    allocate( start_b( sys%nkpts ), stop_b(sys%nkpts ) )

    open(unit=99,file='GWx_GW',form='formatted',status='old')
    read(99,*) gw_nkpt, gw_nspn

!JTV Part of this is dumb because ABINIT won't warn us in the GW file if there is an extra element 
!    at some of the k-points
    do iter = 1, gw_nkpt
      read(99,*) kpt(1)
      read(99,*) num_band
      allocate( re_se( num_band ), im_se( num_band ) )
!JTV
      read(99,*) start_band, e0, re_se( 1 ), im_se( 1 )
! Format is band GW_energy, some_measure_of_conv, imaginary part
!      read(99,*) start_band, re_se( 1 ), e0, im_se( 1 )
      stop_band = start_band + num_band - 1

      do biter = 2, num_band
        read(99,*) ibd, e0, re_se( biter ), im_se( biter )
!        read(99,*) ibd, re_se( biter ), e0, im_se( biter )
      enddo
      re_se( : ) = re_se( : ) * eV2Hartree !/ 27.21138506_DP
!JTV
      im_se( : ) = - im_se( : ) * eV2Hartree !/ 27.21138506_DP
!      im_se( : ) = -( 2.0_DP * im_se( : ) ) / 27.21138506_DP
!      im_se( : ) = abs( im_se( : ) ) / ( 27.21138506_DP * 2.0_DP )
!      im_se( 1 ) = 0.0
!      im_se( 2 ) = 0.0

!JTV TEST
!      re_se( : ) = 0.0_DP
!      im_se( : ) = 0.1_DP / 27.21138506_DP
      
      skip_band = sys%cur_run%start_band - 1
      kiter = 0
      do ix = 1, sys%kmesh(1)
        do iy = 1, sys%kmesh(2)
          do iz = 1, sys%kmesh(3)
            kiter = kiter + 1
            r_kiter = ix + (iy-1) * sys%kmesh(1) + (iz-1)*sys%kmesh(1)*sys%kmesh(2)

            if( kpt_map( r_kiter ) .eq. iter ) then
!JTV need some spin fixing 
!              im_min = im_min + im_se( 1 ) 
!              re_min = re_min + ( re_se( 1 ) - energies( 1, kiter, 1 ) )
              im_max = im_max + im_se( num_band )
!              re_max = re_max + ( re_se( num_band ) - energies( stop_band - sys%cur_run%start_band + 1 , kiter, 1 ) )
              re_max = re_max + re_se( num_band )

              start_b( kiter ) = start_band
              stop_b( kiter ) = stop_band


! biter is actual counter of bands
! se_re/se_im goes from 1 to num_bands
              do biter = max( start_band, sys%cur_run%start_band ), &
                         min( stop_band, sys%cur_run%start_band + sys%num_bands - 1 )
                do ispn = 1, sys%nspn
!                  energies( biter - sys%cur_run%start_band + 1, kiter, ispn ) = re_se( biter - start_band + 1 )
                  energies( biter - sys%cur_run%start_band + 1, kiter, ispn ) =  &
                      energies( biter - sys%cur_run%start_band + 1, kiter, ispn ) + re_se( biter - start_band + 1 )
                  imag_selfenergy( biter - sys%cur_run%start_band + 1, kiter, ispn ) =  &
                        im_se( biter - start_band + 1 )
                  if( kiter .eq. 1 ) then
                    write(6,*) biter, biter - sys%cur_run%start_band + 1, biter - start_band + 1,  &
                               im_se(  biter - start_band + 1 ) * Hartree2eV !* 27.2114_DP
                  endif
                enddo
              enddo


!              do biter = 2 + skip_band - start_band, min( num_band, sys%num_bands )
!!!               if( biter + start_band - 1 .lt. 
!                do ispn = 1, sys%nspn
!                  if( biter == 1+ skip_band ) write(6,*) re_se(biter) * 27.21138506_DP
!                  energies( biter + start_band - 1 - skip_band, kiter, ispn ) = re_se( biter ) &
!                      + energies( biter + start_band - 1 - skip_band, kiter, ispn )
!                  imag_selfenergy( biter + start_band - 1 - skip_band, kiter, ispn ) = im_se( biter )
!                enddo
!              enddo
              
            endif

          enddo
        enddo
      enddo
      deallocate( re_se, im_se )
    enddo
!    return

    im_min = im_min / dble(sys%nkpts)
    re_min = re_min / dble(sys%nkpts)
    im_max = im_max / dble(sys%nkpts)
    re_max = re_max / dble(sys%nkpts)
    
    do ispn = 1, sys%nspn
      do kiter = 1, sys%nkpts
!        do biter = 1, start_b( kiter ) - sys%cur_run%start_band - 1
!!          write(6,*) kiter, start_b( kiter ), sys%cur_run%start_band
!          energies( biter, kiter, ispn ) = energies( biter, kiter, ispn ) + re_min
!          imag_selfenergy( biter, kiter, ispn ) = im_min
!        enddo
        do biter = stop_b( kiter ) + 1 - sys%cur_run%start_band, sys%num_bands
          if( kiter == 1 ) then
            write(6,*) kiter, stop_b( kiter ), sys%num_bands, biter, re_max * Hartree2eV ! 27.2114_DP
          endif
          energies( biter, kiter, ispn ) = energies( biter, kiter, ispn ) + re_max
          imag_selfenergy( biter, kiter, ispn ) = im_max
        enddo
      enddo
    enddo

    deallocate( start_b, stop_b, kpt_map )

    if( .not. keep_imag ) then
      imag_selfenergy(:,:,:) = 0.0_DP
    endif 
!    imag_selfenergy(:,:,:) = -0.014_DP

  end subroutine OCEAN_abinit_fullgw
  


!> @brief Figures out the appropriate Fermi level 
  subroutine core_find_fermi( sys, ener, didGWcorrection, metal, efermi, ierr )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root, comm, MPI_LOGICAL, MPI_DOUBLE_PRECISION
    !
    type(O_system), intent( in ) :: sys
    real(DP), intent( in ) :: ener( sys%num_bands, sys%nkpts, sys%nspn )
    logical, intent( in ) :: didGWcorrection
    logical, intent( out ) :: metal
    real(DP), intent( out ) :: efermi
    integer, intent( inout ) :: ierr
  
    real(DP), allocatable :: sorted_energies(:)
    logical :: useFermiFromFile, doping, overrideFermi
    real(DP) :: per_electron_dope, homo, lumo, fullyOccupiedElectrons
    integer :: n_electron_dope, bandOverlap, tot_electron

    

    if( myid .eq. root ) then
      n_electron_dope = 0

      inquire(file='metal', exist=metal ) 
      if( metal ) then
        open( unit=99, file='metal', form='formatted', status='old' )
        read(99,*) metal
        close(99)
        !
        inquire( file='doping', exist=doping )
        if( doping ) then
          open( unit=99, file='doping', form='formatted', status='old' )
          read( 99, * ) per_electron_dope
          close( 99 )
          n_electron_dope = floor( per_electron_dope * dble( sys%nkpts * sys%nspn ) / 2.d0 )
          write( 6, * ) 'Doping option:'
          write( 6, * ) 'Percent doping: ', per_electron_dope
          write( 6, * ) 'Modifying electron number by ', n_electron_dope
          write( 6, * ) 'Effective doping percent: ', dble( n_electron_dope ) * 2.d0 / ( sys%nkpts * sys%nspn)
        endif
      endif

      bandOverlap = sys%brange(2) - sys%brange(3) + 1
      if( metal .and. bandOverlap .lt. 1 ) then
        write(6,*) 'Metal requested, but no overlap in bands.'
        write(6,*) 'Results may be poorly defined'
      endif

      ! if we have doping, GW+metal, or by request THEN
      !  we need to find the updated Fermi
      ! ELSE we can read from file
      inquire(file='findNewFermi.ipt', exist=overrideFermi )
      if( overrideFermi ) then
        open(unit=99,file='findNewFermi.ipt',form='formatted', status='old' )
        read(99,*) overrideFermi
        close(99)
      endif

      if( doping .or. ( didGWcorrection .and. metal ) ) then
        overrideFermi = .true.
      endif

      if( overrideFermi ) then
        allocate( sorted_energies( sys%num_bands * sys%nkpts * sys%nspn ) )
        sorted_energies = reshape( ener, (/ sys%num_bands * sys%nkpts * sys%nspn /) )
        call do_sort( sorted_energies )
        tot_electron = sys%nelectron * sys%nkpts * sys%nspn
        tot_electron = tot_electron / 2 + n_electron_dope

        fullyOccupiedElectrons = ( sys%brange(3) - 1 ) * sys%nkpts * sys%nspn
        tot_electron = tot_electron - fullyOccupiedElectrons

        homo = sorted_energies( tot_electron )
        lumo = sorted_energies( tot_electron + 1 )
    
        efermi = (homo+lumo) / 2.0_dp

        write( 6, '(A,2(F12.6,1X),E21.14)') 'Found new Fermi:', lumo, homo, efermi

        deallocate( sorted_energies )
      else

        open(unit=99, file='efermiinrydberg.ipt', form='formatted', status='old' )
        read(99,*) efermi
        close( 99 )
        efermi = efermi / 2.0_dp

      endif

    endif

#ifdef MPI
    call MPI_BCAST( metal, 1, MPI_LOGICAL, root, comm, ierr )
    call MPI_BCAST( efermi, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine core_find_fermi

  subroutine do_sort( array )
    real(DP), intent(inout) :: array(:)
    !
    integer :: top, iter, node, node2
    real(DP) :: temp

    top = size( array )
    do iter = size( array )/2, 1, -2
      temp = array(iter)
      node = iter+iter
      node2 = iter

      do
        if( node .gt. top ) goto 10
        if( ( node .lt. top ) .and. ( array(node) .lt. array( node+1 ) ) ) node = node + 1
        if( temp .lt. array( node ) ) then
          array( node2 ) = array( node )
          node2 = node
          node = node + node
        else
          goto 10
        endif
      enddo
10    continue
    enddo

    do iter = size(array), 1, -1
      temp = array( iter )
      array( iter ) = array( 1 )
      node = 2
      node2 = 1

      do
        if( node .gt. iter -1 ) goto 20
        if( ( node .lt. iter - 1 ) .and. ( array( node ) .lt. array( node + 1 ) ) ) node = node+1
        if( temp .lt. array( node ) ) then
          array(node2) = array( node )
          array( node ) = temp
          node2 = node
          node = node + node
        else
          goto 20
        endif
      enddo
20    continue
    enddo

  end subroutine do_sort

end module OCEAN_energies
