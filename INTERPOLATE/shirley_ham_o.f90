!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM shirley_ham 
  !----------------------------------------------------------------------- 

  ! David Prendergast
  ! UCB, Dec 2006

  ! Generates the Hamiltonian for Shirley Brillouin zone interpolation
  ! of electronic wave functions using as basis set an input set of
  ! basis functions which should span this space optimally.

#include "f_defs.h" 
  USE parameters, ONLY : ntypx, npk, lmaxx
  USE io_global,  ONLY : stdout, ionode
  USE kinds,      ONLY : DP 
  USE klist,      ONLY : nkstot
  USE io_files,   ONLY : nd_nmbr
  USE mp_global,  ONLY : npool
  USE control_flags,        ONLY : lscf
  USE basis,                ONLY : starting_pot, starting_wfc

  USE ions_base, ONLY : ntyp=>nsp
  use gvect
  use uspp, only : deeq, dvan, qq
  use uspp_param, only : nh
  use dfunct, only : newd
  use lsda_mod, only : nspin

  use shirley_ham_input, only : get_input, updatepp, pseudo_dir, &
                                nspecies, local_channel, pseudo_file, &
                                calculation, nspin_ham
  use elph_shirley


  IMPLICIT NONE 
  integer :: ios
  integer :: ispin

  CALL start_shirley (nd_nmbr) 

  if( npool /= 1 ) then
    call errore('shirley_ham','number of pools should be 1',abs(npool))
  endif

  call start_clock( 'shirley' )

  IF ( ionode ) THEN

     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')

!     WRITE( unit = stdout, FMT = 9010 ) & 
!         ntypx, npk, lmaxx

9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)


  END IF  

  ! ======================================================================
  ! read input and distribute
  ! ======================================================================
  call get_input

  ! ======================================================================
  ! Setup flags for a nonselfconsistent run
  ! ======================================================================
  lscf = .false.
  starting_pot = 'file'
  starting_wfc = 'file'

  ! ======================================================================
  !   Now allocate space for pwscf variables, read and check them. 
  ! ======================================================================
  call read_file_shirley( nspin_ham )
#ifndef __BUFFER
  call openfil
#endif


  if( nkstot /= 1 ) then
    call errore('shirley_ham','number of k-points should be 1',abs(nkstot))
  endif

  ! I might need to adjust spin here?

  select case (trim(calculation))
  case( 'fullocean')
    ! ======================================================================
    ! Initialize the Hamiltonian from PWSCF
    ! ======================================================================
    ! Should I turn on spin here? Or before this select block?
    CALL hinit0()
    CALL potinit()

    CALL newd()

!    ! are these guys spin-dependent?
!    write(stdout,*) ' dvan: '
!    write(stdout,'(i6,e12.5)') (ios, dvan(ios,ios,1), ios=1,nh(1))
!
!    write(stdout,*) ' qq: '
!    write(stdout,'(i6,e12.5)') (ios, qq(ios,ios,1), ios=1,nh(1))
!
!    do ispin=1,nspin
!      write(stdout,*) ' spin ', ispin
!      write(stdout,*) ' deeq: '
!      write(stdout,'(i6,e12.5)') (ios, deeq(ios,ios,1,ispin), ios=1,nh(1))
!    enddo

    ! ======================================================================
    ! Update pseudopotentials
    ! ======================================================================
    if( updatepp ) then
      call update_pseudo( nspecies, pseudo_dir, pseudo_file )
      write(stdout,*) ' updated deeq: '
      write(stdout,'(i6,e12.5)') (ios, deeq(ios,ios,1,1), ios=1,nh(1))
    else
      ! no specific local-channel - but the array must be allocated anyway
      allocate( local_channel(ntyp) )
      local_channel = -1 ! this is the do-nothing default value
    endif
    ! ======================================================================
    ! make the shirley hamiltonian
    ! ======================================================================
    call hamq()
    !
    ! Run the x-grid for u-functions
    call OCEAN_bofx()
    ! Run the radial grid for screening
    call OCEAN_bofr_multi()
    ! Make OBF <-> PAW projections
    call OCEAN_obf2localbyk()
  case( 'obf2local')
    call OCEAN_obf2localbyk()
  case( 'ocean_bofx')
    call OCEAN_bofx()
  case ('ocean_bofr_multi')
    call OCEAN_bofr_multi()
  case( 'elphmtxel' )
    ! ======================================================================
    ! electron-phonon matrix elements
    ! ======================================================================
    call elphmat()
    !
  case( 'pwmtxel' )
    ! ======================================================================
    ! make the plane-wave matrix elements
    ! ======================================================================
    call pwmat()
    !
  case( 'triples' )
    ! ======================================================================
    ! make the triples < B_k | B_i^* B_j >
    ! ======================================================================
    call triples()
    !
  case( 'interaction_test' )
    call interaction_test()
    !
  case( 'interaction' )
    ! ======================================================================
    ! make the root(interaction) matrix elements
    ! ======================================================================
    call interaction()
    call B0()
    !
  case( 'interaction0' )
    ! ======================================================================
    ! make the interaction matrix elements at the Gamma point 
    ! excluding the divergent term when G=0
    ! ======================================================================
    call interaction0()
    !
  case( 'qinteraction' )
    ! ======================================================================
    ! make the q*root(interaction) matrix elements
    ! ======================================================================
    call qinteraction()
    !
  case( 'B0' )
    ! ======================================================================
    ! dump the zeroth Fourier coefficients of the basis functions
    ! ======================================================================
    call B0()
    !
  case( 'pole_strength' )
    ! ======================================================================
    ! calculate the plasmon-pole strengths in the Shirley basis
    ! ======================================================================
    call pole_strength()
    !
!  case( 'rhoG' )
!    ! ======================================================================
!    ! dump the charge density fourier coefficients
!    ! ======================================================================
!    call rhoG()
!    !
!  case( 'vxc' )
!    ! ======================================================================
!    ! the exchange correlation potential
!    ! ======================================================================
!    ! not sure if hinit0 is necessary apart from core charge density for nlccs
!    CALL hinit0()
!    ! modified version of potinit which only accesses Vxc
!    CALL vxcinit()
!    ! bracket with the optimal basis
!    call vxc()
!    !
  case( 'plot' )
    ! ======================================================================
    ! plot the basis functions on real space grids
    ! ======================================================================
    call plot_basis()
    !
  case( 'proj' )
    ! ======================================================================
    ! Initialize the Hamiltonian from PWSCF
    ! ======================================================================
    CALL hinit0()
    CALL potinit()
    CALL newd()
    ! ======================================================================
    ! Update pseudopotentials
    ! ======================================================================
    if( updatepp ) then
      call update_pseudo( nspecies, pseudo_dir, pseudo_file )
      write(stdout,*) ' updated deeq: '
      write(stdout,'(i6,e12.5)') (ios, deeq(ios,ios,1,1), ios=1,nh(1))
    else
      ! no specific local-channel - but the array must be allocated anyway
      allocate( local_channel(ntyp) )
      local_channel = -1 ! this is the do-nothing default value
    endif
    ! ======================================================================
    ! Calculate the atomic/angular momentum projections of the basis functions
    ! ======================================================================
    ! not implemented - this should be done
    !call proj_basis()
    write(stdout,*) ' done'
    !
  case( 'ham' )
    ! ======================================================================
    ! Initialize the Hamiltonian from PWSCF
    ! ======================================================================
    ! Should I turn on spin here? Or before this select block?
    CALL hinit0()
    CALL potinit()

    CALL newd()

!    ! are these guys spin-dependent?
!    write(stdout,*) ' dvan: '
!    write(stdout,'(i6,e12.5)') (ios, dvan(ios,ios,1), ios=1,nh(1))
!
!    write(stdout,*) ' qq: '
!    write(stdout,'(i6,e12.5)') (ios, qq(ios,ios,1), ios=1,nh(1))
!
!    do ispin=1,nspin
!      write(stdout,*) ' spin ', ispin
!      write(stdout,*) ' deeq: '
!      write(stdout,'(i6,e12.5)') (ios, deeq(ios,ios,1,ispin), ios=1,nh(1))
!    enddo

    ! ======================================================================
    ! Update pseudopotentials
    ! ======================================================================
    if( updatepp ) then
      call update_pseudo( nspecies, pseudo_dir, pseudo_file )
      write(stdout,*) ' updated deeq: '
      write(stdout,'(i6,e12.5)') (ios, deeq(ios,ios,1,1), ios=1,nh(1))
    else
      ! no specific local-channel - but the array must be allocated anyway
      allocate( local_channel(ntyp) )
      local_channel = -1 ! this is the do-nothing default value
    endif
    ! ======================================================================
    ! make the shirley hamiltonian
    ! ======================================================================
    call hamq()
    !
  case default
    call errore('shirley_ham','unrecognized calculation type',1)
  end select
  ! ======================================================================
  ! stop
  ! ======================================================================
  call stop_clock( 'shirley' )
  call stop_shirley

END PROGRAM shirley_ham 
!
