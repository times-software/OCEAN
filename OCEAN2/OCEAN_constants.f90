! Copyright (C) 2016, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_constants
  use AI_kinds

  save

  real(DP), parameter :: eV2Rydberg = 0.073498647554693773_DP
  real(DP), parameter :: Rydberg2eV = 13.60569253_DP  !(30) 2010 CODATA
  real(DP), parameter :: eV2Hartree = 0.036749323777346887_DP
  real(DP), parameter :: Hartree2eV = 27.21138506_DP

  real(DP), parameter :: alphainv = 1.0_DP / 0.0072973525664_DP
  real(DP), parameter :: bohr = 0.52917721067_DP

  ! older data from Rev Mod Phys 80 633 (2008 )
!  real(DP), parameter :: eV2Hartree = 0.0367493253979623245294_DP 
!  real(DP), parameter :: Hartree2eV = 27.21138386_DP
  ! to match ai2nbse
!  real(DP), parameter :: eV2Hartree = 0.036749303600697_DP
!  real(DP), parameter :: Hartree2eV = 27.2114_DP

  real(DP), parameter :: PI_DP = 3.1415926535897932384626433832795028841971693993751_DP
  real(DP), parameter :: PI = 4.0_DP * ATAN( 1.0_DP )
  real(QP), parameter :: PI_QP = 3.1415926535897932384626433832795028841971693993751_QP


  integer, parameter  :: CACHE_DOUBLE = 1


end module OCEAN_constants
  
