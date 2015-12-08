module OCEAN_constants
  use AI_kinds

  save

  real(DP), parameter :: eV2Rydberg = 0.073498647554693773_DP
  real(DP), parameter :: Rydberg2eV = 13.60569253_DP  !(30) 2010 CODATA
  real(DP), parameter :: eV2Hartree = 0.036749323777346887_DP
  real(DP), parameter :: Hartree2eV = 27.21138506_DP

  real(DP), parameter :: PI_DP = 3.1415926535897932384626433832795028841971693993751_DP
  real(DP), parameter :: PI_QP = 3.1415926535897932384626433832795028841971693993751_QP
  real(DP)            :: PI = 4.0_DP * ATAN( 1.0_DP )


end module OCEAN_constants
  
