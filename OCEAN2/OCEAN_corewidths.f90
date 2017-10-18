! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! Copyright on the code only. No claims on data
!
! Data from J.L. Campbell and Tibor Papp,Atomic Data and Nuclear Data Tables 77, 1-56 (2001)
! doi:10.1006/adnd.2000.0848
! Beginning with Z=10 data from Tables IIA and IIB, Z=3,6-8 from Table I
!  WARNING no data for Z=4,Z=5,Z=9 values are made up!!!
! Z  K     L1    L2    L3    M1   M2    M3     M4    M5   N1   N2    N3   N4   N5   N6   N7


module OCEAN_corewidths
  use AI_kinds,only : DP
  implicit none
  private
  

  real(DP),parameter,dimension( 16,3:92 ) :: coreWidths = reshape( &
   [0.03_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.05_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.05_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.1_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.132_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.14_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.16_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.24_dp,0.0_dp,0.01_dp,0.01_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.28_dp,0.28_dp,0.02_dp,0.02_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.33_dp,0.46_dp,0.03_dp,0.03_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.37_dp,0.78_dp,0.04_dp,0.04_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.43_dp,0.9_dp,0.05_dp,0.05_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.47_dp,1.1_dp,0.07_dp,0.07_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.52_dp,1.3_dp,0.09_dp,0.09_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.57_dp,1.5_dp,0.11_dp,0.11_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.66_dp,1.8_dp,0.13_dp,0.13_dp,0.14_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.71_dp,2.1_dp,0.18_dp,0.18_dp,0.7_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.77_dp,2.5_dp,0.21_dp,0.21_dp,1.1_dp,1.2_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.83_dp,3.3_dp,0.36_dp,0.23_dp,1.7_dp,1.2_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.89_dp,3.9_dp,0.52_dp,0.25_dp,2.1_dp,1.2_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    0.96_dp,4.6_dp,0.78_dp,0.28_dp,2.2_dp,1.2_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.02_dp,5.2_dp,0.76_dp,0.32_dp,2.3_dp,1.2_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.11_dp,6.2_dp,0.97_dp,0.36_dp,2.4_dp,1.2_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.19_dp,7.0_dp,1.14_dp,0.41_dp,2.4_dp,1.23_dp,1.23_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.28_dp,7.2_dp,1.13_dp,0.47_dp,2.4_dp,1.25_dp,1.27_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.39_dp,6.4_dp,0.98_dp,0.53_dp,2.3_dp,1.3_dp,1.3_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.49_dp,5.5_dp,1.04_dp,0.61_dp,2.2_dp,1.9_dp,1.8_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.62_dp,4.8_dp,1.06_dp,0.68_dp,2.1_dp,2.1_dp,2.15_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.76_dp,4.1_dp,0.77_dp,0.77_dp,2.0_dp,2.25_dp,2.3_dp,0.012_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    1.92_dp,3.8_dp,0.86_dp,0.86_dp,2.1_dp,2.3_dp,2.3_dp,0.045_dp,0.044_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    2.09_dp,3.8_dp,0.95_dp,0.94_dp,2.4_dp,2.25_dp,2.25_dp,0.058_dp,0.06_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    2.28_dp,3.8_dp,1.05_dp,1.02_dp,2.8_dp,2.2_dp,2.2_dp,0.065_dp,0.066_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    2.49_dp,3.8_dp,1.14_dp,1.11_dp,3.2_dp,2.1_dp,2.15_dp,0.068_dp,0.07_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    2.71_dp,3.75_dp,1.25_dp,1.19_dp,3.5_dp,1.6_dp,1.1_dp,0.07_dp,0.072_dp,0.4_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    2.96_dp,3.75_dp,1.34_dp,1.27_dp,4.0_dp,1.9_dp,1.95_dp,0.067_dp,0.069_dp,1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    3.23_dp,3.75_dp,1.43_dp,1.35_dp,4.4_dp,1.9_dp,1.9_dp,0.061_dp,0.064_dp,1.6_dp,0.4_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    3.51_dp,3.75_dp,1.53_dp,1.43_dp,4.9_dp,1.95_dp,1.95_dp,0.062_dp,0.066_dp,2.0_dp,0.8_dp,0.3_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    3.83_dp,3.75_dp,1.63_dp,1.51_dp,5.4_dp,2.0_dp,2.0_dp,0.07_dp,0.074_dp,2.4_dp,1.2_dp,0.6_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    4.16_dp,3.8_dp,1.73_dp,1.6_dp,5.8_dp,2.05_dp,2.05_dp,0.092_dp,0.095_dp,2.8_dp,1.5_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    4.52_dp,3.8_dp,1.83_dp,1.69_dp,6.3_dp,2.1_dp,2.1_dp,0.22_dp,0.12_dp,3.2_dp,2.2_dp,1.6_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    4.91_dp,3.8_dp,1.93_dp,1.78_dp,6.7_dp,2.15_dp,2.15_dp,0.5_dp,0.14_dp,3.5_dp,2.7_dp,2.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    5.33_dp,3.9_dp,2.03_dp,1.87_dp,7.2_dp,2.2_dp,2.2_dp,0.59_dp,0.17_dp,3.9_dp,3.2_dp,2.8_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    5.77_dp,4.0_dp,2.13_dp,1.96_dp,7.6_dp,2.25_dp,2.25_dp,0.61_dp,0.21_dp,4.2_dp,4.2_dp,3.8_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    6.25_dp,3.9_dp,2.23_dp,2.05_dp,8.0_dp,2.35_dp,2.35_dp,0.26_dp,0.26_dp,4.35_dp,6.4_dp,5.6_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    6.76_dp,3.8_dp,2.32_dp,2.15_dp,8.4_dp,2.45_dp,2.55_dp,0.3_dp,0.31_dp,4.4_dp,8.4_dp,8.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    7.32_dp,3.5_dp,2.42_dp,2.24_dp,8.8_dp,2.55_dp,2.8_dp,0.34_dp,0.35_dp,4.4_dp,10.8_dp,10.5_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    7.9_dp,3.0_dp,2.53_dp,2.34_dp,9.2_dp,2.7_dp,3.05_dp,0.38_dp,0.39_dp,4.2_dp,13.2_dp,14.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
    8.53_dp,2.4_dp,2.64_dp,2.43_dp,9.6_dp,2.85_dp,3.3_dp,0.43_dp,0.44_dp,3.4_dp,17.0_dp,17.0_dp,0.08_dp,0.08_dp,0.0_dp,0.0_dp,&
    9.2_dp,2.3_dp,2.74_dp,2.53_dp,10.0_dp,3.0_dp,3.6_dp,0.47_dp,0.48_dp,2.6_dp,0.0_dp,0.0_dp,0.14_dp,0.14_dp,0.0_dp,0.0_dp,&
    9.91_dp,2.2_dp,2.84_dp,2.62_dp,10.2_dp,3.2_dp,3.9_dp,0.52_dp,0.52_dp,2.4_dp,0.0_dp,0.0_dp,0.17_dp,0.2_dp,0.0_dp,0.0_dp,&
    10.7_dp,2.1_dp,2.95_dp,2.72_dp,10.4_dp,3.35_dp,4.3_dp,0.56_dp,0.56_dp,2.4_dp,0.0_dp,0.0_dp,0.11_dp,0.12_dp,0.0_dp,0.0_dp,&
    11.5_dp,2.0_dp,3.05_dp,2.82_dp,10.6_dp,3.5_dp,4.7_dp,0.6_dp,0.6_dp,2.6_dp,0.0_dp,0.0_dp,0.1_dp,0.08_dp,0.0_dp,0.0_dp,&
    12.3_dp,2.0_dp,3.15_dp,2.92_dp,10.8_dp,3.7_dp,5.0_dp,0.63_dp,0.63_dp,2.8_dp,0.0_dp,0.0_dp,0.08_dp,0.08_dp,0.0_dp,0.0_dp,&
    13.2_dp,2.1_dp,3.25_dp,3.02_dp,11.1_dp,3.9_dp,5.4_dp,0.67_dp,0.67_dp,3.1_dp,5.0_dp,1.3_dp,0.08_dp,0.08_dp,0.0_dp,0.0_dp,&
    14.2_dp,2.2_dp,3.35_dp,3.12_dp,11.4_dp,4.1_dp,5.8_dp,0.7_dp,0.7_dp,3.3_dp,5.03_dp,1.45_dp,0.09_dp,0.1_dp,0.0_dp,0.0_dp,&
    15.2_dp,2.5_dp,3.41_dp,3.19_dp,11.6_dp,4.3_dp,6.2_dp,0.72_dp,0.72_dp,3.5_dp,5.06_dp,1.6_dp,0.61_dp,0.32_dp,0.0_dp,0.0_dp,&
    16.2_dp,2.7_dp,3.48_dp,3.27_dp,11.8_dp,4.5_dp,6.7_dp,0.75_dp,0.75_dp,3.7_dp,5.08_dp,1.75_dp,0.78_dp,0.53_dp,0.0_dp,0.0_dp,&
    17.4_dp,2.9_dp,3.55_dp,3.36_dp,12.0_dp,4.7_dp,7.3_dp,0.78_dp,0.78_dp,4.0_dp,5.1_dp,1.9_dp,1.05_dp,0.8_dp,0.0_dp,0.0_dp,&
    18.5_dp,3.1_dp,3.63_dp,3.44_dp,12.2_dp,5.0_dp,7.8_dp,0.82_dp,0.82_dp,4.2_dp,5.13_dp,2.05_dp,1.38_dp,1.11_dp,0.0_dp,0.0_dp,&
    19.8_dp,3.3_dp,3.7_dp,3.53_dp,12.4_dp,5.2_dp,8.1_dp,0.86_dp,0.86_dp,4.4_dp,5.16_dp,2.2_dp,1.78_dp,1.48_dp,0.0_dp,0.0_dp,&
    21.1_dp,3.6_dp,3.77_dp,3.62_dp,12.6_dp,5.4_dp,8.2_dp,0.9_dp,0.9_dp,4.6_dp,5.2_dp,2.35_dp,2.2_dp,1.9_dp,0.0_dp,0.0_dp,&
    22.4_dp,3.8_dp,3.87_dp,3.72_dp,12.8_dp,5.6_dp,8.3_dp,0.95_dp,0.95_dp,4.9_dp,5.23_dp,2.5_dp,2.45_dp,2.2_dp,0.0_dp,0.0_dp,&
    23.8_dp,4.0_dp,3.93_dp,3.8_dp,13.0_dp,5.8_dp,8.2_dp,1.01_dp,1.01_dp,5.1_dp,5.26_dp,2.65_dp,2.7_dp,2.4_dp,0.0_dp,0.0_dp,&
    25.3_dp,4.3_dp,4.01_dp,3.9_dp,13.2_dp,6.0_dp,8.0_dp,1.065_dp,1.065_dp,5.4_dp,5.3_dp,2.8_dp,2.95_dp,2.6_dp,0.0_dp,0.0_dp,&
    26.9_dp,4.5_dp,4.09_dp,4.0_dp,13.4_dp,6.3_dp,7.8_dp,1.13_dp,1.13_dp,5.6_dp,5.33_dp,2.95_dp,3.15_dp,2.8_dp,0.0_dp,0.0_dp,&
    28.5_dp,4.7_dp,4.18_dp,4.1_dp,13.6_dp,6.6_dp,7.5_dp,1.2_dp,1.2_dp,5.8_dp,5.36_dp,3.15_dp,3.35_dp,2.95_dp,0.0_dp,0.0_dp,&
    30.2_dp,4.9_dp,4.26_dp,4.2_dp,13.8_dp,6.85_dp,7.1_dp,1.27_dp,1.27_dp,6.1_dp,5.4_dp,3.3_dp,3.55_dp,3.1_dp,0.0_dp,0.0_dp,&
    32.0_dp,5.2_dp,4.36_dp,4.31_dp,13.9_dp,7.1_dp,6.7_dp,1.35_dp,1.35_dp,6.3_dp,5.5_dp,3.5_dp,3.7_dp,3.2_dp,0.03_dp,0.03_dp,&
    33.9_dp,5.4_dp,4.46_dp,4.43_dp,14.1_dp,7.3_dp,6.0_dp,1.43_dp,1.43_dp,6.6_dp,5.5_dp,3.65_dp,3.8_dp,3.3_dp,0.03_dp,0.03_dp,&
    35.9_dp,5.7_dp,4.57_dp,4.55_dp,14.2_dp,7.5_dp,5.6_dp,1.52_dp,1.52_dp,6.8_dp,5.6_dp,3.85_dp,3.9_dp,3.5_dp,0.07_dp,0.03_dp,&
    37.9_dp,6.0_dp,4.69_dp,4.68_dp,14.3_dp,7.8_dp,5.7_dp,1.61_dp,1.61_dp,7.0_dp,5.7_dp,4.0_dp,4.0_dp,3.65_dp,0.08_dp,0.04_dp,&
    40.1_dp,6.3_dp,4.82_dp,4.81_dp,14.5_dp,8.1_dp,6.4_dp,1.7_dp,1.7_dp,7.3_dp,5.8_dp,4.2_dp,4.1_dp,3.8_dp,0.1_dp,0.06_dp,&
    42.2_dp,6.7_dp,4.95_dp,4.95_dp,14.6_dp,8.4_dp,6.9_dp,1.79_dp,1.79_dp,7.5_dp,5.9_dp,4.4_dp,4.1_dp,3.9_dp,0.15_dp,0.11_dp,&
    44.6_dp,7.2_dp,5.09_dp,5.09_dp,14.7_dp,8.6_dp,7.5_dp,1.89_dp,1.89_dp,7.7_dp,6.0_dp,4.6_dp,4.1_dp,3.9_dp,0.22_dp,0.18_dp,&
    47.0_dp,7.9_dp,5.23_dp,5.24_dp,14.8_dp,8.9_dp,8.0_dp,1.99_dp,1.99_dp,8.0_dp,6.1_dp,4.75_dp,4.1_dp,4.0_dp,0.31_dp,0.27_dp,&
    49.5_dp,8.8_dp,5.38_dp,5.39_dp,14.9_dp,9.2_dp,8.3_dp,2.08_dp,2.08_dp,8.25_dp,6.25_dp,4.9_dp,4.1_dp,3.95_dp,0.35_dp,0.31_dp,&
    52.1_dp,9.8_dp,5.53_dp,5.54_dp,15.0_dp,9.5_dp,8.5_dp,2.18_dp,2.18_dp,8.5_dp,6.4_dp,5.05_dp,4.1_dp,3.9_dp,0.37_dp,0.33_dp,&
    54.8_dp,10.5_dp,5.69_dp,5.71_dp,15.1_dp,9.8_dp,8.6_dp,2.28_dp,2.28_dp,8.8_dp,6.55_dp,5.3_dp,4.0_dp,3.85_dp,0.33_dp,0.31_dp,&
    57.6_dp,11.1_dp,5.87_dp,5.89_dp,15.1_dp,10.1_dp,8.7_dp,2.38_dp,2.38_dp,9.1_dp,6.7_dp,5.6_dp,3.9_dp,3.8_dp,0.29_dp,0.27_dp,&
    60.6_dp,11.8_dp,5.04_dp,6.07_dp,15.2_dp,10.4_dp,8.7_dp,2.48_dp,2.48_dp,9.35_dp,6.9_dp,5.8_dp,3.8_dp,3.8_dp,0.26_dp,0.23_dp,&
    63.6_dp,12.3_dp,6.22_dp,6.27_dp,15.2_dp,10.7_dp,8.6_dp,2.58_dp,2.58_dp,9.6_dp,7.2_dp,5.95_dp,3.8_dp,3.8_dp,0.22_dp,0.2_dp,&
    66.8_dp,12.7_dp,6.41_dp,6.46_dp,15.3_dp,11.1_dp,8.5_dp,2.68_dp,2.68_dp,9.9_dp,7.35_dp,6.2_dp,3.9_dp,3.8_dp,0.19_dp,0.18_dp,&
    70.0_dp,13.0_dp,6.6_dp,6.66_dp,15.3_dp,11.4_dp,8.4_dp,2.78_dp,2.78_dp,10.1_dp,7.6_dp,6.4_dp,3.9_dp,3.85_dp,0.17_dp,0.16_dp,&
    73.4_dp,13.2_dp,6.81_dp,6.87_dp,15.3_dp,11.7_dp,8.3_dp,2.88_dp,2.88_dp,10.4_dp,7.8_dp,6.6_dp,4.0_dp,3.9_dp,0.16_dp,0.16_dp,&
    76.9_dp,13.5_dp,7.02_dp,7.08_dp,15.4_dp,12.1_dp,8.2_dp,2.98_dp,2.98_dp,10.7_dp,8.1_dp,6.8_dp,4.1_dp,3.95_dp,0.15_dp,0.15_dp,&
    80.6_dp,13.7_dp,7.5_dp,7.29_dp,15.4_dp,12.5_dp,8.2_dp,3.08_dp,3.08_dp,10.95_dp,8.3_dp,7.0_dp,4.15_dp,4.0_dp,0.15_dp,0.16_dp,&
    84.4_dp,14.0_dp,8.0_dp,7.51_dp,15.4_dp,12.9_dp,8.0_dp,3.18_dp,3.18_dp,11.2_dp,8.5_dp,7.25_dp,4.2_dp,4.05_dp,0.15_dp,0.17_dp,&
    88.2_dp,14.3_dp,8.5_dp,7.74_dp,15.5_dp,13.2_dp,8.0_dp,3.28_dp,3.28_dp,11.5_dp,8.75_dp,7.5_dp,4.3_dp,4.1_dp,0.15_dp,0.18_dp,&
    92.1_dp,14.7_dp,9.1_dp,7.97_dp,15.5_dp,13.6_dp,7.9_dp,3.39_dp,3.39_dp,11.6_dp,9.2_dp,7.75_dp,4.4_dp,4.2_dp,0.29_dp,0.25_dp,&
    96.3_dp,16.0_dp,10.0_dp,8.2_dp,15.5_dp,14.1_dp,7.9_dp,3.5_dp,3.5_dp,12.2_dp,9.6_dp,8.0_dp,4.5_dp,4.25_dp,0.37_dp,0.31_dp ]&
    , [16,90] )


  public :: returnLifetime

  contains

  subroutine returnLifetime( ZZ, NN, LL, lifetime )
    integer,intent( in ) :: ZZ, NN, LL
    real(DP),intent(out) :: lifetime
    integer :: indx
    
    lifetime = -1
    if( ZZ .lt. 3 ) then
      lifetime = coreWidths( 1, 3 )
      return
    endif

    if( NN .gt. 4 ) return
    if( LL .ge. NN ) return

    indx = 1
    if( NN .eq. 2 ) indx = 2
    if( NN .eq. 3 ) indx = 5
    if( NN .eq. 4 ) indx = 10
  
    if( LL .eq. 1 ) indx = indx + 2
    if( LL .eq. 2 ) indx = indx + 4
    if( LL .eq. 3 ) indx = indx + 6

    lifetime = coreWidths( indx, ZZ )

  end subroutine returnLifetime

end module OCEAN_corewidths
