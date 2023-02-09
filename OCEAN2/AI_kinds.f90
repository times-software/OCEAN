! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module AI_kinds
#ifdef __HAVE_F03
  use iso_c_binding, only : c_sizeof
#endif
  implicit none
  ! 
  ! The standard numbers should be compatible with MPI and scaLAPACK
  ! Will need some sort of pre-build thingy at some point
  !
  ! SP = MPI_REAL & S for scalapck
  ! DP = MPI_DOUBLE_PRECISION & D for scalapack
  ! QP is unused
!  integer, parameter :: SP = kind(1.0)
  integer, parameter :: SP = selected_real_kind(6,37)
  integer, parameter :: DP = selected_real_kind(15,307)
#ifdef __HAVE_QP
  integer, parameter :: QP = selected_real_kind(30)
#else
  integer, parameter :: QP = DP
#endif
! Running under the assumption that complex numbers can be declared using 
!  the corresponding real definitions.
!  integer, parameter :: CDP = DP
!  integer, parameter :: CSP = SP
!  integer, parameter :: I16 = selected_int_kind(19)
  !  
  ! Standard int is MPI_INTEGER & used for scalapck bounds, ect.
  integer, parameter :: S_INT = kind( 1 )
!  integer, parameter :: I16 = selected_int_kind(32)
  integer, parameter :: I8  = selected_int_kind(18)
  integer, parameter :: I4  = selected_int_kind(9)
  integer, parameter :: I2  = selected_int_kind(4)
  integer, parameter :: I1  = selected_int_kind(1)
  !
!  integer, parameter :: AI_LOG = MPI_LOGICAL
  public :: QP, DP, SP, I4, I8, I2, S_INT !, AI_LOG
  !
  !
contains
  subroutine print_kind_info
    implicit none
#ifdef __HAVE_F03
    write(6,*) 'QP', QP, sizeof(1.0_QP)
    write(6,*) 'DP', DP, c_sizeof(1.0_DP)
    write(6,*) 'SP', SP, c_sizeof(1.0_SP)
    write(6,*) 'I8', I8, c_sizeof(1_I8)
    write(6,*) 'I4', I4, c_sizeof(1_I4)
    write(6,*) 'I2', I2, c_sizeof(1_I2)
    write(6,*) 'S ', S_INT, c_sizeof(1_S_INT)
#else
    write(6,*) 'QP', QP
    write(6,*) 'DP', DP, sizeof(1.0_DP)
    write(6,*) 'SP', SP, sizeof(1.0_SP)
    write(6,*) 'I8', I8, sizeof(1_I8)
    write(6,*) 'I4', I4, sizeof(1_I4)
    write(6,*) 'I2', I2, sizeof(1_I2)  
    write(6,*) 'S ', S_INT, sizeof(1_S_INT)
#endif
  end subroutine print_kind_info

end module AI_kinds
