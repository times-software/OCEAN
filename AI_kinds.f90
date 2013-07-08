module AI_kinds
  use mpi
  implicit none
  ! 
  ! The standard numbers should be compatible with MPI and scaLAPACK
  ! Will need some sort of pre-build thingy at some point
  !
  ! QP is unused
  ! DP = MPI_DOUBLE_PRECISION & D for scalapack
  ! SP = MPI_REAL & S for scalapck
  integer, parameter :: QP = selected_real_kind(33,4900)
  integer, parameter :: DP = selected_real_kind(15,307)
  integer, parameter :: SP = selected_real_kind(6,37)
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
  integer, parameter :: AI_LOG = MPI_LOGICAL
  public :: QP, DP, SP, I4, I8, I2, S_INT, AI_LOG
  !
contains
  subroutine print_kind_info
    implicit none
!    write(6,*) 'QP', QP, sizeof(1.0_QP)
    write(6,*) 'DP', DP, sizeof(1.0_DP)
    write(6,*) 'SP', SP, sizeof(1.0_SP)
    write(6,*) 'I8', I8, sizeof(1_I8)
    write(6,*) 'I4', I4, sizeof(1_I4)
    write(6,*) 'I2', I2, sizeof(1_I2)  
    write(6,*) 'S ', S_INT, sizeof(1_S_INT)
  end subroutine print_kind_info

end module AI_kinds
