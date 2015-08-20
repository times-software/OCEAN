! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module ocean_bloch

  use AI_kinds

  implicit none

  private

  type, public :: bloch_state_holder


!    integer( S_INT ) :: num_bands
!    integer( S_INT ) :: nxpts
!    integer( S_INT ) :: nkpts

    integer( S_INT ) :: my_num_bands
    integer( S_INT ) :: my_nxpts
    integer( S_INT ) :: my_start_nx
    integer( S_INT ) :: my_nkpts


!   bstates are looped by n,k,x
!    first they are summed over n, then fft-ed over k
    complex( DP ), pointer :: bstates( :, :, : )

  end type bloch_state_holder

  public :: create_bloch_state_holder


  contains

  subroutine create_bloch_state_holder( bsh, nb, nx, nk, nproc, iproc, ierr )
    type(bloch_state_holder), intent( inout ) :: bsh
    integer( S_INT ), intent( in ) :: nb, nx, nk, nproc, iproc
    integer, intent( inout ) :: ierr
    !
    integer( S_INT ) :: nx_left, nx_start, nx_tmp

  
    bsh%my_num_bands = nb
    bsh%my_nkpts = nk
  
    
    nx_left = nx
    nx_start = 1
    do i = 0, nproc - 1
      nx_tmp = nx_left / ( nproc - i )
      nx_left = nx_left - nx_tmp
      bsh%my_nxpts = nx_tmp
      bsh%my_start_nx = nx_start 
      nx_start = nx_start + nx_tmp
    enddo



    allocate( bsh%bstates( bsh%my_num_bands, bsh%my_nkpts, bsh%my_nxpts ), STAT = ierr )
    
  end subroutine create_bloch_state_holder

end module ocean_bloch

