! Copyright (C) 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module screen_kxc
  use ai_kinds, only : DP

  implicit none
  private
  save

  public :: dftder3



! Longer-term plan, move everything in here
  ! 0. Each of these has a check to see if we are using Kxc
  ! 1. The density should be read-in ( by one processor ) and stored on all in real-space uniform mesh
  !   probably just kept for pinfo%myid == pinfo%root
  ! 2. For each site in screen_chi_driver_run, then, if pinfo%myid == pinfo%root
  !   A. project density from regular real-space to grid
  !   B. calculate Kxc (overload to handle either diagonal or full-real-space matrix versions


  contains


  ! Original from Eric Shirley
  subroutine dftder3( n, nexc, vxc, kxc, fxc )
    implicit none
    !     
    real( DP ), intent( in ) :: n 
    real( DP ), intent( out ) :: nexc, vxc, kxc, fxc
    !
    integer :: i
    real( DP ) :: dn, n1, ex, ec, xctab( -3 : 3 ), ux1, ux2, uc1, uc2, d1, d2, d3
    !
    real( DP ) :: ecvwn, nec
    !
    nexc = 0
    vxc = 0
    kxc = 0
    fxc = 0
    if ( n .gt. 0.0d0 ) then
       dn = 0.01d0 * n
       do i = -3, 3
          n1 = n + dn * dble( i )
          call cacorr( n1, ex, ec, ux1, ux2, uc1, uc2 )
          call getc2( 0.5d0 * n1, 0.5d0 * n1, ecvwn, nec, uc1, uc2 )
  !       xctab( i ) = n1 * ( ex + ec )
          xctab( i ) = n1 * ( ex + ecvwn )
       end do
       nexc = xctab( 0 )
       d1 = ( xctab( 1 ) - xctab( -1 ) ) / ( 2.0d0 * dn )
       d2 = ( xctab( 2 ) - xctab( -2 ) ) / ( 4.0d0 * dn )
       vxc = ( 4.0d0 * d1 - d2 ) / 3.d0
       d3 = ( xctab( 3 ) - xctab( -3 ) ) / ( 6.0d0 * dn )
       fxc = ( 16.0d0 * d2 - 13.0d0 * d1 - 3.0d0 * d3 ) / ( 4.0d0 * dn ** 2 )
       d1 = ( xctab( 1 ) + xctab( -1 ) - 2.0d0 * xctab( 0 ) ) / dn ** 2
       d2 = ( xctab( 2 ) + xctab( -2 ) - 2.0d0 * xctab( 0 ) ) / ( 2.0d0 * dn ) ** 2
       kxc = ( 4.0d0 * d1 - d2 ) / 3.0d0
    end if
    !
    return
  end subroutine dftder3

  

end module screen_kxc
