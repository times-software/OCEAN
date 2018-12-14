! Copyright (C) 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
!
!
! This module provides for easy access to Lagrange polynomial interpolation
! for uniform grids. 
! By assumption the distances x are rescaled by the spacing between points
! can also pass in a rescaling

! Might want to use polymorphism instead of select to set order
! would be set by the length of the polynomial array?
module ocean_interpolate
  use ai_kinds, only : DP

  implicit none
  private


  public :: makeLagrange, makeSimpleLagrange, evalLagrange

  interface evalLagrange
    module procedure rEvalLagrange, cEvalLagrange, rEvalLagrangeScale, cEvalLagrangeScale
  end interface evalLagrange

  interface makeLagrange
    module procedure rMakeLagrange, cMakeLagrange
  end interface makeLagrange

  interface makeSimpleLagrange
    module procedure rMakeSimpleLagrange, cMakeSimpleLagrange
  end interface makeSimpleLagrange

  contains

  function cEvalLagrangeScale( order, x, delta, p )
    integer, intent( in ) :: order
    real(dp), intent( in ) :: x
    real(dp), intent( in ) :: delta
    complex(dp), intent( in ) :: p(:)
    complex(dp) :: cEvalLagrangeScale
  
    real(dp) :: y

    y = x / delta
    cEvalLagrangeScale = cEvalLagrange( order, y, p )
  end function cEvalLagrangeScale
    
  function rEvalLagrangeScale( order, x, delta, p )
    integer, intent( in ) :: order
    real(dp), intent( in ) :: x
    real(dp), intent( in ) :: delta
    real(dp), intent( in ) :: p(:)
    real(dp) :: rEvalLagrangeScale

    real(dp) :: y

    y = x / delta
    rEvalLagrangeScale = rEvalLagrange( order, y, p )
  end function rEvalLagrangeScale

  function cEvalLagrange( order, x, p )
    integer, intent( in ) :: order
    real(dp), intent( in ) :: x
    complex(dp), intent( in ) :: p(:)
    complex(dp) :: cEvalLagrange
    !

    select case( order )
      case( 4 )
        cEvalLagrange = cEvalLagrange4( x, p )
      case default
        stop

    end select
  end function cEvalLagrange

  function rEvalLagrange( order, x, p )
    integer, intent( in ) :: order
    real(dp), intent( in ) :: x
    real(dp), intent( in ) :: p(:)
    real(dp) :: rEvalLagrange
    !

    select case( order )
      case( 4 )
        rEvalLagrange = rEvalLagrange4( x, p )
      case default
        stop

    end select
  end function rEvalLagrange

  function cEvalLagrange4( x, p )
    real(dp), intent( in ) :: x
    complex(dp), intent( in ) :: p(:)
    complex(dp) :: cEvalLagrange4
    !
    cEvalLagrange4 = p(1) + x * p(2) + x*x*p(3) + x*x*x*p(4)
  end function cEvalLagrange4

  function rEvalLagrange4( x, p )
    real(dp), intent( in ) :: x
    real(dp), intent( in ) :: p(:)
    real(dp) :: rEvalLagrange4
    !
    rEvalLagrange4 = p(1) + x * p(2) + x*x*p(3) + x*x*x*p(4)
  end function rEvalLagrange4

  ! This allows 3d data to be passed in with periodic BC
  subroutine cMakeLagrange( order, ix, iy, iz, d, p )
    integer, intent( in ) :: order, ix, iy,iz
    complex(dp), intent( in ) :: d(:,:,:)
    complex(dp), intent( out ) :: p(:)

    complex(dp) :: contiguousD( order )
    integer :: iix, iiy, iiz, i
    integer :: nx, ny, nz

    nx = size( d, 1 )
    ny = size( d, 2 )
    nz = size( d, 3 )

    iiy = iy
    do while( iiy .lt. 1 )
      iiy = iiy + ny
    enddo
    do while( iiy .gt. ny )
      iiy = iiy - ny
    enddo
    
    iiz = iz
    do while( iiz .lt. 1 ) 
      iiz = iiz + nz
    enddo
    do while( iiz .gt. nz )
      iiz = iiz - nz
    enddo

    iix = ix
    do while( iix .lt. 1 ) 
      iix = iix + nx 
    enddo
    do i = 1, order 
      do while( iix .gt. nx )
        iix = iix - nz
      enddo
      contiguousD( i ) = d(iix,iiy,iiz)
      iix = iix + 1
    enddo

    call MakeSimpleLagrange( order, contiguousD, p )
  end subroutine cMakeLagrange

  ! This allows 3d data to be passed in with periodic BC
  subroutine rMakeLagrange( order, ix, iy, iz, d, p )
    integer, intent( in ) :: order, ix, iy,iz
    real(dp), intent( in ) :: d(:,:,:)
    real(dp), intent( out ) :: p(:)

    real(dp) :: contiguousD( order )
    integer :: iix, iiy, iiz, i
    integer :: nx, ny, nz

    nx = size( d, 1 )
    ny = size( d, 2 )
    nz = size( d, 3 )

    iiy = iy
    do while( iiy .lt. 1 )
      iiy = iiy + ny
    enddo
    do while( iiy .gt. ny )
      iiy = iiy - ny
    enddo

    iiz = iz
    do while( iiz .lt. 1 )
      iiz = iiz + nz
    enddo
    do while( iiz .gt. nz )
      iiz = iiz - nz
    enddo

    iix = ix
    do while( iix .lt. 1 )
      iix = iix + nx
    enddo
    do i = 1, order
      do while( iix .gt. nx )
        iix = iix - nz
      enddo
      contiguousD( i ) = d(iix,iiy,iiz)
      iix = iix + 1
    enddo

    call MakeSimpleLagrange( order, contiguousD, p )
  end subroutine rMakeLagrange

  subroutine cMakeSimpleLagrange( order, d, p )
    integer, intent( in ) :: order
    complex(dp), intent( in ) :: d(:)
    complex(dp), intent( out ) :: p(:)

    select case ( order )
      case( 4 )
        call cSL4( d, p )
      case default
        stop
    end select

  end subroutine cMakeSimpleLagrange

  subroutine rMakeSimpleLagrange( order, d, p )
    integer, intent( in ) :: order
    real(dp), intent( in ) :: d(:)
    real(dp), intent( out ) :: p(:)

    select case ( order )
      case( 4 )
        call rSL4( d, p )
      case default
        stop
    end select

  end subroutine rMakeSimpleLagrange


  subroutine cSL4( d, p )
    complex(dp), intent( in ) :: d(4)
    complex(dp), intent( out ) :: p(4)

    real(dp), parameter :: oneThird = 1.0_dp/3.0_dp
    real(dp), parameter :: oneSixth = 1.0_dp/6.0_dp

    p(1) = d(2)
    p(2) = -oneThird * d(1) - 0.5_dp * d(2) - oneSixth * d(3)
    p(3) = 0.5_dp * d(1) - 0.5_dp * d(3)
    p(4) = - oneSixth * d(1) + 0.5_dp * d(2) - 0.5_dp * d(3) + oneSixth * d(4)
  
  end subroutine cSL4

  subroutine rSL4( d, p )
    real(dp), intent( in ) :: d(4)
    real(dp), intent( out ) :: p(4)

    real(dp), parameter :: oneThird = 1.0_dp/3.0_dp
    real(dp), parameter :: oneSixth = 1.0_dp/6.0_dp

    p(1) = d(2)
    p(2) = -oneThird * d(1) - 0.5_dp * d(2) - oneSixth * d(3)
    p(3) = 0.5_dp * d(1) - 0.5_dp * d(3)
    p(4) = - oneSixth * d(1) + 0.5_dp * d(2) - 0.5_dp * d(3) + oneSixth * d(4)

  end subroutine rSL4



end module ocean_interpolate
