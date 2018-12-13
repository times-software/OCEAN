module ocean_interpolate
  use ai_kinds, only : DP

  implicit none
  private


  public :: makeLagrange, makeSimpleLagrange

  contains

  function evalLagrange( order, x, p )
    integer, intent( in ) :: order
    real(dp), intent( in ) :: x
    complex(dp), intent( in ) :: p(:)
    complex(dp) :: evalLagrange
    !
    select case( order )
      case( 4 )
        call evalLagrange4( x, p )
      case default
        stop

    end select
  end function evalLagrange

  function evalLagrange4( x, p )
    real(dp), intent( in ) :: x
    complex(dp), intent( in ) :: p(:)
    complex(dp) :: evalLagrange4
    !
    evalLagrange4 = p(1) + x * p(2) + x*x*p(3) + x*x*x*p(4)
  end function evalLagrange4

  ! This allows 3d data to be passed in with periodic BC
  subroutine makeLagrange( order, ix, iy, iz, d, p )
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

    call makeSimpleLagrange( order, contiguousD, p )
  end subroutine makeLagrange


  subroutine makeSimpleLagrange( order, d, p )
    integer, intent( in ) :: order
    complex(dp), intent( in ) :: d(:)
    complex(dp), intent( out ) :: p(:)

    select case ( order )
      case( 4 )
        call SL4( d, p )
      case default
        stop
    end select

  end subroutine


  subroutine SL4( d, p )
    complex(dp), intent( in ) :: d(4)
    complex(dp), intent( out ) :: p(4)

    real(dp), parameter :: oneThird = 1.0_dp/3.0_dp
    real(dp), parameter :: oneSixth = 1.0_dp/6.0_dp

    p(1) = d(2)
    p(2) = -oneThird * d(1) - 0.5_dp * d(2) - oneSixth * d(3)
    p(3) = 0.5_dp * d(1) - 0.5_dp * d(3)
    p(4) = - oneSixth * d(1) + 0.5_dp * d(2) - 0.5_dp * d(3) + oneSixth * d(4)
  
  end subroutine SL4



end module ocean_interpolate
