
module ocean_sphericalHarmonics
  use ai_kinds, only : DP
  use ocean_constants, only : PI_DP

  implicit none
  private

  real(DP), parameter :: pref = 1.0_DP / sqrt( 4.0_DP * PI_DP )

  public :: ocean_sphH_getylm

  contains

  pure function ocean_sphH_getylm( rvec, l, m ) result( ylm )
    real(DP), intent( in ) :: rvec( 3 )
    integer, intent( in ) :: l, m
    real(DP) :: ylm

    real(DP) :: x, y, z, rmag

    rmag = 1.0_DP / sqrt( rvec(1)**2 + rvec(2)**2 + rvec(3)**2 )
    x = rvec(1)*rmag
    y = rvec(2)*rmag
    z = rvec(3)*rmag
    ylm = 0.0_DP

    select case ( l )

      case( 0 )
        ylm = pref

      case( 1 )

        select case ( m )
          case( 1 )
            ylm = pref * sqrt(3.0_DP) * x
!            ylm = pref * sqrt(3.0_DP) * y
          case( -1 )
            ylm = pref * sqrt(3.0_DP) * y
!            ylm = pref * sqrt(3.0_DP) * z
          case( 0 )
            ylm = pref * sqrt(3.0_DP) * z
!            ylm = pref * sqrt(3.0_DP) * x
        end select 

      case( 2 )

        select case ( m )
!          case( 0 )
          case( -2 )
            ylm = pref * sqrt(15.0_DP) * x * y
!          case( 2 )
          case( -1 )
            ylm = pref * sqrt(15.0_DP) * y * z
!          case( -2 )
          case( 0 )
            ylm = pref * sqrt(5.0_DP)/2.0_DP * ( 2.0_DP * z**2 - x**2 -y**2 )
!            ylm = pref * sqrt(5.0_DP)/2.0_DP * ( 3.0_DP * z**2 - 1.0_DP )
          case( 1 )
            ylm = pref * sqrt(15.0_DP) * z * x
!          case( -1 )
          case( 2 )
            ylm = pref * sqrt(15.0_DP)/2.0_DP * ( x**2 - y**2 )
        end select

      case( 3 )
        
        select case( m )
          case( -3 )
            ylm = pref * 0.5_dp * sqrt( 35.0_DP * 0.5_DP ) * ( 3.0_DP * x**2 - y**2 ) * y
          case( -2 )
            ylm = pref * sqrt( 105.0_DP ) * x * y * z
          case( -1 ) 
            ylm = pref * 0.5_dp * sqrt( 21.0_DP * 0.5_DP ) * y * ( 4.0_DP * z**2 - x**2 - y**2 )
          case( 0 ) 
            ylm = pref * 0.5_DP * sqrt( 7.0_DP ) * z * ( 2.0_DP*z**2 - 3.0_DP * x**2 - 3.0_DP * y**2 )
          case( 1 )
            ylm = pref * 0.5_DP * sqrt( 21.0_DP * 0.5_DP ) * x * ( 4.0_DP * z**2 - x**2 - y**2 )
          case( 2 ) 
            ylm = pref * 0.5_DP * sqrt( 105.0_DP ) * ( x**2 - y**2 ) * z
          case( 3 )
            ylm = pref * 0.5_DP * sqrt( 35.0_DP * 0.5_DP ) * ( x**2 - 3.0_DP * y**2 ) * x
        end select

    end select

  end function ocean_sphH_getylm

end module ocean_sphericalHarmonics
