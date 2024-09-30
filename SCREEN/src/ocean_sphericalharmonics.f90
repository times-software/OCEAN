
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

    integer :: i
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

      case( 4 )
        select case ( m )
          case( -4 )
            ylm = pref * 1.5_dp * sqrt( 35.0_dp ) * x * y * ( x**2 - y **2 )
          case( -3 )
            ylm = pref * 1.5_dp * sqrt( 35.0_dp * 0.5_dp ) * y * z * ( 3.0_dp * x**2 - y**2 )
          case( -2 )
            ylm = pref * 1.5_dp * sqrt( 5.0_dp ) * x * y * ( 7.0_dp*z**2 - 1.0_dp )
          case( -1 )
            ylm = pref * 1.5_dp * sqrt( 2.5_dp ) * y * z * ( 7.0_dp*z**2 - 3.0_dp )
          case( 0 )
            ylm = pref * (3.0_dp / 8.0_dp ) * ( (35.0_dp * z**2 - 30.0_dp)*z**2 + 3.0_dp )
!            ylm = pref * (3.0_dp / 8.0_dp ) * ( 35.0_dp * z**4 - 30.0_dp*z**2 + 3.0_dp )
          case( 1 )
            ylm = pref * 1.5_dp * sqrt( 2.5_dp ) * x * z * ( 7.0_dp*z**2 - 3.0_dp )
          case( 2 ) 
            ylm = pref * 0.75_dp * sqrt( 5.0_dp ) * ( x**2 - y**2 ) * ( 7.0_dp*z**2 - 1.0_dp )
          case( 3 )
            ylm = pref * 1.5_dp * sqrt( 35.0_dp*0.5_dp ) * x * z * ( x**2 - 3.0_dp*y**2 )
          case( 4 )
            ylm = pref * (3.0_dp / 8.0_dp ) * sqrt( 35.0_dp ) * ( x**2* ( x**2 - 3.0_dp*y**2 ) &
                                                              - y**2 * ( 3.0_dp*x**2 - y**2 ) ) 

        end select

      case default

        ylm = pref * sqrt(real(2*l+1,DP)) * AssocLegendrePlm( l, abs(m), z )
        do i= l-abs(m)+1, l+abs(m)
          ylm = ylm * sqrt( 1.0_DP/real(i,DP) )
        enddo
        if( mod(abs(m),2) .eq. 1 ) ylm = -ylm
!        if( m .ne. 0 ) ylm = ylm * (-1.0_DP)**(m)
        if( m .lt. 0 ) then
          ylm = ylm * sqrt(2.0_DP) * sin( real(-m,DP) * atan2(y,x))
        elseif( m .gt. 0 ) then
          ylm = ylm * sqrt(2.0_DP) * cos( real(m,DP) * atan2(y,x))
        endif


        

    end select

  end function ocean_sphH_getylm

!!!!
  pure function AssocLegendrePlm( l, m, x ) result( plm )
    integer, intent( in ) ::  l, m
    real(DP), intent( in ) :: x
    real(DP) :: plm
    integer :: i, ll, absm
    real(DP) :: factorial, pll, pmm, pmmp1, sqrtOneMinusX2, pmp1, pmm1
  
    plm = 0.0_DP
    absm = abs( m )
    if( absm .gt. l ) return
    pmm = 1.0_DP
    
    ! P^{m=l}_{l} = ( -1 )^l (2l-1)!! (1-x^2)^{l/2}  !! note l/2 not 1/2
    !             = ( -1 )^m (2m-1)!! (1-x^2)^{1/2}^m
    ! sqrt( 1 - x**2 )
    sqrtOneMinusX2=sqrt((1.0_DP-x)*(1.0_DP+x)) 
    factorial = 1.0_DP
    do i = 1, absm
      pmm = -pmm * factorial * sqrtOneMinusX2
      factorial = factorial + 2.0_DP 
    enddo
   
    ! 
    
    ! (l-m+1)P^m_{l+1} = (2l+1)xP^m_l - (l+m)P^m_{l-1}
    ! for first iteration m=l so P^m_{l-1}=0
    ! in what follows, ll = l+1 (from the recursion above)
    !   (l-m+1) -> ll-m
    !   (2l+1)  -> ll-1
    !   (l+m)   -> ll-1+m
    pmm1 = 0.0_DP
    do ll = absm+1, l
      pmp1 = (real(2*ll-1,DP)*x*pmm - real(ll-1+absm,DP)*pmm1 ) / real(ll-absm,DP)
      pmm1 = pmm
      pmm = pmp1
    enddo
    plm = pmm

  return
  end function

end module ocean_sphericalHarmonics
