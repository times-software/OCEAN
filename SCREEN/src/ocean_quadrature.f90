module ocean_quadrature
  use ai_kinds, only : DP 

  implicit none
  private
  save

!  type quadHolder
!    real(WP), allocatable, dimension( : ) :: abscissae
!    real(WP), allocatable, dimension( : ) :: weights
!  end type quadHolder

!  type( quadHolder ), allocatable, protected :: LegendreQuad( : )

!  public :: quadHolder
  public :: ocean_quadrature_loadLegendre

  contains

  subroutine ocean_quadrature_loadLegendre( order, abscissae, weights, ierr )
    integer, intent( in ) :: order
    real(DP), intent( out ) :: abscissae( : )
    real(DP), intent( out ) :: weights( : )
    integer, intent( inout ) :: ierr

    integer :: i, j

    if( mod( order, 2 ) .ne. 0 ) then
      ierr = 1
      return
    endif

    open( unit=99, file='EvenQuadHalf.txt', form='formatted', status='old', iostat=ierr, ERR=10 )
    do
      read(99,*) i
      if( i .eq. order ) then
        do j = 1, i/2
          read(99,*) abscissae(j), weights(j)
          abscissae(i-j+1) = -abscissae(j)
          weights(i-j+1) = weights(j)
        enddo
        
        goto 11
      else
        do j = 1, i/2
          read(99,*)
        enddo
      endif
    enddo
    
11 continue
    close( 99 )

10  continue
    
  end subroutine ocean_quadrature_loadLegendre

end module
