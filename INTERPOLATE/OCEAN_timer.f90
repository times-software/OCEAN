module OCEAN_timer
  private
  save

  integer(8) :: start_time
  integer(8) :: max_time
  integer(8) :: tics_per_second

  public :: OCEAN_t_reset, OCEAN_t_time, OCEAN_t_printtime

  contains

  subroutine OCEAN_t_reset
    implicit none

    call system_clock( start_time, tics_per_second, max_time )
  end subroutine OCEAN_t_reset

  subroutine OCEAN_t_printtime( stringname, fh )
    implicit none
    character(*), intent(in)  :: stringname
    integer, intent(in ) :: fh
    !
    integer(8) :: tics
    real(kind=kind(1.d0)) :: secs

    call OCEAN_t_time( tics, secs )
    write(fh,'(a,a,i19,a,f15.7,a)') stringname, ': ', tics, ' tics, ', secs, ' secs'

    
  end subroutine OCEAN_t_printtime

  subroutine OCEAN_t_time( tics, secs )
    implicit none
    integer(8), intent(out) :: tics
    real(kind=kind(1.d0)), intent(out) :: secs
    !
    integer(8) :: tmp_time

    call system_clock( tmp_time )

    if( tmp_time .lt. start_time ) then
    ! roll over
      tics = tmp_time + ( max_time - start_time )
    else
      tics = tmp_time - start_time
    endif

    secs = dble( tics ) / dble( tics_per_second )
  end subroutine OCEAN_t_time

end module OCEAN_timer
