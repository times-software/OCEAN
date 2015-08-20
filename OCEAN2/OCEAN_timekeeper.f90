module OCEAN_timekeeper
  implicit none
  save
  integer, private, parameter :: ndivs = 5
  integer(8), private :: max_count, count_rate
  integer, public, parameter :: tk_lr = 1, &
                                tk_mult = 2, &
                                tk_e0 = 3, &
                                tk_psisum = 4, &
                                tk_inv = 5

  character(LEN=10), private, parameter :: tk_label( ndivs ) = (/ 'long range', 'multiplet ', &
    'energies  ', 'sum vector', 'inversion ' /)

  integer(8), private :: total( ndivs )
  integer(8), private :: prev( ndivs )

  contains

  subroutine OCEAN_tk_init
    implicit none
    integer(8) :: cl
    call SYSTEM_CLOCK( cl, count_rate, max_count )
    total( : ) = 0
  end subroutine OCEAN_tk_init

  subroutine OCEAN_tk_start( id )
    implicit none
    integer, intent(in) :: id
    !
    call SYSTEM_CLOCK( prev( id ) )
    !
  end subroutine OCEAN_tk_start

  subroutine OCEAN_tk_stop( id )
    implicit none
    integer, intent(in) :: id
    integer(8) :: cl
    !
    call SYSTEM_CLOCK( cl )
    if( prev( id ) .gt. cl ) then
      total( id ) = total( id ) + ( cl - ( prev( id ) - max_count ) )
    else
      total( id ) = total( id ) + ( cl - prev( id ) )
    endif
    !
  end subroutine OCEAN_tk_stop


  subroutine OCEAN_tk_printtimes( myid )
    implicit none
    integer, intent( in ) :: myid

    character(LEN=16) :: filnam
    integer :: iter

!   This can be changed if for some reason we run on all the computers
    if( myid .gt. 999999 ) return 

    write(filnam,'(A6,I6.6,A4)' ) 'timing', myid, '.txt'
    open(unit=100+myid,file=filnam,form='formatted',status='unknown')
    rewind(100+myid)

    do iter = 1, ndivs
      write(100+myid,*) tk_label( iter ), total( iter ), ' tics', &
                (dble( total( iter ) )/dble(count_rate)), 'secs'
    enddo

    close(100+myid)

  end subroutine OCEAN_tk_printtimes


end module OCEAN_timekeeper
