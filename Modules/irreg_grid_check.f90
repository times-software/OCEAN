module irreg_grid_check
  use AI_kinds
  use OCEAN_mpi
  save
  ! variable and constant declarations
  type irreg_grid
          integer num_coord
          real(DP), allocatable :: curvi_coord(:,:), shifted_curvi(:,:), weights(:)
          logical :: have_curvi, have_weights
  end type

  type(irreg_grid), public, protected :: grid

contains 

  !> Checks for a non-uniform grid and its integration weights
  !! then broadcasts the results to all processes.
  subroutine check_for_grid( ierr ) 
    implicit none
    integer, intent(inout) :: ierr
    integer :: i, j

    if (myid .eq. root ) then
            inquire(file='reduced_nonuniform.txt', exist=grid%have_curvi)
            if ( grid%have_curvi ) then
                    open(unit=99, file='reduced_nonuniform.txt', form='formatted', status='old', action='read')
                    ! read number of coordinates from first line
                    read(99, *) grid%num_coord
                    allocate( grid%curvi_coord(3, grid%num_coord) )
                    allocate( grid%shifted_curvi(3, grid%num_coord) )
                    allocate( grid%weights(grid%num_coord) )
                    do i = 1, grid%num_coord
                      read(99, *) grid%curvi_coord(1, i), grid%curvi_coord(2, i), grid%curvi_coord(3, i)
                    enddo
                    close(99)
                    grid%shifted_curvi = grid%curvi_coord
            endif
            
            inquire(file='integration_weight.txt', exist=grid%have_weights)
            if (grid%have_curvi .and. grid%have_weights) then
                    open(unit=99, file='integration_weight.txt', form='formatted', status='old', action='read')
                    do j = 1, grid%num_coord
                      read(99, *) grid%weights(j)
                    enddo
                    close(99)
            else
                    ! assume all weights = 1
                    do j = 1, grid%num_coord
                      grid%weights(j) = 1
                    enddo
            endif
    endif

#ifdef MPI
    call MPI_BCAST(grid%have_curvi, 1, MPI_LOGICAL, 0, comm, ierr)
    if (ierr /= 0) goto 111
    if (grid%have_curvi) then
            call MPI_BCAST(grid%num_coord, 1, MPI_INTEGER, 0, comm, ierr)
            if (ierr /= 0) goto 111
            ! if not root, allocate array of size curvi_coord
            if (myid /= 0) then
                    allocate( grid%curvi_coord(3, grid%num_coord) )
                    allocate( grid%shifted_curvi(3, grid%num_coord) )
                    allocate( grid%weights(grid%num_coord) )
            endif
            call MPI_BCAST(grid%curvi_coord, grid%num_coord*3, MPI_DOUBLE_PRECISION, 0, comm, ierr)
            if (ierr /= 0) goto 111
            call MPI_BCAST(grid%shifted_curvi, grid%num_coord*3, MPI_DOUBLE_PRECISION, 0, comm, ierr)
            if (ierr /= 0) goto 111
    endif
    
    call MPI_BCAST(grid%have_weights, 1, MPI_LOGICAL, 0, comm, ierr)
    if (ierr /= 0) goto 111
    if (grid%have_curvi .and. grid%have_weights) then
            call MPI_BCAST(grid%weights, grid%num_coord, MPI_DOUBLE_PRECISION, 0, comm, ierr)
            if (ierr /= 0) goto 111
    endif
#endif
111 continue
  end subroutine check_for_grid

  !> Shifts the curvilinear grid by xshift and stores the result in grid%shifted_curvi 
  subroutine do_shift(i, xshift)
    integer, intent(in) :: i
    real(DP), intent(in) :: xshift(3)
    integer :: j

    do j = 1, 3
      ! shift the coordinate
      grid%shifted_curvi(j, i) = grid%curvi_coord(j, i) - xshift(j)
      ! move it back into the unit cell
      if (grid%curvi_coord(j, i) < 0) then 
              grid%shifted_curvi(j, i) = grid%curvi_coord(j, i) + 1
      elseif (grid%curvi_coord(j, i) > 1) then
              grid%shifted_curvi(j, i) = grid%curvi_coord(j, i) - 1
      endif
    enddo

  end subroutine do_shift

  !> Deallocates every array
  subroutine irreg_grid_final(this)
    type(irreg_grid), intent(inout) :: this
    if (allocated(this%curvi_coord)) then
      deallocate(this%curvi_coord)
    endif
    if (allocated(this%shifted_curvi)) then
      deallocate(this%shifted_curvi)
    endif
    if (allocated(this%weights)) then
      deallocate(this%weights)
    endif
  end subroutine irreg_grid_final

  subroutine finalize_grid()
    call irreg_grid_final(grid)
  end subroutine finalize_grid

end module irreg_grid_check
