! Reads through ximat and checks against a reference ximat

program builder_check
  implicit none

  integer, parameter :: DP = kind(1.0d0)
  integer :: npt, ipt, jpt
  real(DP), allocatable :: xi_ref(:,:), xi(:,:)

  real(DP) :: max_diff, diff
  real(DP), parameter :: tol = 1.0d-14



  open(unit=99,file='rbfile.bin',form='unformatted',status='old')
  rewind(99)
  read(99) npt
  close(99)
  write(6,*) npt

  allocate( xi_ref( npt, npt ), xi( npt, npt ) )

  write(6,*) 'Reading in ximat_reference'
  open(unit=99,file='ximat_reference',form='unformatted',status='old')
  rewind(99)
  do ipt = 1, npt
    read( 99 ) xi_ref( :, ipt )
  enddo
  close(99)
  
  write(6,*) 'Reading in ximat'
  open(unit=99,file='ximat',form='unformatted',status='old')
  rewind(99)
  do ipt = 1, npt
    read( 99 ) xi( :, ipt )
  enddo
  close(99)

  max_diff = 0.0_DP
  do jpt = 1, npt
    do ipt = 1, npt
      diff = abs( xi_ref( ipt, jpt ) - xi( ipt, jpt ) )
      if( diff .gt. max_diff ) max_diff = diff
      if( diff .gt. tol ) write(6,*) ipt, jpt, diff
    enddo
  enddo

  write(6,*) 'Max difference: ', max_diff

  deallocate( xi, xi_ref)

end program builder_check
