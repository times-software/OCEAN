! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module schi_sinqr
  use ai_kinds, only : DP

  implicit none
  private
  save

  public :: 

  contains


  subroutine schi_sinqr_printSite( grid, FullChi, ierr )
    use screen_grid, only : sgrid
    use ocean_constants, only : PI_DP
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullChi(:,:,:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: vipt( : ), basfcn(:,:), qtab(:), rhs(:)
    real(DP) :: q, pref, arg, su
    integer :: i, j ,ii
    integer :: nbasis, npt, nr

    nbasis = size( FullChi, 1 )
    npt = grid%npt
    nr = grid%nr

    write(6,*) 'mkvipt'
    allocate( vipt( npt ), basfcn( npt, nbasis ), qtab( nbasis ), rhs( nbasis ) )
    call mkvipt( npt, grid%drel, vipt )

    write(6,*) 'basis'
    do i = 1, nbasis
      q = PI_DP * real( i, DP ) / grid%rmax
      pref = 2.0_DP * PI_DP * grid%rmax / q**2
      pref = 1.0_DP / sqrt( pref )
      qtab( i ) = q

      do j = 1, npt
        arg = q * grid%drel( j )
        if( arg .gt. 0.0002_DP ) then
          basfcn( j, i ) = grid%wpt(j) * pref * sin( arg ) / arg
        else
          basfcn( j, i ) = grid%wpt(j) * pref * (1.0_DP - arg**2/4.0_DP )
        endif
      enddo
    enddo

    write(6,*) 'rhs'
    open( unit=99, file='rhs', form='formatted', status='unknown' )
    rewind 99
    do i = 1, nbasis
       su = 0
       do ii = 1, npt
          su = su + basfcn( ii, i ) * vipt( ii ) !* grid%wpt( ii )
       end do
       rhs( i ) = su
       write ( 99, '(1i5,2(1x,1e15.8))' ) i, qtab( i ), rhs( i )
    end do
    close( unit=99 )



!    open(unit=99, file='xifull', form='formatted', 

    deallocate( vipt, basfcn, qtab )

  end subroutine schi_sinqr_printSite


end module schi_sinqr
