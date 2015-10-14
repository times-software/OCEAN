! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine longit( nr, r, dl, ph1, ph2, mel, lprobe, q )
  implicit none
  !
  integer :: nr, lprobe
  double precision :: dl, mel, q
  double precision :: r( nr ), ph1( nr ), ph2( nr )
  !
  integer :: i, j, ii
  double precision :: f( 0 : 4 ), arg, fcn
  !
  mel = 0.d0
  do i = 1, nr - 4, 4
     do j = 0, 4
        ii = i + j
        arg = q * r( ii )
        select case( lprobe )
        case( 0 )
           ! perfect orthogonality simulated by dropping first term ... 
           fcn = -( arg ** 2 / 6.0d0 ) + arg ** 4 / 120.0d0
           if ( arg .gt. 0.001d0 ) then
              fcn = sin( arg ) / arg - 1.0d0
           end if
        case( 1 )
           fcn = arg / 3.0d0 - arg ** 3 / 30.0d0
           if ( arg .gt. 0.001d0 ) then
              fcn = sin( arg ) / arg ** 2 - cos( arg ) / arg
           end if
        case( 2 )
           fcn = arg ** 2 / 15.0d0 - arg ** 4 / 210.0d0
           if ( arg .gt. 0.001d0 ) then
              fcn = ( 3.0d0 / arg ** 3 - 1.0d0 / arg ) * sin( arg ) - 3.0d0 / arg ** 2 * cos( arg )
           end if
        case( 3 )
           fcn = arg ** 3 / 105.0d0 - arg ** 5 / 1890.0d0
           if ( arg .gt. 0.001d0 ) then
              fcn = ( 15.0d0 / arg ** 4 - 6.0d0 / arg ** 2 ) * sin( arg ) - ( 15.0d0 / arg ** 3 - 1.0d0 / arg ) * cos( arg )
           end if
        end select
        f( j ) = ph1( ii ) * ph2( ii ) * r( ii ) * fcn * dble( 2 * lprobe + 1 ) / q
     end do
     mel = mel + 14.d0 * ( f( 0 ) + f( 4 ) )
     mel = mel + 64.d0 * ( f( 1 ) + f( 3 ) )
     mel = mel + 24.d0 * f( 2 )
  end do
  !
  mel = mel * dl / 45.d0
  !
  return
end subroutine longit
