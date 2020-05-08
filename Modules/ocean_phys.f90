! Copyright (C) 2020 OCEAN collaboration
! 
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! This module will contain some of the needed physical cell operations for reuse between sections
! At the moment it just has a way to transform between crystal coordinates and real-space
! There is also a one-way conversion from angstrom to Bohr
!
! The module upscales things to quad precision to avoid numerical problems with the 
! inversion of the lattice vectors and avoid round-off errors in this 
module ocean_phys
  use ai_kinds, only : DP, QP
  implicit none

  private

  public :: ophys_fixCoords, ophys_getOmega, ophys_getBvecs

  contains

!> @brief returns the unit cell volume
  subroutine ophys_getOmega( avecs, omega )
    real(DP), intent( in ) :: avecs(3,3)
    real(DP), intent( out ) :: omega

    real(QP) :: a(3,3), o

    a(:,:) = avecs(:,:)
    call getOmega( a, o )
    omega = o
  end subroutine ophys_getOmega

!> @brief calculates the unit cell volume in QP
  subroutine getOmega( a, o )
    !
    real(QP), intent( in ) :: a( 3, 3 )
    real(QP), intent( out ) :: o

    integer :: i, j, k
    !
    o = 0.0_QP
    do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     o = o + a( i, 1 ) * a( j, 2 ) * a( k, 3 ) &
       - a( i, 1 ) * a( k, 2 ) * a( j, 3 )
    end do
    o = abs( o )
    !
    return
  end subroutine getomega


!> @brief calculates the b-vectors, optionally returns the cell volume
  subroutine ophys_getBvecs( avec, bvec, omega )
    use ocean_constants, only : PI_QP
    !
    real(DP), intent( in ) :: avec( 3, 3 )
    real(DP), intent( out ) :: bvec( 3, 3 )
    real(DP), intent( out ), optional :: omega
    !
    real(QP) :: a(3,3), b(3,3), o
    !

    a(:,:) = avec(:,:)
    call getomega( a, o )
    if( present( omega ) ) omega = o

    o = 2.0_QP * PI_QP / o
    call crossProduct( a(:,2), a(:,3), b(:,1) )
    call crossProduct( a(:,3), a(:,1), b(:,2) )
    call crossProduct( a(:,1), a(:,2), b(:,3) )

    b(:,:) = b(:,:) * o
    
    bvec(:,:) = b(:,:) 

  end subroutine ophys_getBvecs

  subroutine crossProduct( a, b, c )
    real(QP), intent( in ) :: a(3), b(3)
    real(QP), intent( out ) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = - a(1) * b(3) + a(3) * b(1)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end subroutine crossProduct


!> @brief takes the a-vectors and atom coordinates in ssome format and returns both 
!> the real-space coordinates in Bohr and the crystal coordinates
  subroutine ophys_fixCoords( avecsInBohr, inputFormat, xinput, xred, xcoord, ierr )

    use ocean_constants, only : angstrom2Bohr
    real(DP), intent( in ) :: avecsInBohr( 3, 3 )
    character(len=*) :: inputFormat
    real(DP), intent( in ) :: xinput(:,:)
    real(DP), intent( out ) :: xred(:,:)
    real(DP), intent( out ) :: xcoord(:,:)
    integer, intent( inout ) :: ierr

    ! Here we test in/out array sizes

    select case( inputFormat )
      case( 'angstrom' , 'xangst' )
        xcoord(:,:) = xinput(:,:) * angstrom2bohr
        call coord2red( avecsInBohr, xcoord, xred, ierr )
      case( 'bohr' , 'xcart' )
        xcoord(:,:) = xinput(:,:)
        call coord2red( avecsInBohr, xcoord, xred, ierr )
      case( 'xred' )
        xred(:,:) = xinput(:,:)
        call red2coord( avecsInBohr, xred, xcoord, ierr )

      case default
        ierr = -2259
    end select

  end subroutine ophys_fixCoords
    


  subroutine red2coord( avecsInBohr, xred, xcoord, ierr )
    real(DP), intent( in ) :: avecsInBohr( 3, 3 )
    real(DP), intent( in ) :: xred(:,:)
    real(DP), intent( out ) :: xcoord(:,:)
    integer, intent( inout ) :: ierr

    real(QP) :: x(3), y(3)
    integer :: i, j

    do i = 1, size( xred, 2 )
      x(:) = xred(:,i)
      y(:) = 0.0_QP

      do j = 1, 3
        y( : ) = y( : ) + avecsInBohr( :, j ) * x( j )
      enddo

      xcoord(:,i) = y(:)
    enddo
  end subroutine red2coord


  subroutine coord2red( avecsInBohr, xcoord, xred, ierr )
    real(DP), intent( in ) :: avecsInBohr( 3, 3 )
    real(DP), intent( in ) :: xcoord(:,:)
    real(DP), intent( out ) :: xred(:,:)
    integer, intent( inout ) :: ierr

    real(QP) :: invA(3,3), det, x(3), y(3), a(3,3)
!    real(QP) :: b(3,3)
    integer :: i, j

    a(:,:) = avecsInBohr(:,:)
    det = 0.0_QP
    do i = 1, 3
      det = det + a(i,1) * ( a(mod(i,3)+1,2) * a(mod(i+1,3)+1,3) - a(mod(i,3)+1,3) * a(mod(i+1,3)+1,2) )
    enddo

    det = 1.0_QP / det
    do i = 1, 3
      do j = 1, 3
        invA(j,i) = det * &
                  ( a(mod(i,3)+1,mod(j,3)+1)   * a(mod(i+1,3)+1,mod(j+1,3)+1) &
                  - a(mod(i,3)+1,mod(j+1,3)+1) * a(mod(i+1,3)+1,mod(j,3)+1) )
      enddo
    enddo

!    b = matmul( invA, a )
!    do i = 1, 3
!      write(6,*) b(:,i)
!    enddo

    do i = 1, size( xcoord, 2 )
      x(:) = xcoord(:,i)
      y(:) = 0.0_QP

      do j = 1, 3
        y( : ) = y( : ) + invA( :, j ) * x( j )
      enddo

      xred(:,i) = y(:)
    enddo


  end subroutine

end module
