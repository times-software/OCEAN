! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! Intended to be the new home of methods for generating W(r,r';R) 
module OCEAN_WRR
  use AI_kinds, only : DP
  
  implicit none
  private


  public :: OCEAN_WRR_generate

  contains

  subroutine OCEAN_WRR_generate( sys, screeningMethod, nkpts_pad, nxpts_pad, nypts, ladder, & 
                                 nx, nx_start, nkret, kret, ierr, ladcap, kk )
    use OCEAN_system, only : O_system
    use OCEAN_hyb_louie_levine, only : OS_hyb_louie_levine
    !
    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: screeningMethod, nkpts_pad, nxpts_pad, nypts, nx, nx_start
    real(dp), intent( out ) :: ladder( nkpts_pad, nxpts_pad, nypts )
    integer, intent( out ) :: nkret, kret( sys%nkpts )
    integer, intent( inout ) :: ierr
    integer, intent( out ) :: ladcap(2,3)
    integer, intent( out ) :: kk( sys%nkpts, 3 )
    !

    select case( screeningMethod )
      case( 1 )
        call OS_hyb_louie_levine( sys, nkpts_pad, nxpts_pad, nypts, ladder, nx, nx_start, &
                                  nkret, kret, ierr, ladcap, kk )

      case( 2 )
        ! LAZY!!
        !TODO fold OS_hyb_louie_levine generator in here?
        call OS_hyb_louie_levine( sys, nkpts_pad, nxpts_pad, nypts, ladder, nx, nx_start, &
                                  nkret, kret, ierr, ladcap, kk )
        if( ierr .ne. 0 ) return
        call OCEAN_WRR_LocalRPA( sys, nkpts_pad, nxpts_pad, nypts, ladder, nx, nx_start, &
                                  nkret, kret, ierr, ladcap, kk )


      case default
        ierr = -1
        
    end select

    return
  end subroutine OCEAN_WRR_generate
  


  subroutine OCEAN_WRR_LocalRPA( sys, nkpts_pad, nxpts_pad, nypts, ladder, nx, nx_start, &
                                  nkret, kret, ierr, ladcap, kk )
    use OCEAN_system, only : O_system
    use OCEAN_mpi, only : myid, root, MPI_INTEGER, comm, MPI_DOUBLE_PRECISION
    use OCEAN_constants, only : PI_DP, eV2Hartree, bohr
    use OCEAN_sphericalHarmonics, only : ocean_sphH_getylm
    !
    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: nkpts_pad, nxpts_pad, nypts, nx, nx_start
    !TODO ladder is only in/out because of the dumb above where we first calc all of the HLL first
    real(dp), intent( inout ) :: ladder( nkpts_pad, nxpts_pad, nypts )
    !TODO  this stuff is only in because of prev. 
    integer, intent( in ) :: nkret, kret( sys%nkpts )
    integer, intent( inout ) :: ladcap(2,3)
    integer, intent( in ) :: kk( sys%nkpts, 3 )
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: valencePots(:,:,:), tempPot(:,:), r_array(:,:), x_array(:,:), LMvalencePots(:,:,:)
    real(DP) :: decut, valRadius, rmax, dr, wx, wy, w0, cutoff, fcn, dx, f1, f2, qde, q, xred(3), ww, de, ff1, wold, woldx, woldy
    integer :: err, clip, lmax, nLM, nr, ir, iLM, ix, iy, l, m, iter1, vindLen

    character(len=128) :: filnam


    if( myid .eq. root ) then
      open( unit=99, file='decut', form='formatted', status='old', IOSTAT=err )
      if ( err .ne. 0 ) then
        write(6,*) 'problem with file decut, iostat=', err
        goto 111
      endif
      read( 99, * ) decut
      close( 99 )

      open( unit=99, file='screen.lmax', form='formatted', status='old', IOSTAT=err )
      if ( err .ne. 0 ) then
        write(6,*) 'problem with file screen.lmax, iostat=', err
        goto 111 
      endif
      read( 99, * ) lmax
      close( 99 )

      open(unit=99, file='screen.final.rmax', form='formatted', status='old', IOSTAT=err )
      if ( err .ne. 0 ) then
        write(6,*) 'problem with file screen.final.rmax, iostat=', err
        goto 111
      endif
      read( 99, * ) rmax
      close( 99 )

      open(unit=99, file='screen.final.dr', form='formatted', status='old', IOSTAT=err )
      if ( err .ne. 0 ) then
        write(6,*) 'problem with file screen.final.dr, iostat=', err
        goto 111
      endif
      read( 99, * ) dr
      close( 99 )

      open(unit=99, file='cnbse.rad', form='formatted', status='old' )
      read( 99, * ) cutoff
      close( 99 )

      if( lmax > 0 ) then
        open( unit=99, file='vind.length', form='formatted', status='old', IOSTAT=err )
        if( err .ne. 0 ) then
          vindLen = 0
        else
          read( 99, * ) vindLen
          close( 99 )
        endif
      else
        vindLen = 0
      endif

111 continue
    endif

#ifdef MPI
    call MPI_BCAST( err, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( err .ne. 0 ) then
      ierr = err
      return
    endif
    call MPI_BCAST( decut, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( rmax, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( dr, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( cutoff, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( lmax, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( vindLen, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    q = ( 6.0d0 * pi_DP ** 2 / ( sys%celvol * dble( sys%nkpts ) ) ) ** ( 1.0d0 / 3.0d0 )

    nLM = (lmax + 1 ) * ( lmax + 1 )
    nr = floor( rmax / dr )
    valRadius = ( 0.75_DP * sys%celvol / ( PI_DP * real( sys%nxpts, DP ) ) ) ** (1.0_DP/3.0_DP )

    ! To start, we will simply store all the potentials on all the processors. 
    ! In the future, each processor should only store its own, and then some sharing must be done to symmetrize

    if( myid .eq. 0 ) write(6,*) 'VAL RADIUS:', valRadius
    allocate( valencePots( 2, nr, sys%nxpts ), LMvalencePots( vindLen, nLM+1, sys%nxpts ) )
    if( myid .eq. root ) then
      do ix = 1, sys%nxpts
        write(filnam, '(A5,I6.6)' ) 'rpotx', ix
        open( unit=99, file=trim(filnam), form='formatted', status='old' )
        do ir = 1, nr
!          read(99,*) valencePots( :, ir, ix )
           read(99,*) valencePots( 2, ir, ix ), valencePots( 1, ir, ix )
        enddo
        close(99)
        if( vindLen > 0 ) then
          write(filnam, '(A5,I6.6)' ) 'vindx', ix
          open( unit=99, file=trim(filnam), form='formatted', status='old' )
          do ir = 1, vindLen
            read(99,*) LMvalencePots( ir, :, ix )
          enddo
          close(99)
        endif
      enddo
    endif
#ifdef MPI
    call MPI_BCAST( valencePots, 2*nr*sys%nxpts, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( vindLen > 0 ) then
      call MPI_BCAST( LMvalencePots, (nLM+1)*vindLen*sys%nxpts, MPI_DOUBLE_PRECISION, root, comm, ierr )
      if( ierr .ne. 0 ) return
    endif
#endif

    allocate( x_array( 3, sys%nxpts ), r_array( 3, sys%nkpts ) )
    call formx( sys%xmesh(1), sys%xmesh(2), sys%xmesh(3), x_array ) 
    call formr( sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), r_array, ladcap )

    allocate( tempPot( nr, nLM+1 ) )

    ! Loop over my sites
    do ix = nx_start, nx_start+nx-1

!      tempPot = transpose( valencePots( :, :, ix ) )

!      call integratePot( tempPot, valRadius, W0 )
      call integratePot( valencePots( :, :, ix ), valRadius, W0 )
!      W0 = valencePots( 2, 1, ix )

      do iy = 1, sys%nxpts
!        tempPot2 = transpose( valencePots( :, :, iy ) )
      
        do iter1 = 1, nkret

          call distAndAngle( x_array( :, ix ), x_array( :, iy ), r_array( :, kret( iter1 ) ), sys%avec, &
                     sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), clip, de, xred )

          if( de .le. cutoff ) then
!          if( de .lt. ( decut / bohr) .and. de .lt. tempPot( nr, 1 ) ) then

            if( de .lt. 0.01_DP ) then
              ww = w0
              ! this is for some checks
              wx = w0
              wy = w0
            else
              ir = min( floor( de / dr ) + 1, nr - 1 )
            
              if( valencePots( 1, ir, ix ) .gt. de .or. valencePots( 1, ir+1, ix ) .lt. de ) & 
                  write(1000+myid,'(A,3(1X,F20.8))') '!!', valencePots( 1, ir, ix ), de, &
                                                     valencePots( 1, ir+1, ix )
!              dx = valencePots( 1, ir+1, ix ) - valencePots( 1, ir, ix )
              
              ! Do x -> y
              f1 = 0.0_DP
              f2 = 0.0_DP

              f1 = valencePots( 2, ir, ix ) / (2.0_DP * sqrt( pi_dp ) )
              f2 = valencePots( 2, ir+1, ix ) / (2.0_DP * sqrt( pi_dp ) )

              ! valencePots start with radius
              wx = 0.0_DP
              if( de .le. LMvalencePots( vindLen, 1, ix ) ) then
                iLM = 2
                do l = 1, lmax
                  do m = -l, l
                    iLM = iLM + 1
                    call intval( vindLen, LMvalencePots( :, 1, ix ), LMvalencePots( :, iLM, ix ), &
                                 de, ff1, 'cap', 'cap' )
                    wx = wx - ocean_sphH_getylm( xred, l, m ) * ff1
                  enddo
                enddo
              endif

              wx = wx + f2 * ( de - valencePots( 1, ir, ix ) ) / dr &
                      + f1 * ( valencePots( 1, ir+1, ix ) - de ) / dr
              woldx = f2 * ( de - valencePots( 1, ir, ix ) ) / dr &
                      + f1 * ( valencePots( 1, ir+1, ix ) - de ) / dr

              ! now do y -> x
              xred(:) = -xred(:)
              f1 = 0.0_DP
              f2 = 0.0_DP
              ! valencePots start with radius

              f1 = valencePots( 2, ir, iy ) / (2.0_DP * sqrt( pi_dp ) )
              f2 = valencePots( 2, ir+1, iy ) / (2.0_DP * sqrt( pi_dp ) )

              wy = 0.0_DP
              if( de .le. LMvalencePots( vindLen, 1, iy ) ) then
                iLM = 2
                do l = 1, lmax
                  do m = -l, l
                    iLM = iLM + 1
                    call intval( vindLen, LMvalencePots( :, 1, iy ), LMvalencePots( :, iLM, iy ), &
                                 de, ff1, 'cap', 'cap' )
                    wy = wy - ocean_sphH_getylm( xred, l, m ) * ff1
                  enddo
                enddo
              endif

              
              wy = wy + f2 * ( de - valencePots( 1, ir, iy ) ) / dr &
                      + f1 * ( valencePots( 1, ir+1, iy ) - de ) / dr
              woldy = f2 * ( de - valencePots( 1, ir, iy ) ) / dr &
                      + f1 * ( valencePots( 1, ir+1, iy ) - de ) / dr

!              ww = 0.5_DP * ( wx + wy )
              ww = sqrt( pi_dp ) * ( wx + wy )
              wold = sqrt( pi_dp ) * (woldx + woldy )

            endif
            
            qde = q * de
            if( qde .gt. 0.000001d0 ) then
              fcn = 3.0d0 * ( sin( qde ) - qde * cos( qde ) ) / qde ** 3
            else
              fcn = 1.0d0 - 0.1d0 * qde ** 2
            end if
            write(5000+myid,'(3(1X,I6),4(1X,E24.16))') iter1, ix, iy, de, &
                ladder( iter1, ix - nx_start + 1, iy ), ww * fcn**2, wold * fcn**2

            ! already in Ha.
            ladder( iter1, ix - nx_start + 1, iy ) = ww * fcn**2 
  
            ! unlike in HLL, everything here should be in Bohr
          endif
        
        enddo
      enddo
    enddo

    flush( 5000+myid )
    deallocate( valencePots, tempPot, r_array, x_array, LMvalencePots )

    
    ! NEED some scalable way of symmetrizing W(r,r') = W(r'r)
    ! Possibly, do the thing, and then 1 by 1, exchange W(x,y) with one proc, average and save
    !   continue through all of them
    ! Can't save to the same location

    ! Possibly just send the needed ones now?
    ! IE, if my x and y are within range, but I don't have y I know whoever has y needs my x
    

  end subroutine OCEAN_WRR_LocalRPA


  subroutine distAndAngle( x, y, r, avec, nkx, nky, nkz, clip, d, xred )
    implicit none
    !
    integer, intent( in ) :: nkx, nky, nkz, clip
    real( DP ), intent( in ) :: x( 3 ), y( 3 ), r( 3 ), avec( 3, 3 )
    real( DP ), intent( out ) :: d, xred(3)
    !
    integer :: i, iix, iiy, iiz
    real( DP ) :: dsqd, r0( 3 ), r1( 3 ), r2( 3 ), r3( 3 ), tempD
    ! 
    r0( : ) = 0
    do i = 1, 3
       r0( : ) = r0( : ) + avec( :, i ) * ( x( i ) - y( i ) - r( i ) )
    end do
    dsqd = r0( 1 ) ** 2 + r0( 2 ) ** 2 + r0( 3 ) ** 2
    if( dsqd .lt. 0.000001_DP ) then
      xred(:) = 0.0_DP
    endif
    xred(:) = r0(:)
    do iix = -clip, clip
       r1( : ) = r0( : ) + iix * nkx * avec( :, 1 )
       do iiy = -clip, clip
          r2( : ) = r1( : ) + iiy * nky * avec( :, 2 )
          do iiz = -clip, clip
             r3( : ) = r2( : ) + iiz * nkz * avec( :, 3 )
             tempD = r3( 1 ) ** 2 + r3( 2 ) ** 2 + r3( 3 ) ** 2
             if( tempD .lt. d ) then
                d = tempD
                xred(:) = r3(:)
             endif
!             dsqd = min( dsqd, r3( 1 ) ** 2 + r3( 2 ) ** 2 + r3( 3 ) ** 2 )
          end do
       end do
    end do
    d = sqrt( dsqd ) 
    xred(:) = xred(:) / d
    !
    return
  end subroutine distAndAngle

  subroutine integratePot( tempPot, valRadius, W0 )
    use OCEAN_constants, only : PI_DP
    use ocean_mpi
    real(DP), intent( in ) :: tempPot( :, : )
    real(DP), intent( in ) :: valRadius
    real(DP), intent( out ) :: W0
    !
    real(DP) :: m, b, x1, x2, xstop
    integer :: nr, ir, i

    nr = size( tempPot, 2 )

    ! r^2 dr of a line segment
    W0 = 0.0_DP
    do ir = 1, nr-1
      if( tempPot( 1, ir ) .gt. valRadius ) exit
      xstop = min( tempPot( 1, ir+1 ), valRadius )

      x1 = tempPot( 1, ir )**3
      x2 = xstop**3
!      m = ( tempPot( 2, ir+1 ) - tempPot( 2, ir ) ) / ( xstop - tempPot( 1, ir ) )
      m = ( tempPot( 2, ir+1 ) - tempPot( 2, ir ) ) / ( tempPot( 1, ir+1 ) - tempPot( 1, ir ) )
      b = tempPot( 2, ir ) - tempPot( 1, ir ) * m

      
      W0 = W0 + (1.0_DP/3.0_DP) * b * ( x2 - x1) & 
         + 0.25_DP * ( x2*xstop - x1*tempPot( 1, ir ) ) * m
!      write(1000+myid, '(8(E24.16,X))' ) tempPot( 1, ir ), tempPot( 1, ir+1 ), tempPot( 2, ir ), tempPot( 2, ir+1 ), m, b, 4*pi_dp*W0, 4*pi_dp*(W0 * 0.75_DP /( PI_DP *x2 ))
    enddo
!    W0 = W0 * 0.75_DP / ( PI_DP * valRadius**3 )


    ! factor of 4pi from the integral (over angles not actually done)
    ! then divide by volume = 4/3 pi r^3
    W0 = W0 * 3.0_DP / ( valRadius**3 )

    return

    W0 = 0.0_DP
    i = 0
    do ir = 1, nr-1
      if( tempPot( 1, ir ) .gt. valRadius ) exit
      W0 = W0 + tempPot( 2, ir )
      i = i + 1
    enddo
    if( i .gt. 0 ) then
      W0 = W0 / real( i, DP )
    endif
    


  end subroutine integratePot


  subroutine formx( ngx, ngy, ngz, x )
    implicit none
    !
    integer, intent( in ) :: ngx, ngy, ngz
    real(dp), intent( out ) :: x( 3, ngx * ngy * ngz )
    !
    integer :: i, ix, iy, iz
    !
    !
    ! do not change this convention for x( 1 ), x( 2 ) ...
    ! note that the z coordinate loops innermost.
    !
    i = 0
    do ix = 1, ngx
       do iy = 1, ngy
          do iz = 1, ngz
             i = i + 1
             x( 1, i ) = dble( ix - 1 ) / dble( ngx )
             x( 2, i ) = dble( iy - 1 ) / dble( ngy )
             x( 3, i ) = dble( iz - 1 ) / dble( ngz )
          end do
       end do
    end do
    !
    return
  end subroutine formx
  !
  subroutine formr( ngx, ngy, ngz, x, ladcap )
    implicit none
    !
    integer, intent(in) :: ngx, ngy, ngz
    integer, intent( out ) :: ladcap(2,3)
    real(dp), intent(out) :: x( 3, ngx * ngy * ngz )
    !
    integer :: i, ix, iy, iz, xmesh(3)
    !
    !
    ! the code as such makes x innermost, because x is the first 
    ! index in the fft mesh when iterating the BSE.
    !
    xmesh( 1 ) = ngx
    xmesh(2) = ngy
    xmesh(3) = ngz
    do i = 1, 3
      ladcap( 2, i ) = xmesh( i ) / 2
      ladcap( 1, i ) = 1 + ladcap( 2, i ) - xmesh( i )
    enddo
    !
    i = 0
    do iz = 1+ngz/2-ngz, ngz/2
       do iy = 1+ngy/2-ngy, ngy/2
          do ix = 1+ngx/2-ngx, ngx/2
             i = i + 1
             x( 1, i ) = dble( ix )
             x( 2, i ) = dble( iy )
             x( 3, i ) = dble( iz )
          end do
       end do
    end do
    !
    return
  end subroutine formr

end module OCEAN_WRR
