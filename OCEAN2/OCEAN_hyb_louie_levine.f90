! Copyright (C) 2016, 2017, 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_hyb_louie_levine
  use AI_kinds
  use OCEAN_constants, only : eV2Rydberg, eV2Hartree, Hartree2eV

  private

  integer, parameter :: nrad = 30      ! 30 rads for x=y sphere
  integer, parameter :: nrs = 100
  integer, parameter :: nr = 400
  integer, parameter :: nstep = 2000
  real(dp), parameter :: rsphh = 1.d0
  real(dp), parameter :: rsphl = 0.1d0
!
  integer, parameter :: nq = 2000
  real(dp), parameter :: dqq = 0.005d0
  real(dp), parameter :: rh = 80.d0
  real(dp), parameter :: rl = 0.01d0

  public :: OS_hyb_louie_levine

  contains

  subroutine OS_hyb_louie_levine( sys, nkpts_pad, nxpts_pad, nypts, ladder, nx, nx_start, nkret, & 
                                  kret, ierr, ladcap, kk )
    use OCEAN_system
    use OCEAN_mpi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: nkpts_pad, nxpts_pad, nypts, nx, nx_start
    real(dp), intent( out ) :: ladder( nkpts_pad, nxpts_pad, nypts )
    integer, intent( out ) :: nkret, kret( sys%nkpts )
    integer, intent( inout ) :: ierr
    integer, intent( out ) :: ladcap(2,3)
    integer, intent( out ) :: kk( sys%nkpts, 3 )
    !
    real(dp), allocatable :: rho(:)

    real(dp) :: whomdat( nr, nrs ), w0dat( nrad, nrs ), d_array( nr ), rs_array( nrs ), rad_array( nrad )
    real(dp) :: x_array( 3, sys%nxpts ), r_array( 3, sys%nkpts ), ftab( sys%nxpts )
    real(dp) :: vol
    real(dp) :: pi, rsx, decut, rad0, smear_vol, rsmin, rsmax, fy, gy, de, fff, fx, wx, &
              wy, gx, ww, fcn, qde, dspc, avec( 3, 3 ), vtest( 3 ), gap( 3 ), mds, rmag, maxxy
    
    integer :: iter1, iter2, rad_floor, rad_ceil, clip, rs_floor, rs_ceil, &
              rs_cur, rad_cur, iyr, ixr, irad0, jd, irtab( sys%nxpts ), ix, iy
    integer :: err
    real(dp) :: xk( sys%nkpts, 3 )

#ifdef DEBUG
    complex(dp), allocatable :: compl_rho(:)
    real(dp), allocatable :: ladder2(:,:,:) 
    integer, allocatable :: kret_( : )
    integer :: nkret_, i
    logical :: override_rho, override_ladder
#endif
    
!    real( DP ), external :: AI_max_length


    pi = 4.0d0 * atan( 1.0d0 )


    allocate( rho( sys%nxpts ), STAT=ierr ) 
    if( ierr .ne. 0 ) return

    call OCEAN_get_rho( sys%xmesh, sys%celvol, rho, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) then
      open( unit=99, file='decut', form='formatted', status='old', IOSTAT=err )
      if ( err .ne. 0 ) then
        write(6,*) 'problem with file decut, iostat=', err
        goto 111
      endif
      read( 99, * ) decut
      close( 99 )
111 continue
    endif

!   This allows grabbing in the rho written by the legacy ai2nbse routines
#ifdef DEBUG
    if( myid .eq. root ) then
      inquire(file='rho.xpts', exist=override_rho)
      if( override_rho ) then
        allocate( compl_rho( sys%nxpts ) )
        write(6,*) 'Overriding rho with rho.xpts!'
        open( unit=99, file='rho.xpts', form='unformatted', status='old' )
        rewind 99
        read ( 99 ) compl_rho( : )
        close( unit=99 )
        open(unit=99,file='rho_compare.txt',form='formatted',status='unknown')
        do i = 1, sys%nxpts
          write(99,*) rho(i), real( compl_rho(i), dp )
        enddo
        close(99)
        rho(:) = real( compl_rho(:), dp )
        deallocate( compl_rho )
      endif
    endif
#endif
! end debug flag


#ifdef MPI
    call MPI_BCAST( err, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( err .ne. 0 ) then
      ierr = err
      return
    endif
    call MPI_BCAST( decut, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
!    call MPI_BCAST( clip, 1, MPI_INTEGER, root, comm, ierr )
!    if( ierr .ne. 0 ) return
#endif
    
    call brcapper( sys%nxpts, rho, rsmin, rsmax, ierr )
    if( ierr .ne. 0 ) return 
    call llft( sys%epsilon0, rsmin, rsmax, whomdat, d_array, rs_array, ierr )
    if( ierr .ne. 0 ) return
    call whom0( whomdat, d_array, w0dat, rad_array, ierr )
    if( ierr .ne. 0 ) return
      
    vol = sys%celvol * 0.529177249d0 ** 3
    rad0 = ( 3.0d0 * vol / ( 4.0d0 * pi * dble( sys%nxpts ) ) ) ** ( 1.0d0 / 3.0d0 )
    smear_vol = ( 6.0d0 * pi ** 2 / ( vol * dble( sys%nkpts ) ) ) ** ( 1.0d0 / 3.0d0 )
    if( myid .eq. root ) then
      write(6,*)  " vol =", vol,       "  cubic angstroms."
      write(6,*)  "rad0 =", rad0,      "        angstroms."
      write(6,*)  "   q =", smear_vol, "inverse angstroms."
    endif
    write(1000+myid,*)  " vol =", vol,       "  cubic angstroms."
    write(1000+myid,*)  "rad0 =", rad0,      "        angstroms."
    write(1000+myid,*)  "   q =", smear_vol, "inverse angstroms."

    ! convert avecs to angstroms
    call AI_ladder_formx( sys%xmesh(1), sys%xmesh(2), sys%xmesh(3), x_array )
    call AI_ladder_formr( sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), r_array, ladcap )
    call AI_kmapr( sys%kmesh, xk, kk, ladcap )
   !
    avec( :, : ) = sys%avec( :, : ) * 0.529177d0
    !
    ! find largest distance between any two x-points
    maxxy = AI_max_length( avec, x_array, sys%nxpts )
    !
    ! figure which r-points count for anything
    ! tabulate the ones to retain in kret( 1 : nkret )
    ! mds is a rigorous upper bound on the maximum distance sampled
    mds = 0.d0
    nkret = 0 
    kret( : ) = 0 
    do iter1 = 1, sys%nkpts
      vtest( : ) = 0.d0
      do iter2 = 1, 3
        vtest( : ) = vtest( : ) + avec( :, iter2 ) * r_array( iter2, iter1  )
      enddo
      rmag = sqrt( dot_product( vtest, vtest ) )
      if ( rmag .lt. decut + maxxy ) then
        nkret = nkret + 1
        kret( nkret ) = iter1
        mds = max( mds, rmag + maxxy )
      end if
      if( myid .eq. 0 ) write ( 73, '(2i8,1f20.10)' ) iter1, nkret, rmag

    enddo
    if( myid .eq. 0 ) then
      write ( 6, '(1a21,4f20.10)' ) 'maxxy decut rmag mds ', maxxy, decut, rmag, mds
      write ( 6, '(1a9,2i8)' ) 'nk nkret ', sys%nkpts, nkret
    endif
    !
    call getspacings( avec, gap )
    clip = 0
    do iter1 = 1, 3
      do
        if ( ( 1 + clip ) * sys%kmesh( iter1 ) * gap( iter1 ) - mds .gt. mds ) exit
        clip = clip + 1
      end do
     if( myid .eq. 0 ) write ( 6, '(1a4,1i1,5x,1a7,1i3)' ) 'i = ', iter1, 'clip = ', clip
    end do
    if( myid .eq. 0 ) write ( 6, '(2(1a7,1x,1e15.8))' ) 'rad1 = ', rad_array( 1 ), ' radn = ', rad_array( nrad )

    !
    ! check that rad0 is enclosed in the W0 grid
    if ( rad0 .lt. rad_array( 1 ) ) then
      write(6,*) 'rad0 too small, adjust xmesh', rad0
      ierr = 20
      return
    endif
    if ( rad0 .gt. rad_array( nrad ) )then
      write(6,*) 'rad0 too large, adjust xmesh', rad0
      ierr = 21
      return
    endif

    ! find rad0 location
    rad_floor = 1
    rad_ceil = nrad - 1
    rad_cur = ( rad_ceil + rad_floor ) / 2
    do
      if ( ( rad0 .ge. rad_array( rad_cur ) ) .and. ( rad0 .lt. rad_array( rad_cur + 1 ) ) ) goto 10
      if ( rad_ceil - rad_floor .eq. 1 ) then
        write(6,*) 'bad rad0 and/or irad0'
        ierr = 22
        return
      endif
      if ( rad0 .gt. rad_array( rad_cur ) ) then
        rad_floor = rad_cur
      else
        rad_ceil = rad_cur
      endif
      rad_cur = ( rad_ceil + rad_floor ) / 2
    enddo
  10 continue
    irad0 = rad_cur
    !
    !
    do iter1 = 1, sys%nxpts
      rsx = ( 3.0d0 / ( 4.0d0 * pi * rho( iter1 ) ) ) ** ( 1.0d0 / 3.0d0 )
      rs_floor = 1
      rs_ceil = nrs - 1
      rs_cur = ( rs_floor + rs_ceil ) / 2
      do
        if ( ( rsx .ge. rs_array( rs_cur ) ) .and. ( rsx .le. rs_array( rs_cur + 1 ) ) ) goto 20
        if ( rs_ceil - rs_floor .eq. 1 ) then
          write(6,*) 'bad ixr'
          ierr = 23
          return
        endif
        if ( rsx .gt. rs_array( rs_cur ) ) then
          rs_floor = rs_cur
        else
          rs_ceil = rs_cur
        endif
        rs_cur = ( rs_ceil + rs_floor ) / 2
      enddo
  20  continue
      irtab( iter1 ) = rs_cur
      ftab( iter1 ) = ( rsx - rs_array( rs_cur ) ) / ( rs_array( rs_cur + 1 ) - rs_array( rs_cur ) )
    enddo
    !
    if( myid .eq. 0 ) write(6,*)  'Done with loop'
    dspc = d_array(2) - d_array(1)
    !
    !
    ladder( :, :, : ) = 0.0_dp

    ! This allows grabbing in the ladder.dat file created by the legacy ai2nbse routines
#ifdef DEBUG
    inquire( file='ladder.dat', exist=override_ladder ) 
    if( override_ladder .and. nproc .eq. 1 ) then
      write(6,*) 'Override ladder!!!'
      open(unit=99,file='ladder.dat',form='unformatted', status='old')
      read(99) nkret_
      allocate( kret_( nkret_ ) )!, tempor( 1 : nkret_ )
      read(99) kret_( 1 : nkret_ )

      if( nkret_ .ne. nkret ) then
        override_ladder = .false.
        write(6,*) 'Override failed!!!'
      else
        allocate( ladder2( nkret, sys%nxpts, sys%nxpts ) )
        do ix = 1, sys%nxpts
          do iy = 1, sys%nxpts
            read( 99) ladder2( 1 : nkret, ix, iy ) 
            ladder2( 1 : nkret, ix, iy ) = ladder2( 1 : nkret, ix, iy ) * eV2Hartree
          enddo
        enddo
      endif      
      close( 99 )

!      deallocate( kret_) !, tempor )
    endif
#endif
!   end debug




!    if( .not. override_ladder ) then
!$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( NONE ) SCHEDULE( STATIC ) &
#ifdef DEBUG
!$OMP SHARED( ladder2 )
#endif
!$OMP PRIVATE( ix, ixr, fx, gx, iy, iyr, fy, gy, iter1, de, fff, jd, ww, wx, wy, qde, fcn ) &
!$OMP SHARED( nkret, x_array, r_array, kret, sys, clip, avec, decut, d_array, irad0, dspc ) & 
!$OMP SHARED( rad_array, whomdat, smear_vol, ladder, nx_start, nx, irtab, ftab, w0dat, rad0 )
    do ix = nx_start, nx_start+nx-1
      do iy = 1, sys%nxpts
        ixr = irtab( ix )
        fx = ftab( ix )
        gx = 1.0d0 - fx
        iyr = irtab( iy )
        fy = ftab( iy )
        gy = 1.0d0 - fy
        !
        do iter1 = 1, nkret
          call dist( x_array( :, ix ), x_array( :, iy ), r_array( :, kret( iter1 ) ), de, avec, &
                     sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), clip )
          if( ieee_is_nan( de ) ) then
            write(6,*) 'de', ix, iy, iter1
          endif
          if ( de .le. decut ) then
            if ( de .le. d_array( 1 ) ) then
                   fff = ( rad0 - rad_array( irad0 ) ) / ( rad_array( irad0 + 1 ) - rad_array( irad0 ) )
                   call dublinterp( fx, fff, w0dat( irad0, ixr ), w0dat( irad0, ixr + 1 ), wx )
                   call dublinterp( fy, fff, w0dat( irad0, iyr ), w0dat( irad0, iyr + 1 ), wy )
                   ww = 0.5d0 * ( wx + wy )
            elseif ( de .le. d_array( nr ) ) then
                   jd = 1 + ( de - d_array( 1 ) ) / dspc 
                   if ( jd .lt. 1 ) jd = 1
                   if ( jd .ge. nr ) jd = nr - 1
                   fff = ( de - d_array( jd ) ) / ( d_array( jd + 1 ) - d_array( jd ) )
                   call dublinterp( fx, fff, whomdat( jd, ixr ), whomdat( jd, ixr + 1 ), wx )
                   call dublinterp( fy, fff, whomdat( jd, iyr ), whomdat( jd, iyr + 1 ), wy )
                   ww = 0.5d0 * ( wx + wy ) / de
            else
                   wx = fx * whomdat( jd, ixr + 1 ) + gx * whomdat( jd, ixr )
                   wy = fy * whomdat( jd, iyr + 1 ) + gy * whomdat( jd, iyr )
                   ww = 0.5d0 * ( wx + wy ) / de
            endif

            qde = smear_vol * de
            if ( qde .gt. 0.000001d0 ) then
              fcn = 3.0d0 * ( sin( qde ) - qde * cos( qde ) ) / qde ** 3
            else
              fcn = 1.0d0 - 0.1d0 * qde ** 2
            end if
            ladder( iter1, ix - nx_start + 1, iy ) = 14.400d0 * ww * fcn ** 2 * eV2Hartree ! e^2/Angstrom = 14.4 eV
          endif
#if 0
          if( ieee_is_nan( ladder( iter1, ix - nx_start + 1, iy ) ) ) then
            ierr = 1001
            write(6,*) iter1, ix, iy, ladder( iter1, ix - nx_start + 1, iy )
            write(6,*) de, smear_vol, qde
            write(6,*) ww, wx, wy
            write(6,*) gx, gy, fx, fy
            write(6,*) whomdat( jd, ixr + 1 ), whomdat( jd, ixr )
            write(6,*) whomdat( jd, iyr + 1 ), whomdat( jd, iyr )
            write(6,*) jd, fff, fcn      
            return
          endif
#endif
#ifdef DEBUG
          if( override_ladder ) then
            if( de .le. decut) then
            if( abs(ladder( iter1, ix - nx_start + 1, iy ) - ladder2(iter1, ix - nx_start + 1, iy ) ) &
                .gt. 0.0001 ) then
              write(6,*) iter1, kret( iter1 ), ix - nx_start + 1, iy, ladder( iter1, ix - nx_start + 1, iy ) , & 
                         ladder2(iter1, ix - nx_start + 1, iy ), de, decut
            endif
            elseif( abs(ladder( iter1, ix - nx_start + 1, iy ) - ladder2(iter1, ix - nx_start + 1, iy ) ) &
                .gt. 0.0000001 ) then
              write(6,*) 'decut', iter1, kret( iter1 ), ix - nx_start + 1, iy, ladder( iter1, ix - nx_start + 1, iy ) , &
                         ladder2(iter1, ix - nx_start + 1, iy ), de, decut
            endif
            if( kret( iter1 ) .ne. kret_(iter1 ) ) then
              write(6,*) 'kret', iter1, kret( iter1 ), kret_(iter1 )
            endif
!            ladder( iter1, ix - nx_start + 1, iy ) = ladder2(iter1, ix - nx_start + 1, iy ) / (27.2114_dp*eV2Hartree)
            ladder( iter1, iy, ix - nx_start + 1 ) = ladder2(iter1, iy, ix - nx_start + 1 ) / (27.2114_dp*eV2Hartree)
          endif
#endif
        enddo ! iter1 = 1, nkret
      enddo ! iy = 1, num_xpoints
    enddo ! ix = 1, num_xpoints
!$OMP END PARALLEL DO
!  endif

#ifdef DEBUG
    if( allocated( ladder2 ) ) deallocate( ladder2 )
    if( allocated( kret_ ) ) deallocate( kret_ )
#endif
    if( myid .eq. root ) write( 6, * ) 'Hybertsen-Louie-Levine model W constructed'


  end subroutine OS_hyb_louie_levine



  ! walks through the density and finds the minimum and maximum
  subroutine brcapper( nxprod, rho, rsmin, rsmax, ierr )
    use OCEAN_mpi, only : myid
    implicit none
    !
    integer, intent(in) :: nxprod
    real( DP ), intent(in) :: rho( nxprod )
    real( DP ), intent(out) :: rsmin, rsmax
    integer, intent(inout) :: ierr
    !
    integer :: i
    real( DP ) :: rhomin, rhomax, pi
    !
  !  write(6,*) nxprod, rho(1)
    pi = 4.0d0 * atan( 1.0d0 )
    rhomin = .000001
    do i = 1, nxprod
       if ( i .eq. 1 ) then
          rhomin = rho( 1 )
          if ( rho(1) .gt. 0 ) then
            rhomax = rho( 1 )
          else
            write (6, '(1a26,1e15.8)') 'Warning: negative density ', rho(1)
            ierr = 10
            goto 111
          endif
       end if
       if ( rho( i ) .gt. 0 )  then
         rhomin = min( rhomin, rho( i ) )
       else
         write (6, '(1a26,1e15.8,1i10)') 'Warning: negative density ', rho(i), i
         ierr = 11
         goto 111
       endif
       rhomax = max( rhomax, rho( i ) )
    end do
    write ( 1000 + myid, '(1a9,1x,1e15.8)') 'rhomax = ', rhomax
    write ( 1000 + myid, '(1a9,1x,1e15.8)') 'rhomin = ', rhomin
    rsmin = ( 3.0d0 / ( 4.0d0 * pi * rhomax ) ) ** ( 1.0d0 / 3.0d0 )
    rsmax = ( 3.0d0 / ( 4.0d0 * pi * rhomin ) ) ** ( 1.0d0 / 3.0d0 )
    write ( 1000 + myid, '(2(1a8,1e15.8))' ) 'rsmax = ', rsmax, 'rsmin = ', rsmin
    !
  111 continue
    return
  end subroutine brcapper

  subroutine llft( epsilon0, rsmin, rsmax,  whomdat, d_array, rs_array, ierr )
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none
    !
    real( DP ), intent(in) :: epsilon0, rsmin, rsmax
    real( DP ), intent(out) :: whomdat( nr, nrs ), d_array( nr ), rs_array( nrs )
    integer, intent(inout) :: ierr
    !
    real( DP ) :: remd( nq ), rsh, rsl
    !
    integer :: irs
    real( DP ) :: n, qf, wf, wpsqd, lam
    !
    integer :: ir
    real( DP ) :: r !, rang
    !
    integer j
    real( DP ) :: intr, q, qq, dq
    real( DP ) :: cold, sold, cnew, snew, cs, ss
    !
    complex( DP ) :: ra, rb, c1, a, b, aa
    real( DP ) :: eill, eift, beta, gamma0, c0, c2, c4, u, w, raw
!    real( DP ), external :: levlou
    !
    real( DP ) :: pi
    pi = 4.0d0 * datan( 1.0d0 )
    !
    rsl = rsmin * 0.95d0
    rsh = rsmax * 1.05d0
    !
  !  open( unit=whomdat, file='W.els', form='formatted',
  ! &      status='unknown' )
  !  rewind whomdat 
    !
    do ir = 1, nr
      d_array( ir ) = rl + ( rh - rl ) * dble( ir - 1 ) / dble( nr - 1 )
    enddo
    !
    do irs = 1, nrs
    !
      rs_array( irs ) = rsl + ( rsh - rsl ) * dble( irs - 1 ) / dble( nrs - 1 )

  ! move back from r to the density
  ! The following mimics Hybertsen and Louie Phys Rev B 37. 2733
      n = 1.d0 / ( 4.d0 * pi * rs_array( irs ) ** 3 / 3.d0 )
  ! calculate fermi momentum associated with n
      qf = ( 3.d0 * pi ** 2 * n ) ** ( 1.d0 / 3.d0 )
  ! calculate the plasma frequency (squared) asspciated with n
      wpsqd = 4.d0 * pi * n
  ! fermi energy for n
      wf = 0.5d0 * qf ** 2
  ! lambda(rs) eq. 6 of above
      lam = dsqrt( wpsqd / ( wf ** 2 * ( epsilon0 - 1 ) ) )
  !
      c0 =   1.d0 / epsilon0
      c2 =   c0 ** 2 * 64.d0 / ( 5.d0 * pi * lam ** 4 * qf )
      c4 = - 16.d0 / ( 3.d0 * pi * qf )
  !
      beta = - c2 * dabs( c4 ) / ( 1.d0 - c0 ) ** 2
      gamma0 = c4 / ( c0 - 1.d0 )
  !
      u = - 0.5d0 * beta
      w = - 0.25d0 * ( beta ** 2 - 4.d0 * gamma0 )
  !
      raw = dsqrt( dabs( w ) )
      if ( w .ge. 0.d0 ) then
        a = cmplx( u, - raw, DP )
        b = cmplx( u, + raw, DP )
        aa = cmplx( 0.d0, - 0.5d0 * c4 / raw, DP )
      else
        a = cmplx( u + raw, 0.d0, DP )
        b = cmplx( u - raw, 0.d0, DP )
        aa = cmplx( - 0.5d0 * c4 / raw, 0.d0, DP )
      end if
      ra = sqrt( a )
      if ( dble( ra ) .lt. 0.d0 ) ra = - ra
      rb = sqrt( b )
      if ( dble( rb ) .lt. 0.d0 ) rb = - rb
      c1 = cmplx( 1.d0, 0.d0, DP )
      dq = dqq * qf
  !
      do j = 1, nq
        qq = dqq * ( dble( j ) - 0.5d0 )
        q = qq * qf
        eill = levlou( qq, qf, lam )
        eift = 1.d0 + c4 / ( ( u + qq ** 2 ) ** 2 + w )
        remd( j ) = ( eill - eift ) / q
      end do
  !
      do ir = 1, nr
  !      rang = rl + ( rh - rl ) * dble( ir - 1 ) / dble( nr - 1 )
  !      r = rang / 0.529177d0 ! get r into bohrs
         r = d_array( ir ) / 0.529177d0
  !
  !      open( unit=99, file='prog', form='formatted',
  ! &          status='unknown' )
  !      rewind 99
  !      write ( 99, '(2i5,2x,2i5)' ) ir, nr, irs, nrs
  !      close( unit=99 )
  !
        cs = dcos( dq * r )
        ss = dsin( dq * r )
        q = - dq * 0.5d0
        cold = dcos( q * r )
        sold = dsin( q * r )
        intr = 0.d0
        do j = 1, nq
          q = q + dq
          cnew = cold * cs - sold * ss
          snew = sold * cs + cold * ss
          intr = intr + dq * remd( j ) * snew
          sold = snew
          cold = cnew
        end do ! j = 1, nq
        intr = intr * ( 2.d0 / pi )
        intr = intr + c1 + aa / a * ( c1 - exp( - ra * qf * r ) )  &
                       - aa / b * ( c1 - exp( - rb * qf * r ) )
        if( ieee_is_nan( intr ) ) then
          write( 6, * ) c1, aa, a
          write(6,*) ra, rb, qf, r
          ierr = 1002
          return
        endif
        whomdat( ir, irs ) = intr
  !
  !      write ( whomdat, '(3e16.8)' ) rs, rang, intr
  !
      end do ! ir = 1, nr
  !    write ( whomdat, * )
  !
    end do ! irs = 1, nrs
  !  close( unit=whomdat )
  !
  end subroutine llft

  subroutine whom0( whomdat, d_array, w0dat, rad_array, ierr )
    implicit none
  !
    real( DP ), intent( in ) :: whomdat( nr, nrs )
    real( DP ), intent( in ) :: d_array( nr )
    real( DP ), intent( out ) :: w0dat( nrad, nrs ), rad_array( nrad )
    integer, intent( inout ) :: ierr
  !  include 'whom.h'
  !  integer, parameter :: wdat = 30
  !  integer, parameter :: w0dat = 31
  !
  !  double precision whom( nr, nrs ), d( nr ), rs( nrs )
  !
    integer :: irs, jl, jh, i, irad
    real( DP ) :: whomdat2( nr, nrs )
  !
    real( DP ) :: radl, radh, num, den, rho, drho, w, brack, r
  !
  !  open( unit=wdat, file='W.els', form='formatted', status='old' )
  !  rewind wdat
  !  do irs = 1, nrs
  !    do ir = 1, nr
  !      read ( wdat, * ) rs( irs ), d( ir ),  whom( ir, irs )
  !      whom( ir, irs ) = whom( ir, irs ) - 1.d0
  !    end do
  !  end do
  !  close(unit=wdat)
  !
    whomdat2 = whomdat - 1.d0
    !
  !  open( unit=w0dat, file='W0.els', form='formatted',
  ! &      status='unknown' )
  !  rewind w0dat
  !
    do irad = 1, nrad
      rad_array( irad ) = rsphl + ( rsphh - rsphl ) * dble( irad - 1 ) / dble( nrad - 1 )
    enddo
    do irs = 1, nrs
  !
  !    open( unit=99, file='prog', form='formatted',
  ! &        status='unknown' )
  !    rewind 99
  !    write ( 99, '(2i5)' ) irs, nrs
  !    close( unit=99 )
  !
      do irad = 1, nrad
        r = rsphl + ( rsphh - rsphl ) * dble( irad - 1 ) / dble( nrad - 1 )
        num= 0.d0
        den= 0.d0
        drho = 2.d0 * r / dble( nstep )
        jl = 1
        radl = d_array( 1 )
        jh = 2
        radh = d_array( 2 )
        do i = 1, nstep
          rho = drho * ( dble( i ) - 0.5d0 )
          do while ( radh .lt. rho )
            jl = jh
            radl = radh
            jh = jh + 1
            radh = d_array( jh )
          end do
          if ( radl .gt. rho ) then
            w = whomdat2( jl, irs ) * rho / radl
          else
            w = whomdat2( jl, irs ) +  &
                ( whomdat2( jh, irs ) - whomdat2( jl, irs ) ) *  &
                ( rho - radl ) / ( radh - radl )
          end if
          brack = ( r ** 3 - 0.125d0 * rho ** 3 ) / 3.d0 -  &
                  rho * ( r ** 2 - 0.25d0 * rho ** 2 ) / 4.d0
          num = num + drho * brack * rho * w
          den = den + drho * brack * rho ** 2
        end do ! i = 1, nstep
  !      write ( w0dat, '(3f16.10)' )
  ! &      rs( irs ), r , num / den + 1.2d0 / r
        w0dat( irad, irs ) = num / den + 1.2d0 / r
      end do  ! irad = 1, nrad
  !    write ( w0dat, * )
    end do ! irs = 1, nrs
  !
  !  close( unit=w0dat )
  !
  end subroutine whom0

  subroutine AI_ladder_formx( ngx, ngy, ngz, x )
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
  end subroutine AI_ladder_formx
  !
  subroutine AI_ladder_formr( ngx, ngy, ngz, x, ladcap )
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
    write(6,*) 'ladcap'
    do i = 1, 3
      ladcap( 2, i ) = xmesh( i ) / 2
      ladcap( 1, i ) = 1 + ladcap( 2, i ) - xmesh( i )
      write(6,*) ladcap(:,i)
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
  end subroutine AI_ladder_formr

  function levlou( q, qf, lam )
      implicit none
      double precision q, qf, lam, levlou
      double precision xp, xm, term1, term2, term3
      double precision tmp, c0, c2, c4
      real( DP ), parameter :: pi =3.1415926535897932384d0
!     include 'general.h'
      if ( q .gt. 0.01d0 ) then
        xp = q * ( 2.d0 + q ) / lam
        xm = q * ( 2.d0 - q ) / lam
        term1 =   1.d0 / q ** 2
        term2 = - lam / ( 2.d0 * q ** 3 ) * &
                   ( datan( xp ) + datan( xm ) )
        term3 =   ( lam ** 2 / ( 8.d0 * q ** 5 ) + &
                    1.d0 / ( 2.d0 * q ** 3 ) - &
                    1.d0 / ( 8.d0 * q ) ) * &
                  dlog( ( 1.d0 + xp ** 2 ) / ( 1.d0 + xm ** 2 ) )
        tmp = term1 + term2 + term3
      else
        c0 =   2.6666667d0 / lam ** 2
        c2 = - 6.4d0 / lam ** 4
        c4 =   ( 384.d0 / lam ** 6 - 56.d0 / lam ** 4 ) / 21.d0
        tmp = c0 + c2 * q ** 2 + c4 * q ** 4
      end if
      levlou = 1.d0 / ( 1.d0 + 2.d0 / ( pi * qf ) * tmp )
      return
    end function

  function AI_max_length( avec, x_array, nx ) result (maxxy)
  implicit none
  !
  integer :: nx, i, j, k
  real( DP ) :: avec(3,3), x_array(3,nx), vtest(3), maxxy
  !
  maxxy = 0
  do i = 1, nx
     do j = 1, nx
        vtest( : ) = 0
        do k = 1, 3
           vtest( : ) = vtest( : ) + avec( :, k ) * ( x_array( k, i ) - x_array( k, j ) )
        end do
        maxxy = max( maxxy, dot_product( vtest, vtest ) )
     end do
  end do
  maxxy = sqrt( maxxy )
  return
  end function AI_max_length

subroutine dist( x, y, r, d, avec, nkx, nky, nkz, clip )
  implicit none
  !
  integer :: nkx, nky, nkz, clip
  real( DP ) :: x( 3 ), y( 3 ), r( 3 )
  real( DP ) :: d, avec( 3, 3 )
  !
  integer :: i, iix, iiy, iiz
  real( DP ) :: dsqd, r0( 3 ), r1( 3 ), r2( 3 ), r3( 3 )
  ! 
  r0( : ) = 0
  do i = 1, 3
     r0( : ) = r0( : ) + avec( :, i ) * ( x( i ) - y( i ) - r( i ) )
  end do
  dsqd = r0( 1 ) ** 2 + r0( 2 ) ** 2 + r0( 3 ) ** 2
  do iix = -clip, clip
     r1( : ) = r0( : ) + iix * nkx * avec( :, 1 )
     do iiy = -clip, clip
        r2( : ) = r1( : ) + iiy * nky * avec( :, 2 )
        do iiz = -clip, clip
           r3( : ) = r2( : ) + iiz * nkz * avec( :, 3 )
           dsqd = min( dsqd, r3( 1 ) ** 2 + r3( 2 ) ** 2 + r3( 3 ) ** 2 )
        end do
     end do
  end do
  d = sqrt( abs( dsqd ) )
  !
  return
end subroutine dist



! a function f(x,y) is interpolated in 2-d and the result is returned in [rslt]
! one has 
!
!     vecx1( 1 ) = f(x1,y1), vecx1( 2 )=f(x1,y2)
!     vecx2( 1 ) = f(x2,y1), vecx2( 2 )=f(x1,y2)
! 
!     fx = ( x - x1 ) / ( x2 - x1 )
!     fy = ( y - y1 ) / ( y2 - y1 )
!
subroutine dublinterp( fx, fy, vecx1, vecx2, rslt )
  implicit none
  !
  real( DP ) :: fx, fy, rslt, vecx1( 2 ), vecx2( 2 )
  !
  real( DP ) :: gx, gy
  !
  gx = 1.0d0 - fx
  gy = 1.0d0 - fy
  rslt = gy * ( gx * vecx1( 1 ) + fx * vecx2( 1 ) ) + fy * ( gx * vecx1( 2 ) + fx * vecx2( 2 ) )
  !
  return
end subroutine dublinterp


subroutine AI_kmapr( kmesh, xk, kk, ladcap )
  implicit none
  !
  integer, intent( in ) :: kmesh(3)
  integer, intent( out ) :: kk( kmesh(3), kmesh(2), kmesh(1), 3 )
  integer, intent( in ) ::  ladcap(2,3)
  real( DP ), intent( out ) :: xk( kmesh(3), kmesh(2), kmesh(1), 3 )
  !
  integer :: ikx, iky, ikz, j
  real( DP ) :: c( 3 ), pi
  !
  pi = 4.0d0 * datan( 1.0d0 )
  !
  c( : ) = 2.0d0 * pi / dble( kmesh( : ) )
  do ikx = 1, kmesh(1)
    do iky = 1, kmesh(2)
      do ikz = 1, kmesh(3)
        xk( ikz, iky, ikx, : ) = c( : ) * dble( (/ ikx-1, iky-1, ikz-1 /) ) 
        kk( ikz, iky, ikx, : ) = (/ ikx-1, iky-1, ikz-1 /)
        ! mapping k points onto (approximately) zero-centered grid
        do j = 1, 3
          if( kk( ikz, iky, ikx, j ) .gt. ladcap( 2, j ) ) &
              kk( ikz, iky, ikx, j ) = kk( ikz, iky, ikx, j ) - kmesh( j )
        enddo
!        write(74,*) xk( ikz, iky, ikx, : )
      enddo
    enddo
  enddo
  !

end subroutine AI_kmapr

!  open( unit=99, file='ladcap', form='formatted', status='unknown' )
!  rewind 99
!  read ( 99, * ) kl( : ), kh( : )
!  close ( unit=99 )
!  open( unit=99, file='kpts.dat', form='formatted', status='unknown' )
!  rewind 99
!  nkv( 1 ) = nkx; nkv( 2 ) = nky; nkv( 3 ) = nkz
!  c( : ) = 2.0d0 * pi / dble( nkv( : ) )
!  ik = 0
!  do ikx = 1, nkx
!     do iky = 1, nky
!        do ikz = 1, nkz
!           ikv( 1 ) = ikx - 1; ikv( 2 ) = iky - 1; ikv( 3 ) = ikz - 1
!           ik = ik + 1
!           xk( ik, : ) = c( : ) * dble( ikv( : ) )
!           kk( ik, : ) = ikv( : )
!           read ( 99, * ) ick, xck( : )
!           if ( ick .ne. ik ) stop 'bad ick'
           ! mapping k points onto (approximately) zero-centered grid
!           do j = 1, 3
!              if ( abs( xck( j ) - xk( ik, j ) ) .gt. 0.0001d0 ) stop 'bad kpt'
!              if ( ikv( j ) .gt. ladcap( 2, j ) ) kk( ik, j ) = ikv( j ) - nkv( j )
!           end do
!        end do
!     end do
!  end do
!  close( unit=99 )
!  !
!  return
!end subroutine AI_kmapr

  subroutine getspacings( u, gap )
    implicit none
    !
    real( DP ) :: u( 3, 3 ), gap( 3 )
    !
    integer :: i, j, k
    real( DP ) :: x, perp( 3 )
    !
    do i = 1, 3
       j = i + 1
       if ( j .eq. 4 ) j = 1
       k = j + 1
       if ( k .eq. 4 ) k = 1
       call cross( u( :, j ), u( :, k ), perp( : ) )
       x = sqrt( dot_product( perp( : ), perp( : ) ) )
       perp( : ) = perp( : ) / x
       gap( i ) = abs( dot_product( perp( : ), u( :, i ) ) )
    end do
    !
    return
  end subroutine getspacings

  subroutine cross( v1, v2, c )
    implicit none
    !
    real( DP ) :: v1( 3 ), v2( 3 ), c( 3 )
    !
    integer :: i, j, k
    !
    do i = 1, 3
       j = i + 1; if ( j .eq. 4 ) j = 1
       k = j + 1; if ( k .eq. 4 ) k = 1
       c( i ) = v1( j ) * v2( k ) - v1( k ) * v2( j )
    end do
    !
    return
  end subroutine cross


end module
