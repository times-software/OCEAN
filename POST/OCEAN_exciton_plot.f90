! Copyright (C) 2017 - 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
program OCEAN_exciton_plot
  use periodic, only : get_atom_number
  use ocean_interpolate
  implicit none

  integer, parameter :: DP = kind(1.0d0 )
  complex(DP), allocatable :: exciton(:,:,:), cond_exciton(:,:,:), Rspace_exciton(:,:,:,:), u2(:,:), rk_exciton(:,:,:)
  complex(DP) :: cphs, Pgrid(4), P(4,4), Qgrid(4), Q(4), R

  real(DP), allocatable :: z_stripe( : ), xyz(:,:), atom_loc(:,:), cubeExciton(:,:,:)
  real(DP) :: qinb(3), avecs(3,3), su, k0(3), qvec(3), Rvec(3), xphs, yphs, zphs, twopi, tau(3), ur, ui, suTarg
  real(DP) :: realSpaceBox(3), rsDelta, inverseA(3,3), temp1(3), distance(3), isoTarg, isoMin, isoMax, su2, isos(2,19)

  integer, allocatable :: ibeg(:,:)
  integer :: Rmesh(3), kmesh(3), nband, nalpha, nkpts, NR, Riter, kiter, xmesh(3), nspn, ispin, ivh2
  integer :: RmeshMin(3), RmeshMax(3)
  integer :: ikx, iky, ikz, iRx, iRy, iRz, NX, i, ix, x_count, xiter, iy, iz, izz, bloch_selector, icms, ivms, icml
  integer :: brange(4), u2size, u2start, Rshift(3), natom, kiter_break, Rstart(3), idum(3), ZNL(3)
  integer :: xmeshOut(3), xtarg(3), natom2, atomTarg, ixx, iyy, ixxx, iyyy, iter, j, k
  character(len=25) :: filname
  character(len=128) :: outname
  character(len=2), allocatable :: elname(:)

  logical :: metal, legacy_ibeg, oldStyle, foundMax

  real(DP), external :: DZNRM2
  complex(DP), parameter :: one = 1.0_dp
  complex(DP), parameter :: zero = 0.0_dp

  twopi = 2.0_dp * 4.0_dp * atan( 1.0_dp )

!  oldStyle = .false.
!  atomTarg = 3
!  realSpaceBox(1) = 16.0_DP
!  realSpaceBox(2) = 16.0_DP
!  realSpaceBox(3) = 16.0_DP
!  rsDelta = 0.12_DP

  open(unit=99,file='exciton_plot.ipt',form='formatted',status='old')
  read(99,*) filname
  read(99,*) outname
  read(99,*) oldStyle
  if( oldStyle ) then
    read(99,*) Rmesh(:)
    read(99,*) Rstart(:)
  else
    read(99,*) atomTarg
    read(99,*) realSpaceBox(:)
    read(99,*) rsDelta
    do i = 1, 3
      xmeshOut(i) = 1+ realSpaceBox(i) / rsDelta
    enddo
  endif
!  read(99,*) tau(:)
  close(99)
  tau(:) = 0.0_dp

  open(unit=99,file='nbuse.ipt',form='formatted',status='old')
  read(99,*) nband
  close(99)

  open(unit=99,file='brange.ipt',form='formatted',status='old')
  read(99,*) brange(1:4)
  close(99)

  open(unit=99,file='kmesh.ipt',form='formatted',status='old')
  read(99,*) kmesh(:)
  close(99)
  nkpts = product( kmesh(:) )

  open(unit=99,file='nspin',form='formatted',status='old')
  read(99,*) nspn
  close(99)
!  if( nspn .ne. 1 ) then
!    write(6,*) 'WARNING! Spin not yet supported!'
!    goto 111
!  endif
!  nspn = 1

  open(unit=99,file='qinunitsofbvectors.ipt',form='formatted',status='old')
  read(99,*) qinb(:)
  close(99)

  open(unit=99,file='k0.ipt',form='formatted',status='old')
  read(99,*) k0(:)
  close(99)

  open(unit=99,file='xmesh.ipt',form='formatted',status='old')
  read(99,*) xmesh(:)
  close(99)
  NX = product( xmesh(:) )

  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  read(99,*) avecs(:,:)
  close(99)

  call invA( avecs, inverseA )

  open(unit=99,file='xyz.wyck',form='formatted',status='old')
  read(99,*) natom
  allocate( xyz(3,natom), atom_loc(3,natom), elname( natom ) ) 
  do i = 1, natom
    read(99,*) elname( i ), xyz(:,i)
  enddo
  close(99)

  open(unit=99,file='bloch_selector',form='formatted',status='old')
  read(99,*) bloch_selector
  close(99)


  open(unit=99,file='ZNL',form='formatted',status='old')
  read(99,*) ZNL(:)
  close(99)
  !(2l+1)
  nalpha = (2*ZNL(3)+1) * 4

  inquire(file='force_legacy_ibeg.ipt', exist=legacy_ibeg )
  if( legacy_ibeg ) then
    open( unit=99, file='force_legacy_ibeg.ipt', form='formatted',status='old')
    read( 99, * ) legacy_ibeg
    close( 99 )
  endif

  atom_loc = 0.0_DP
  do i = 1, natom
    do ix = 1, 3
      atom_loc(:,i) = atom_loc(:,i) + avecs(:,ix) * xyz(ix,i)
    enddo
  enddo


  inquire(file='metal', exist=metal )
  if( metal ) then
    open(unit=99,file='metal',form='formatted',status='old')
    read(99,*) metal
    close(99)
  endif

  ! For our purposes metal means play games with ibeg
  ! This means skipping over bands and the like
  ! But the newer version of OCEAN doesn't bother
  if( legacy_ibeg .eqv. .false. ) metal = .false.

  if( metal ) then
    allocate( ibeg( nkpts, nspn ) )
    open( unit=99,file='ibeg.h',form='formatted',status='old')
    do ispin = 1, nspn
      do kiter = 1, nkpts
        read(99,*) idum(1), ibeg( kiter, ispin )
      enddo
    enddo
    close( 99 )
  else
    allocate( ibeg( 1, 1 ) )
    ivh2 = brange( 2 )
  endif


  ! Surely there is a smatter way, but this will be fast enough
  if( .not. oldStyle ) then
    RmeshMin(:) = 0
    RmeshMax(:) = 0

!    do ix = 1, xmeshOut(1)
!      rvec(1) = dble(ix-1) * rsDelta - realSpaceBox(1)/2.0_DP + atom_loc(1, atomTarg )
!      do iy = 1, xmeshOut(2)
!        rvec(2) = dble(iy-1) * rsDelta - realSpaceBox(2)/2.0_DP + atom_loc(2, atomTarg )
    do ix = 0, 1
      rvec(1) = dble(ix) * rsDelta * dble( xmeshOut(1)-2) - realSpaceBox(1)/2.0_DP + atom_loc(1, atomTarg )
      do iy = 0,1
        rvec(2) = dble(iy) * rsDelta * dble( xmeshOut(2)-2) - realSpaceBox(2)/2.0_DP + atom_loc(2, atomTarg )
        do iz = 0, 1
          rvec(3) = dble(iz) * rsDelta * dble( xmeshOut(3) -2 ) - realSpaceBox(3)/2.0_DP + atom_loc(3, atomTarg )

          temp1 = matmul( inverseA, rvec )
          write(6,*) temp1
          do i = 1, 3
            RmeshMin(i) = min( RmeshMin(i), floor( temp1(i) ) )
            RmeshMax(i) = max( RmeshMax(i), ceiling( temp1(i) ) )
          enddo
        enddo
      enddo
    enddo

    write(6,*) RmeshMin(:)
    write(6,*) RmeshMax(:)
    write(6,*) RmeshMax(:)-RmeshMin(:)+1

    Rmesh(:) = RmeshMax(:)-RmeshMin(:)+1
    Rstart(:) = RmeshMin(:)
    
  endif


!JTV!
!  nalpha = 4
!  Rmesh(:) = 3
!  Rstart(:) = -1
  write(6,*) Rmesh(:)
  write(6,*) Rstart(:)
  NR = product( Rmesh(:) )


  write(6,*) nband, nkpts, nalpha

  allocate( exciton( nband, nkpts, nalpha ), cond_exciton( nband, nkpts, nspn ) )
  write(6,*) filname
  open(unit=99,file=filname,form='unformatted',status='old')
  read(99) exciton
  close(99)

  ! Want to normalize excitonic wvfn   
  su = DZNRM2( nband*nkpts*nalpha, exciton, 1 )
  su = 1.0_DP / su
  write(6,*) su

  cond_exciton(:,:,:) = 0.0_DP
  
!  do i = 1, nalpha
  i = 0
  do icms = 1, 2
    do icml = -ZNL(3), ZNL(3)
      do ivms = 1, 2
        i = i + 1
          cond_exciton(:,:,min(ivms,nspn)) = cond_exciton(:,:,min(ivms,nspn)) + su * exciton(:,:,i)
      enddo
    enddo
  enddo

 
  u2size = brange(4)-brange(3)+brange(2)-brange(1)+2
  u2start = brange(2)-brange(1)+2
  write(6,*) u2size, u2start, nband
  allocate( u2( NX, u2size ) )
  allocate( rk_exciton( NX, nkpts, nspn ) )
  rk_exciton(:,:,:) = 0.0d0

  write(6,*) 'Opening u2'
  kiter_break = nkpts / 20

  select case( bloch_selector )

  case( 2, 3 )
    deallocate( u2 )
    allocate( u2( NX, nband) )
    open( unit=99,file='con.u2.dat',access='stream',status='old',form='unformatted' )
    do ispin = 1, nspn
      do kiter = 1, nkpts
        read(99) u2
        call ZGEMV( 'N', NX, nband, one, u2, NX, cond_exciton( 1, kiter, ispin ), 1, &
                     one, rk_exciton( 1, kiter, ispin ), 1 )
      enddo
    enddo
    close( 99 )
  
  case( 1 )
    open(unit=99,file='u2par.dat',access='stream',status='old',form='unformatted' )

    if( nspn .ne. 1 ) then  
      write(6,*) 'Spin not implemented for u2par.dat'
      goto 111
    endif
    ispin = 1
    do kiter = 1, nkpts
      if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
      read(99) u2
      if( metal ) ivh2 = ibeg( kiter, 1 ) - 1
      u2start = ivh2 + 2
      call ZGEMV( 'N', NX, nband, one, u2(1,u2start), NX, cond_exciton( 1, kiter, ispin ), 1, &
                   one, rk_exciton( 1, kiter, ispin ), 1 )
    enddo

    close( 99 )

  case( 0 )
    open(unit=99,file='u2.dat',form='unformatted',status='old')
    do ispin = 1, nspn
      do kiter = 1, nkpts
        if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
        if( metal ) ivh2 = ibeg( kiter, 1 ) - 1
        do i = 1, ivh2 - brange( 1 ) + 1 !brange(2)-brange(1)+1
          do ix = 1, NX
            read(99) 
          enddo
        enddo
        do i = 1, nband
          do ix =1, NX
            read(99) idum(1:3), ur, ui
            u2( ix, i ) = cmplx( ur, ui, DP )
          enddo
        enddo
        do i = ivh2 + nband + 1, brange(2)-brange(1)+brange(4)-brange(3) + 2
          do ix = 1, NX
             read ( 99 )
          end do
        enddo

        call ZGEMV( 'N', NX, nband, one, u2, NX, cond_exciton( 1, kiter, ispin ), 1, &
                     one, rk_exciton( 1, kiter, ispin ), 1 )
      enddo
    enddo
    close( 99 )
  case default
    stop
  end select

  write(6,*) 'Done loading u2'

  
  deallocate( u2, cond_exciton )


  select case( bloch_selector )

  case( 0 , 1 )
! Add in phases for rk_exciton
  kiter = 0
  do ikx = 0, kmesh(1)-1
    qvec(1) = qinb(1) + (k0(1) + dble(ikx))/dble(kmesh(1))
    do iky = 0, kmesh(2)-1
      qvec(2) = qinb(2) + (k0(2) + dble(iky))/dble(kmesh(2))
      do ikz = 0, kmesh(3)-1
        kiter = kiter + 1
        qvec(3) = qinb(3) + (k0(3) + dble(ikz))/dble(kmesh(3))

        xiter = 0
        do ix = 0, xmesh(1)-1
          Rvec(1) = twopi * (dble(ix)/dble(xmesh(1)) - tau(1))
          xphs = Rvec(1) * qvec(1)
          do iy = 0, xmesh(2)-1
            Rvec(2) = twopi * (dble(iy)/dble(xmesh(2)) - tau(2))
            yphs = Rvec(2) * qvec(2) + xphs
            do iz = 0, xmesh(3)-1
              xiter = xiter + 1
              Rvec(3) = twopi * (dble(iz)/dble(xmesh(3)) - tau(3))
              zphs = Rvec(3) * qvec(3) + yphs
              cphs = cmplx( cos( zphs ), sin( zphs ) )
              ! Not sure how bad the cost of out-of-order memory will be
              rk_exciton( xiter, kiter, : ) = rk_exciton( xiter, kiter, : ) * cphs
              
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! new con.u2.dat and val.u2.dat are (x,y,z), not (z,y,x)
  case( 2, 3 )
  kiter = 0
  do ikx = 0, kmesh(1)-1
    qvec(1) = qinb(1) + (k0(1) + dble(ikx))/dble(kmesh(1))
    do iky = 0, kmesh(2)-1
      qvec(2) = qinb(2) + (k0(2) + dble(iky))/dble(kmesh(2))
      do ikz = 0, kmesh(3)-1
        kiter = kiter + 1
        qvec(3) = qinb(3) + (k0(3) + dble(ikz))/dble(kmesh(3))

        xiter = 0
        do iz = 0, xmesh(3)-1
          Rvec(3) = twopi * (dble(iz)/dble(xmesh(3)) - tau(3))
          zphs = Rvec(3) * qvec(3)
          do iy = 0, xmesh(2)-1
            Rvec(2) = twopi * (dble(iy)/dble(xmesh(2)) - tau(2))
            yphs = Rvec(2) * qvec(2) + zphs
            do ix = 0, xmesh(1)-1
              xiter = xiter + 1
              Rvec(1) = twopi * (dble(ix)/dble(xmesh(1)) - tau(1))
              xphs = Rvec(1) * qvec(1) + yphs
              cphs = cmplx( cos( xphs ), sin( xphs ) )
              ! Not sure how bad the cost of out-of-order memory will be
              rk_exciton( xiter, kiter, : ) = rk_exciton( xiter, kiter, : ) * cphs

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo



  end select
  write(6,*) 'Done with rk exciton'


! Do Fourier transform into Rspace_exciton
!  allocate( Rspace_exciton( NX, NR, nspn ) )
  allocate( Rspace_exciton( xmesh(3)*Rmesh(3), xmesh(2)*Rmesh(2), xmesh(1)*Rmesh(1), nspn ) )
  Rspace_exciton = 0.0_DP

! Using slow FT and 6 loops for clarity
  
!  qvec(1) = qinb(1) + k0(1)/dble(kmesh(1))
  do ispin = 1, nspn
    kiter = 0
    do ikx = 0, kmesh(1)-1
      qvec(1) = qinb(1) + (k0(1) + dble(ikx))/dble(kmesh(1))
  !    qvec(1) = qinb(1) + dble(ikx)/dble(kmesh(1))
  !    qvec(2) = qinb(2) + k0(2)/dble(kmesh(2))
      do iky = 0, kmesh(2)-1
        qvec(2) = qinb(2) + (k0(2) + dble(iky))/dble(kmesh(2))
  !      qvec(2) = qinb(2) + dble(iky)/dble(kmesh(2))
  !      qvec(3) = qinb(3) + k0(3)/dble(kmesh(3))
        do ikz = 0, kmesh(3)-1
          qvec(3) = qinb(3) + (k0(3) + dble(ikz))/dble(kmesh(3))
  !        qvec(3) = qinb(3) + dble(ikz)/dble(kmesh(3))
          
          kiter = kiter + 1
          Riter = 0

          do iRx = Rstart(1), Rstart(1) + Rmesh(1) - 1
            xphs = dble( iRx ) * qvec(1) 
            do iRy = Rstart(2), Rstart(2) + Rmesh(2) - 1
              yphs = xphs + dble(iRy) * qvec(2) 
              do iRz = Rstart(3), Rstart(3) + Rmesh(3) - 1
                Riter = Riter + 1
                zphs = yphs + dble( iRz ) * qvec(3)
                cphs = cmplx( cos( twopi * zphs ), sin( twopi * zphs ) )

  !JTV PHASE!!
                ! phase info is wrong here
                xiter = 0
                do ix = 1, xmesh(1)
                  do iy = 1, xmesh(2)
                    do iz = 1, xmesh(3)
                      if( bloch_selector .eq. 2 .or. bloch_selector .eq. 3 ) then
                        xiter = ix + (iy-1)*xmesh(1) + (iz-1)*xmesh(1)*xmesh(2)
                      else
                        xiter = xiter + 1
                      endif
                      Rspace_exciton( iz + (iRz-Rstart(3))*xmesh(3), &
                                      iy + (iRy-Rstart(2))*xmesh(2), &
                                      ix + (iRx-Rstart(1))*xmesh(1), ispin ) = &
                      Rspace_exciton( iz + (iRz-Rstart(3))*xmesh(3), &
                                      iy + (iRy-Rstart(2))*xmesh(2), &
                                      ix + (iRx-Rstart(1))*xmesh(1), ispin ) + &
                                    cphs * rk_exciton( xiter, kiter, ispin )
                      enddo
                    enddo
                  enddo
!                do xiter = 1, NX
!                  Rspace_exciton( xiter, Riter, ispin ) = Rspace_exciton( xiter, Riter, ispin ) & 
!                                                        + cphs * rk_exciton( xiter, kiter, ispin )
!                enddo

              enddo
            enddo
          enddo

  !        qvec(3) = qvec(3) + 1.0_dp / dble( kmesh(3) )
        enddo
  !      qvec(2) = qvec(2) + 1.0_dp / dble( kmesh(2) )
      enddo
  !    qvec(1) = qvec(1) + 1.0_dp / dble( kmesh(1) )
    enddo
  enddo


  write(6,*) 'Done with Rspace exciton'
  deallocate( rk_exciton )


  if( oldStyle ) then
    do ix = 1, 3
      if( Rmesh(ix) .gt. 1 ) then
        Rshift(ix) = Rmesh(ix)/2
      else
        Rshift(ix) = 0
      endif
      Rshift( ix ) = - Rstart( ix )
    enddo

    do i = 1, natom
      do ix = 1, 3
        atom_loc( :, i ) = atom_loc( :, i ) + Rshift( ix ) * avecs(:, ix )
      enddo
    enddo
  endif

  
  ! Now ready to write out 
!  open(unit=99,file='out.cube',form='formatted')
  
!  outname = trim(filname)//'.cube'
  open(unit=99,file=outname,form='formatted')
  write(99,*) "OCEAN exciton plot"
  write(99,*) "---"
  if( .not. oldStyle ) then
    natom2 = 0
    do iRx = Rstart(1), Rstart(1) + Rmesh(1) - 1
      do iRy = Rstart(2), Rstart(2) + Rmesh(2) - 1
        do iRz = Rstart(3), Rstart(3) + Rmesh(3) - 1
          do i = 1, natom
!            call get_atom_number( ix, elname( i ) )
            tau(:) = atom_loc(:,i) - atom_loc(:, atomTarg )
            tau( : ) = tau( : ) + iRx * avecs(:, 1 )
            tau( : ) = tau( : ) + iRy * avecs(:, 2 )
            tau( : ) = tau( : ) + iRz * avecs(:, 3 )

            if( abs( tau( 1 ) ) .gt. realSpaceBox(1)/2.0_DP ) cycle
            if( abs( tau( 2 ) ) .gt. realSpaceBox(2)/2.0_DP ) cycle
            if( abs( tau( 3 ) ) .gt. realSpaceBox(3)/2.0_DP ) cycle
            natom2 = natom2 + 1
          enddo
        enddo
      enddo
    enddo

    allocate(z_stripe(xmeshOut(3)),cubeExciton(xmeshOut(3),xmeshOut(2),xmeshOut(1)))
    write(99,'(I5,3(F12.6))') natom2, 0.0_dp, 0.0_dp, 0.0_dp
    write(99,'(I5,3(F12.6))') xmeshOut(1), rsDelta, 0.0_DP, 0.0_DP
    write(99,'(I5,3(F12.6))') xmeshOut(2), 0.0_DP, rsDelta, 0.0_DP
    write(99,'(I5,3(F12.6))') xmeshOut(3), 0.0_DP, 0.0_DP, rsDelta
    do iRx = Rstart(1), Rstart(1) + Rmesh(1) - 1
      do iRy = Rstart(2), Rstart(2) + Rmesh(2) - 1 
        do iRz = Rstart(3), Rstart(3) + Rmesh(3) - 1 
          do i = 1, natom
            call get_atom_number( ix, elname( i ) ) 
            tau(:) = atom_loc(:,i) - atom_loc(:, atomTarg ) 
            tau( : ) = tau( : ) + iRx * avecs(:, 1 ) 
            tau( : ) = tau( : ) + iRy * avecs(:, 2 ) 
            tau( : ) = tau( : ) + iRz * avecs(:, 3 ) 

            if( abs( tau( 1 ) ) .gt. realSpaceBox(1)/2.0_DP ) cycle
            if( abs( tau( 2 ) ) .gt. realSpaceBox(2)/2.0_DP ) cycle
            if( abs( tau( 3 ) ) .gt. realSpaceBox(3)/2.0_DP ) cycle
            tau( 1 ) = tau( 1 ) + realSpaceBox(1)/2.0_DP
            tau( 2 ) = tau( 2 ) + realSpaceBox(2)/2.0_DP
            tau( 3 ) = tau( 3 ) + realSpaceBox(3)/2.0_DP
            write(99,'(I5,4(F12.6))') ix, 0.0, tau(:)
          enddo
        enddo
      enddo
    enddo
    do ix = 1, xmeshOut(1)
      rvec(1) = dble(ix-1) * rsDelta - realSpaceBox(1)/2.0_DP + atom_loc(1, atomTarg )
      do iy = 1, xmeshOut(2)
        rvec(2) = dble(iy-1) * rsDelta - realSpaceBox(2)/2.0_DP + atom_loc(2, atomTarg )
        do iz = 1, xmeshOut(3)
          rvec(3) = dble(iz-1) * rsDelta - realSpaceBox(3)/2.0_DP + atom_loc(3, atomTarg )

!          temp1 = matmul( inverseA, rvec ) - dble(Rstart(:))
!          temp1(:) = temp1(:) * dble(xmesh(:))
          temp1 = matmul( inverseA, rvec )
          temp1(:) = temp1(:) * dble(xmesh(:)) - Rstart(:)*xmesh(:) + 1
          xtarg(:) = floor( temp1(:) )
!          xtarg(:) = floor(temp1(:)) + 1 - Rstart(:)*xmesh(:)
          distance(:) = temp1(:) - xtarg(:)

          do ixx = xtarg(1)-1, xtarg(1)+2
            ixxx = ixx
            if( ixx .lt. 1 ) ixxx = 1
            if( ixx .gt. xmesh(1)*Rmesh(1) ) ixxx = xmesh(1)*Rmesh(1)
            do iyy = xtarg(2)-1, xtarg(2)+2
              iyyy = iyy
              if( iyy .lt. 1 ) iyyy = 1
              if( iyy .gt. xmesh(2)*Rmesh(2) ) iyyy = xmesh(2)*Rmesh(2)
              call MakeLagrange( 4, xtarg(3), iyyy, ixxx, Rspace_exciton(:,:,:,1), &
                                 Pgrid(:) )
              P( iyy -xtarg(2) + 2, ixx-xtarg(1)+2 ) = evalLagrange( 4, distance(3), 1.0_DP, Pgrid(:) )
            enddo
          enddo

          do ixx = 1, 4
            call MakeLagrange( 4, P(:,ixx), Qgrid(:) )
            Q(ixx) = evalLagrange( 4, distance(2), 1.0_DP, Qgrid(:) )
          enddo
          call makeLagrange( 4, Q, Qgrid )
          R = evalLagrange( 4, distance(1), 1.0_DP, Qgrid )
            

!          write(6,'(3(F12.6,X),3(I5,X))') rvec(:), xtarg(:)

          z_stripe(iz) = (R * conjg(R) )
          cubeExciton(iz,iy,ix) = z_stripe(iz)
!          do ispin = 1, nspn
!            z_stripe(iz) = z_stripe(iz) &
!                        + (Rspace_exciton( xtarg(3), xtarg(2), xtarg(1), ispin ) &
!                   * conjg(Rspace_exciton( xtarg(3), xtarg(2), xtarg(1), ispin )))
!          enddo
        enddo
        write(99,'(5(E13.6,X))') z_stripe(:)
      enddo
    enddo

    su = 0.0_DP
    do ix = 1, xmeshOut(1)
      do iy = 1, xmeshOut(2)
        do iz = 1, xmeshOut(3)
          su = su + cubeExciton(iz,iy,ix)
        enddo
      enddo
    enddo
    write(6,*) 'Total density: ', su 

    isoTarg = su / dble( product( xmeshOut(:) ) )
    isos(:,:) = -1
    do j = 1, 19
      isoMin = 0.0_DP
      isoMax = maxval( cubeExciton )
      foundMax = .false.
      suTarg = dble(j)/20.0_DP * su
      do k = 1, j - 1
        if( isos(1,k) .gt. 0.0_DP ) then
          isoMax = isos(1,k)
          foundMax = .true.
        endif
      enddo
      do k = 19, j + 1
        if( isos(1,k) .gt. 0 ) then
          isoMin = isos(1,k)
        endif
      enddo
      if( isos(1,j) .gt. 0 ) then
        isoTarg = isos(1,j)
      else
        isoTarg = 0.5_DP * ( isoMin + isoMax )
      endif
      write(6,'(I4,X,4E14.6,X,L1)') j, isoTarg, suTarg, isoMin, isoMax, foundMax

      do i = 1, 50
        su2 = 0.0_DP
        do ix = 1, xmeshOut(1)-1
          do iy = 1, xmeshOut(2)-1
            do iz = 1, xmeshOut(3)-1
              iter = 0
              su2 = su2 + interpHeavi( iter, isoTarg, cubeExciton(iz,iy,ix), cubeExciton(iz+1,iy,ix), &
                                     cubeExciton(iz,iy+1,ix), cubeExciton(iz+1,iy+1,ix), &
                                     cubeExciton(iz,iy,ix+1), cubeExciton(iz+1,iy,ix+1), &
                                     cubeExciton(iz,iy+1,ix+1), cubeExciton(iz+1,iy+1,ix+1) )
            enddo
          enddo
        enddo

        k = nint( 20.0_DP * su2 / su )
        write(6,*) 'Sub density: ', su2, isoTarg
        if( k .ne. j .and. k .ge. 1 .and. k .le. 19 ) then
          if( abs( su2 / su - dble(k)/20.0_DP ) .lt. abs( isos(2,k) - dble(k)/20.0_DP ) ) then
            isos(1,k) = isoTarg
            isos(2,k) = su2 / su
            write(6,*) 'Stored:', k, (su2 / su)
          endif
        endif
      ! exit early if close enough
        if( abs( su2 / su - dble(j)/20.0_DP ) .lt. 0.00001_DP ) then
          write(6,*) '   Early', j, i
!          isos(1,j) = isoTarg
!          isos(2,j) = su2 / su
          exit
        endif
        if( su2 .gt. suTarg ) then
         isoMin = isoTarg
        else
          foundMax = .true.
          isoMax = isoTarg
        endif
        if( foundMax ) then
          isoTarg = 0.5_DP * ( isoMin + isoMax )
        else
          isoTarg = isoTarg * 2.0_DP
        endif
      enddo
      isos(1,j) = isoTarg
      isos(2,j) = su2 / su
      write(6,*) j, isos(:,j)
    enddo

    do j = 1, 19
      write(6,*) dble(j)/20.0_DP, isos(:,j)
    enddo
    deallocate( cubeExciton )

  else
  write(99,'(I5,3(F12.6))') natom*product(Rmesh(:)), 0.0_dp, 0.0_dp, 0.0_dp
  do ix = 1, 3
    x_count = Rmesh( ix ) * xmesh( ix )
    write(99,'(I5,3(F12.6))') x_count, (Rmesh(:) * avecs(:,ix))/dble(x_count)
  enddo

  do iRx = Rstart(1), Rstart(1) + Rmesh(1) - 1
    do iRy = Rstart(2), Rstart(2) + Rmesh(2) - 1
      do iRz = Rstart(3), Rstart(3) + Rmesh(3) - 1
        do i = 1, natom 
!          call get_atom_number( elname( i ), ix )
          call get_atom_number( ix, elname( i ) )
          tau(:) = atom_loc(:,i)
          tau( : ) = tau( : ) + iRx * avecs(:, 1 )
          tau( : ) = tau( : ) + iRy * avecs(:, 2 )
          tau( : ) = tau( : ) + iRz * avecs(:, 3 )

          write(99,'(I5,4(F12.6))') ix, 0.0, tau(:)
!          write(99,'(I5,4(F12.6))') ix, 0.0, atom_loc(:,i)
        enddo
      enddo
    enddo
  enddo


  ! Need Z to be fast axis

  allocate( z_stripe( xmesh(3) * Rmesh(3) ) )

  select case( bloch_selector )

  case( 0, 1 )
  Riter = 0
  do iRx = 1, Rmesh(1)
    do ix = 1, xmesh(1)

      do iRy = 1, Rmesh(2)
        do iy = 1, xmesh(2)
          Xiter = ( iy - 1 ) * xmesh(3) + ( ix - 1 ) * xmesh(3) * xmesh(2)


          izz = 0
          do iRz = 1, Rmesh(3)
            Riter = iRz + ( iRy - 1 ) * Rmesh(3) + ( iRx - 1 ) * Rmesh(3) * Rmesh(2)

            do iz = 1, xmesh(3)
              izz = izz + 1
              z_stripe( izz ) = 0.0_DP
              do ispin = 1, nspn
!                z_stripe( izz ) = z_stripe( izz ) &
!                                + (Rspace_exciton( Xiter + iz, Riter, ispin )) &
!                                 *conjg(Rspace_exciton( Xiter + iz, Riter, ispin  ))
              enddo
            enddo
          enddo

          write(99,'(5(E13.6,X))') z_stripe(:)
            
        enddo
      enddo

    enddo
  enddo              

  case( 2, 3 )
    Riter = 0
!    do iRx = 1, Rmesh(1)
      do ix = 1, xmesh(1)*Rmesh(1)

!        do iRy = 1, Rmesh(2)
          do iy = 1, xmesh(2)*Rmesh(2)
!            Xiter = ( iy - 1 ) * xmesh(3) + ( ix - 1 ) * xmesh(3) * xmesh(2)


!            izz = 0
!            do iRz = 1, Rmesh(3)
!              Riter = iRz + ( iRy - 1 ) * Rmesh(3) + ( iRx - 1 ) * Rmesh(3) * Rmesh(2)

              do iz = 1, xmesh(3)*Rmesh(3)
                z_stripe( iz ) = 0.0_DP
                do ispin = 1, nspn
                  z_stripe( iz ) = z_stripe( iz ) &
                                 + (Rspace_exciton( iz, iy, ix, ispin )) &
                                   *conjg(Rspace_exciton( iz, iy, ix, ispin  ))
                enddo

!                Xiter = ix + (iy - 1 ) * xmesh(1) + (iz-1) * xmesh(1)*xmesh(2)
!                izz = izz + 1
!                z_stripe( izz ) = 0.0_DP
!                do ispin = 1, nspn
!                  z_stripe( izz ) = z_stripe( izz ) &
!                                  + (Rspace_exciton( Xiter, Riter, ispin )) &
!                                   *conjg(Rspace_exciton( Xiter, Riter, ispin  ))
!                enddo
              enddo
!            enddo

            write(99,'(5(E13.6,X))') z_stripe(:)

          enddo
        enddo

!      enddo
!    enddo


  end select
  
  endif
  close( 99 )

  deallocate( Rspace_exciton, z_stripe, ibeg )


  111 continue

  contains
  subroutine invA( a, b )
    real(DP), intent( in ) :: a(3,3)
    real(DP), intent( out ) :: b(3,3)
    real(DP) :: det

    det = a(1,1) * ( a(2,2) * a(3,3) - a(3,2) * a(2,3) ) &
        - a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

    b(1,1) = ( a(2,2) * a(3,3) - a(3,2) * a(2,3) ) / det
    b(1,2) = ( a(1,3) * a(3,2) - a(1,2) * a(3,3) ) / det
    b(1,3) = ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

    b(2,1) = ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) / det
    b(2,2) = ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
    b(2,3) = ( a(2,1) * a(1,3) - a(1,1) * a(2,3) ) / det

    b(3,1) = ( a(2,1) * a(3,2) - a(3,1) * a(2,2) ) / det
    b(3,2) = ( a(3,1) * a(1,2) - a(1,1) * a(3,2) ) / det
    b(3,3) = ( a(1,1) * a(2,2) - a(2,1) * a(1,2) ) / det

  end subroutine invA

  recursive function interpHeavi( iii, a, c000, c001, c010, c011, c100, c101, c110, c111 ) result( val )
    implicit none
    integer, intent( inout ) :: iii
    real(DP), intent( in ) :: a
    real(DP), intent( in ) :: c000, c001, c010, c011, c100, c101, c110, c111
    real(DP) :: val

!    real(DP) :: d000, d001, d010, d011, d100, d101, d110, d111
    real(DP) :: d(3,3,3)

    integer :: iix, iiy, iiz

    if( c000 .ge. a .and. c001 .ge. a .and. c010 .ge. a .and. c011 .ge. a .and. &
        c100 .ge. a .and. c101 .ge. a .and. c110 .ge. a .and. c111 .ge. a ) then
      val = 0.125_DP * ( c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111 )
    elseif( c000 .lt. a .and. c001 .lt. a .and. c010 .lt. a .and. c011 .lt. a .and. &
            c100 .lt. a .and. c101 .lt. a .and. c110 .lt. a .and. c111 .lt. a ) then
      val = 0.0_DP
    else
      if( iii .gt. 8 ) then
        val  = 0.0_DP
      else
        iii = iii + 1
        d(1,1,1) = c000
        d(2,1,1) = 0.5_DP * ( c000 + c100 )
        d(3,1,1) = c100
  
        d(1,2,1) = 0.5_DP * ( c000 + c010 )
        d(2,2,1) = 0.25_DP * ( c000 + c010 + c100 + c110 )
        d(3,2,1) = 0.5_DP * (c100+c110)

        d(1,3,1) = c010
        d(2,3,1) = 0.5_DP * (c010 + c110)
        d(3,3,1) = c110

        d(1,1,2) = 0.5_DP * ( c000 + c001 )
        d(2,1,2) = 0.25_DP * (c000 + c001 + c100 + c101 )
        d(3,1,2) = 0.5_DP * (c100 + c101 )

        d(1,2,2) = 0.25_DP * (c000 + c001 + c010 +c011 )
        d(2,2,2) = 0.125_DP * ( c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111 )
        d(3,2,2) = 0.25_DP * ( c100 + c101 + c110 + c111 )

        d(1,3,2) = 0.5_DP * ( c010 + c011 )
        d(2,3,2) = 0.25_DP * ( c010 + c011 + c110 + c111 )
        d(3,3,2) = 0.5_DP * ( c110 + c111 )
    
        d(1,1,3) = c001
        d(2,1,3) = 0.5_DP * ( c001 + c101 )
        d(3,1,3) = c101

        d(1,2,3) = 0.5_DP * ( c001 + c011 )
        d(2,2,3) = 0.25_DP * ( c001 + c011 + c101 + c111 )
        d(3,2,3) = 0.5_DP * ( c101 + c111 )

        d(1,3,3) = c011
        d(2,3,3) = 0.5_DP * ( c011 + c111 )
        d(3,3,3) = c111

        val = 0.0_DP
        do iiz = 1, 2
          do iiy = 1, 2
            do iix = 1,2
              val = val + 0.125_DP * interpHeavi( iii, a, d(iix,iiy,iiz), d(iix,iiy,iiz+1), d(iix,iiy+1,iiz), &
                                                  d(iix,iiy+1,iiz+1), d(iix+1,iiy,iiz), d(iix+1,iiy,iiz+1), &
                                                  d(iix+1,iiy+1,iiz), d(iix+1,iiy+1,iiz+1) )
            enddo
          enddo
        enddo
        
      endif

    endif

  end function interpHeavi
    

end program OCEAN_exciton_plot
