! Copyright (C) 2016 - 2017, 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
program OCEAN_val_exciton_plot
  implicit none
  integer, parameter :: DP = kind(1.0d0 )
  complex(DP), allocatable :: exciton(:,:,:,:), plot_exciton(:,:), Rspace_exciton(:,:), &
                              u2(:,:), rk_exciton(:,:), cv_exciton(:,:,:),cvkex(:,:),pointwf(:)
  complex(DP) :: cphs
  real(DP), allocatable :: z_stripe( : ), xyz(:,:), atom_loc(:,:)
  real(DP) :: qinb(3), avecs(3,3), su, k0(3), qvec(3), Rvec(3), xphs, yphs, zphs, twopi, tau(3), ur, ui, ehcoor(3)
  integer :: Rmesh(3), kmesh(3), nband, nalpha, nkpts, NR, Riter, kiter, xmesh(3), nspn, iehcoor(3)
  integer :: ikx, iky, ikz, iRx, iRy, iRz, NX, i, ix, x_count, xiter, iy, iz, izz, bloch_selector
  integer :: brange(4), u2size, u2start, Rshift(3), natom, kiter_break, Rstart(3), idum(3)
  integer :: nvb, ncb, ehflag, ixctr
  character(len=25) :: filname
  character(len=128) :: outname
  character(len=2), allocatable :: elname(:)
  real(DP), external :: DZNRM2
  complex(DP), parameter :: one = 1.0_dp
  complex(DP), parameter :: zero = 0.0_dp
  twopi = 2.0_dp * 4.0_dp * atan( 1.0_dp )
  open(unit=99,file='exciton_plot.ipt',form='formatted',status='old')
  read(99,*) filname
  read(99,*) outname
  read(99,*) Rmesh(:)
  read(99,*) Rstart(:)
!  read(99,*) tau(:)
  tau(:) = 0.0_dp
  close(99)
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
  nspn = 1
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
  open(unit=99,file='xyz.wyck',form='formatted',status='old')
  read(99,*) natom
  allocate( xyz(3,natom), atom_loc(3,0:natom), elname( natom ) ) 
  do i = 1, natom
    read(99,*) elname( i ), xyz(:,i)
  enddo
  close(99)
  open(unit=99,file='bloch_selector',form='formatted',status='old')
  read(99,*) bloch_selector
  close(99)
  open(unit=99,file='ehflag.ipt',form='formatted',status='old')
  read(99,*) ehflag
  close(99)
  if(ehflag < 1) ehflag=1
  if(ehflag > 2) ehflag=2
  open(unit=99,file='ehcoor.ipt',form='formatted',status='old')
  read(99,*) ehcoor
  close(99)
  
  atom_loc = 0.0_DP
  do i = 1, natom
    do ix = 1, 3
      atom_loc(:,i) = atom_loc(:,i) + avecs(:,ix) * xyz(ix,i)
    enddo
  enddo
!JTV!
  nalpha = 1 !non-spin valence calcs
!  Rmesh(:) = 3
!  Rstart(:) = -1
  write(6,*) Rmesh(:)
  write(6,*) Rstart(:)
  NR = product( Rmesh(:) )
  nvb=brange(2)-brange(1)+1
  ncb=brange(4)-brange(3)+1
  allocate( exciton( ncb, nvb, nkpts, nalpha ), cv_exciton( ncb, nvb, nkpts ) )
  write(6,*) filname
  open(unit=99,file=filname,form='unformatted',status='old')
  read(99) exciton
  close(99)
  
  ! Want to normalize excitonic wvfn   
  write(6,*) "ncb, nvb, nkpts, nalpha:", ncb, nvb, nkpts, nalpha
  su = DZNRM2( ncb*nvb*nkpts*nalpha, exciton, 1 )
  su = 1.0_DP / su
  write(6,*) su
  cv_exciton(:,:,:) = 0.0_DP
  do i = 1, nalpha
    cv_exciton(:,:,:) = cv_exciton(:,:,:) + su * exciton(:,:,:,i)
  enddo
  deallocate(exciton) !Dont need this anymore
 
  u2size = brange(4)-brange(3)+brange(2)-brange(1)+2
  allocate( u2( NX, u2size ) )
  allocate( rk_exciton( NX, nkpts ) )
  rk_exciton(:,:) = 0.0d0
  !select point near centroid
  write (6,'(A,3f12.6)') "Reference point selected:", ehcoor
 ! do i=1,3
 !    ehcoor(i) = fraction( 1 + ehcoor(i)) ! ged rid of negative frac coor
 ! enddo
  
  do i=1,3
    iehcoor(i) = nint( xmesh(i)*ehcoor(i) ) + 1
    write(6,*) ehcoor(i), iehcoor(i)
    if( iehcoor(i) .ge. xmesh(i) ) iehcoor(i) = iehcoor(i) - xmesh(i)
    if( iehcoor(i) .lt. 0        ) iehcoor(i) = iehcoor(i) + xmesh(i)
     
    ehcoor(i) = dble(iehcoor(i)-1)/dble(xmesh(i))
    write(6,*) ehcoor(i), iehcoor(i)
  enddo
  do ix = 1, 3
    atom_loc(:,0) = atom_loc(:,0) + avecs(:,ix) * ehcoor(ix)
  enddo
  
!  ixctr = floor(xmesh(1)*ehcoor(1)) + ( floor(xmesh(2)*ehcoor(2)) - 1 ) * xmesh(1) & 
!        + ( floor(xmesh(3)*ehcoor(3)) - 1 ) * xmesh(1) * xmesh(2)
!  write(6,*) 'Selected nearest real space point index:', ixctr
  ixctr = iehcoor(1) + (iehcoor(2) - 1 ) * xmesh(1) &
        + ( iehcoor(3) - 1 ) * xmesh(1) * xmesh(2)
  write(6,*) 'Selected nearest real space point index:', ixctr
  ixctr = iehcoor(3) + (iehcoor(2) - 1 ) * xmesh(3) &
        + ( iehcoor(1) - 1 ) * xmesh(3) * xmesh(2)
  write(6,*) 'Selected nearest real space point index:', ixctr

  ixctr = iehcoor(1) + (iehcoor(2) - 1 ) * xmesh(1) &
        + ( iehcoor(3) - 1 ) * xmesh(1) * xmesh(2)
  write(6,*) 'Selected nearest real space point index:', ixctr

  ixctr = iehcoor(3) + (iehcoor(2) - 1 ) * xmesh(3) &
        + ( iehcoor(1) - 1 ) * xmesh(3) * xmesh(2)
  write(6,*) 'Selected nearest real space point index:', ixctr
  
  write(6,*) 'Opening u2'
  kiter_break = nkpts / 20
  select case( ehflag)
  case (1)   ! plot hole part
     
     allocate(plot_exciton(nvb, nkpts))
     plot_exciton(:,:)=0.0d0
     allocate(pointwf(ncb))
     
     u2start = brange(1)
     write(6,*) "u2size, u2start, nvb:", u2size, u2start, nvb
     
     select case ( bloch_selector )
     case(1)
!        write(6, *) "u2par.dat Not implemented"
!        STOP
        open(unit=99,file='u2par.dat',access='stream',status='old',form='unformatted' )
        do kiter = 1, nkpts
          if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
          read(99) u2
          pointwf(1:ncb)=u2( ixctr, u2size-ncb+1:u2size )  !get value near center
          !Contract electron part to specified point
          call ZGEMV('T', ncb, nvb, one, cv_exciton(1,1,kiter), ncb, pointwf(1), 1, one, plot_exciton(1, kiter), 1)
          !Convert hole part to real space           
          call ZGEMV( 'N', NX, nvb, one, u2(1,1), NX, plot_exciton( 1, kiter ), 1, &
                      one, rk_exciton( 1, kiter ), 1 )
        enddo
      close( 99 )
     case(0)
        open(unit=99,file='u2.dat',form='unformatted',status='old')

        if( .true. ) then
        write(6,*) iehcoor(:)
        kiter = 0
        do ikx = 0, kmesh(1)-1
          qvec(1) = qinb(1) + (k0(1) + dble(ikx))/dble(kmesh(1))
          xphs = twopi * dble( iehcoor(1)-1)/dble( xmesh(1) )  * qvec( 1 )
          do iky = 0, kmesh(2)-1
            qvec(2) = qinb(2) + (k0(2) + dble(iky))/dble(kmesh(2))
            yphs = twopi * dble( iehcoor(2)-1)/dble( xmesh(2) )  * qvec( 2 ) + xphs
            do ikz = 0, kmesh(3)-1
              kiter = kiter + 1
              qvec(3) = qinb(3) + (k0(3) + dble(ikz))/dble(kmesh(3))
              zphs = twopi * dble( iehcoor(3)-1)/dble( xmesh(3) ) * qvec( 3 ) + yphs
              if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
              do i = u2start, u2size
                 do ix =1, NX
                    read(99) idum(1:3), ur, ui
                    u2( ix, i ) = cmplx( ur, -ui, DP )
                 enddo
              enddo
              cphs = cmplx( cos( -zphs ), sin(-zphs) )
              ! The hole gets phase and then gets conjugated
              pointwf(1:ncb)= conjg(u2( ixctr, u2size-ncb+1:u2size ) * cphs)  !get value near center
              !Contract hole part to specified point
              call ZGEMV('T', ncb, nvb, one, cv_exciton(1,1,kiter), ncb, pointwf(1), 1, one, plot_exciton(1, kiter), 1)
              !Convert electron part to real space           
              call ZGEMV( 'N', NX, nvb, one, u2(1,1), NX, plot_exciton( 1, kiter ), 1, &
                      one, rk_exciton( 1, kiter ), 1 )
            enddo
          enddo
        enddo
        else
          do kiter = 1, nkpts
             if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
             do i = u2start, u2size
                do ix =1, NX
                   read(99) idum(1:3), ur, ui
                   u2( ix, i ) = cmplx( ur, ui, DP )
                enddo
             enddo
             pointwf(1:ncb)=u2( ixctr, u2size-ncb+1:u2size )  !get value near center
             
             !Contract electron part to specified point
             
             call ZGEMV('T', ncb, nvb, one, cv_exciton(1,1,kiter), ncb, pointwf(1), 1, one, plot_exciton(1, kiter), 1)
             
             !Convert hole part to real space           
             call ZGEMV( 'N', NX, nvb, one, u2(1,1), NX, plot_exciton( 1, kiter ), 1, &
                     one, rk_exciton( 1, kiter ), 1 )
             
          enddo
        endif
        close( 99 )
     case default
        stop
     end select
     qinb(:) = 0.0_DP
     write(6,*) 'Loaded'
  case(2)  !plot electron part
     
     allocate(plot_exciton(ncb, nkpts))
     plot_exciton(:,:)=0.0d0
     allocate(pointwf(nvb))
     
     u2start = brange(1)
     write(6,*) "u2size, u2start, ncb:", u2size, u2start, ncb
     
     select case ( bloch_selector )
     case(1)
        write(6, *) "u2par.dat Not implemented"
        STOP
     case(0)
        open(unit=99,file='u2.dat',form='unformatted',status='old')
        !!! This is the new section !!!
        if( .true. ) then 
        kiter = 0
        do ikx = 0, kmesh(1)-1
          qvec(1) = (k0(1) + dble(ikx))/dble(kmesh(1))
          xphs = twopi * dble( iehcoor(1)-1)/dble( xmesh(1) )  * qvec( 1 )
          do iky = 0, kmesh(2)-1
            qvec(2) = (k0(2) + dble(iky))/dble(kmesh(2))
            yphs = twopi * dble( iehcoor(2)-1)/dble( xmesh(2) )  * qvec( 2 ) + xphs
            do ikz = 0, kmesh(3)-1
              kiter = kiter + 1
              qvec(3) = (k0(3) + dble(ikz))/dble(kmesh(3))
              zphs = twopi * dble( iehcoor(3)-1)/dble( xmesh(3) ) * qvec( 3 ) + yphs
              if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
              do i = u2start, u2size
                 do ix =1, NX
                    read(99) idum(1:3), ur, ui
                    u2( ix, i ) = cmplx( ur, ui, DP )
                 enddo
              enddo
              cphs = cmplx( cos( zphs ), sin(zphs) )
              ! The hole gets phase and then gets conjugated
              pointwf(1:nvb)= conjg(u2( ixctr, 1:nvb ) * cphs)  !get value near center
              !Contract hole part to specified point
              call ZGEMV('N', ncb, nvb, one, cv_exciton(1,1,kiter), ncb, pointwf(1), 1, one, plot_exciton(1, kiter), 1)
              !Convert electron part to real space           
              call ZGEMV( 'N', NX, ncb, one, u2(1,u2size-ncb+1), NX, plot_exciton( 1, kiter ), 1, &
                      one, rk_exciton( 1, kiter ), 1 )
            enddo
          enddo
        enddo
        else
        do kiter = 1, nkpts
           if( mod( kiter, kiter_break ) .eq. 0 ) write(6,*) kiter
           do i = u2start, u2size
              do ix =1, NX
                 read(99) idum(1:3), ur, ui
                 u2( ix, i ) = cmplx( ur, ui, DP )
              enddo
           enddo
           pointwf(1:nvb)=u2( ixctr, 1:nvb )  !get value near center
           
           !Contract hole part to specified point
           call ZGEMV('N', ncb, nvb, one, cv_exciton(1,1,kiter), ncb, pointwf(1), 1, one, plot_exciton(1, kiter), 1) 
           !Convert electron part to real space           
           call ZGEMV( 'N', NX, ncb, one, u2(1,u2size-ncb+1), NX, plot_exciton( 1, kiter ), 1, &
                   one, rk_exciton( 1, kiter ), 1 )
           
        enddo
        endif
        close( 99 )
        
     case default
        stop
     end select
     
  case default
     stop
  end select
  write(6,*) 'Done loading u2'
  
  deallocate( pointwf, u2, plot_exciton, cv_exciton )
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
              cphs = cmplx( cos( -zphs ), sin( -zphs ) )
              rk_exciton( xiter, kiter ) = rk_exciton( xiter, kiter ) * cphs
              
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
! Do Fourier transform into Rspace_exciton
  allocate( Rspace_exciton( NX, NR ) )
  Rspace_exciton = 0.0_DP
! Using slow FT and 6 loops for clarity
  
  kiter = 0
!  qvec(1) = qinb(1) + k0(1)/dble(kmesh(1))
  do ikx = 0, kmesh(1)-1
    qvec(1) = qinb(1) + (k0(1) + dble(ikx))/dble(kmesh(1))
!    qvec(1) = (k0(1) + dble(ikx))/dble(kmesh(1))
!    qvec(1) = qinb(1) + dble(ikx)/dble(kmesh(1))
!    qvec(2) = qinb(2) + k0(2)/dble(kmesh(2))
    do iky = 0, kmesh(2)-1
      qvec(2) = qinb(2) + (k0(2) + dble(iky))/dble(kmesh(2))
!      qvec(2) = (k0(2) + dble(iky))/dble(kmesh(2))
!      qvec(2) = qinb(2) + dble(iky)/dble(kmesh(2))
!      qvec(3) = qinb(3) + k0(3)/dble(kmesh(3))
      do ikz = 0, kmesh(3)-1
        qvec(3) = qinb(3) + (k0(3) + dble(ikz))/dble(kmesh(3))
!        qvec(3) = (k0(3) + dble(ikz))/dble(kmesh(3))
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
              cphs = cmplx( cos( twopi * -zphs ), sin( twopi * -zphs ) )
!JTV PHASE!!
              ! phase info is wrong here
              do xiter = 1, NX
                Rspace_exciton( xiter, Riter ) = Rspace_exciton( xiter, Riter ) + cphs * rk_exciton( xiter, kiter )
              enddo
            enddo
          enddo
        enddo
!        qvec(3) = qvec(3) + 1.0_dp / dble( kmesh(3) )
      enddo
!      qvec(2) = qvec(2) + 1.0_dp / dble( kmesh(2) )
    enddo
!    qvec(1) = qvec(1) + 1.0_dp / dble( kmesh(1) )
  enddo
  deallocate( rk_exciton )

  su = DZNRM2( NX*Rmesh(1)*Rmesh(2)*Rmesh(3), Rspace_exciton, 1 )
  write(6,*) su
  su = 1.0_DP / su
  call ZDSCAL( NX*Rmesh(1)*Rmesh(2)*Rmesh(3), su, Rspace_exciton, 1 )



  do ix = 1, 3
    if( Rmesh(ix) .gt. 1 ) then
      Rshift(ix) = Rmesh(ix)/2
    else
      Rshift(ix) = 0
    endif
    Rshift( ix ) = - Rstart( ix )
  enddo
  do i = 0, natom
    do ix = 1, 3
      atom_loc( :, i ) = atom_loc( :, i ) + Rshift( ix ) * avecs(:, ix )
    enddo
  enddo
  
  ! Now ready to write out 
!  open(unit=99,file='out.cube',form='formatted')
  
!  outname = trim(filname)//'.cube'
  open(unit=99,file=outname,form='formatted')
  write(99,*) "OCEAN exciton plot"
  write(99,*) "---"
  write(99,'(I5,3(F12.6))') 1+natom*product(Rmesh(:)), 0.0_dp, 0.0_dp, 0.0_dp
  do ix = 1, 3
    x_count = Rmesh( ix ) * xmesh( ix )
    write(99,'(I5,3(F12.6))') x_count, (Rmesh(:) * avecs(:,ix))/dble(x_count)
  enddo
  
  ! This is where I write out the location of the electron/hole that is frozen out
  ! The "atom" type is hardwired as Z=99, could be other things
  write(99,'(I5,4(F12.6))') 99, 0.0, atom_loc(:,0)
  do iRx = Rstart(1), Rstart(1) + Rmesh(1) - 1
    do iRy = Rstart(2), Rstart(2) + Rmesh(2) - 1
      do iRz = Rstart(3), Rstart(3) + Rmesh(3) - 1
        do i = 1, natom 
          call get_atom_number( elname( i ), ix )
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
              z_stripe( izz ) = (Rspace_exciton( Xiter + iz, Riter ))*conjg(Rspace_exciton( Xiter + iz, Riter ))
            enddo
          enddo
          write(99,'(5(E13.6,X))') z_stripe(:)
            
        enddo
      enddo
    enddo
  enddo              
  
  close( 99 )
  deallocate( Rspace_exciton, z_stripe )
  contains
  subroutine get_atom_number( elnam, elnum )
    implicit none
    character(len=2), intent( in ) :: elnam
    integer, intent( out ) :: elnum
    elnum = 1
    if( elnam .eq. 'Li' ) elnum = 3
    if( elnam .eq. 'N_' ) elnum = 7
    if( elnam .eq. 'O_' ) elnum = 8
    if( elnam .eq. 'S_' ) elnum = 16
    if( elnam .eq. 'Cd' ) elnum = 38
  end subroutine
end program OCEAN_val_exciton_plot
