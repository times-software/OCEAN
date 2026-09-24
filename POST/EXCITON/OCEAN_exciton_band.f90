! Copyright (C) 2025 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
program OCEAN_exciton_band
!  use periodic, only : get_atom_number
!  use ocean_interpolate
  use poly_interp_mod, only : interp_poly_cut
  
  implicit none
  integer, parameter :: DP = kind(1.0d0 )

  integer :: nband, brange(4), kmesh(3), nspn, nalpha, ZNL(3), nkpts, kpathLength, nbval, ierr
  
  integer :: icms, icml, ivms, i, j, ix, iix, iy, iiy, iz, iiz, x0, y0, z0, iik, iband, ik, &
             ispin, ialpha

  real(DP) :: su, k0(3), klen, bandE, fff(8), ff(4), f(2), deltax, deltay, deltaz, kpoint(3), &
              photonq(3)
  real(DP) :: ggg(8), gg(4), g(2)
  real(DP) :: efermi, qinb(3)

  complex(DP), allocatable :: exciton(:,:,:,:)
  real(DP), allocatable :: cond_exciton(:,:,:), kpathExciton(:,:,:), kpath(:,:), &
                           val_exciton(:,:,:), kpathValExciton(:,:,:)

  logical :: is_core, actualZero, have_fermi, is_xas

  character(len=25) :: filname
  character(len=25) :: inbandfile
  character(len=128) :: outname
  character(len=3) :: calc
  real(DP), external :: DZNRM2
  real(DP), parameter :: eps = 1.0d-12

  have_fermi = .true.

  open(unit=99,file='exciton_band.ipt',form='formatted',status='old')
  read(99,*) filname
  read(99,*) inbandfile
  read(99,*) outname
  read(99,*,IOSTAT=ierr) calc
  if( ierr .ne. 0 ) calc = '---'
  read(99,*,IOSTAT=ierr) efermi
  if( ierr .ne. 0 ) then
    have_fermi = .false.
    efermi = 0.0_DP
  endif
  close(99)


  select case (calc)
    case ( 'xes' )
      is_core = .true.
      is_xas = .false.
    case ( 'xas' )
      is_core = .true.
      is_xas = .true.
    case ( 'rxs', 'val' )
      is_core = .false.
      is_xas = .false.
    case default
      write(6,*) 'Unrecognized calc type, default to xas'
      is_core = .true.
      is_xas = .true.
  end select

  open(unit=99,file='nbuse.ipt',form='formatted',status='old')
  read(99,*) nband
  close(99)

  open(unit=99,file='brange.ipt',form='formatted',status='old')
  read(99,*) brange(1:4)
  close(99)
  nbval = brange(2)-brange(1)+1

  open(unit=99,file='kmesh.ipt',form='formatted',status='old')
  read(99,*) kmesh(:)
  close(99)
  nkpts = product( kmesh(:) )

  open(unit=99,file='k0.ipt',form='formatted',status='old')
  read(99,*) k0(:)
  close(99)

  photonq(:) = 0.0_DP
  if( is_core ) then
!    open(unit=99,file='photon_q',form='formatted',status='old')
!    read(99,*) photonq(:)
!    close(99)
  endif

  open(unit=99,file='nspin',form='formatted',status='old')
  read(99,*) nspn
  close(99)

  if( is_core ) then
    open(unit=99,file='ZNL',form='formatted',status='old')
    read(99,*) ZNL(:)
    close(99)
    !(2l+1)
    nalpha = (2*ZNL(3)+1) * 4
  elseif( calc .eq. 'rxs' ) then
    if( nspn .eq. 1 ) then
      open(unit=99,file='ZNL',form='formatted',status='old')
      read(99,*) ZNL(:)
      close(99)
      if( ZNL(3) .gt. 0 ) nspn = 2
    endif
    nalpha = nspn**2
  else
    nalpha = nspn**2
  endif

  if( .not. is_core ) then
    open(unit=99,file='qinunitsofbvectors.ipt',form='formatted',status='old')
    read(99,*) qinb(:)
    close( 99 )
  endif
    

  open(unit=99,file='kpath.inp',form='formatted',status='old')
  read(99,*) kpathLength
  allocate( kpath(3,kpathLength), kpathExciton(nband,kpathLength,nspn) )
  if( .not. is_core ) then
    allocate( kpathValExciton(nbval,kpathLength,nspn) )
    kpathValExciton(:,:,:) = 0.0_DP
  endif
  kpathExciton(:,:,:) = 0.0_DP
  do i = 1, kpathLength
    read(99,*) kpath(:,i)
  enddo
  close(99)
  
  write(6,*) nband, nkpts, nalpha

  if( is_core ) then
    allocate( exciton( nband, nkpts, nalpha, 1 ), cond_exciton( nband, nkpts, nspn ), &
              val_exciton(1,1,1) )
  else
    allocate( exciton( brange(4)-brange(3)+1, brange(2)-brange(1)+1, nkpts, nalpha ) )
    allocate( cond_exciton(1:brange(4)-brange(3)+1, nkpts, nspn ) )
    allocate( val_exciton(brange(2)-brange(1)+1, nkpts, nspn ) )
    write(6,*) (brange(4)-brange(3)+1), (brange(2)-brange(1)+1),nkpts,nalpha
  endif
  write(6,*) filname
  open(unit=99,file=filname,form='unformatted',status='old')
  read(99) exciton
  close(99)

  ! Want to normalize excitonic wvfn   
  if( is_core ) then
    su = DZNRM2( nband*nkpts*nalpha, exciton, 1 )
  else
    su = DZNRM2( (brange(4)-brange(3)+1)*(brange(2)-brange(1)+1)*nkpts*nalpha, exciton, 1 )
  endif
!  su = 1.0_DP / su
  su = real(nkpts,DP) / su
  write(6,*) su

  cond_exciton(:,:,:) = 0.0_DP
  val_exciton(:,:,:) = 0.0_DP
  
  if( is_core ) then
    i = 0
    do icms = 1, 2
      do icml = -ZNL(3), ZNL(3)
        do ivms = 1, 2
          i = i + 1
            cond_exciton(:,:,min(ivms,nspn)) = cond_exciton(:,:,min(ivms,nspn)) &
                                             + su * real( exciton(:,:,i,1) * conjg( exciton(:,:,i,1)),DP )
        enddo
      enddo 
    enddo
  else
!    if( nspn .ne. 1 ) then
!      write(6,*) 'Finish checking valence spin first!'
!      stop
!    endif
    ialpha = 0
    do icms = 1, nspn
      do ivms = 1, nspn
        ialpha = ialpha + 1
        do ik = 1, nkpts
          do i = 1, brange(2)-brange(1)+1
            do j = 1, brange(4)-brange(3)+1
              cond_exciton(j,ik,min(icms,nspn)) = cond_exciton(j,ik,min(icms,nspn)) &
                                                + su * real( exciton(j,i,ik,ialpha)* &
                                                       conjg(exciton(j,i,ik,ialpha)),DP)
              val_exciton(i,ik,min(ivms,nspn)) = val_exciton(i,ik,min(ivms,nspn)) &
                                               + su * real( exciton(j,i,ik,ialpha) &
                                                    *conjg(exciton(j,i,ik,ialpha)),DP)
!                                                 + su * sum( real(exciton(:,i,ik,ialpha) &
!                                                           *conjg(exciton(:,i,ik,ialpha)),DP))
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  ispin = 1
  do ik = 1, kpathLength
    kpoint(:) = kpath(:,ik)
    do j = 1, 3
      do while( kpoint(j) .lt. 0.0_DP )
        kpoint(j) = kpoint(j)+1.0_DP
      enddo
      do while( kpoint(j) .ge. 1.0_DP )
        kpoint(j) = kpoint(j)-1.0_DP
      enddo
    enddo
    
    x0 = floor( kmesh(1)*kpoint(1)-k0(1) )
    deltax = kpoint(1) - (k0(1)+dble(x0))/dble(kmesh(1))
    deltax = deltax*dble(kmesh(1))
    ! Do this after delta!
    if( x0 .lt. 0 ) then
      x0 = x0 + kmesh(1)
    endif
    y0 = floor( kmesh(2)*kpoint(2)-k0(2) )
    deltay = kpoint(2) - (k0(2)+dble(y0))/dble(kmesh(2))
    deltay = deltay*dble(kmesh(2))
    if( y0 .lt. 0 ) then
      y0 = y0 + kmesh(2)
    endif
    z0 = floor( kmesh(3)*kpoint(3)-k0(3) )
    deltaz = kpoint(3) - (k0(3)+dble(z0))/dble(kmesh(3))
    deltaz = deltaz*dble(kmesh(3))
    if( z0 .lt. 0 ) then
      z0 = z0 + kmesh(3)
    endif
    write(6,*) kpoint(:)
    write(6,*) x0, y0, z0, deltax, deltay, deltaz
!    write(6,*) (k0(1)+dble(x0))/dble(kmesh(1)), (k0(2)+dble(y0))/dble(kmesh(2)), &
!                (k0(3)+dble(z0))/dble(kmesh(3))
!    write(6,*) (k0(1)+dble(x0+1))/dble(kmesh(1)), (k0(2)+dble(y0+1))/dble(kmesh(2)), &
!                (k0(3)+dble(z0+1))/dble(kmesh(3))

    do iband = 1, nband
      i = 0
      actualZero = .false.
      do ix = 0, 1
        iix = x0 + ix
        if( iix .ge. kmesh(1) ) then
!          if( iix .ne. kmesh(1)) write(6,*) '---', iix, kmesh(1)
          iix = iix - kmesh(1)
        endif
        do iy = 0, 1
          iiy = y0 + iy
          if( iiy .ge. kmesh(2) ) then
!            if( iiy .ne. kmesh(2)) write(6,*) '---', iiy, kmesh(2)
            iiy = iiy - kmesh(2)
          endif
          do iz = 0, 1
            iiz = z0 + iz
            if( iiz .ge. kmesh(3) ) then
!              if( iiz .ne. kmesh(3)) write(6,*) '---', iiz, kmesh(3)
              iiz = iiz - kmesh(3)
            endif
            i = i + 1
            iik = 1 + iix*kmesh(2)*kmesh(3) + iiy*kmesh(3) + iiz
            fff(i) = cond_exciton(iband,iik,ispin)
            if( fff(i) .lt. eps ) then 
              actualZero = .true.
!              write(6,*) ix, iy, iz, iband, fff(i)
            endif
          enddo
        enddo
      enddo
      ff(1) = fff(1) + (fff(2)-fff(1))*deltaz
      ff(2) = fff(3) + (fff(4)-fff(3))*deltaz
      ff(3) = fff(5) + (fff(6)-fff(5))*deltaz
      ff(4) = fff(7) + (fff(8)-fff(7))*deltaz
      f(1) = ff(1) + (ff(2)-ff(1))*deltay
      f(2) = ff(3) + (ff(4)-ff(3))*deltay
      kpathExciton(iband,ik,ispin) = f(1) + (f(2)-f(1))*deltax
!      write(6,*) kpathExciton(iband,ik,ispin), interp_poly_cut(fff, u=deltax, v=deltay, w=deltaz, cutoff=-eps)
      if( actualZero ) then !kpathExciton(iband,ik,ispin) = 0.0_DP
        write(6,*) deltax, deltay, deltaz
        write(6,*) fff(:)
        write(6,*) kpathExciton(iband,ik,ispin), interp_poly_cut(fff, u=deltax, v=deltay, w=deltaz, cutoff=eps)
        kpathExciton(iband,ik,ispin) = interp_poly_cut(fff, u=deltax, v=deltay, w=deltaz, cutoff=eps)
      endif
      if( iband .eq. 1 ) then
!        write(6,*) fff(:)
!        write(6,*) ff(:)
!        write(6,*) f(:)
!        write(6,*) kpathExciton(iband,ik,ispin), sum(fff(:))*0.125_DP
      endif
    enddo
  enddo

  if( .not. is_core) then
    k0(:) = k0(:) - real(kmesh(:),DP)*qinb(:)
    do ik = 1, kpathLength
      kpoint(:) = kpath(:,ik)
      do j = 1, 3
        do while( kpoint(j) .lt. 0.0_DP )
          kpoint(j) = kpoint(j)+1.0_DP
        enddo
        do while( kpoint(j) .ge. 1.0_DP )
          kpoint(j) = kpoint(j)-1.0_DP
        enddo
      enddo
    x0 = floor( kmesh(1)*kpoint(1)-k0(1) )
    deltax = kpoint(1) - (k0(1)+dble(x0))/dble(kmesh(1))
    deltax = deltax*dble(kmesh(1))
    ! Do this after delta!
    do while( x0 .lt. 0 ) 
      x0 = x0 + kmesh(1)
    enddo
    do while( x0 .ge. kmesh(1) )
      x0 = x0 - kmesh(1)
    enddo
    y0 = floor( kmesh(2)*kpoint(2)-k0(2) )
    deltay = kpoint(2) - (k0(2)+dble(y0))/dble(kmesh(2))
    deltay = deltay*dble(kmesh(2))
    do while( y0 .lt. 0 ) 
      y0 = y0 + kmesh(2)
    enddo
    do while( y0 .ge. kmesh(2) )
      y0 = y0 - kmesh(2)
    enddo
    z0 = floor( kmesh(3)*kpoint(3)-k0(3) )
    deltaz = kpoint(3) - (k0(3)+dble(z0))/dble(kmesh(3))
    deltaz = deltaz*dble(kmesh(3))
    do while( z0 .lt. 0 ) 
      z0 = z0 + kmesh(3)
    enddo
    do while( z0 .ge. kmesh(3) )
      z0 = z0 - kmesh(3)
    enddo

      do iband = 1, nbval
        i = 0
        actualZero = .false.
        do ix = 0, 1
          iix = x0 + ix
          if( iix .ge. kmesh(1) ) then
            iix = iix - kmesh(1)
          endif
          do iy = 0, 1
            iiy = y0 + iy
            if( iiy .ge. kmesh(2) ) then
              iiy = iiy - kmesh(2)
            endif
            do iz = 0, 1
              iiz = z0 + iz
              if( iiz .ge. kmesh(3) ) then
                iiz = iiz - kmesh(3)
              endif
              i = i + 1
              iik = 1 + iix*kmesh(2)*kmesh(3) + iiy*kmesh(3) + iiz
              fff(i) = val_exciton(iband,iik,ispin)
              if( fff(i) .lt. eps ) then
                actualZero = .true.
              endif
            enddo 
          enddo 
        enddo 
        ff(1) = fff(1) + (fff(2)-fff(1))*deltaz
        ff(2) = fff(3) + (fff(4)-fff(3))*deltaz
        ff(3) = fff(5) + (fff(6)-fff(5))*deltaz
        ff(4) = fff(7) + (fff(8)-fff(7))*deltaz
        f(1) = ff(1) + (ff(2)-ff(1))*deltay
        f(2) = ff(3) + (ff(4)-ff(3))*deltay
        kpathValExciton(iband,ik,ispin) = f(1) + (f(2)-f(1))*deltax
        if( actualZero ) then
          kpathValExciton(iband,ik,ispin) = interp_poly_cut(fff, u=deltax, v=deltay, w=deltaz, cutoff=eps)
        endif
      enddo
  
    enddo
  endif

  if( is_core ) then
    open(unit=98,file=inbandfile,form='formatted',status='old')
    open(unit=99,file=outname,form='formatted',status='unknown')
    do iband = 1, brange(3)-1
      do ik = 1, kpathLength
        read(98,*) klen, bandE
        write(99,*) klen, bandE-eFermi, 0.0_DP
      enddo
      read(98,*)
      write(99,*) ''
    enddo
    do iband = 1, nband
      do ik = 1, kpathLength
        read(98,*) klen, bandE
        if( have_fermi ) then 
          if( ( is_xas .and. bandE .lt. efermi ) .or. &
              ( ( .not. is_xas ) .and. bandE .gt. efermi ) ) kpathExciton(iband,ik,ispin) = 0.0_DP
        endif
        write(99,*) klen, bandE-eFermi, kpathExciton(iband,ik,ispin)
      enddo
      read(98,*)
      write(99,*) ''
    enddo
    close(98)
    close(99)
  else
    open(unit=98,file=inbandfile,form='formatted',status='old')
    open(unit=99,file=outname,form='formatted',status='unknown')
    do iband = 1, brange(1)-1
      do ik = 1, kpathLength
        read(98,*) klen, bandE
        write(99,*) klen, bandE-eFermi, 0.0_DP, 0.0_DP
      enddo
      read(98,*)
      write(99,*) ''
    enddo
    do iband = brange(1), brange(3)-1
      do ik = 1, kpathLength
        read(98,*) klen, bandE
        if( have_fermi .and. ( bandE .gt. efermi ) ) kpathValExciton(iband-brange(1)+1,ik,ispin) = 0.0_DP
        write(99,*) klen, bandE-eFermi, kpathValExciton(iband-brange(1)+1,ik,ispin), 0.0_DP
      enddo
      read(98,*)
      write(99,*) ''
    enddo
    do iband = brange(3), brange(2)
      do ik = 1, kpathLength
        read(98,*) klen, bandE
        if( have_fermi ) then
          if( bandE .gt. efermi ) kpathValExciton(iband-brange(1)+1,ik,ispin) = 0.0_DP
          if( bandE .lt. efermi ) kpathExciton(iband-brange(3)+1,ik,ispin) = 0.0_DP
        endif
        write(99,*) klen, bandE-eFermi, kpathValExciton(iband-brange(1)+1,ik,ispin), &
                                 kpathExciton(iband-brange(3)+1,ik,ispin)
      enddo
      read(98,*)
      write(99,*) ''
    enddo
    do iband = brange(2)+1,brange(4) !brange(3)+nband-1
      do ik = 1, kpathLength
        read(98,*) klen, bandE
        if( have_fermi ) then
          if( bandE .lt. efermi ) kpathExciton(iband-brange(3)+1,ik,ispin) = 0.0_DP
        endif
        write(99,*) klen, bandE-eFermi, 0.0_DP, kpathExciton(iband-brange(3)+1,ik,ispin)
      enddo
      read(98,*)
      write(99,*) ''
    enddo
    close(98)
    close(99)
  endif
  deallocate( kpathExciton, cond_exciton)

end program OCEAN_exciton_band
