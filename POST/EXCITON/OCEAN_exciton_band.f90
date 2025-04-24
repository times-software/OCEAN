! Copyright (C) 2025 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
program OCEAN_exciton_band
!  use periodic, only : get_atom_number
!  use ocean_interpolate
  
  implicit none
  integer, parameter :: DP = kind(1.0d0 )

  integer :: nband, brange(4), kmesh(3), nspn, nalpha, ZNL(3), nkpts, kpathLength
  integer :: icms, icml, ivms, i, j, ix, iix, iy, iiy, iz, iiz, x0, y0, z0, iik, iband, ik, ispin

  real(DP) :: su, k0(3), klen, bandE, fff(8), ff(4), f(2), deltax, deltay, deltaz, kpoint(3)

  complex(DP), allocatable :: exciton(:,:,:)
  real(DP), allocatable :: cond_exciton(:,:,:), kpathExciton(:,:,:), kpath(:,:)

  character(len=25) :: filname
  character(len=25) :: inbandfile
  character(len=128) :: outname
  real(DP), external :: DZNRM2


  open(unit=99,file='exciton_band.ipt',form='formatted',status='old')
  read(99,*) filname
  read(99,*) inbandfile
  read(99,*) outname
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

  open(unit=99,file='k0.ipt',form='formatted',status='old')
  read(99,*) k0(:)
  close(99)


  open(unit=99,file='nspin',form='formatted',status='old')
  read(99,*) nspn
  close(99)

  open(unit=99,file='ZNL',form='formatted',status='old')
  read(99,*) ZNL(:)
  close(99)
  !(2l+1)
  nalpha = (2*ZNL(3)+1) * 4

  open(unit=99,file='kpath.inp',form='formatted',status='old')
  read(99,*) kpathLength
  allocate( kpath(3,kpathLength), kpathExciton(nband,kpathLength,nspn) )
  kpathExciton(:,:,:) = 0.0_DP
  do i = 1, kpathLength
    read(99,*) kpath(:,i)
  enddo
  close(99)
  
  write(6,*) nband, nkpts, nalpha

  allocate( exciton( nband, nkpts, nalpha ), cond_exciton( nband, nkpts, nspn ) )
  write(6,*) filname
  open(unit=99,file=filname,form='unformatted',status='old')
  read(99) exciton
  close(99)

  ! Want to normalize excitonic wvfn   
  su = DZNRM2( nband*nkpts*nalpha, exciton, 1 )
!  su = 1.0_DP / su
  su = real(nkpts,DP) / su
  write(6,*) su

  cond_exciton(:,:,:) = 0.0_DP
  
!  do i = 1, nalpha
  i = 0
  do icms = 1, 2
    do icml = -ZNL(3), ZNL(3)
      do ivms = 1, 2
        i = i + 1
          cond_exciton(:,:,min(ivms,nspn)) = cond_exciton(:,:,min(ivms,nspn)) &
                                           + su * real( exciton(:,:,i) * conjg( exciton(:,:,i)),DP )
      enddo
    enddo 
  enddo

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
    ! Do this after delta!
    if( x0 .lt. 0 ) then
      x0 = x0 + kmesh(1)
    endif
    y0 = floor( kmesh(2)*kpoint(2)-k0(2) )
    deltay = kpoint(2) - (k0(2)+dble(y0))/dble(kmesh(2))
    if( y0 .lt. 0 ) then
      y0 = y0 + kmesh(2)
    endif
    z0 = floor( kmesh(3)*kpoint(3)-k0(3) )
    deltaz = kpoint(3) - (k0(3)+dble(z0))/dble(kmesh(3))
    if( z0 .lt. 0 ) then
      z0 = z0 + kmesh(3)
    endif
!    write(6,*) kpoint(:)
!    write(6,*) x0, y0, z0, deltax, deltay, deltaz
!    write(6,*) (k0(1)+dble(x0))/dble(kmesh(1)), (k0(2)+dble(y0))/dble(kmesh(2)), &
!                (k0(3)+dble(z0))/dble(kmesh(3))

    do iband = 1, nband
      i = 0
      do ix = 0, 1
        iix = x0 + ix
        if( iix .eq. kmesh(1) ) then
          iix = 0
        endif
        do iy = 0, 1
          iiy = y0 + iy
          if( iiy .eq. kmesh(2) ) then
            iiy = 0
          endif
          do iz = 0, 1
            iiz = z0 + iz
            if( iiz .eq. kmesh(3) ) then
              iiz = 0
            endif
            i = i + 1
            iik = 1 + iix*kmesh(2)*kmesh(3) + iiy*kmesh(3) + iiz
            fff(i) = cond_exciton(iband,iik,ispin)
          enddo
        enddo
      enddo
      ff(1) = fff(1) + (fff(2)-fff(1))*real(kmesh(3),DP)*deltaz
      ff(2) = fff(3) + (fff(4)-fff(3))*real(kmesh(3),DP)*deltaz
      ff(3) = fff(5) + (fff(6)-fff(5))*real(kmesh(3),DP)*deltaz
      ff(4) = fff(7) + (fff(8)-fff(7))*real(kmesh(3),DP)*deltaz
      f(1) = ff(1) + (ff(2)-ff(1))*real(kmesh(2),DP)*deltay
      f(2) = ff(3) + (ff(4)-ff(3))*real(kmesh(2),DP)*deltay
      kpathExciton(iband,ik,ispin) = f(1) + (f(2)-f(1))*real(kmesh(1),DP)*deltax
      if( iband .eq. 1 ) then
!        write(6,*) fff(:)
!        write(6,*) ff(:)
!        write(6,*) f(:)
        write(6,*) kpathExciton(iband,ik,ispin), sum(fff(:))*0.125_DP
      endif
    enddo
      

  enddo

  open(unit=98,file=inbandfile,form='formatted',status='old')
  open(unit=99,file=outname,form='formatted',status='unknown')
  do iband = 1, brange(3)-1
    do ik = 1, kpathLength
      read(98,*) klen, bandE
      write(99,*) klen, bandE, 0.0_DP
    enddo
    read(98,*)
    write(99,*) ''
  enddo
  do iband = 1, nband
    do ik = 1, kpathLength
      read(98,*) klen, bandE
      write(99,*) klen, bandE, kpathExciton(iband,ik,ispin)
    enddo
    read(98,*)
    write(99,*) ''
  enddo
  close(98)
  close(99)
  deallocate( kpathExciton, cond_exciton)

end program OCEAN_exciton_band
