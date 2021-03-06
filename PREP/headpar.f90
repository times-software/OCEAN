! Copyright (C) 2015, 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine headerpar( filenum, fermie, maxnpw, maxband, nsppol, nspinor, nkpt, kout, ierr )
  implicit none
  
  integer, intent( in )  :: filenum, kout
  integer, intent( out ) ::  maxnpw, maxband, nsppol, nspinor, nkpt
  integer, intent( inout ) :: ierr
  real(kind=kind(1.0d0)), intent( out ) :: fermie

  character(len=6) :: codvsn
  integer :: fform, headform

  ierr = 0
  read(filenum) codvsn,headform,fform

  if( headform < 80 ) then
    call headerpar_53(filenum, fermie, maxnpw, maxband, nsppol, nspinor, nkpt,kout)
  else
    call headerpar_80(filenum, fermie, maxnpw, maxband, nsppol, nspinor, nkpt, ierr )
  endif
end subroutine headerpar

subroutine headerpar_80(filenum, fermie, maxnpw, maxband, nsppol, nspinor, nkpt, ierr )
  implicit none
  integer, intent( in )  :: filenum
  integer, intent( out ) ::  maxnpw, maxband, nsppol, nspinor, nkpt
  integer, intent( inout ) :: ierr
  real(kind=kind(1.0d0)), intent( out ) :: fermie
  !
  integer :: bndtot, date, intxc, ixc, natom, ngfft(3), nsym, npsp, ntypat, occopt, pertcase, & 
             usepaw, nspden, i
  real(kind=kind(1.0d0)) :: residm, etot
  real(kind=kind(1.0d0)), allocatable :: xred(:,:)
  integer, allocatable :: istwfk(:), nband(:), npwarr(:)

  read(filenum) bndtot, date, intxc, ixc, natom, ngfft, nkpt, nspden, nspinor, nsppol, nsym, & 
                npsp, ntypat, occopt, pertcase, usepaw

  allocate( istwfk(nkpt), nband(nkpt*nsppol), npwarr(nkpt) )
  ! istwfk etc
  read(filenum) istwfk, nband, npwarr
  maxband = maxval( nband )
  maxnpw = maxval( npwarr )
  
  deallocate( istwfk, nband, npwarr )
  
  
  allocate( xred( 3, natom ) )
  read(filenum) residm, xred, etot, fermie
  deallocate( xred )

  ! kptop ... shiftl
  read(filenum)

  !psps
  do i = 1, npsp
    read(filenum)
  enddo

  if( usepaw .eq. 1 ) then
    ierr = -1
  endif

  return

end subroutine headerpar_80

subroutine headerpar_53(filenum, fermie, maxnpw, maxband, nsppol, nspinor, nkpt,kout)
  implicit none

  integer :: filenum, maxnpw, maxband,i,kout

  character(len=6) :: codvsn
  integer :: headform,fform
  integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,             &
  nspden,nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw
  double precision :: ecut,ecutdg,ecutsm,ecut_eff,qptn(3),          &
  rprimd(3,3),stmbias,tphysel,tsmear
  integer, allocatable :: istwfk(:),nband(:),npwarr(:),so_typat(:), &
  symafm(:),symrel(:,:,:),typat(:)
  double precision, allocatable :: kpt(:,:),occ(:),tnons(:,:),      &
  znucltypat(:),xred(:,:) 
!      double precision, allocatable, intent (inout) :: kptlist(:,:)
  character(len=132) :: title
  double precision :: znuclpsp,zionpsp
  integer :: pspso,pspdat,pspcod,pspxc,npsp,ipsp
  double precision :: residm,etotal,fermie


!  read(filenum) codvsn,headform,fform
  read(filenum) bantot,date,intxc,ixc,natom,ngfft(1:3),             &
 & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,     &
 & usepaw,ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),    &
 & stmbias,tphysel,tsmear

  allocate(istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),            &
 & so_typat(ntypat),symafm(nsym),symrel(3,3,nsym),typat(natom))

  allocate(kpt(3,nkpt))

  allocate(occ(bantot),tnons(3,nsym),znucltypat(ntypat),            &
 &  xred(3,natom))

  read(filenum) istwfk(1:nkpt),nband(1:nkpt*nsppol),                &
 & npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym),                &
 & symrel(1:3,1:3,1:nsym),typat(1:natom),kpt(1:3,1:nkpt),       &
 & occ(1:bantot),tnons(1:3,1:nsym),znucltypat(1:ntypat)
  do ipsp=1,npsp
! (npsp lines, 1 far each pseudopotential ; npsp=ntypat, 
!  except if alchemical pseudo-atoms)
   read(filenum) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc
  enddo
!(final record: residm, coordinates, total energy, Fermi energy)
  read(filenum) residm,xred(1:3,1:natom),etotal,fermie

!      write(6,*) dummy
!      if ( density .eq. 1 ) then
!       fermie= dummy
!      elseif (dummy .gt. fermie ) then
!       fermie=dummy
!      endif 

  maxnpw = maxval(npwarr)
  maxband = maxval(nband)

  if (kout .ne. 0 ) then
   do i=1,nkpt
    write(kout,*)kpt(1,i),kpt(2,i),kpt(3,i)
   enddo
  endif

  deallocate(istwfk,nband,npwarr,so_typat,symafm,symrel,typat)
  deallocate(occ,tnons,znucltypat,xred,kpt)


end subroutine headerpar_53
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
