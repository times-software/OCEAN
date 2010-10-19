      subroutine headerpar(filenum, fermie, maxnpw, maxband, nsppol,    &
     & nspinor, nkpt,kout)
      implicit none

      integer :: filenum, maxnpw, maxband,i,kout,density

      character*6 :: codvsn
      integer :: headform,fform
      integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,             &
     & nspden,nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw
      double precision :: ecut,ecutdg,ecutsm,ecut_eff,qptn(3),          &
     & rprimd(3,3),stmbias,tphysel,tsmear
      integer, allocatable :: istwfk(:),nband(:),npwarr(:),so_typat(:), &
     & symafm(:),symrel(:,:,:),typat(:)
      double precision, allocatable :: kpt(:,:),occ(:),tnons(:,:),      &
     & znucltypat(:),xred(:,:) 
!      double precision, allocatable, intent (inout) :: kptlist(:,:)
      character*132 :: title
      double precision :: znuclpsp,zionpsp
      integer :: pspso,pspdat,pspcod,pspxc,npsp,ipsp
      double precision :: residm,etotal,fermie,dummy


      read(filenum) codvsn,headform,fform
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
 

      end 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
