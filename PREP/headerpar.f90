      subroutine headerpar(filenum, fermie, maxnpw, maxband, nsppol,    &
     & nspinor, nkpt)
      implicit none

      integer :: filenum, maxnpw, maxband

      character*6 :: codvsn
      integer :: headform,fform
      integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,             &
     & nspden,nspinor,nsppol,isppol,nsym,ntypat,occopt,pertcase,usepaw
      double precision :: ecut,ecutdg,ecutsm,ecut_eff,qptn(3),          &
     & rprimd(3,3),stmbias,tphysel,tsmear
      integer, allocatable :: istwfk(:),nband(:),npwarr(:),so_typat(:), &
     & symafm(:),symrel(:,:,:),typat(:)
      double precision, allocatable :: kpt(:,:),occ(:),tnons(:,:),      &
     & znucltypat(:),xred(:,:)
      character*132 :: title
      double precision :: znuclpsp,zionpsp
      integer :: pspso,pspdat,pspcod,pspxc
      double precision :: residm,etotal,fermie,dummy

      read(filenum) codvsn,headform,fform
      read(filenum) bantot,date,intxc,ixc,natom,ngfft(1:3),             &
     & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,     &
     & usepaw,ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),    &
     & stmbias,tphysel,tsmear

      allocate(istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),            &
     & so_typat(ntypat),symafm(nsym),symrel(3,3,nsym),typat(natom))
      allocate(kpt(3,nkpt),occ(bantot),tnons(3,nsym),znucltypat(ntypat),&
     & xred(3,natom))

      read(filenum) istwfk(1:nkpt),nband(1:nkpt*nsppol),                &
     & npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym),                &
     & symrel(1:3,1:3,1:nsym),typat(1:natom),kpt(1:3,1:nkpt),           &
     & occ(1:bantot),tnons(1:3,1:nsym),znucltypat(1:ntypat)
      do ipsp=1,npsp
! (npsp lines, 1 far each pseudopotential ; npsp=ntypat, 
!  except if alchemical pseudo-atoms)
       read(filenum) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc
      enddo
!(final record: residm, coordinates, total energy, Fermi energy)
      read(filenum) residm,xred(1:3,1:natom),etotal,fermie

      deallocate(istwfk,nband,npwarr,so_typat,symafm,symrel,typat)
      deallocate(kpt,occ,tnons,znucltypat,xred)
      
      maxnpw = maxval(npwarr)
      maxband = maxval(nband)

      end 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine grabwf(filenum, maxband, maxnpw, kg_unshift, kg_shift, &
     &  eigen_un, eigen_sh, occ_un, occ_sh, cg_un, cg_sh, occ_max,      &
     &  unocc_max, nband(2))

      integer :: filenum, maxband, maxnpw, master_iter, nband(2),       &
     &  iband,ii,ikpt,maxnpw,un_npw,sh_npw,nspinor, occ_max, unocc_max
      double precision :: eigen_un(maxband), eigen_sh(maxband),         &
     &  occ_un(maxband), occ_sh(maxband), cg_un(maxband,2*maxnpw),       &
     &  cg_sh(maxband,2*maxnpw)
      integer :: kg_unshift(3,maxnpw), kg_shift(3,maxnpw), band_sh
      integer :: unocc_min,unocc_max,occ_min,occ_max,bandtot,maxband


      read(filenum) un_npw,nspinor,nband(1) 
      read(filenum) kg_unshift(1:3,1:un_npw)
      read(filenum) eigen_un(1:nband(1)),occ(1:nband(1))
      do iband=1,nband(1)
       read(filenum) (cg_un(iband,ii),ii=1,2*un_npw)
      enddo
      if ( nband(1) .lt. occ_max ) then
       stop 'occ_max is less than nbands from abinit run'
      endif

      read(filenum) sh_npw,nspinor,nband(2)
!     write(6,*) sh_npw,nspinor,nband(ikpt*2)
      read(filenum) kg_shift(1:3,1:sh_npw)
      read(filenum) eigen_sh(1:nband(2)),occ(1:nband(2))
      do iband=1,nband(2)
       read(filenum) (cg_sh(iband,ii),ii=1,2*sh_npw)
      enddo
      if ( nband(2) .lt. unocc_max ) then
       stop 'unocc_max is less than nbands from abinit run'
      endif

      end 
      
