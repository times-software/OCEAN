! program makecoords
!
      program makecoords

      use periodic
!
      implicit none
      character*2, allocatable :: satom(:),zsymb(:)
      integer :: i,ityp
      integer :: natoms, numtyp,iostatus
      integer, allocatable    :: typat(:), znum(:), zatom(:)
      real(kind=kind(1.d0))              :: rscale(3)
      real(kind=kind(1.d0)), allocatable :: pos(:,:)

      logical :: have_zsymb
!
!
      write(6,*) " in makecoords"
!
      open(unit=99,file='natoms',form='formatted',status='old')
      read(99,*) natoms
      close(99)
!
      open(unit=99,file='rscale',form='formatted',status='old')
      read(99,*) rscale(:)
      close(99)
!
      allocate( pos(3,natoms) )
      open(unit=99,file='taulist',form='formatted',status='old')
      do i = 1, natoms
         read(99,*) pos(:,i)
!         pos(:,i) = pos(:,i) * rscale(:) * 0.529177249
      end do
      close(99)
!
!
      allocate( typat(natoms) )
      open(unit=99,file='typat',form='formatted',status='old')
      read(99,*) typat(:)
      close(99)
!
      numtyp = 0
      do i = 1, natoms
         if ( typat(i) .gt. numtyp ) numtyp = typat(i)
      end do
!
      allocate( znum(numtyp) )
      open(unit=99,file='znucl',form='formatted',status='old')
      read(99,*) znum(:)
!      do i = 1, numtyp
!         read(99,*) znum(i)
!         write(6,*) i, znum(i)
!      end do
      close(99)
!
      allocate( zatom(natoms), satom(natoms) )
      do i = 1, natoms
         zatom(i) = znum( typat(i) )
      end do
!
      allocate( zsymb(numtyp) )
      open(unit=99,file='zsymb',form='formatted',status='old')
      read(99,*,IOSTAT=iostatus) zsymb(:)
      close(99)
      if( iostatus .eq. 0 ) then 
        have_zsymb = .true.
        write(6,*) 'User has set ZSYMB'
      else
        have_zsymb = .false.
        write(6,*) 'Using default ZSYMB'
      endif
!
!
      open(unit=98,file='coords',form='formatted',status='unknown')
      do i = 1, natoms
         call getsymbol( zatom(i), satom(i) )
         ityp=typat(i)
!         if(trim(zsymb(ityp)) .ne. '') satom(i)=trim(zsymb(ityp))
          if( have_zsymb ) satom(i)=trim(zsymb(ityp))
         write( 98, '(a,x,f16.10,x,f16.10,x,f16.10)') satom(i), pos(:,i)
      enddo
      close(98)

      deallocate( satom )
      deallocate( typat, zatom, znum, zsymb )
      deallocate( pos )

      end program makecoords


subroutine getsymbol_old(zatom,satom)
  integer, intent(in) :: zatom
  character*2, intent(out) :: satom

    select case( zatom )

      case ( 1) 
        satom = "H "
      case ( 2) 
        satom = "He"
      case ( 3) 
        satom = "Li"
      case ( 4) 
        satom = "Be"
      case ( 5) 
        satom = "B "
      case ( 6) 
        satom = "C "
      case ( 7) 
        satom = "N "
      case ( 8) 
        satom = "O "
      case ( 9) 
        satom = "F "
      case (10) 
        satom = "Ne"
      case (11)
        satom = "Na"
      case (12)
        satom = "Mg"
      case (13)
        satom = "Al"
      case (14)
        satom = "Si"
      case (15)
        satom = "P "
      case (16)
        satom = "S "
      case (17)
        satom = "Cl"
      case (18)
        satom = "Ar"
      case (19)
        satom = "K "
      case (20)
        satom = "Ca"
      case (21)
        satom = "Sc"
      case (22)
        satom = "Ti"
      case (23)
        satom = "V "
      case (24)
        satom = "Cr"
      case (25)
        satom = "Mn"
      case (26)
        satom = "Fe"
      case (27)
        satom = "Co"
      case (28) 
        satom = "Ni"
      case (29)
        satom = "Cu"
      case (30)
        satom = "Zn"
      case (38)
        satom = "Sr"

    end select


end subroutine getsymbol_old



