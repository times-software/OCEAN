! program makecoords
!
      program makecoords
!
      implicit none
      character(len=2), allocatable :: satom(:)
      integer :: i
      integer :: ibrav
      integer :: natoms, numtyp
      integer, allocatable    :: typat(:), znum(:), zatom(:)
      real(kind=kind(1.d0))              :: rscale(3)
      real(kind=kind(1.d0)), allocatable :: pos(:,:)
!
!
      open(unit=99,file='natoms',form='formatted',status='old')
      read(99,*) natoms
      close(99)
!
      open(unit=99,file='rscale',form='formatted',status='old')
      read(99,*) rscale(:)
      close(99)
!
      open(unit=99,file='ibrav',form='formatted',status='old')
      read(99,*) ibrav
      close(99)
!
      allocate( pos(3,natoms) )
      open(unit=99,file='taulist',form='formatted',status='old')
      if ( ibrav .eq. 0 ) then
         do i = 1, natoms
            read(99,*) pos(:,i)
         end do
      else
         do i = 1, natoms
            read(99,*) pos(:,i)
         end do
      endif
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
      close(99)
!
      allocate( zatom(natoms), satom(natoms) )
      do i = 1, natoms
         zatom(i) = znum( typat(i) )
      end do
!
!
      open(unit=98,file='coords',form='formatted',status='unknown')
      do i = 1, natoms
         call getsymbol( zatom(i), satom(i) )
         write( 98, *) satom(i), pos(:,i)
      enddo
      close(98)

      deallocate( satom )
      deallocate( typat, zatom, znum )
      deallocate( pos )

      end program makecoords


!!! move this to a header file

subroutine getsymbol(zatom,satom)
  integer, intent(in) :: zatom
  character(len=2), intent(out) :: satom

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
      case (31)
        satom = "Ga"
      case (32)
        satom = "Ge"
      case (33)
        satom = "As"
      case (34)
        satom = "Se"
      case (35)
        satom = "Br"
      case (36)
        satom = "Kr"
      case (37)
        satom = "Rb"
      case (38)
        satom = "Sr"
      case (39)
        satom = "Y "
      case (40)
        satom = "Zr"
      case (41)
        satom = "Nb"
      case (42)
        satom = "Mo"
      case (43)
        satom = "Tc"
      case (44)
        satom = "Ru"
      case (45)
        satom = "Rh"
      case (46)
        satom = "Pd"
      case (47)
        satom = "Ag"
      case (48)
        satom = "Cd"
      case (82)
        satom = "Pb"

    end select


end subroutine getsymbol
