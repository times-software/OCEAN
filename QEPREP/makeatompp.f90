! program makeatompp
!
      program makeatompp
!
      implicit none
      character(len=2) , allocatable :: satom(:)
      character(len=7) , allocatable :: mass(:)
      character(len=99), allocatable :: ppname(:), ppline(:)
      integer :: i
      integer :: ntype
      integer, allocatable    :: znucl(:)
!
!
      open(unit=99,file='ntype',form='formatted',status='old')
      read(99,*) ntype
      close(99)
!
      allocate( znucl(ntype) )
      open(unit=99,file='znucl',form='formatted',status='old')
      read(99,*) znucl(:)
      close(99)
!
      allocate( ppname(ntype) )
      open(unit=99,file='pplist',form='formatted',status='old')
      do i = 1, ntype
         read(99,*) ppname(i)
      end do
      close(99)
!
! get symbol & mass, concatenate
!
      allocate( satom(ntype), mass(ntype), ppline(ntype) )
      do i = 1, ntype
         call getsymbol( znucl(i), satom(i) )
         call getmass  ( znucl(i), mass (i) )
         ppline(i) = satom(i) // '   ' // mass(i) // '   ' // trim(ppname(i)) &
            &        // '.UPF'
      enddo
!
      open(unit=98,file='atompp',form='formatted',status='unknown')
      do i = 1, ntype
         write(98,*) trim( ppline(i) )
      end do
      close(98)
!
      deallocate( satom, mass, ppline )
      deallocate( znucl, ppname )
!
      end program makeatompp


!!! we should move this to a periodic table header file

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


subroutine getmass(zatom,mass)
  integer, intent(in) :: zatom
  character(len=7), intent(out) :: mass

    select case( zatom )

      case ( 1)
        mass = " 1.0079"
      case ( 2)
        mass = " 4.0026"
      case ( 3)
        mass = " 6.9410"
      case ( 4)
        mass = " 9.0122"
      case ( 5)
        mass = "10.8110"
      case ( 6)
        mass = "12.0107"
      case ( 7)
        mass = "14.0067"
      case ( 8)
        mass = "15.9994"
      case ( 9)
        mass = "18.9984"
      case (10)
        mass = "20.1797"
      case (11)
        mass = "22.9900"
      case (12)
        mass = "24.3050"
      case (13)
        mass = "26.9820"
      case (14)
        mass = "28.0850"
      case (15)
        mass = "30.9740"
      case (16)
        mass = "32.0600"
      case (17)
        mass = "35.4500"
      case (18)
        mass = "39.9480"
      case (19)
        mass = "39.0980"
      case (20)
        mass = "40.0780"
      case (21)
        mass = "44.9560"
      case (22)
        mass = "47.8670"
      case (23)
        mass = "50.9420"
      case (24)
        mass = "51.9960"
      case (25)
        mass = "54.9380"
      case (26)
        mass = "55.8450"
      case (27)
        mass = "58.9330"
      case (28)
        mass = "58.6930"  
      case (29)
        mass = "63.5460"
      case (30)
        mass = "65.3800"
      case (31)
        mass = "69.7230"
      case (32)
        mass = "72.6300"
      case (33)
        mass = "74.9220"
      case (34)
        mass = "78.9600"
      case (35)
        mass = "79.9040"
      case (36)
        mass = "83.7980"
      case (37)
        mass = "85.4680"
      case (38)
        mass = "87.6200"
      case (39)
        mass = "88.9060"
      case (40)
        mass = "91.2240"
      case (41)
        mass = "92.906"
      case (42)
        mass = "95.96"
      case (43)
        mass = "97.91"
      case (44)
        mass = "101.07"
      case (45)
        mass = "102.91"
      case (46)
        mass = "106.42"
      case (47)
        mass = "107.87"
      case (48)
        mass = "112.41"
      case (82)
        mass = "207.2"

    end select


end subroutine getmass

