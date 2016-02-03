! program makeatompp
!
      program makeatompp

      use periodic
!
      implicit none
      character(len=3) , allocatable :: satom(:), zsymb(:)
      character(len=7) , allocatable :: mass(:)
      character(len=99), allocatable :: ppname(:), ppline(:)
      integer :: i, iostatus
      integer :: ntype
      integer, allocatable    :: znucl(:)
      logical :: have_zsymb = .false.
!
!
      write(6,*) " in makeatompp"
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
      allocate( zsymb(ntype) )
      open(unit=99,file='zsymb',form='formatted',status='old')
      read(99,*,IOSTAT=iostatus) zsymb(:)
      close(99)
!
      allocate( ppname(ntype) )
      open(unit=99,file='pplist',form='formatted',status='old')
      do i = 1, ntype
         read(99,'(a)') ppname(i)
      end do
      close(99)
      if( iostatus .eq. 0 ) have_zsymb = .true.
!
! get symbol & mass, concatenate
!
      allocate( satom(ntype), mass(ntype), ppline(ntype) )
      do i = 1, ntype
         call getsymbol( znucl(i), satom(i) )
         if( have_zsymb ) satom(i)=trim(zsymb(i))
!         if(trim(zsymb(i)) .ne. '') satom(i)=trim(zsymb(i))
         call getmass  ( znucl(i), mass (i) )
         ppline(i) = trim(satom(i)) // '   ' // mass(i) // '   ' // trim(ppname(i)) &
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


subroutine getsymbol_old(zatom,satom)
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
      case (38)
        satom = "Sr"
      case (74)
        satom = "W "
      case (78)
        satom = "Pt"

      case default
        satom = "--"

    end select


end subroutine getsymbol_old


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
        mass = "74.9210"
      case (34)
        mass = "78.9710"
      case (35)
        mass = "79.9040"
      case (36)
        mass = "83.7980"
      case (37)
        mass = "85.4678"
      case (38)
        mass = "87.6200"
      case (39)
        mass = "88.9058"
      case (40)
        mass = "91.2240"
      case (41)
        mass = "92.9064"
      case (42)
        mass = "95.9500"
      case (43)
        mass = "98.0000"
      case (44)
        mass = "101.070"
      case (45)
        mass = "102.900"
      case (46)
        mass = "106.420"
      case (47)
        mass = "107.868"
      case (48)
        mass = "112.414"
      case (49)
        mass = "114.818"
      case (50)
        mass = "118.710"
      case (51)
        mass = "121.760"
      case (52)
        mass = "127.600"
      case (53)
        mass = "126.900"
      case (54)
        mass = "131.293"
      case default
        mass = "120"
    end select

end subroutine getmass

