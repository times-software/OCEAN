! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program spec_average
 implicit none
  !
  integer :: Nvecs, Nedges, Nspecs, Npnts
  character*14 :: avename
  character*5  :: mode
  !
  integer, dimension(:), allocatable :: nval, lval, site
  character*2, dimension(:), allocatable :: elm
  character*20, dimension(:), allocatable :: polavename
  character*21, dimension(:), allocatable :: specname
  character*22, dimension(:), allocatable :: siteavename
  real( kind = kind(1.0d0) ), dimension(:), allocatable :: energy, dummyI
  real( kind = kind(1.0d0) ), dimension(:,:), allocatable :: specvecave, specpolave
  !
  logical, external :: checkmode
   !
   ! get the dimensions of the spectra
   call readinfo(mode,Nvecs,Nedges,Npnts)
   if ( checkmode(mode) .eqv. .false. ) stop
   Nspecs = Nedges * Nvecs
   !
   ! get the edge information
   allocate( nval(Nedges), lval(Nedges), elm(Nedges), site(Nedges) )
   call readhfin(Nedges,nval,lval,elm,site)
   !
   ! build the list of spectral file names
   allocate( specname(Nspecs), polavename(Nvecs), siteavename(Nedges) )
   call specnamelist(Nedges,Nvecs,Nspecs,mode,elm,nval,lval,site,specname,siteavename,polavename,avename)
   !
   ! get the energy grid
   allocate( energy(Npnts), dummyI(Npnts) )
   call readspec(specname(1),Npnts,energy,dummyI)
   !
   ! average over all polarizations by site
   allocate( specvecave(Nedges,Npnts), specpolave(Nvecs,Npnts) )
   call polarization_ave(Nedges,Nvecs,Npnts,Nspecs,specname,siteavename,specvecave)
   !
   ! average over all sites by polarization
   call site_ave(Nedges,Nvecs,Npnts,Nspecs,specname,polavename,specpolave)
   !
   ! total average; all sites and all polarizations
   call total_ave(Nvecs,Npnts,energy,specpolave,avename)
   !
   !
   deallocate( nval, lval, elm, site )
   deallocate( specname, polavename, siteavename )
   deallocate( specvecave, specpolave )
   deallocate( energy, dummyI )
   !
  !
end program spec_average


subroutine readinfo(mode,Nvecs,Nedges,Npnts)
 implicit none
  !
  integer, intent(out) :: Nvecs, Nedges, Npnts
  character*5, intent(out) :: mode
  !
  real( kind = kind(1.0d0) ) :: Emin, Emax
   !
   open( unit=99, file='cnbse.mode', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) mode
!   write(6,*) "Mode : ", mode, ':'
   close( 99 )
   !
   open( unit=99, file='cnbse.ways', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) Nvecs
!   write(6,*) "Nvecs : ", Nvecs
   close( 99 )
   !
   open( unit=99, file='nedges', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) Nedges
!   write(6,*) "Nedges : ", Nedges
   close( 99 )
   !
   open( unit=99, file='cnbse.spect_range', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) Npnts, Emin, Emax
   close( 99 )
  !
end subroutine readinfo


logical function checkmode( mode )
 implicit none
  !
  character*3, intent(in) :: mode
   !
   checkmode = .false.
   if ( mode .eq. 'xas' ) checkmode = .true.
   if ( mode .eq. 'xes' ) checkmode = .true.
  ! 
 return 
end function checkmode


subroutine readhfin(Nedges,nval,lval,elm,site)
 implicit none
  !
  integer, intent(in) :: Nedges
  !
  integer, dimension(Nedges), intent(out) :: nval, lval, site
  character*2, dimension(Nedges), intent(out) :: elm
  !
  integer :: i, Idum
  character*51 :: Adum51
   !
   open( unit=99, file='hfinlist', form='formatted', status='unknown' )
   rewind 99
   do i = 1, Nedges
      read( 99, '(A51,I3,I4,I4,A3,I4)' ) Adum51, Idum, nval(i), lval(i), elm(i), site(i)
   end do
   close( 99 )
  !
end subroutine readhfin


subroutine specnamelist(Nedges,Nvecs,Nspecs,mode,elm,nval,lval,site,specname,siteavename,polavename,avename)
 implicit none
  !
  integer, intent(in) :: Nedges, Nvecs, Nspecs
  character*3, intent(in) :: mode
  integer, dimension(Nedges), intent(in) :: nval, lval, site
  character*2, dimension(Nedges), intent(in) :: elm
  !
  character*14, intent(out) :: avename
  character*20, dimension(Nvecs), intent(out) :: polavename
  character*21, dimension(Nspecs), intent(out) :: specname
  character*22, dimension(Nedges), intent(out) :: siteavename
  !
  integer :: i,j,k
  character*4 :: coresym
  character*8 :: basename
   !
   if ( mode .eq. 'xas' ) then
      basename = 'absspct_'
   else if ( mode .eq. 'xes' ) then
      basename = 'xesspct_'
   end if
!   alternative to above
!   write(basename,'(A3,A5)') mode, 'spct_'
   !
   k = 1
   do i = 1, Nedges
      !
      select case ( lval(i) )
        case (0)
            write(coresym,'(A1,I1,A2)') '_', nval(i), 's_'
        case (1)
            write(coresym,'(A1,I1,A2)') '_', nval(i), 'p_'
        case (2)
            write(coresym,'(A1,I1,A2)') '_', nval(i), 'd_'
        case (3)
            write(coresym,'(A1,I1,A2)') '_', nval(i), 'f_'
      end select
      !
      if ( i == 1 ) then
         do j = 1, Nvecs
            write(polavename(j),'(A8,A2,A4,A4,I0.2)') basename, elm(i), '.ave', coresym, j
!            write(6,*) polavename(j)
         end do
      end if
      !
      do j = 1, Nvecs
         write(specname(k),'(A8,A2,A1,I0.4,A4,I0.2)') basename, elm(i), '.', site(i), coresym, j
         k = k + 1
      end do
      write(siteavename(i),'(A8,A2,A1,I0.4,A4,A3)') basename, elm(i), '.', site(i), coresym, 'ave'
   end do
   !
   write(avename,'(A8,A2,A4)') basename, elm(1), '.ave'
  !
end subroutine specnamelist


subroutine readspec(specname,Npnts,energy,intensity)
 implicit none
  !
  integer, intent(in) :: Npnts
  character*21, intent(in) :: specname
  !
  real( kind = kind(1.0d0) ), dimension(Npnts), intent(out) :: energy, intensity
  !
  integer :: l, dummyI
  real( kind = kind(1.0d0) ) :: dummyR
   !
   open( unit=99, file=specname, form='formatted', status='unknown' )
   rewind 99
   do l = 1,Npnts
      read ( 99, * ) energy(l), intensity(l), dummyR, dummyR, dummyI, dummyR, dummyR, dummyI
   end do
   close( 99 )
  !
end subroutine readspec


subroutine polarization_ave(Nedges,Nvecs,Npnts,Nspecs,specname,specavename,specint)
 implicit none
  !
  integer, intent(in) :: Nedges, Nvecs, Npnts, Nspecs
  character*21, dimension(Nspecs), intent(in) :: specname
  character*22, dimension(Nedges), intent(in) :: specavename
  !
  real( kind = kind(1.0d0) ), dimension(Nedges,Npnts), intent(out) :: specint
  !
  integer :: i,j,k,l
  character*26 :: tmpspecavename
  real( kind = kind(1.0d0) ), dimension(Npnts) :: energy, intensity
   !
   k = 1
   do i = 1, Nedges
      !
      specint(i,:) = 0.0d0
      do j = 1, Nvecs
         !   
         call readspec(specname(k),Npnts,energy,intensity)
         !
         do l = 1, Npnts
            specint(i,l) = specint(i,l) + ( intensity(l) / (1.0*Nvecs) )
         end do
         !
         k = k + 1
      end do
      !
      open( unit=99, file=specavename(i), form='formatted', status='unknown' )
      do l = 1, Npnts
         write ( 99, '(E15.8,E16.8)' ) energy(l), specint(i,l)
      end do
      close( 99 )
      !
   end do
  !
end subroutine polarization_ave


subroutine site_ave(Nedges,Nvecs,Npnts,Nspecs,specname,specavename,specint)
 implicit none
  !
  integer, intent(in) :: Nedges, Nvecs, Npnts, Nspecs
  character*21, dimension(Nspecs), intent(in) :: specname
  character*20, dimension(Nvecs), intent(in) :: specavename
  !
  real( kind = kind(1.0d0) ), dimension(Nvecs,Npnts), intent(out) :: specint
  !
  integer :: i,j,k,l
  character*26 :: tmpspecavename
  real( kind = kind(1.0d0) ), dimension(Npnts) :: energy, intensity
   !
   do j = 1, Nvecs
      !
      specint(j,:) = 0.0d0
      do i = 1, Nedges
         k = (i-1)*Nvecs + j
         call readspec(specname(k),Npnts,energy,intensity)
         do l = 1, Npnts
            specint(j,l) = specint(j,l) + ( intensity(l) / (1.0*Nedges) )
         end do
      end do
      !
      open( unit=99, file=specavename(j), form='formatted', status='unknown' )
      do l = 1, Npnts
         write ( 99, '(E15.8,E16.8)' ) energy(l), specint(j,l)
      end do
      close( 99 )
      !
   end do
  !
end subroutine site_ave


subroutine total_ave(Nvecs,Npnts,energy,specpolave,avename)
 implicit none
  !
  integer, intent(in) :: Nvecs, Npnts
  character*14, intent(in) :: avename
  real( kind = kind(1.0d0) ), dimension(Npnts), intent(in) :: energy
  real( kind = kind(1.0d0) ), dimension(Nvecs,Npnts), intent(in) :: specpolave
  !
  integer :: j,l
  real( kind = kind(1.0d0) ), dimension(Npnts) :: totave
   !
   totave(:) = 0.0d0
   do j = 1, Nvecs
      do l = 1, Npnts
         totave(l) = totave(l) + ( specpolave(j,l) / Nvecs )
      end do
   end do
   !
   open( unit=99, file=avename, form='formatted', status='unknown' )
   do l = 1, Npnts
      write ( 99, '(E15.8,E16.8)' ) energy(l), totave(l)
   end do
   close( 99 )
  !
end subroutine

