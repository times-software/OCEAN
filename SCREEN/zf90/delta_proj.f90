! Copyright (C) 2014 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program delta_proj
  implicit none

  integer, parameter :: DP = kind(1.d0)  

  integer :: ZNL(3), lmin, lmax, nr, nrtot
  integer, allocatable :: nprj( : )

  real(DP) :: rmax, prefs(0:1000), ae_val, ps_val
  real(DP), allocatable :: ae_prj(:,:), ps_prj(:,:), drel(:), wpt(:), posn(:,:),ang(:,:)

  character(len=11) :: prj_filename, rad_filename
  character(len=7) :: ae_file, ps_file


  integer :: il, im, ii, iang, nang, npt, tot_nproj, ipt, ir, iproj, ierr

  complex(DP) :: ylm
  complex(DP), allocatable :: delta_prj(:,:), ps_alone_prj(:,:)


  open( unit=99, file='ZNL', status='old', form='formatted' )
  rewind( 99 )
  read( 99, * ) ZNL(:)
  close(99)

  write(prj_filename,'(a8,i3.3)') 'prjfilez', ZNL(1)
  open( unit=99, file=prj_filename, status='old', form='formatted' )
  rewind(99)
  read(99,*) lmin, lmax 
  allocate( nprj( lmin : lmax ) )
  do il = lmin, lmax
    read( 99, * ) nprj( il )
  enddo
  close( 99 )



  write(rad_filename,'(a8,i3.3)') 'radfilez', ZNL(1)
  open( unit=99, file=rad_filename, status='old', form='formatted' )
  rewind(99)
  read(99,*) rmax, nr, nrtot
  close(99)
  


  open(unit=99,file='specpnt',form='formatted',status='old')
  read(99,*) nang
  allocate( ang(4,nang) )
  do il = 1, nang
    read(99,*) ang(:,il)
  enddo
  close(99)


  write(6,*) lmin, lmax
  write(6,*) nprj(:)
  write(6,*) rmax, nr, nrtot

  call newgetprefs( prefs, 4, nang, ang(4,:), ang(1,:), ang(2,:), ang(3,:) )


  open(unit=99,file='rbfile.bin',form='unformatted',status='old')
  read(99) npt
  allocate( posn(3,npt), wpt(npt), drel(npt) )
  read(99) posn
  read(99) wpt
  read(99) drel
  close(99)

  ipt = 1
  do while( drel( ipt ) .le. rmax )
    ipt = ipt + 1
    write(6,*) drel(ipt)
  enddo
  ipt = ipt -1

  write(6,*) 'new npt', ipt
  npt = ipt
   


  ii = 0
  do il = lmin, lmax
    do iproj = 1, nprj( il )
      do im = -il, il
        ii = ii + 1
      enddo
    enddo
  enddo

  tot_nproj = ii
  
  allocate( delta_prj( npt, tot_nproj ), ps_alone_prj( npt, tot_nproj ) )

  ii = 0 
  do il = lmin, lmax
    ! Read in projectors
    if( il .eq. lmin ) then
      allocate( ae_prj( 0 : nprj( il ), nr ), ps_prj( 0 : nprj( il ), nr ) )
    else
      deallocate( ae_prj, ps_prj ) 
      allocate( ae_prj( 0 : nprj( il ), nr ), ps_prj( 0 : nprj( il ), nr ) )
    endif

    write(ps_file,'(a2,i1.1,a1,i3.3)') 'ps', il, 'z', ZNL(1)
    write(6,*) ps_file
    open(unit=99,file=ps_file,form='formatted',status='old')
    do ir = 1, nr
      read(99,*) ps_prj(:,ir)
    enddo
    close(99)
    write(ae_file,'(a2,i1.1,a1,i3.3)') 'ae', il, 'z', ZNL(1)
    write(6,*) ae_file
    open(unit=99,file=ae_file,form='formatted',status='old')
    do ir = 1, nr
      read(99,*) ae_prj(:,ir)
    enddo
    close(99)
    

    do iproj = 1, nprj( il )
!      prev_rad = -1

      do im = -il, il
        ii = ii + 1
        write(6,*) ii
        iang = 0
        do ipt = 1, npt

!        ! cache radial interpolation
!        if( drel( ipt ) .ne. prev_rad ) then
!          prev_rad = drel( ipt )
          call radial_interpolate( drel(ipt), ae_val, nprj(il), nr, iproj, ae_prj, ierr )
          if( ierr .ne. 0 ) goto 111
          call radial_interpolate( drel(ipt), ps_val, nprj(il), nr, iproj, ps_prj, ierr )
          if( ierr .ne. 0 ) goto 111
!        endif

          iang = iang + 1
          if( iang .gt. nang ) iang = 1

          call newgetylm( il, im, ang(1,iang), ang(2,iang), ang(3,iang), ylm, prefs )
!          ylm = 1.0_DP
          delta_prj( ipt, ii ) = ylm * ( ae_val - ps_val )
          ps_alone_prj( ipt, ii ) = ylm * ps_val

        enddo
      enddo

!      ii = ii + ( 2*il + 1 )
    enddo
  

  enddo



  open(unit=99,file='delta_nrpoj.dat',form='unformatted',status='unknown')
!  do ipt = 1, npt
!    write(99,*) drel(ipt), delta_prj( ipt, : )
!  enddo
  write(99) npt, rmax
  write(99) delta_prj
  close(99)

  open(unit=99,file='ps_nrpoj.dat',form='unformatted',status='unknown')
  write(99) npt, rmax
  write(99) ps_alone_prj
  close(99)






  deallocate( nprj )


111 continue

contains


end program delta_proj

subroutine radial_interpolate( rad, val, nproj, nr, il, proj, ierr )
  implicit none
  
  integer, intent(in) :: nproj, nr, il
  integer, intent(out) :: ierr
  real(kind=kind(1.d0)), intent(out) :: val
  real(kind=kind(1.d0)), intent(in) :: rad, proj(0:nproj,nr)
  !
  integer :: low, high, mid

  ierr = 0
  low = 1
  high = nr

  if( rad .le. proj(0,1) ) then
    val = proj(il,1)
    return
  endif

  if( rad .ge. proj(0,nr) ) then
    ierr = -1
    return
  endif

!  do while ( high - low .gt. 1 )
!    mid = ( high - low + 1 ) / 2
!    if( rad .lt. proj(0,mid) ) then
!      high = mid
!    else
!      low = mid
!    endif
!  end do
  do high = 1, nr
    if( proj(0,high) .gt. rad ) goto 11
  enddo
11 continue
  low = high - 1

  val = proj(il,low) + ( rad - proj(0,low) ) * ( proj(il,low+1) - proj(il,low) ) / (proj(0,low+1) - proj(0,low) )

end subroutine radial_interpolate



