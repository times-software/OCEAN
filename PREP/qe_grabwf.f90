! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
      subroutine qe_grabwf(ikpt, maxband, maxnpw, kg_unshift, kg_shift, &
     &  eigen_un, eigen_sh, occ_un, occ_sh, cg_un, cg_sh, occ_max,      &
     &  unocc_max, nband, un_npw, sh_npw, noshift)
      use iotk_module

      integer :: ikpt, maxband, maxnpw, nband(2),                       &
     &  iband,ii,un_npw,sh_npw,nspinor, occ_max, unocc_max
      double precision :: eigen_un(maxband), eigen_sh(maxband),         &
     &  occ_un(maxband), occ_sh(maxband), cg_un(maxband,2*maxnpw),      &
     &  cg_sh(maxband,2*maxnpw)
      integer :: kg_unshift(3,maxnpw), kg_shift(3,maxnpw)
      logical :: noshift
      character(len=22) :: dirname
      character(len=16) :: prefix = 'Out/system.save/'
      character(len=128) :: filename
      integer :: ierr, npw

      write( dirname, '(16a,a,i5.5)' prefix, 'K', ikpt
      write(6,*) dirname

      ! Open eigval.xml
      filename = trim( dirname ) // 'eigval.xml'
      call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr )
      if( ierr .ne. 0 ) return
      call iotk_scan_dat( 99, "EIGENVALUES", eigen_un, IERR=ierr )
      if( ierr .ne. 0 ) return
      call iotk_scan_dat( 99, "OCCUPATIONS", occ_un, IERR=ierr )
      if( ierr .ne. 0 ) return
      call iotk_close_read( 99 )

      ! gkvectors.dat
      filename = trim( dirname ) // 'gkvectors.dat'
      call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr )
      if( ierr .ne. 0 ) return
      call iotk_scan_dat( 99, "NUMBER_OF_GK-VECTORS", npw, IERR=ierr )
      if( ierr .ne. 0 ) return
      call iotk_scan_dat( 99, "GRID", kg_unshift, IERR=ierr )
      if( ierr .ne. 0 ) return
      if( ierr .ne. 0 ) return
      call iotk_close_read( 99 )  
    
      ! evc.dat
      allocate( tbuffer( npw ) )
      filename = trim( dirname ) // 'evc.dat'
      call iotk_open_read ( 99, FILE = trim(filename), IERR=ierr ) 
      if( ierr .ne. 0 ) return
      do i = nband(1), nband(2)
        call iotk_scan_dat( 99, "evc" // trim(iotk_index( i ) ), tbuffer, IERR=ierr )
        if( ierr .ne. 0 ) return
        do j = 1, npw
          cg_un( i, 2*(j-1)+1 ) = real( tbuffer( j ) )
          cg_un( i, 2*(j-1)+2 ) = aimag( tbuffer( j ) )
        enddo
      enddo
      call iotk_close_read( 99 )
      deallocate( tbuffer )

      end 
      
