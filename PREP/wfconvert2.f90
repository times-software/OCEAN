! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
!      program wfconvert
!    This is to convert wf binaries from abinit to nbse formats
!    First step is to read in the abinit files
!    Written by JTV Nov 08

program wfconvert
      IMPLICIT NONE

      integer :: nkpt(2),nspinor(2),nsppol(2),isppol
      integer :: nband(2), fh
      double precision :: dummy
      integer :: iband,ikpt,maxnpw(2),un_npw,sh_npw 
      double precision, allocatable :: eigen_un(:),eigen_sh(:),         &
     &  cg_un(:,:),cg_sh(:,:), cg(:,:,:),cg_imag(:,:,:),kr(:,:),ki(:,:),&
     &  occ_sh(:),occ_un(:)
      integer, allocatable :: kg_shift(:,:),kg_unshift(:,:),            &
     &   g_occ(:,:,:),gordered(:,:)
      integer :: brange(4),bandtot,maxband(2), nspin
      integer :: xstart,xend,ystart,yend,zstart,zend
      integer :: xrange,yrange,zrange,np_counter
      integer :: g_un_min(3),g_un_max(3),g_sh_min(3),g_sh_max(3), umk(3)
      integer ::i,j,k,hkpt,gtot,kr_iter,ki_iter,nfiles, ierr
      integer :: master_iter,g_iter,nkpts,kpt_counter, qe_kpt
      character(len=11) :: wfkin!,wfkout
      character(len=12) :: wfkout
      double precision, allocatable :: enklisto(:,:,:),enklistu(:,:,:)
      double precision :: lda_low, lda_high,qval
      double precision :: orthcr,orthci,q1,q2,q3,bv1(3),bv2(3),bv3(3)
      logical :: noshift, split_dft
      character(len=3) :: dft_flavor
#ifdef __HAVE_IOTK
      character(len=128) :: qe_filename
#endif
      character(len=9), parameter :: f9 = 'formatted'
!
      integer, parameter :: enkfile =   40
      integer, parameter :: tmels   =   41
      integer, parameter :: listwfile = 42
      integer, parameter :: mastwfile = 43
      integer, parameter :: enk_un  =   45
      integer, parameter :: enk_sh  =   46
      integer, parameter :: wfkinfile = 47
      integer, parameter :: wfkinfile_shift = 50
      integer, parameter :: wfkoutfile= 48
      integer, parameter :: umklapp =   49
! 
      open(unit=umklapp,file='umklapp', form=f9,status='unknown')

      open(unit=36,file='klist',form=f9,status='unknown')
      open(unit=99,file='brange.ipt',form=f9,status='old')
        read(99,*)brange(:)
      close(99)
      write(6,'(4I5)')brange(1),brange(2),brange(3),brange(4) 
      bandtot=brange(4)-brange(3)+brange(2)-brange(1)+2
   
      open(unit=99,file='Nfiles',form=f9,status='old')
        read(99,*) nfiles
      close(99)

      open(unit=99,file='nkpts',form=f9,status='old')
        read(99,*)nkpts
      close(99) 
      kpt_counter = 0

      open(unit=99,file='qinunitsofbvectors.ipt',form=f9,      &
     &   status='unknown')
      read(99,*) q1,q2,q3
      close(99)

      open(unit=99,file='bvecs',form=f9,status='old')
        read(99,*) bv1(1),bv1(2),bv1(3)
        read(99,*) bv2(1),bv2(2),bv2(3)
        read(99,*) bv3(1),bv3(2),bv3(3)
      close(99)

      open(unit=99,file='dft',form=f9,status='old')
      read(99,'(a)') dft_flavor
      close(99)

      open(unit=99,file='nspin',form=f9,status='old')
      read(99,*) nspin
      close(99)


      noshift = .false.
      if ( (abs(q1) + abs(q2) + abs(q3) ) .eq. 0 ) noshift = .true.

      if( .not. noshift ) then
        inquire( file='dft.split', exist=split_dft )
        if( split_dft ) then
          open( unit=99, file='dft.split', form='formatted', status='old' )
          read( 99, * ) split_dft
          close( 99 )
        endif
      else
        split_dft = .false.
      endif

      qval =dsqrt( (q1*bv1(1)+q2*bv2(1)+q3*bv3(1))**2                   &
     &            +(q1*bv1(2)+q2*bv2(2)+q3*bv3(2))**2                   &
     &            +(q1*bv1(3)+q2*bv2(3)+q3*bv3(3))**2)

      if ( (.not. noshift) .and. (qval .le. 0) ) then
        write(6,*) 'Un-physical qval. Quitting...'
        stop
      endif
      if (noshift) then
        qval = 1.d0
        allocate(enklisto(brange(4)-brange(1)+1,nkpts,nspin),                 &
                 enklistu(brange(4)-brange(3)+1,nkpts,nspin) )
      else
        allocate(enklisto(brange(2)-brange(1)+1,nkpts,nspin),                 &
                 enklistu(brange(4)-brange(3)+1,nkpts,nspin ) )
      endif

      enklisto(:,:,:) = 0.0
      enklistu(:,:,:) = 0.0
      open(unit=enk_un, file='enk_un',form=f9,status='unknown')
      open(unit=enk_sh, file='enk_sh',form=f9,status='unknown')

      lda_low = 0.0
      lda_high = 0.0

      if (nfiles .gt. 9999) then
!     Only giving four places for file names. 9999 Is an arbitrary limit 
         write(6,*) 'Impending DOOM!'
         stop 'nfiles too great'
      endif

      open(unit=tmels,file='tmels',form=f9,status='unknown')
      open(unit=listwfile,file='listwfile',form=f9,status='unknown')
      master_iter = 0

!   For each run file grab the matching k, k+q pts and write for NBSE
      select case( dft_flavor )
        case( 'qe' )
#ifdef __HAVE_IOTK
          qe_filename = 'Out/system.save/data-file.xml'
          call qehead( qe_filename, maxband(1), maxnpw(1), nsppol(1), nspinor(1), nkpt(1), ierr )
          if( split_dft ) then
            qe_filename = 'Out_shift/system.save/data-file.xml'
            call qehead( qe_filename, maxband(2), maxnpw(2), nsppol(2), nspinor(2), nkpt(2), ierr )
          endif
#else
          open(unit=99,file='qe_data.txt',form='formatted',status='old')
          read(99,*) maxband(1), maxnpw(1), nsppol(1), nspinor(1), nkpt(1)
          if( split_dft ) then
            read(99,*) maxband(2), maxnpw(2), nsppol(2), nspinor(2), nkpt(2)
          endif
          close(99)
#endif
          if( nsppol(1) .ne. nspin ) then
            write(6,'(a,i1.1,a,i1.1)' ) 'Was expecting spin=',nspin, ' but found spin=', nsppol
            stop
          endif
          
        case default
          
          call getwfkin(wfkin,1,wfkinfile)
          write(6,*) wfkin
          call headerpar(wfkinfile, dummy, maxnpw(1), maxband(1), nsppol(1), nspinor(1),nkpt(1), 36, ierr)

          if( split_dft) then
            call getwfkin(wfkin,2,wfkinfile_shift)
            write(6,*) wfkin
            call headerpar(wfkinfile_shift, dummy, maxnpw(2), maxband(2), nsppol(2), nspinor(2),nkpt(2), 37, ierr)
          endif
      end select

      if( .not. split_dft ) then
        maxband(2) = maxband(1)
        maxnpw(2) = maxnpw(1)
        nsppol(2) = nsppol(1)
        nspinor(2) = nspinor(1)
        nkpt(2) = nkpt(1)
      else ! test
        if( nkpt(1) .ne. nkpt(2) ) then
          write(6,*) "Mismatch between unshifted/shifted kpts:", nkpt(:)
          stop
        endif
        if( nspinor(1) .ne. nspinor(2) ) then
          write(6,*) "Mismatch between unshifted/shifted spinors:", nspinor(:)
          stop
        endif
        if( nsppol(1) .ne. nsppol(2) ) then
          write(6,*) "Mismatch between unshifted/shifted spin:", nsppol(:)
          stop
        endif
      endif

      allocate(kg_unshift(3,maxnpw(1)), eigen_un(maxband(1)) )
      !if ( .not. noshift ) 
      allocate(kg_shift(3,maxnpw(2)), eigen_sh(maxband(2)) )

      if (nspinor(1) .ne. 1 ) then
       write(6,*) "Spinor stuff is currently not supported by AI2NBSE"
       stop 'problem with nspinor ne 1'
       ! quit, this is because I'm assuming npsinor = 1 for reading cg
      endif
    
      ! test nkpt is even when shift and unshifted are together
      if( .not. split_dft ) then
        if ( (.not. noshift) .and. (mod(nkpt(1),2) .ne. 0) ) then
          write(6,*) 'Error: nkpt must be even. k and k+q together.'
          stop 'Cannot continue'
        endif
      endif

      if (noshift) then
         hkpt = nkpt(1)
      elseif( .not. split_dft ) then
        hkpt = nkpt(1) / 2
      else
        hkpt = nkpt(1)
      endif

      allocate(cg_un(maxband(1),2*maxnpw(1)),cg_sh(maxband(2),2*maxnpw(2)))
      allocate(occ_un(maxband(1)),occ_sh(maxband(2)))

      isppol=1
      do isppol=1,nsppol(1)   ! this is currently unsupported
        do ikpt=1,hkpt
          kpt_counter = kpt_counter + 1

          !if (.not. noshift) read(umklapp, * ) umk(:)

          call getwfkout( wfkout, isppol, ikpt, wfkoutfile)
          if (.not. noshift)  kg_shift(:,:) = 0
          kg_unshift(:,:) = 0

          master_iter = master_iter + 1

          write(listwfile,'(I6,2X,A12)') master_iter,wfkout


!     Read in the wavefunctions
          select case( dft_flavor )

          case( 'qe' )

            if( noshift .or. split_dft ) then
              qe_kpt = ikpt
            else
              qe_kpt = ( ikpt -1 ) * 2 + 1
            endif

#ifdef __HAVE_IOTK
            call qe_grabwf(qe_kpt, isppol, nsppol(1), maxband(1), maxnpw(1), kg_unshift,  &
     &                     eigen_un, occ_un, cg_un, un_npw, ierr)
#else
            call qe_grabwf_noiotk(qe_kpt, isppol, nsppol(1), maxband(1), maxnpw(1), kg_unshift,  &
     &                     eigen_un, cg_un, un_npw, ierr)
#endif
            if( ierr .ne. 0 ) then
              write(6,*) ikpt, ierr
              stop
            endif

            if( .not. noshift ) then
              if( .not. split_dft ) then
                qe_kpt = qe_kpt + 1
              endif
#ifdef __HAVE_IOTK
              call qe_grabwf(qe_kpt, isppol, nsppol(2), maxband(2), maxnpw(2), kg_shift,  &
     &                       eigen_sh, occ_sh, cg_sh, sh_npw, ierr)
#else
              call qe_grabwf_noiotk(qe_kpt, isppol, nsppol(2), maxband(2), maxnpw(2), kg_shift,  &
     &                     eigen_sh, cg_sh, sh_npw, ierr)
#endif
            endif

          case default

            if( noshift ) then
              ierr = 0
              call grabwf_abi(wfkinfile, maxband(1), maxnpw(1), brange(4), un_npw, kg_unshift, eigen_un, &
                          cg_un, ierr )
              if( ierr .ne. 0 ) then
                write(6,*) "Failed to read unshifted at kpt: ", ikpt
                stop
              endif

            else
              ierr = 0
              call grabwf_abi(wfkinfile, maxband(1), maxnpw(1), brange(2), un_npw, kg_unshift, eigen_un, &
                          cg_un, ierr )
              if( ierr .ne. 0 ) then
                write(6,*) "Failed to read unshifted at kpt: ", ikpt
                stop
              endif

              if( split_dft ) then
                fh = wfkinfile_shift
              else
                fh = wfkinfile
              endif
              call grabwf_abi(fh, maxband(2), maxnpw(2), brange(4), sh_npw, kg_shift, eigen_sh, &
                          cg_sh, ierr )
              if( ierr .ne. 0 ) then
                write(6,*) "Failed to read unshifted at kpt: ", ikpt
                stop
              endif
            endif


          end select

          if (.not. noshift) then
            if( split_dft ) then
              write(enk_un,*) 2.0*eigen_un(brange(1):brange(2))
              write(enk_sh,*) 2.0*eigen_sh(brange(1):brange(4))
            else
              write(enk_un,*) 2.0*eigen_un(brange(1):brange(4))
              write(enk_sh,*) 2.0*eigen_sh(brange(1):brange(4))
            endif
          else
            write(enk_un,*) 2.0*eigen_un(brange(1):brange(4))
          endif

          if (kpt_counter .eq. 1) then
            lda_low  = eigen_un(1)
            lda_high = lda_low !eigen_sh(nband(2))
          endif

!        This can be changed to use fewer vars probably
          if ( noshift ) then
            g_un_min = minval(kg_unshift, DIM=2) 
            g_un_max = maxval(kg_unshift, DIM = 2)
            xstart = g_un_min(1)
            xend = g_un_max(1)
            ystart = g_un_min(2)
            yend = g_un_max(2)
            zstart = g_un_min(3)
            zend = g_un_max(3)
          else
            ! account for umklapp by shifting mapping
            read(umklapp, * ) umk(:)
            do i = 1, sh_npw
              kg_shift(:,i) = kg_shift(:,i) - umk(:)
            enddo
            !
            g_sh_min = minval(kg_shift, DIM=2)
            g_un_min = minval(kg_unshift, DIM=2)
            g_sh_max = maxval(kg_shift, DIM = 2)
            g_un_max = maxval(kg_unshift, DIM = 2)
        
            !! account for umklapp
            !!g_sh_min(:) = g_sh_min(:) = umk(:)
            !!g_sh_max(:) = g_sh_max(:) = umk(:)
            !! 
            xstart = min(g_sh_min(1),g_un_min(1))
            xend   = max(g_sh_max(1),g_un_max(1))
            ystart = min(g_sh_min(2),g_un_min(2))
            yend   = max(g_sh_max(2),g_un_max(2))
            zstart = min(g_sh_min(3),g_un_min(3))
            zend   = max(g_sh_max(3),g_un_max(3))         
          endif
!         write(6,*)xstart,xend,ystart,yend,zstart,zend

         xrange = xend-xstart+1
         yrange = yend-ystart+1
         zrange = zend-zstart+1
         allocate(cg(xrange,yrange,zrange))
         allocate(cg_imag(xrange,yrange,zrange))
         allocate(g_occ(xrange,yrange,zrange))

!         write(31,'(I3)')bandtot


!        write out gvectors
          g_occ(:,:,:) = 0
          do i=1,un_npw
            g_occ(kg_unshift(1,i)-xstart+1,kg_unshift(2,i)-ystart+1,    &
     &           kg_unshift(3,i)-zstart+1) = 1
          enddo
          if (.not. noshift) then
            do i=1,sh_npw
              g_occ(kg_shift(1,i)-xstart+1,kg_shift(2,i)-ystart+1,      &
     &              kg_shift(3,i)-zstart+1) = 1
            enddo
          endif
!
          gtot = sum(g_occ)
!         write(31,*) gtot

         allocate(kr(gtot,bandtot),ki(gtot,bandtot))
         allocate(gordered(gtot,3))
         kr_iter = 0
         ki_iter = 0
         g_iter  = 0
         do i=1,xrange
           do j=1,yrange
             do k=1,zrange
              if (g_occ(i,j,k) .eq. 1 ) then
                 g_iter = g_iter + 1
                 gordered(g_iter,1) = i+xstart-1
                 gordered(g_iter,2) = j+ystart-1
                 gordered(g_iter,3) = k+zstart-1
!                 write(31,'(3(i5))')i+xstart-1,j+ystart-1,k+zstart-1
               endif
             enddo
           enddo
         enddo

        if ( noshift) then
         do iband=brange(1),brange(2)
           enklisto(iband,ikpt,isppol) = eigen_un(iband)
           cg(:,:,:) = 0.0
           cg_imag(:,:,:) = 0.0
           do np_counter=1,un_npw
!   Need to shift since our gvectors are -x to x and we have them as 
!   0 to 2x+1
             cg(kg_unshift(1,np_counter)-xstart+1,                      &
     &          kg_unshift(2,np_counter)-ystart+1,                      &
     &          kg_unshift(3,np_counter)-zstart+1) =                    &
     &       cg_un(iband,np_counter*2-1)
             cg_imag(kg_unshift(1,np_counter)-xstart+1,                 &
     &          kg_unshift(2,np_counter)-ystart+1,                      &
     &          kg_unshift(3,np_counter)-zstart+1) =                    &
     &       cg_un(iband,np_counter*2)
           enddo
           g_iter = 0
           do i=1,xrange
             do j=1,yrange
               do k=1,zrange
                 if (g_occ(i,j,k) .eq. 1 ) then
                   g_iter = g_iter + 1
                   kr(g_iter,iband) = cg(i,j,k)
                   ki(g_iter,iband) = cg_imag(i,j,k)
!                   write(31,*)iband,cg(i,j,k),cg_imag(i,j,k)
                 endif
               enddo
             enddo
           enddo
         enddo

         do iband=brange(3),brange(4)
           enklisto(iband,ikpt,isppol) = eigen_un(iband)
           cg(:,:,:) = 0.0
           cg_imag(:,:,:) = 0.0
           do np_counter=1,un_npw
!   Need to shift since our gvectors are -x to x and we have them as 
!   0 to 2x+1
             cg(kg_unshift(1,np_counter)-xstart+1,                      &
     &          kg_unshift(2,np_counter)-ystart+1,                      &
     &          kg_unshift(3,np_counter)-zstart+1) =                    &
     &       cg_un(iband,np_counter*2-1)
             cg_imag(kg_unshift(1,np_counter)-xstart+1,                 &
     &          kg_unshift(2,np_counter)-ystart+1,                      &
     &          kg_unshift(3,np_counter)-zstart+1) =                    &
     &       cg_un(iband,np_counter*2)
           enddo
           g_iter = 0
           do i=1,xrange
             do j=1,yrange
               do k=1,zrange
                 if (g_occ(i,j,k) .eq. 1 ) then
                   g_iter = g_iter + 1
                   kr(g_iter,iband+brange(2)-brange(3)+1) = cg(i,j,k)
                   ki(g_iter,iband+brange(2)-brange(3)+1) = cg_imag(i,j,k)
!                   write(31,*)iband,cg(i,j,k),cg_imag(i,j,k)
                 endif
               enddo
             enddo
           enddo
         enddo

        else
!    write out all the real components, band by band
!    first the occ bands, then the unocc bands
!
!  Also, this is a fine time to store the eigenvalues
         do iband=brange(1),brange(2)
           enklisto(iband,ikpt,isppol) = eigen_un(iband)
           cg(:,:,:) = 0.0
           cg_imag(:,:,:) = 0.0
           do np_counter=1,un_npw
!   Need to shift since our gvectors are -x to x and we have them as 
!   0 to 2x+1
             cg(kg_unshift(1,np_counter)-xstart+1,                      &
     &          kg_unshift(2,np_counter)-ystart+1,                      &
     &          kg_unshift(3,np_counter)-zstart+1) =                    &
     &       cg_un(iband,np_counter*2-1)
             cg_imag(kg_unshift(1,np_counter)-xstart+1,                 &
     &          kg_unshift(2,np_counter)-ystart+1,                      &
     &          kg_unshift(3,np_counter)-zstart+1) =                    &
     &       cg_un(iband,np_counter*2)
           enddo
           g_iter = 0
           do i=1,xrange
             do j=1,yrange
               do k=1,zrange
                 if (g_occ(i,j,k) .eq. 1 ) then
                   g_iter = g_iter + 1
                   kr(g_iter,iband) = cg(i,j,k)
                   ki(g_iter,iband) = cg_imag(i,j,k)
!                   write(31,*)iband,cg(i,j,k),cg_imag(i,j,k)
                 endif
               enddo 
             enddo
           enddo
         enddo


         do iband=brange(3),brange(4)
           enklistu(iband-brange(3)+1,ikpt,isppol) = eigen_sh(iband)
           cg(:,:,:) = 0.0
           cg_imag(:,:,:) = 0.0
           do i=1,sh_npw
             cg(kg_shift(1,i)-xstart+1,kg_shift(2,i)-ystart+1,          &
     &          kg_shift(3,i)-zstart+1) = cg_sh(iband,i*2-1)
             cg_imag(kg_shift(1,i)-xstart+1,kg_shift(2,i)-ystart+1,     &
     &          kg_shift(3,i)-zstart+1) = cg_sh(iband,i*2)

           enddo
           g_iter = 0
           do i=1,xrange
             do j=1,yrange
               do k=1,zrange
                 if (g_occ(i,j,k) .eq. 1 ) then
                   g_iter = g_iter + 1
                   kr(g_iter,iband+brange(2)-brange(3)+1) = cg(i,j,k)
                   ki(g_iter,iband+brange(2)-brange(3)+1)=cg_imag(i,j,k)
!                   write(31,*)iband,cg(i,j,k),cg_imag(i,j,k)
                 endif
               enddo
             enddo
           enddo 

         enddo

!   Write out the tmels for this set
!   This is written to be safe for all manner of brange
! 
       endif !! switch for noshift/shift

         do i=brange(2)+1,brange(2)+1+brange(4)-brange(3)
           do j=1,brange(2)-brange(1)+1
             orthcr = 0.d0
             orthci = 0.d0
             do k=1,gtot
               orthcr = orthcr + kr(k,j)*kr(k,i) + ki(k,j)*ki(k,i)
               orthci = orthci + kr(k,j)*ki(k,i) - ki(k,j)*kr(k,i)
             enddo
             write(tmels,'(8(1x,1e22.15))')orthcr/qval,orthci/qval,0.0,0.0,&
     &                            0.0,0.0,0.0,0.0
!              write(32,'(2(1x,1e22.15),6(1x,I2))')orthcr/qval,       &
!     &               orthci/qval,0,0,0,0,0,0
           enddo
         enddo

         write(wfkoutfile) gtot
         write(wfkoutfile) gordered
         write(wfkoutfile) kr
         write(wfkoutfile) ki

         deallocate(cg)
         deallocate(cg_imag)
         deallocate(kr,ki)
         deallocate(g_occ)
         deallocate(gordered)
         close(wfkoutfile)
!         close(31)
       enddo 

!      deallocate(istwfk,nband,npwarr,so_typat,symafm,symrel,typat)
!      deallocate(kpt,occ,tnons,znucltypat,xred)
      if (.not. noshift) deallocate(kg_shift,eigen_sh)
      deallocate(kg_unshift,eigen_un)      
      deallocate(cg_un, cg_sh)
      deallocate(occ_un, occ_sh)

      close(wfkinfile)


!      close(32)
      enddo  ! loop over run files.

      close(tmels)
      close(listwfile)
      open(unit=mastwfile,file='masterwfile',form=f9,status='unknown')
      write(mastwfile,*) master_iter
      close(mastwfile)

      open(unit=enkfile,file='enkfile',form=f9,status='unknown')
!      write(6,*)enklisto
      do isppol = 1, nsppol(1)
        do i=1,nkpts
          if (noshift) then
            write(enkfile,*) (2*enklisto(j,i,isppol),j=brange(1),brange(2))
            write(enkfile,*) (2*enklisto(j,i,isppol),j=brange(3),brange(4))
          else
            write(enkfile,*)(2*enklisto(j,i,isppol),j=brange(1),brange(2))
            write(enkfile,*)(2*enklistu(j,i,isppol),j=1,brange(4)-brange(3)+1)
          endif
        enddo
      enddo
      close(enkfile)
      close(36)

      deallocate(enklisto,enklistu)

!      open(unit=26,file='efermiinrydberg.ipt',form=f9,         &
!     & status='unknown')
!       fermie = fermie*2
!        write(26,*)fermie
!      close(26) 

!      open(unit=27,file='ldagap',form=f9,status='unknown')
!        ldagap=2*13.605698*(lda_high-lda_low)
!        write(27,*)ldagap
!      close(27)

       close(enk_un)
       close(enk_sh)
       close(umklapp)
end program wfconvert
