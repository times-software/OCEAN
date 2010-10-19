!      program wfconvert
!    This is to convert wf binaries from abinit to nbse formats
!    First step is to read in the abinit files
!    Written by JTV Nov 08

      program wfconvert
      IMPLICIT NONE

      integer :: nkpt,nspinor,nsppol,isppol,kout
      integer :: nband(2)
      double precision :: fermie,dummy
      integer :: iband,ikpt,maxnpw,un_npw,sh_npw 
      double precision, allocatable :: eigen_un(:),eigen_sh(:),         &
     &  cg_un(:,:),cg_sh(:,:), cg(:,:,:),cg_imag(:,:,:),kr(:,:),ki(:,:),&
     &  occ_sh(:),occ_un(:),kptlist(:,:)
      integer, allocatable :: kg_shift(:,:),kg_unshift(:,:),            &
     &   g_occ(:,:,:),gordered(:,:)
      integer :: brange(4),bandtot,maxband
      integer :: xstart,xend,ystart,yend,zstart,zend
      integer :: xrange,yrange,zrange,np_counter
      integer :: g_un_min(3),g_un_max(3),g_sh_min(3),g_sh_max(3)
      integer ::i,j,k,hkpt,gtot,kr_iter,ki_iter,nfiles
      integer :: files_iter,master_iter,g_iter,nkpts,kpt_counter
      character*11 :: wfkin!,wfkout
      character*12 :: wfkout
      character*4 :: wfknum
      character*4 :: kptnum
      double precision, allocatable :: enklisto(:,:),enklistu(:,:)
      double precision :: lda_low, lda_high, ldagap,qval
      double precision :: orthcr,orthci,q1,q2,q3,bv1(3),bv2(3),bv3(3)

      open(unit=45,file='enk_un',form='formatted',status='unknown')
      open(unit=46,file='enk_sh',form='formatted',status='unknown')
 
      open(unit=21,file='brange.ipt',form='formatted',status='old')
        read(21,*)brange(:)
      close(21)
      write(6,'(4I5)')brange(1),brange(2),brange(3),brange(4) 
      bandtot=brange(4)-brange(3)+brange(2)-brange(1)+2
   
      open(unit=22,file='Nfiles',form='formatted',status='old')
        read(22,*) nfiles
      close(22)

      open(unit=22,file='nkpts',form='formatted',status='old')
        read(22,*)nkpts
      close(22)
      nkpts = nkpts / 2
      kpt_counter = 0

      open(unit=75,file='qinunitsofbvectors.ipt',form='formatted',      &
     &   status='unknown')
      read(75,*) q1,q2,q3
      close(75)

      open(unit=22,file='bvecs',form='formatted',status='old')
        read(22,*) bv1(1),bv1(2),bv1(3)
        read(22,*) bv2(1),bv2(2),bv2(3)
        read(22,*) bv3(1),bv3(2),bv3(3)
      close(22)
      qval =dsqrt( (q1*bv1(1)+q2*bv2(1)+q3*bv3(1))**2                   &
     &            +(q1*bv1(2)+q2*bv2(2)+q3*bv3(2))**2                   &
     &            +(q1*bv1(3)+q2*bv2(3)+q3*bv3(3))**2)

      if (qval .le. 0 ) then
        write(6,*) 'Un-physical qval. Quitting...'
        stop
      endif

      allocate(enklisto(nkpts,brange(2)-brange(1)+1),                   &
     &         enklistu(nkpts,brange(4)-brange(3)+1))
      enklisto(:,:) = 0.0
      enklistu(:,:) = 0.0
      lda_low = 0.0
      lda_high = 0.0

      if (nfiles .gt. 9999) then
!     Only giving four places for file names. 9999 Is an arbitrary limit 
         write(6,*) 'Impending DOOM!'
         stop 'nfiles too great'
      endif

      open(unit=32,file='tmels',form='formatted',status='unknown')
      open(unit=23,file='listwfile',form='formatted',status='unknown')
      master_iter = 0

! open up klist, aren't currently using this
      open(unit=36,file='klist',form='formatted',status='unknown')
      rewind(36)

!   Now iterate over all the run files.
!   For each run file grab the matching k, k+q pts and write for NBSe


      do files_iter=1,nfiles
      
        write(wfknum,'(I4)') files_iter
        do i=1,3
          if ( wfknum(i:i) .eq. ' ' ) then
            wfknum(i:i) = '0'
          endif
        enddo

        wfkin = 'RUN' // wfknum // '_WFK'
        write(6,*) wfkin
        open(unit=20,file=wfkin,form='unformatted',status='old')
        call headerpar(20, dummy, maxnpw, maxband, nsppol, nspinor,nkpt,&
     &                 36)
!        do kout=1,nkpt
!         write(36,*)kptlist(1,kout),kptlist(2,kout),kptlist(3,kout)
!        enddo
!        deallocate(kptlist)

        allocate(kg_shift(3,maxnpw),kg_unshift(3,maxnpw),               &
     &      eigen_un(maxband), eigen_sh(maxband))

        if (nsppol .ne. 1) then
         write(6,*) "Spin stuff is currently not supported by AI2NBSE"
         ! quit
         stop 'problem with nsppol ne 1'
        endif
        if (nspinor .ne. 1 ) then
         write(6,*) "Spin stuff is currently not supported by AI2NBSE"
         stop 'problem with nspinor ne 1'
         ! quit, this is because I'm assuming npsinor = 1 for reading cg
        endif
      
        ! test nkpt is even,
        if ( mod(nkpt,2) .ne. 0 ) then
         write(6,*) 'Error: nkpt must be even. k and k+q together.'
         stop 'Cannot continue'
        endif
        hkpt = nkpt/2
      
        allocate(cg_un(maxband,2*maxnpw),cg_sh(maxband,2*maxnpw))
        allocate(occ_un(maxband),occ_sh(maxband))
        isppol=1
!       do isppol=1,nsppol   ! this is currently unsupported
        do ikpt=1,hkpt
          kpt_counter = kpt_counter + 1

         !  Set up the output files
          wfkout(1:8) = '.Psi   .' 
          wfkout(5:7) = wfkin(5:7)
          write(kptnum,'(I4.4)') ikpt
!          if ( kptnum(1:1) .eq. ' ' ) then
!            kptnum(1:1) = '0'
!          endif
          wfkout(9:12)= kptnum(1:4)
!         write(6,*) '     ',wfkout

          kg_shift(:,:) = 0
          kg_unshift(:,:) = 0

          open(unit=30,file=wfkout,form='unformatted',status='unknown')
          rewind(30)
          master_iter = master_iter + 1
!         write(23,'(I6,2X,A10,A11)') master_iter,'../ABINIT/',wfkout
          write(23,'(I6,2X,A11)') master_iter,wfkout
!         wfkout(1:4) = 'text'
!         open(unit=31,file=wfkout,form='formatted',status='unknown')
!     Currently the occ array is just thrown away info. Could/Should use 
!     populate the occ.dat written out now by NBSE (and maybe not used).


!     Read in the wavefunctions
          call grabwf(20, maxband, maxnpw, kg_unshift,                  &
     &     kg_shift, eigen_un, eigen_sh, occ_un, occ_sh, cg_un, cg_sh,  &
     &     brange(2), brange(4), nband, un_npw, sh_npw)


         write(45,*) 2.0*eigen_un(brange(1):brange(4))
         write(46,*) 2.0*eigen_sh(brange(1):brange(4))

         if (kpt_counter .eq. 1) then
           lda_low  = eigen_un(1)
           lda_high = eigen_sh(nband(2))
         endif

!        This can be changed to use fewer vars probably
         g_sh_min = minval(kg_shift, DIM=2)
         g_un_min = minval(kg_unshift, DIM=2)
         g_sh_max = maxval(kg_shift, DIM = 2)
         g_un_max = maxval(kg_unshift, DIM = 2)
         
         xstart = min(g_sh_min(1),g_un_min(1))
         xend   = max(g_sh_max(1),g_un_max(1))
         ystart = min(g_sh_min(2),g_un_min(2))
         yend   = max(g_sh_max(2),g_un_max(2))
         zstart = min(g_sh_min(3),g_un_min(3))
         zend   = max(g_sh_max(3),g_un_max(3))         
       
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
           g_occ(kg_unshift(1,i)-xstart+1,kg_unshift(2,i)-ystart+1,     &
     &           kg_unshift(3,i)-zstart+1) = 1
         enddo
         do i=1,sh_npw
           g_occ(kg_shift(1,i)-xstart+1,kg_shift(2,i)-ystart+1,         &
     &           kg_shift(3,i)-zstart+1) = 1
         enddo
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

!    write out all the real components, band by band
!    first the occ bands, then the unocc bands
!
!  Also, this is a fine time to store the eigenvalues
         do iband=brange(1),brange(2)
           enklisto(kpt_counter,iband) = eigen_un(iband)
!           if ((eigen_un(iband) .gt. lda_low) .and. (eigen_un(iband)    &
!      &         .lt. fermie)) then
!             lda_low = eigen_un(iband)
!           endif
!           write(6,*)kpt_counter,iband, eigen_un(iband)
!           write(6,*)enklisto(kpt_counter,iband)
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

! Currently this will break down for overlap. Need to five kr/ki allocation
! and storage here in the unocc bands.
         do iband=brange(3),brange(4)
           enklistu(kpt_counter,iband-brange(3)+1) = eigen_sh(iband)
!           if ((eigen_sh(iband) .lt. lda_high) .and. (eigen_sh(iband)   &
!     &         .gt. fermie)) then
!             lda_high = eigen_un(iband)
!           endif

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
                   ki(g_iter,iband+brange(2)-brange(3)+1) = cg_imag(i,j,k)
!                   write(31,*)iband,cg(i,j,k),cg_imag(i,j,k)
                 endif
               enddo
             enddo
           enddo 

         enddo

!   Write out the tmels for this set
!   This is written to be safe for all manner of brange
! 
         do i=brange(2)+1,brange(2)+1+brange(4)-brange(3)
           do j=1,brange(2)-brange(1)+1
             orthcr = 0.d0
             orthci = 0.d0
             do k=1,gtot
               orthcr = orthcr + kr(k,j)*kr(k,i) + ki(k,j)*ki(k,i)
               orthci = orthci + kr(k,j)*ki(k,i) - ki(k,j)*kr(k,i)
             enddo
             write(32,'(8(1x,1e22.15))')orthcr/qval,orthci/qval,0.0,0.0,&
     &                            0.0,0.0,0.0,0.0
!              write(32,'(2(1x,1e22.15),6(1x,I2))')orthcr/qval,       &
!     &               orthci/qval,0,0,0,0,0,0
           enddo
         enddo

         write(30) gtot
         write(30) gordered
         write(30) kr
         write(30) ki

         deallocate(cg)
         deallocate(cg_imag)
         deallocate(kr,ki)
         deallocate(g_occ)
         deallocate(gordered)
         close(30)

!         close(31)
       enddo 
!      enddo
!      deallocate(istwfk,nband,npwarr,so_typat,symafm,symrel,typat)
!      deallocate(kpt,occ,tnons,znucltypat,xred)
      deallocate(kg_shift,kg_unshift,eigen_un,eigen_sh)      
      deallocate(cg_un,cg_sh)
      deallocate(occ_un,occ_sh)

      close(20)



!      close(32)
      enddo  ! loop over run files.
      close(32)
      close(23)
      open(unit=24,file='masterwfile',form='formatted',status='unknown')
      write(24,*) master_iter
      close(24)

      open(unit=25,file='enkfile',form='formatted',status='unknown')
!      write(6,*)enklisto
      do i=1,nkpts
!        do j=brange(1),brange(2)
        write(25,*)(2*enklisto(i,j),j=brange(1),brange(2))
!        enddo
!        do j=1,brange(4)-brange(3)+1
        write(25,*)(2*enklistu(i,j),j=1,brange(4)-brange(3)+1)
!        enddo
      enddo
      close(25)
      close(36)

      deallocate(enklisto,enklistu)

!      open(unit=26,file='efermiinrydberg.ipt',form='formatted',         &
!     & status='unknown')
!       fermie = fermie*2
!        write(26,*)fermie
!      close(26) 

!      open(unit=27,file='ldagap',form='formatted',status='unknown')
!        ldagap=2*13.605698*(lda_high-lda_low)
!        write(27,*)ldagap
!      close(27)

       close(45)
       close(46)
      end program wfconvert
