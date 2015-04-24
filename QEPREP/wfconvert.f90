!      program wfconvert
!    This is to convert wf binaries from espresso to nbse formats
!    First step is to read in the espresso files
!    Written by JTV Nov 08; modified by KG July 2012
!     All energies should be in Rydberg

      program wfconvert
      IMPLICIT NONE

!      integer   :: nargs
!      character :: prefix*80

      REAL(kind=kind(1.d0)), PARAMETER :: ha2ryd = 2.0d0
      REAL(kind=kind(1.d0)), PARAMETER :: ev2ha  = 0.036749309
      REAL(kind=kind(1.d0)), PARAMETER :: ha2ev  = 27.211396132

      integer :: nkpt,nspinor,nsppol,isppol
      integer :: nband(2)
      integer :: iband,ikpt,maxnpw,un_npw,sh_npw 
      double precision, allocatable :: eigen_un(:),eigen_sh(:),         &
     &  cg_un(:,:),cg_sh(:,:), cg(:,:,:),cg_imag(:,:,:),kr(:,:),ki(:,:),&
     &  occ_sh(:),occ_un(:)
      integer, allocatable :: kg_shift(:,:),kg_unshift(:,:),            &
     &   g_occ(:,:,:),gordered(:,:)
      integer :: brange(4),bandtot,maxband
      integer :: xstart,xend,ystart,yend,zstart,zend
      integer :: xrange,yrange,zrange,np_counter
      integer :: g_un_min(3),g_un_max(3),g_sh_min(3),g_sh_max(3), umk(3)
      integer ::i,j,k,hkpt,gtot,kr_iter,ki_iter,nfiles
      integer :: files_iter,master_iter,g_iter,nkpts,kpt_counter
      character(len=11) :: wfkin
      character(len=12) :: wfkout
      double precision, allocatable :: enklisto(:,:),enklistu(:,:)
      double precision :: lda_low, lda_high,qval
      double precision :: orthcr,orthci,q1,q2,q3,bv1(3),bv2(3),bv3(3)
      logical :: noshift
      character(len=9), parameter :: f9 = 'formatted'
      !
      character(len=256) :: prefix, outdir
      !
      integer, parameter :: enkfile =   40
      integer, parameter :: tmels   =   41
      integer, parameter :: listwfile = 42
      integer, parameter :: mastwfile = 43
      integer, parameter :: enk_un  =   45
      integer, parameter :: enk_sh  =   46
      integer, parameter :: wfkinfile = 47
      integer, parameter :: wfkoutfile= 48
      integer, parameter :: umklapp =   49
      !
#ifdef __HAVE_F03 
      call get_command_argument(1,VALUE=prefix)
#else
      call getarg(1,prefix)
#endif
      !
      open(unit=umklapp,file='umklapp', form=f9,status='unknown')

      open(unit=36,file='klist',form=f9,status='unknown')
      open(unit=99,file='brange.ipt',form=f9,status='old')
        read(99,*)brange(:)
      close(99)
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

      noshift = .false.
      if ( (abs(q1) + abs(q2) + abs(q3) ) .eq. 0 ) noshift = .true.

      qval =dsqrt( (q1*bv1(1)+q2*bv2(1)+q3*bv3(1))**2                   &
     &            +(q1*bv1(2)+q2*bv2(2)+q3*bv3(2))**2                   &
     &            +(q1*bv1(3)+q2*bv2(3)+q3*bv3(3))**2)

      if ( (.not. noshift) .and. (qval .le. 0) ) then
        write(6,*) 'Un-physical qval. Quitting...'
        stop
      endif
      if (noshift) then
        qval = 1.d0
        allocate(enklisto(nkpts,brange(4)-brange(1)+1),                 &
     &         enklistu(nkpts,brange(4)-brange(3)+1))
      else
        allocate(enklisto(nkpts,brange(2)-brange(1)+1),                 &
     &         enklistu(nkpts,brange(4)-brange(3)+1))
      endif
        enklisto(:,:) = 0.0
        enklistu(:,:) = 0.0
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

!   Now iterate over all the run files.
!   For each run file grab the matching k, k+q pts and write for NBSE
      files_iter = 1

        call getwfkin(wfkin,files_iter,wfkinfile)
        write(6,*) wfkin

        write(6,*) " PREFIX = ", prefix

        outdir = "./Out"
        call qehead2(prefix, outdir, maxband, nsppol, nspinor, nkpt)
!        call qehead2("system", outdir, maxband, nsppol, nspinor, nkpt)

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
        if ( (.not. noshift) .and. (mod(nkpt,2) .ne. 0) ) then
         write(6,*) 'Error: nkpt must be even. k and k+q together.'
         stop 'Cannot continue'
        endif
        if (noshift) then
           hkpt = nkpt
        else
          hkpt = nkpt / 2
        endif

!        allocate(cg_un(maxband,2*maxnpw),cg_sh(maxband,2*maxnpw))
        allocate(occ_un(maxband),occ_sh(maxband))
        isppol=1
!       do isppol=1,nsppol   ! this is currently unsupported
        do ikpt=1,hkpt

          kpt_counter = kpt_counter + 1

          call qenpw2(prefix, outdir, ikpt, maxnpw )

          ! get the number of pw used for the current k-point

          allocate(cg_un(maxband,2*maxnpw),cg_sh(maxband,2*maxnpw))
          allocate(kg_unshift(3,maxnpw), eigen_un(maxband) )
          allocate(kg_shift(3,maxnpw), eigen_sh(maxband) )

          !if (.not. noshift) read(umklapp, * ) umk(:)

          call getwfkout( wfkout, files_iter, ikpt, wfkoutfile)
          if (.not. noshift)  kg_shift(:,:) = 0
          kg_unshift(:,:) = 0

          master_iter = master_iter + 1

          write(listwfile,'(I6,2X,A12)') master_iter,wfkout


!     Read in the wavefunctions
          call grabwf(prefix,outdir,ikpt,maxband,maxnpw,kg_unshift,kg_shift, &
     &  eigen_un,eigen_sh,occ_un,occ_sh,cg_un,cg_sh,brange(2),brange(4),  &
     &  nband,un_npw,sh_npw,noshift)

          kg_shift = kg_unshift

!KG!          ! eigenvalues are in eV
!KG!          ! eV to Hartree conversion
          eigen_un(:) = eigen_un(:) * ev2ha
          eigen_sh(:) = eigen_sh(:) * ev2ha

          if (.not. noshift) then
            write(enk_un,*) ha2ryd*eigen_un(brange(1):brange(4))
            write(enk_sh,*) ha2ryd*eigen_sh(brange(1):brange(4))
          else
            write(enk_un,*) ha2ryd*eigen_un(brange(1):brange(4))
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

         xrange = xend-xstart+1
         yrange = yend-ystart+1
         zrange = zend-zstart+1
         allocate(cg(xrange,yrange,zrange))
         allocate(cg_imag(xrange,yrange,zrange))
         allocate(g_occ(xrange,yrange,zrange))


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
               endif
             enddo
           enddo
         enddo


        if ( noshift) then
         do iband=brange(1),brange(2)
           enklisto(kpt_counter,iband) = eigen_un(iband)
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
           enklisto(kpt_counter,iband) = eigen_un(iband)
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
           enklisto(kpt_counter,iband) = eigen_un(iband)
           cg(:,:,:) = 0.0
           cg_imag(:,:,:) = 0.0
           do np_counter=1,un_npw  !? - un_npw - ?!
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
           enklistu(kpt_counter,iband-brange(3)+1) = eigen_sh(iband)
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
         deallocate(cg_un, cg_sh)
         deallocate(eigen_un)
         deallocate(kg_unshift)
         deallocate(kg_shift,eigen_sh)
         close(wfkoutfile)
         close(61)
         close(62)
         close(63)
         close(64)
!         close(31)
       enddo  ! end ink 


      deallocate(occ_un, occ_sh)


      close(wfkinfile)


      close(tmels)
      close(listwfile)
      open(unit=mastwfile,file='masterwfile',form=f9,status='unknown')
      write(mastwfile,*) master_iter
      close(mastwfile)

      open(unit=enkfile,file='enkfile',form=f9,status='unknown')
      do i=1,nkpts
        if (noshift) then
          write(enkfile,*) (ha2ev*ha2ryd*enklisto(i,j),j=brange(1),brange(2))
          write(enkfile,*) (ha2ev*ha2ryd*enklisto(i,j),j=brange(3),brange(4))
        else
          write(enkfile,*)(ha2ev*ha2ryd*enklisto(i,j),j=brange(1),brange(2))
          write(enkfile,*)(ha2ev*ha2ryd*enklistu(i,j),j=1,brange(4)-brange(3)+1)
        endif
      enddo
      close(enkfile)
      close(36)

      deallocate(enklisto,enklistu)


       close(enk_un)
       close(enk_sh)
       close(umklapp)
      end program wfconvert
