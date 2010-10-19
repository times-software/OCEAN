      subroutine vcbder2( ibl, ibh, vlev, vhev, clev, chev )
!      program vcdber
      implicit none
      integer :: ibl, ibh
      double precision :: vlev, vhev, clev, chev

      integer :: brange(4), nkpts,numocc,numunocc,bandcount,kcount
      integer :: kmesh(3)
      double precision :: efermi, ry
      double precision, allocatable :: enklist(:,:)
    
      ry = 13.6057d0

!      ibl = 1
!      ibh = 12

      open(unit=99,file='brange.ipt',form='formatted',status='old')
      read(99,*)brange(:)
!      write(6,*) brange(:)
      close( 99 )

      open(unit=99,file='kmesh.ipt',form='formatted',status='old')
      read(99,*)kmesh(:)
!      write(6,*) "nkpts = ",nkpts
      close(99)
      nkpts = kmesh(1)*kmesh(2)*kmesh(3)
      
      numocc = brange(2)-brange(1)+1
      numunocc = brange(4)-brange(3)+1
      if ( (ibl .gt. brange(2)) .or. (ibh .gt. brange(4)) ) then
       write(6,*) "ibl or ibh too big"
       stop 
      endif

      allocate(enklist(numocc+numunocc,nkpts))
      open(unit=99,file='enkfile',form='formatted',status='old')
      read(99,*) enklist(:,:)
!      write(6,*) enklist(:,1)
      close(99)

      open(unit=99,file='efermiinrydberg.ipt',form='formatted',         &
     &     status='old')
      read(99,*)efermi
      close(99)


! these energies are not in eV until later
      vlev = enklist(ibl,1)
      vhev = enklist(ibl,1)
      clev = enklist(ibh+numocc-brange(3)+1,1)
      chev = enklist(ibh+numocc-brange(3)+1,1)
      write(6,*) vlev,vhev,clev,chev 
      do kcount=1,nkpts
        do bandcount=ibl,brange(2)
          if (enklist(bandcount,kcount) .lt. vlev ) then
            vlev = enklist(bandcount,kcount)
          elseif ((enklist(bandcount,kcount) .gt. vhev) .and.           &
     &            (enklist(bandcount,kcount) .lt. efermi ) )  then
            vhev = enklist(bandcount,kcount)
          endif
        enddo
        do bandcount=numocc+1,numocc+ibh-brange(3)+1
          if (enklist(bandcount,kcount) .gt. chev ) then
            chev = enklist(bandcount,kcount)
          elseif ((enklist(bandcount,kcount) .lt. clev ).and.           &
                   enklist(bandcount,kcount) .gt. efermi ) then
            clev = enklist(bandcount,kcount)
          endif
        enddo
      enddo !kpts

      deallocate(enklist)
!      write(6,*) vlev,vhev,clev,chev

! convert energies to eV
      vlev = vlev * ry
      vhev = vhev * ry
      clev = clev * ry
      chev = chev * ry

      return
      end subroutine vcbder2
