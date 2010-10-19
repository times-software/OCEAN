! program mkrbfile
!
! includes: mkcmesh.f90, snatch.f90
! constructs the rbfile.bin
!
!
      program mkrbfile
      implicit none
!
      integer, parameter :: stdin = 5
      character *2 :: element
      integer :: nang, nr, indx
      real(kind=kind(1.d0)) :: rmax, avec(3,3)
      real(kind=kind(1.d0)), allocatable :: posn(:,:), wpt(:), drel(:)
!
!
      read(stdin,*) rmax, nr, element, indx
!
      open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
      read(99,*) avec(:,:)
      close(99)
!
      open(unit=99,file='specpnt',form='formatted',status='old')
      read(99,*) nang
      close(99)
!
      allocate( posn(3,nr*nang), wpt(nr*nang), drel(nr*nang) )
!
      call mkcmesh(nang, nr, rmax, element, indx, posn, wpt, drel, avec)
!
      open(unit=99,file='rbfile.bin',form='unformatted')
        rewind 99
        write(99) nr*nang, rmax
        write(99) posn
        write(99) wpt
        write(99) drel
      close(99)
!
      deallocate(posn, wpt, drel )
!
      end program mkrbfile    
