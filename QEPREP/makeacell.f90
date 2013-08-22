! program makeacell
!
      program makeacell
!
      implicit none
      integer :: i
      real(kind=kind(1.d0))   :: rscale(3), rprim(3,3)
!
!
      open(unit=99,file='rscale',form='formatted',status='old')
      read(99,*) rscale(:)
      close(99)
!
      open(unit=99,file='rprim',form='formatted',status='old')
      open(unit=98,file='acell',form='formatted',status='unknown')
      do i = 1, 3
         read(99,*) rprim(1:3,i)
         rprim(:,i) = rprim(:,i) * rscale(:)
         write( 98, *) rprim(:,i)
      end do
      close(98)
      close(99)
!

      end program makeacell
