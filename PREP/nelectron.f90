! program nelectron
! gets the number of electrons per unit cell (via the realspace density)
!
      program nelectron
!
      implicit none
      character*20 :: lineburn
      integer :: xpnt(3), counter
      real(kind=kind(1.d0)) :: rho, nelect, avecs(3,3), vol
!
      open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
      read(99,*) avecs(:,:)
      close(99)
      vol = avecs(1,1) * (avecs(2,2)*avecs(3,3)-avecs(3,2)*avecs(2,3))  &
     &    - avecs(2,1) * (avecs(1,2)*avecs(3,3)-avecs(3,2)*avecs(1,3))  &
     &    + avecs(3,1) * (avecs(1,2)*avecs(2,3)-avecs(2,2)*avecs(1,3))
!      write(6,*) vol

      nelect = 0.d0;
      counter = 0
      open(unit=99,file='rhoofr',form='formatted',status='old')
      read(99,*) lineburn
      do
        read(99,*, END = 10) xpnt(:),rho
        nelect = nelect + rho
        counter = counter + 1
      enddo

 10   close(99)
      nelect = nelect / dble(counter) * vol
      write(6,*) nelect, anint(nelect)
      open(unit=99,file='nelectron',form='formatted',status='unknown')
      write(99,*) int(anint(nelect))
      close(99)

      end program nelectron
