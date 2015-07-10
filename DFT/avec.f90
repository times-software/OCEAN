program avec
  implicit none
  double precision :: ascale(3),bvec(3,3),pi
  character(len=5) :: gtxt
  integer :: i,j
  pi = 4.d0* datan(1.d0)
  open(unit=17,file='rcell',form='formatted',status='unknown')
   read(17,*)gtxt,ascale(1),ascale(2),ascale(3)
   read(17,*)gtxt,bvec(1,1),bvec(1,2),bvec(1,3)
  do i=2,3
   read(17,*)bvec(i,1),bvec(i,2),bvec(i,3)
  enddo
  close(17)
  open(unit=17,file='AVECS',form='formatted',status='unknown')
  do i=1,3
  write(17,*) (ascale(i)*bvec(i,j)/(2*pi),j=1,3)
  enddo
  close(17)
  open(unit=17,file='avecsinbohr.ipt',form='formatted',status='unknown')
  do i=1,3
  write(17,*) (ascale(i)*bvec(i,j),j=1,3)
  enddo
  close(17)
end

