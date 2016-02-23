program smag_qe
!
   implicit none
!
    character(len=50) :: magline
    integer :: io, itype
    real(kind=kind(1.0d0)) :: smag
!
!
      open(unit=98,file='qesmag',form='formatted',status='unknown')
!
      open(unit=99,file='smag',form='formatted',status='old')
      do
         read(99,*,IOSTAT=io)  itype, smag
         if (io > 0) then
            write(*,*) 'Error reading smag'
            exit
         else if (io < 0) then
            exit
         else
            write(6,*) itype, smag
            write(magline,'(A23,I0,A4,F8.5)') 'starting_magnetization(', itype, ') = ', smag
            write(98,'(A2,A50)') '  ', magline
         end if
      end do
!
20 continue
!
      close(99)
      close(98)
!
end program smag_qe
