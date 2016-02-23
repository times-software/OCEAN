program hubbard_qe
!
   implicit none
!
    logical :: ldau
    character(len=50) :: magline
    integer :: io, itype
    real(kind=kind(1.0d0)) :: hub_u
!
!
      open(unit=98,file='hubbard_u',form='formatted',status='unknown')
!
      open(unit=99,file='ldau',form='formatted',status='old')
      read(99,*) ldau
!      write(98,'(A13,A7)') 'lda_plus_u = ', ldau

      if ( ldau ) then
         write(98,'(A21)') '  lda_plus_u = .true.'
         do
            read(99,*,IOSTAT=io)  itype, hub_u
            if (io > 0) then
               write(*,*) 'Error reading smag'
               exit
            else if (io < 0) then
               exit
            else
               write(6,*) itype, hub_u
               write(magline,'(A10,I0,A4,F8.5)') 'Hubbard_U(', itype, ') = ', hub_u
               write(98,'(A2,A50)') '  ', magline
            end if
         end do
      else
         write(98,'(A22)') '  lda_plus_u = .false.'
      end if
!
20 continue
!
      close(99)
      close(98)
!
end program hubbard_qe
