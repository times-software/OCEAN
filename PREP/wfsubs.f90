! Opens the wavefunction file as specified to the file handle
      subroutine getwfkin( wfkin, files_iter, handle)
      implicit none
!
      character(len=11), intent(out) :: wfkin
      integer, intent(in) :: files_iter, handle
!
      write(wfkin(1:11),'(A3,I4.4,A4)') 'RUN', files_iter, '_WFK'
      open(unit=handle, file=wfkin, form='unformatted', status='old')
      rewind(handle)
!
      end subroutine getwfkin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Opens up the wavefunction file to be written out
      subroutine getwfkout( wfkout, files_iter, kptnum, handle)
      implicit none
!
      character(len=12), intent(out) :: wfkout
      integer, intent(in) :: kptnum, handle, files_iter
!
      write(wfkout,'(A4,I3.3,A1,I4.4)') '.Psi', files_iter, '.', kptnum
      open(unit=handle,file=wfkout,form='unformatted',status='unknown')
      rewind(handle)
!
      end subroutine getwfkout
