! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
      subroutine snatch( element, indx, tau )
      implicit none
!
      integer :: indx
      real( kind = kind( 1.0d0 ) ) :: tau( 3 )
      character(len=2) :: element
!
      integer :: ntot, nmatch, i
      real( kind = kind( 1.0d0 ) ) :: tmp( 3 )
      logical :: found
      character(len=2) :: ein
!
      open( unit=99, file='xyz.wyck', form='formatted', status='old' )
      rewind 99
      read ( 99, * ) ntot
      nmatch = 0
      found = .false.
      do i = 1, ntot
        read ( 99, * ) ein, tmp( : )
        if ( ein .eq. element ) then
          nmatch = nmatch + 1
          if ( nmatch .eq. indx ) then
            tau( : ) = tmp( : )
            found = .true.
          end if
        end if
        if ( found ) exit
      end do
      if ( .not. found ) then 
        write(6,*) element, indx
        stop 'atom coord not found!'
      endif
      write ( 6, '(1a15,3f20.15)' ) 'snatched alpha=', tau( : )
!
      return
      end subroutine snatch

