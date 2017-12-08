! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
 program shifted_conugtoux
   implicit none
   !
   integer, parameter :: u1dat = 22, iwf = 23
   integer :: nk, nb, nx( 3 ), ng, i, nwfile, iwfile, j, k, nspin, ispin, nbv, nbc, nbands
   integer, allocatable :: g( :, : ), ikl( :, : ), flip_g( :, : )
   real( kind = kind( 1.0d0 ) ), allocatable :: zr( :, : ), zi( :, : )
   real(kind=kind(1.0d0)) :: su
   character(len=20), allocatable :: wfnam( : )
   character(len=40) :: fileName
   character(len=3) :: DFT
   logical :: is_jdftx
   !
   open( unit=99, file='kandb.h', form='formatted', status='unknown' )
   call igetval( nk, 'nk   ' )
   call igetval( nb, 'nb   '  )
   call igetval( nbv, 'nbv  ' )
   call igetval( nbc, 'nbc  ' )
   call igetval( nx( 1 ), 'ngx  ' )
   call igetval( nx( 2 ), 'ngy  ' )
   call igetval( nx( 3 ), 'ngz  ' )
   call igetval( nspin,   'nspin' )
   close( unit=99 )
   !
   inquire( file='dft', exist=is_jdftx )
   if( is_jdftx ) then
     open( unit=99, file='dft', form='formatted', status='old' )
     read(99,*) DFT
     close(99)
     if( DFT .eq. 'jdf' ) then
       is_jdftx = .true. 
     else 
       is_jdftx = .false.
     endif
   endif
   !
   open( unit=99, file='masterwfile', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) nwfile
   close( unit=99 )
   allocate( wfnam( nwfile ), ikl( 2, nwfile ) )
   open( unit=99, file='listwfile', form='formatted', status='unknown' )
   rewind 99
   do i = 1, nwfile
      read ( 99, * ) ikl( 1, i ), wfnam( i )
      if ( i .gt. 1 ) ikl( 2, i - 1 ) = ikl( 1, i ) - 1
   end do
   close( unit=99 )
   ikl( 2, nwfile ) = nk
   !
   open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown' )
   rewind u1dat
   iwfile = 1
   k = 0
   do ispin = 1, nspin
     do i = 1, nk
        k = k + 1
        open( unit=99, file='prog', form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(2x,1a6,3x,4i5)' ) 'conv03g', i,nk, ispin, nspin
        close( unit=99 )

        ! here we double-up, first read shifted, then non-shifted
        filePrefix = 'unshifted/'
        nbands = nbv
        do ishift = 1, 2

           write(fileName,'(A,A)') filePrefix, wfnam( iwfile )
           write(6,*) fileName
           if( is_jdftx ) then
             open ( unit=iwf, file=fileName, form='stream',status='old'  )
             read( iwf ) ng
           else
             open ( unit=iwf, file=fileName, form='unformatted',status='old' )
             rewind iwf
             read ( iwf ) ng
           endif
           write(6,*) i, ng, nbands
           allocate( g( ng, 3 ), zr( ng, nbands ), zi( ng, nbands ) )
           if( is_jdftx ) then
             allocate( flip_g( 3, ng ) )
             read( iwf ) flip_g
             g = transpose( flip_g )
             deallocate( flip_g )
           else
             read ( iwf ) g
           endif

          read ( iwf ) zr
          read ( iwf ) zi

          if( i .eq. 1 .and. is_jdftx ) then
            do j = 1, ng
              write(31,'(3I6,4X,E23.16,4X,E23.16)') g(j,1), g(j,2), g(j,3), zr(j,1), zi(j,1)
            enddo
          endif
          if( i .eq. 1 ) then
            su = 0
            do j = 1, ng
              su = su + zr(j,1)**2 + zi(j,i)**2
            enddo
            write(6,*) 'Norm:', su
          endif

          call gentoreal( nx, nb, zr, zi, ng, g, u1dat, ( ( i .eq. 1) .and. ( ispin .eq. 1 ) ) )
      
          prefix = 'shifted/'
          nbands = nb

          deallocate( g, zr, zi )
        enddo
        iwfile = iwfile + 1
     end do
   end do
   close( unit=u1dat )
   !
   deallocate( wfnam, ikl )
   !
 end program conugtoux
