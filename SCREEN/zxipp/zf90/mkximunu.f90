! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program mkximunu
  implicit none
  !
  integer :: nr, nang, nang2, i, j, ir, jr, ia, ja, chan1, chan2
  complex( kind = kind( 1.0d0 ) ) :: csu, rm1
  character(len=8) :: fnam
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: rad( : ), ymu( :, : ), xpt( : ), ypt( : ), zpt( : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: rerow( :, : ), imrow( :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: xi( :, :, :, : )
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
!  read ( 5, * ) chan1, chan2
  chan1 = 1
  chan2 = 1
  !
  open( unit=99, file='projsupp', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nr, nang
  allocate( rad( nr ), xi( nang, nr, nang, nr ), rerow( nang, nr ), imrow( nang, nr ), &
            xpt( nang ), ypt( nang ), zpt( nang ), wpt( nang ) )
  read ( 99, * ) rad( : )
  close( unit=99 )
  !
  open( unit=99, file='specpnt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nang2
  if ( nang .ne. nang2 ) stop 'nang nang2 mismatch in mkximunu'
  do i =1, nang
    read ( 99, * ) xpt( i ), ypt( i ), zpt( i ), wpt( i )
  enddo
  close( unit=99 )
  !
  allocate( ymu( nang, 9 ) )
  call formreytab( nang, xpt, ypt, zpt, ymu, 9 )
!  open( unit=99, file='ytab', form='formatted', status='unknown' )
!  rewind 99
!  do i = 1, 9
!     do j = 1, nang
!        read ( 99, * ) ymu( j, i )
!     end do
!  end do
!  close( unit=99 )
  !
  write ( 6, * ) 'reading ximat'
  open( unit=99, file='ximat', form='unformatted', status='unknown' )
  rewind 99
  do i = 1, nr
     do j = 1, nang
        read ( 99 ) rerow( :, : )
        xi( :, :, j, i ) = rerow( :, : ) 
     end do
  end do 
  close( unit=99 )
  write ( 6, * ) 'done reading ximat'
  !
! do i = 1, 9
!    do j = 1, 9
  i = chan1; j = chan2
        write ( fnam, '(1a6,2i1)' ) '.ymunu', i, j
        open( unit=99, file=fnam, form='formatted', status='unknown' )
        rewind 99
        do ir = 1, nr
           do jr = 1, nr
              csu = 0.0d0
              do ia = 1, nang
                 do ja = 1, nang
                    csu = csu + xi( ia, ir, ja, jr ) * ymu( ia, i ) * wpt( ia ) * ymu( ja, j ) * wpt( ja )
                 end do
              end do
              write ( 99, '(4(1x,1e15.8))' ) rad( ir ), rad( jr ), csu
           end do
           write ( 99, * )
        end do
        close( unit=99 )
!    end do
! end do
  !
end program mkximunu
