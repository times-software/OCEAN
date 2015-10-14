! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program rscombine
  implicit none
  !
  integer :: rc, vc, mc, n1, n2, n3, nr, i
  !
  real( kind = kind( 1.0d0 ) ) :: dr, r, val1, val2, val3
  real( kind = kind( 1.0d0 ) ), allocatable :: dat( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: r1( : ), v1( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: r2( : ), v2( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: r3( : ), v3( : )
  !
  integer :: nvps, nvae
  real( kind = kind( 1.0d0 ) ) :: valps, valae
  logical :: amend
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rvps, vvps, rvae, vvae
  !
  read ( 5, * ) rc, vc
  read ( 5, * ) n1
  mc = max( rc, vc )
  allocate( dat( mc ), r1( n1 ), v1( n1 ) )
  do i = 1, n1
     read ( 5, * ) dat( : )
     r1( i ) = dat( rc )
     v1( i ) = dat( vc ) 
  end do
  deallocate( dat )
  !
  read ( 5, * ) rc, vc
  read ( 5, * ) n2
  mc = max( rc, vc )
  allocate( dat( mc ), r2( n2 ), v2( n2 ) )
  do i = 1, n2
     read ( 5, * ) dat( : )
     r2( i ) = dat( rc )
     v2( i ) = dat( vc ) 
  end do
  deallocate( dat )
  !
  read ( 5, * ) rc, vc
  read ( 5, * ) n3
  mc = max( rc, vc )
  allocate( dat( mc ), r3( n3 ), v3( n3 ) )
  do i = 1, n3
     read ( 5, * ) dat( : )
     r3( i ) = dat( rc )
     v3( i ) = dat( vc ) 
  end do
  deallocate( dat )
  !
  read ( 5, * ) amend
  if ( amend ) then
     read ( 5, * ) nvps
     allocate( rvps( nvps ), vvps( nvps ) )
     do i = 1, nvps
        read ( 5, * ) rvps( i ), vvps( i )
     end do 
     read ( 5, * ) nvae
     allocate( rvae( nvae ), vvae( nvae ) )
     do i = 1, nvae
        read ( 5, * ) rvae( i ), vvae( i )
     end do
  end if
  !
  read ( 5, * ) dr, nr
  !
  open( unit=99, file='rpot', form='formatted', status='unknown' )
  rewind 99
  do i = 0, nr - 1
     r = dr * dble( i )
     call intval( n1, r1, v1, r, val1, 'cap', 'cap' )
     call intval( n2, r2, v2, r, val2, 'cap', 'cap' )
     call intval( n3, r3, v3, r, val3, 'cap', 'cap' )
     if ( amend ) then
        call intval( nvps, rvps, vvps, r, valps, 'cap', 'cap' )
        call intval( nvae, rvae, vvae, r, valae, 'cap', 'cap' )
        val2 = val2 + ( valae - valps )
     end if 
     write ( 6, '(2x,4(1x,1e15.8))' ) r, val1, val2, val1 + val2
     write ( 99, '(2x,2(1x,1e15.8))' ) -( val1 + val2 + val3 ), r
  end do
  close( unit=99 )
  !
  open( unit=99, file='rpothires', form='formatted', status='unknown' )
  rewind 99
  do i = 0, 10 * nr - 1
     r = 0.1d0 * dr * dble( i )
     call intval( n1, r1, v1, r, val1, 'cap', 'cap' )
     call intval( n2, r2, v2, r, val2, 'cap', 'cap' )
     call intval( n3, r3, v3, r, val3, 'cap', 'cap' )
     if ( amend ) then
        call intval( nvps, rvps, vvps, r, valps, 'cap', 'cap' )
        call intval( nvae, rvae, vvae, r, valae, 'cap', 'cap' )
        val2 = val2 + ( valae - valps )
     end if 
     write ( 99, '(2x,2(1x,1e15.8))' ) -( val1 + val2 + val3 ), r
  end do
  close( unit=99 )
  !
  deallocate( r1, v1, r2, v2 )
  !
end program rscombine
