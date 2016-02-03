! Copyright (C) 2014 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program dotter
  implicit none
  !
  integer :: nptot, ntot, mc, i, lc, idum, nspn, ivms, iufmt, ifmt
  real( kind = kind( 1.0d0 ) ) :: rr, ri, ir, ii, br, bi, tau( 3 )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: mer, mei
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: pcr, pci
  character * 8 :: str
  character * 9 :: str2
  character * 9 :: filnam
  logical :: have_spin
  !
  read ( 5, * ) filnam
  inquire(file='nspin',exist=have_spin)
  open( unit=99, file=filnam, form='unformatted', status='unknown' )
  rewind 99
  if( have_spin ) then
    read ( 99 ) nptot, ntot, nspn
  else
    read ( 99 ) nptot, ntot
    nspn = 1
  endif
  read ( 99 ) tau( : )
  allocate( pcr( nptot, ntot, nspn ), pci( nptot, ntot, nspn ) )
  read ( 99 ) pcr
  read ( 99 ) pci
  close( unit=99 )
  !
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) idum, idum, lc
  close( unit=99 )
  !
  allocate( mer( nptot ), mei( nptot ) )
  open( unit=99, file='mels', form='formatted', status='unknown' )
  rewind 99
  if(nspn .eq. 1) then
     do mc = -lc, lc
        do i = 1, nptot
           read ( 99, * ) mer( i ), mei( i )
        end do
        write ( str, '(1a4,1i2.2)' ) 'beff', 1 + mc + lc
        write ( str2, '(1a5,1i2.2)' ) 'abeff', 1 + mc + lc
        open( unit=98, file=str, form='unformatted', status='unknown' )
        open( unit=97, file=str2, form='formatted', status='unknown' )
        rewind 98
        write ( 98 ) tau( : )
        rewind 97
        write ( 97, '(3(1x,1e15.8))' ) tau( : )
        do i = 1, ntot
           rr = dot_product( mer( : ), pcr( :, i, 1 ) )
           ri = dot_product( mer( : ), pci( :, i, 1 ) )
           ir = dot_product( mei( : ), pcr( :, i, 1 ) )
           ii = dot_product( mei( : ), pci( :, i, 1 ) )
           br = rr - ii
           bi = -ri - ir
           write ( 98 ) br, bi
           write ( 97, '(2(1x,1e15.8))' ) br, bi
        end do
     end do
  else
     do mc = -lc, lc
        do i = 1, nptot
           read ( 99, * ) mer( i ), mei( i )
        end do
        do ivms = 1 , 2
           iufmt=98
           ifmt=97
           write ( str, '(1a4,1i2.2,1a1,1i1.1)' ) 'beff', 1 + mc + lc, '.', ivms
           write ( str2, '(1a5,1i2.2,1a1,1i1.1)' ) 'abeff', 1 + mc + lc, '.', ivms
           open( unit=iufmt, file=str, form='unformatted', status='unknown' )
           open( unit=ifmt, file=str2, form='formatted', status='unknown' )
           rewind iufmt
           write ( iufmt ) tau( : )
           rewind ifmt
           write ( ifmt, '(3(1x,1e15.8))' ) tau( : )
           do i = 1, ntot
              rr = dot_product( mer( : ), pcr( :, i, ivms ) )
              ri = dot_product( mer( : ), pci( :, i, ivms ) )
              ir = dot_product( mei( : ), pcr( :, i, ivms ) )
              ii = dot_product( mei( : ), pci( :, i, ivms ) )
              br = rr - ii
              bi = -ri - ir
              write ( iufmt ) br, bi
              write ( ifmt, '(2(1x,1e15.8))' ) br, bi
           end do
           close(iufmt)
           close(ifmt)
        enddo
     end do     
  endif
  !
end program dotter
