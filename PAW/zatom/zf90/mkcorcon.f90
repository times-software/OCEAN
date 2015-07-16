! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine mkcorcon( alfa, rel, z, zc, njrc, rcocc, rsocc, fdirac )
  implicit none
  !
  integer :: njrc( 4 ) ! njrc( l + 1 ) > 0 iff ang mom l was pseudized
  real ( kind = kind( 1.0d0 ) ) :: alfa, rel, z, zc, rcocc( 0 : 3 ) ! z = at no; zc = core chg see notes below for alfa, rel
  real ( kind = kind( 1.0d0 ) ) :: rsocc( 0 : 3 ) ! 
  logical :: fdirac
  !
  integer :: fs( 0 : 3 ), nco, l, n
  character * 6 :: functionale, relat, dirsca
  logical :: df
  !
  open( unit=99, file='atomoptions', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) z, fs, dirsca, relat, functionale
  read ( 99, * ) rcocc( : )
  read ( 99, * ) rsocc( : )
  close( unit=99 )
  !
  open( unit=99, file='valcon', form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) 4
  do l = 0, 3
     write ( 99, '(2i5,1f10.6)' ) fs( l ) + l + 1, l, rcocc( l )
  end do
  close( unit=99 )
  !
  select case( functionale ) ! alfa = 1 implies LDA, alfa = 0 implies HF 
  case( 'lda' )
     alfa = 1.0d0
  case( 'hf' )
     alfa = 0.0d0
  end select
  !
  select case( relat ) ! for relativistic, rel = 1; for nonrel, rel = 0
  case( 'rel' ) 
     rel = 1.0d0
  case( 'nonrel' )
     rel = 0.0d0  
  end select
  !
  select case( dirsca )
  case( 'dirac' )
     df = .true.
  case( 'scalar' )
     df = .false.  
  end select 
  if ( fdirac ) df = .true.
  !
  ! df is diracify; .true. means split j= l +/- 0.5; subshells, .false. does scalar rel or nonrel
  if ( .not. df ) then
     nco = sum( fs( : ) )
  else
     nco = fs( 0 ) + 2 * ( fs( 1 ) + fs( 2 ) + fs( 3 ) )
  end if
  zc = z - 2 * fs( 0 ) - 6 * fs( 1 ) - 10 * fs( 2 ) - 14 * fs( 3 )
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5)' ) nco
  do l = 0, 3
     if ( fs( l ) .gt. 0 ) then
        do n = l + 1, l + fs( l )
           if ( l .eq. 0 ) then
              call dssh( n, l, 0.5d0, 2 )
           else
              if ( df ) then
                 call dssh( n, l, dble( l ) - 0.5d0, 2 * l )
                 call dssh( n, l, dble( l ) + 0.5d0, 2 * l + 2 )
              else
                 call dssh( n, l, dble( l ), 4 * l + 2 )
              end if
           end if
        end do
     end if
  end do
  close( unit=99 )
  !
  open( unit=99, file='skip', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5)' ) 3
  do l = 0, 3
     write ( 99, '(2i5)' ) l, fs( l )
  end do
  close( unit=99 )
  !
  return
end subroutine mkcorcon
!
subroutine dssh( n, l, j, occ )
  implicit none
  !
  integer :: n, l, occ
  real( kind = kind( 1.0d0 ) ) :: j
  !
  write ( 99, '(3i3,1f6.2,2i4)' ) n, l, 0, -j, 1, occ
  !
  return
end subroutine dssh
