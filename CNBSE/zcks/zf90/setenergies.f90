! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine setenergies( energies, nbands, style )
  implicit none
  !
  integer, intent( in ) :: nbands, style
  real( kind = kind( 1.0d0 ) ), intent( inout ) :: energies( nbands )
  !
  !
  !
  select case( style )
    case( 0 )
      call setenergy_shiftscale( energies, nbands )
    case( 1 )
      call setenergy_mpse( energies, nbands )
    case default
      stop
  end select 
  !
  return
end subroutine setenergies
!
!
subroutine setenergies_shiftscale( energies, nbands )
  implicit none
  !
  integer, intent( in ) :: nbands
  real( kind = kind( 1.0d0 ) ), intent( inout ) :: energies( nbands )
  !
  real( kind = kind( 1.0d0 ) ) :: edge, sc, eshift
  !
  open( unit=99, file='cksshift', form='formatted', status='old')
  read( 99, * ) edge
  close( 99 )
  !
  open( unit=99, file='cksstretch', form='formatted', status='old')
  read( 99, * ) sc
  close( 99 )
  !
  open( unit=99, file='eshift.ipt', form='formatted', status='old' )
  read ( 99, * ) eshift
  close( unit=99 ) 
  !
  w( : ) = ( edge + sc * ( w( : ) * 13.6057d0 + eshift - edge ) ) / 27.2114d0
  !
  return
end subroutine setenergies_shiftscale
!
!
subroutine setenergies_mpse( energies, nbands )
  implicit none
  !
  integer, intent( in ) :: nbands
  real( kind = kind( 1.0d0 ) ), intent( inout ) :: energies( nbands )
  !
  real( kind = kind( 1.0d0 ) ) :: efermi, ldagap, egrid( 2 ), re_se( 2 ), im_se( 2 ), newcbm
  integer :: iter, file_err
  logical :: SEfile_fault
  !
  !
  open( unit=99, file='efermiinrydberg.ipt', form='formatted', status='old' )
  read( 99, * ) efermi
  close( 99 )
  !
  open( unit=99, file='ldagap', form='formatted', status='old' )
  read( 99, * ) ldagap
  close( 99 )
  !
!  open( unit=99, file='cksshift', form='formatted', status='old')
!  read( 99, * ) edge
!  close( 99 )
  !
  open( unit=99, file='eshift.ipt', form='formatted', status='old' )
  read ( 99, * ) eshift
  close( unit=99 )
  eshift = eshift * -1.0d0
  !
  open( unit=99, file='SelfEnergy.dat', form='formatted', status='old' ) !, ERR=20, iostat=file_err )
  read( 99, * ) egrid( 1 ), re_se( 1 ), im_se( 1 )
  read( 99, * ) egrid( 2 ), re_se( 2 ), im_se( 2 )
  !
  do iter = 1, nbands
    if( eshift .gt. egrid( 1 ) .and. eshift .lt. egrid( 2 ) ) then
      newcbm = eshift +  re_se( 1 ) + ( ( re_se( 2 ) - re_se( 1 ) ) / ( egrid( 2 ) - egrid( 1 ) ) &
          * ( tmp - egrid( 1 ) )
    endif
    !
    if( energies( iter ) .gt. re_se( 2 ) ) then
      egrid( 1 ) = egird( 2 ); re_se( 1 ) = re_se( 2 ); im_se( 1 ) = im_se( 2 )
    else
      tmp = energies( nband )
      ! energies come in in Ryd, but ldagap and SE are in eV
      tmp = ( tmp - efermi ) * 13.6057d0 + 0.5d0 * ldagap
      tmp = tmp + re_se( 1 ) + ( ( re_se( 2 ) - re_se( 1 ) ) / ( egrid( 2 ) - egrid( 1 ) ) &
          * ( tmp - egrid( 1 ) )
      energies( nband ) = tmp
    endif
  enddo
  close( 99 )
  !
  energies( : ) = ( energies( : ) - newcbm ) / 27.2114d0
  !
  return
end subroutine setenergies_mpse
