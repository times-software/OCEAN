! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! J. Vinson
!
! Simple/dumb routine to read in an echamp file (assuming currently it is from rixs, ie valence )
!  and then it dumps the energy histograms of the underlying DFT states that make up that exciton.
program rixs_energy_scatter
  implicit none
  complex(kind=kind(1.0d0)), allocatable :: vec(:,:,:)
  real(kind=kind(1.0d0)), allocatable :: val_energies(:,:), con_energies(:,:)
  real(kind=kind(1.0d0)), allocatable :: val_vec(:,:), con_vec(:,:), val_en_out(:), con_en_out(:)
  

  real(kind=kind(1.0d0)) :: deltaE, min_E, max_E

  integer :: brange(4), kpts(3), nbv, nbc, nkpts, E_count
  integer :: ikpt, ib, ibc, ie, etarg

  character(len=128) :: echamp_file



  write(6,*) 'Enter in file name of echamp'
  read(5,*) echamp_file

  open(unit=98,file=echamp_file, form='unformatted',status='old')

  write(6,*) 'Enter in spacing of energy grid for output'
  read(5,*) deltaE

  open( 99, file="brange.ipt", status="old", form='formatted' )
  read(99,*) brange(1:2)
  read(99,*) brange(3:4)
  close(99)
  nbv = brange(2) - brange(1) + 1
  nbc = brange(4) - brange(3) + 1

  open( unit=99,file='kmesh.ipt', status="old", form='formatted' )
  read(99,*) kpts(:)
  close(99)
  nkpts = product(kpts(:))


  allocate( val_energies( nbv, nkpts ), con_energies( nbc, nkpts ) )
  open( unit=99, file='enkfile',form='formatted',status='old')
  do ikpt = 1, nkpts
    read(99,*) val_energies( :, ikpt )
    read(99,*) con_energies( :, ikpt )
  enddo

  
  ! Energies are in Ryd
  val_energies( :, : ) = val_energies( :, : ) * 13.60569253d0
  con_energies( :, : ) = con_energies( :, : ) * 13.60569253d0
  ! now energies are in eV

  ! go ahead and load up the ehcamp
  allocate( vec( nbc, nbv, nkpts ) )
  read(98) vec
  close(98)
  ! close the echamp file


  ! figure out the energy range of the DFT states
  min_E = minval( val_energies )
  max_E = maxval( con_energies )
  E_count = ( max_E - min_E + deltaE ) / deltaE

  write(6,*) 'Energy grid: '
  write(6,*) min_E, max_E, E_count



  ! Real arrays that will hold the exciton vector**2 and summed over valence or conduction bands
  allocate( val_vec( nbv, nkpts ), con_vec( nbc, nkpts ) )
  val_vec = 0.0d0
  con_vec = 0.0d0

  do ikpt = 1, nkpts
    do ib = 1, nbv
      val_vec( ib, ikpt ) = sum( vec( :, ib, ikpt ) * conjg( vec( :, ib, ikpt ) ) )
      do ibc = 1, nbc
        con_vec( ibc, ikpt ) = con_vec( ibc, ikpt ) + vec( ibc, ib, ikpt ) * conjg( vec( ibc, ib, ikpt ) )
      enddo
    enddo
  enddo



  ! Energy histogram arrays
  allocate( val_en_out(0:E_count), con_en_out(0:E_count) )
  val_en_out = 0.0d0
  con_en_out = 0.0d0

  do ikpt = 1, nkpts
    do ib = 1, nbv
      etarg = ( val_energies( ib, ikpt ) - min_E ) / deltaE
      val_en_out( etarg ) = val_en_out( etarg ) + val_vec( ib, ikpt )
    enddo
    do ib = 1, nbc
      etarg = ( con_energies( ib, ikpt ) - min_E ) / deltaE
      con_en_out( etarg ) = con_en_out( etarg ) + con_vec( ib, ikpt )
    enddo
  enddo



  ! write out both valence and conduction to file
  open(unit=99, file='test.txt', form='formatted', status='unknown' )
  do ie = 0, E_count
    write(99,*) (min_E + dble(ie)*deltaE), val_en_out(ie), con_en_out(ie)
  enddo
  close(99)


end program rixs_energy_scatter
