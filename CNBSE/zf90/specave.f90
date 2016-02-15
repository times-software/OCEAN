! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program specave
  implicit none
  !
  character*50 :: pseudo
  character*2  :: elm
  character*4  :: edgename
  character*11 :: filebase
  character*21 :: filename
  character*22 :: outname1
  character*21 :: outname2
  character*18 :: outname3
  !
  integer :: nedges, npol, nlines
  integer :: znucl, ncore, lcore
  integer :: i,j,l
  integer, allocatable :: edgenum(:)
  !
  real(kind=kind(1.0d0)) :: absint
  real(kind=kind(1.0d0)), allocatable :: energy(:), tot_abs(:)
  real(kind=kind(1.0d0)), allocatable :: tot_pol(:,:), tot_atom(:,:)
  real(kind=kind(1.0d0)) :: dummy1, dummy2, dummy3, dummy4, dummy5, dummy6


  open( unit=99, file='nedges', form='formatted', status='old' )
  rewind 99
  read( 99, *) nedges
  close( 99 )

  allocate( edgenum(nedges) )

  open( unit=99, file='cnbse.ways', form='formatted', status='old' )
  rewind 99
  read( 99, *) npol
  close( 99 )

  open( unit=99, file='cnbse.spect_range', form='formatted', status='old' )
  rewind 99
  read( 99, *) nlines, dummy1, dummy2
  close( 99 )

  open( unit=99, file='hfinlist', form='formatted', status='old' )
  rewind 99
  do i = 1, nedges
    read( 99, *) pseudo, znucl, ncore, lcore, elm, edgenum(i)
  end do
  close( 99 )

  write(filebase,"(A8,A2,A1)") 'absspct_', elm, '.' 

  if ( lcore .eq. 0 ) then
     write(edgename,"(A1,I1,A2)") '_', ncore, 's_'
  else if ( lcore .eq. 1 ) then
     write(edgename,"(A1,I1,A2)") '_', ncore, 'p_'
  else if ( lcore .eq. 2 ) then
     write(edgename,"(A1,I1,A2)") '_', ncore, 'd_'
  end if


  allocate( energy(nlines), tot_abs(nlines) )
  allocate( tot_pol(npol,nlines), tot_atom(nedges,nlines) )
  energy(:) = 0.0d0
  tot_abs(:) = 0.0d0
  tot_pol(:,:) = 0.0d0
  tot_atom(:,:) = 0.0d0

  do i = 1, nedges

     do j = 1, npol

        write(filename,"(A11,I0.4,A4,I0.2)") filebase, edgenum(i), edgename, j
        open( unit=99, file=filename, form='formatted', status='old' )
        rewind 99
        do l = 1, nlines
          read( 99, *) energy(l), absint, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
          tot_abs(l) = tot_abs(l) + ( absint / (npol*nedges) )
          tot_atom(i,l) = tot_atom(i,l) + ( absint / npol )
          tot_pol(j,l) = tot_pol(j,l) + ( absint / nedges )
        end do
        close( 99 )

     end do

  end do

  ! output polarization average for each absorption site
  do i = 1, nedges
     write(outname1,"(A11,I0.4,A4,A3)") filebase, edgenum(i), edgename, 'ave'
     open( unit=98, file=outname1, form='formatted', status='unknown' )
     do l = 1, nlines
        write( 98, *) energy(l), tot_atom(i,l)
     end do
     close( 98 )
  end do

  ! output site average for each polarization
  do j = 1, npol
     write(outname2,"(A11,A4,A4,I0.2)") filebase, edgename, 'ave_', j
     open( unit=98, file=outname2, form='formatted', status='unknown' )
     do l = 1, nlines
        write( 98, *) energy(l), tot_pol(j,l)
     end do
     close( 98 )
  end do

  ! output the total average
   write(outname3,"(A11,A4,A3)") filebase, edgename, 'ave'
   open( unit=98, file=outname3, form='formatted', status='unknown' )
   do l = 1, nlines
      write( 98, *) energy(l), tot_abs(l)
   end do
   close( 98 )

   deallocate( edgenum )
   deallocate( energy, tot_abs )
   deallocate( tot_pol, tot_atom )
  !
end program specave
