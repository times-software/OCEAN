! Copyright (C) 2014, 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! Modification of the dotter program to make modified mels
! 
! The output is a new mels as opposed to the normal output of dotter (which is semi-deprecated for ocean.x)
program tpa_dotter
  implicit none
  !
  integer, parameter :: DP = kind(1.0d0)
  integer :: nptot, ntot, mc, i, lc, idum, nspn, ivms, ispn, k, j
  real( DP ) :: rr, ri, ir, ii, br, bi, tau( 3 ), norm
  real( DP ), allocatable, dimension( : ) :: mer, mei
  real( DP ), allocatable, dimension( :, : ) :: reSecondPhoton, imSecondPhoton
  real( DP ), allocatable, dimension( :, :, : ) :: pcr, pci, reSecondary, imSecondary, reOuter, imOuter
  character(len=11) :: cksnam
  character(len=7) :: secondPhotonName
  character(len=18) :: melInName, melOutName
  logical :: have_spin
  !
  read ( 5, * ) cksnam
  read ( 5, * ) melInName
  read ( 5, * ) secondPhotonName
  read ( 5, * ) melOutName
  !
  inquire(file='nspin',exist=have_spin)
  open( unit=99, file=cksnam, form='unformatted', status='old' )
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

  ! First step is the < v'l'm' | band states > < bands states | v l m >
  ! The use of the cks output (pcr,pci) is that it adds a Fermi factor
  ! Otherwise, the | v l m > states are approximately complete for the 
  ! localized basis and finite energy range that we've specified
  
  allocate( reOuter( nptot, nptot, nspn ), imOuter( nptot, nptot, nspn ) )
  reOuter(:,:,:) = 0.0_dp
  imOuter(:,:,:) = 0.0_dp

  norm = 1.0_dp / ( real( ntot, dp ) )
  ! Testing option here
  if( .false. ) then
    do ispn = 1, nspn
      do i = 1, nptot
        reOuter( i, i, ispn ) = 1.0_dp
      enddo
    enddo
  else
    do ispn = 1, nspn
      do k = 1, ntot
        do j = 1, nptot
          do i = 1, nptot
            ! THIS IS ALSO NO GOOD. CLEARLY SHOULD BE RR + II , +/-(RI +/- IR)
            reOuter( i, j, ispn ) = reOuter( i, j, ispn ) + pcr( j, k, ispn ) * pcr( i, k, ispn ) &
                                  + pci( j, k, ispn ) * pci( i, k, ispn )
            imOuter( i, j, ispn ) = imOuter( i, j, ispn ) + pcr( j, k, ispn ) * pci( i, k, ispn ) &
                                  - pci( j, k, ispn ) * pcr( i, k, ispn )
          enddo
        enddo
      enddo
    enddo
    reOuter(:,:,:) = reOuter(:,:,:) * norm
    imOuter(:,:,:) = imOuter(:,:,:) * norm
  endif

  deallocate( pcr, pci )

    


  
  !
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) idum, idum, lc
  close( unit=99 )
  !
  allocate( reSecondary( nptot, nspn, -lc:lc ), imSecondary( nptot, nspn, -lc:lc ) )
  !
  allocate( mer( nptot ), mei( nptot ) )
  open( unit=99, file=melInName, form='formatted', status='old' )
  rewind 99
  if(nspn .eq. 1) then
     do mc = -lc, lc
        do i = 1, nptot
           read ( 99, * ) mer( i ), mei( i )
        end do
        do i = 1, nptot
           rr = dot_product( mer( : ), reOuter( :, i, 1 ) )
           ri = dot_product( mer( : ), imOuter( :, i, 1 ) )
           ir = dot_product( mei( : ), reOuter( :, i, 1 ) )
           ii = dot_product( mei( : ), imOuter( :, i, 1 ) )
           br = rr + ii
           bi = ri - ir
           reSecondary( i, 1, mc ) = br
           imSecondary( i, 1, mc ) = bi
        end do
     end do
  else
    write(6,*) 'Spin not supported yet!!!!!!!!'
    stop
     do mc = -lc, lc
        do i = 1, nptot
           read ( 99, * ) mer( i ), mei( i )
        end do
        do ivms = 1 , 2
           do i = 1, ntot
!              rr = dot_product( mer( : ), pcr( :, i, ivms ) )
!              ri = dot_product( mer( : ), pci( :, i, ivms ) )
!              ir = dot_product( mei( : ), pcr( :, i, ivms ) )
!              ii = dot_product( mei( : ), pci( :, i, ivms ) )
              br = rr - ii
              bi = -ri - ir
           end do
        enddo
     end do     
  endif
  close(99)

  deallocate( reOuter, imOuter, mer, mei )
  allocate( reSecondPhoton( nptot, nptot ), imSecondPhoton( nptot, nptot ) )
  open(unit=99, name=secondPhotonName, form='formatted', status='old' )
  do i = 1, nptot
    do j = 1, nptot
      read( 99, * ) reSecondPhoton( j, i ), imSecondPhoton( j, i )
    enddo
  enddo
  close( 99 )

  open( unit=99, file=melOutName, form='formatted', status='unknown' )
  if( nspn .eq. 1 ) then
    do mc = -lc, lc
      do i = 1, nptot
        rr = dot_product( reSecondary( :, 1, mc ), reSecondPhoton( :, i ) )
        ri = dot_product( reSecondary( :, 1, mc ), imSecondPhoton( :, i ) )
        ir = dot_product( imSecondary( :, 1, mc ), reSecondPhoton( :, i ) )
        ii = dot_product( imSecondary( :, 1, mc ), imSecondPhoton( :, i ) )
          
        ! secondPhoton is starred
        br = rr + ii
        bi = - ri + ir
  
        write( 99, * ) br, bi
      enddo
    enddo
  else
    stop
  endif
  close( 99 )
  deallocate( reSecondPhoton, imSecondPhoton, reSecondary, imSecondary )
  !
  !
end program tpa_dotter
