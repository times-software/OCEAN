! sorts energies

subroutine dump_energies( band_subset, nbands, nkpts, nspin, nshift, e0, lumo_shift, start_band, ierr )
  use kinds, only : dp
  implicit none

  integer, intent( in ) :: band_subset(2), nbands, nkpts, nspin, nshift, start_band( nkpts, nspin, nshift )
  real(dp), intent( in ) :: e0( nbands, nkpts, nspin, nshift ), lumo_shift
  integer, intent( inout ) :: ierr
  integer, external :: freeunit
  !
  integer :: fh, ispin, ik, ibd, ishift, brange(4), nbuse
  real(dp), allocatable :: temp_energy(:,:,:)


  write(6,*) lumo_shift

  fh = freeunit()
  open(unit=fh,file='ibeg.h',form='formatted',status='unknown')!,buffered='yes')
  do ishift = 1, nshift
    do ispin = 1, nspin
      write(fh,'(2I8)') (ik,start_band(ik,ispin,ishift), ik=1,nkpts)
    enddo
  enddo
  close(fh)

  brange(1) = band_subset(1)
  brange(4) = band_subset(2)
  brange(2) = maxval( start_band )
  brange(2) = brange(2) - 1
  brange(3) = minval( start_band )

  write(6,*) maxval( start_band ), minval( start_band )

  open(unit=fh,file='brange.ipt',form='formatted',status='unknown')!,buffered='yes')
  write(fh,*) brange(1), brange(2)
  write(fh,*) brange(3), brange(4)
  close(fh)

  nbuse = brange( 4 ) - brange( 2 ) 
  open(unit=fh,file='nbuse.ipt',form='formatted',status='unknown')
  write(fh,*) nbuse
  close(fh)

  allocate( temp_energy( nbuse, nkpts, nspin ) )
  do ispin = 1, nspin
    do ik = 1, nkpts
      temp_energy( :, ik, ispin ) = &
          e0( start_band( ik, ispin, nshift ) : start_band( ik, ispin, nshift ) + nbuse - 1, &
              ik, ispin, nshift ) - lumo_shift
    enddo
  enddo
  open(unit=fh,file='wvfcninfo',form='unformatted',status='unknown')!,buffered='yes')
  write(fh) nbuse, nkpts, nspin
  write(fh) temp_energy
  close(fh)


  nbuse = brange( 2 ) - brange( 1 ) + 1
  open(unit=fh,file='nbuse_xes.ipt',form='formatted',status='unknown')
  write(fh,*) nbuse
  close(fh) 

  deallocate( temp_energy )
  allocate( temp_energy( nbuse, nkpts, nspin ) )
  temp_energy = 0.0_dp
  do ispin = 1, nspin 
    do ik = 1, nkpts
      temp_energy( :, ik, ispin ) = e0( 1 : start_band( ik, ispin, 1 ) -1, ik, ispin, 1 ) - lumo_shift
    enddo
  enddo
  open(unit=fh,file='wvfvainfo',form='unformatted',status='unknown')!,buffered='yes')
  write(fh) nbuse, nkpts, nspin
  write(fh) temp_energy
  close(fh)

  deallocate( temp_energy )

end subroutine dump_energies



subroutine find_startband( nbands, nkpts, nspin, nshift, nelect, fermi, band_style, e0, start_band )
  use kinds, only : dp
  implicit none

  integer, intent( in ) :: nbands, nkpts, nspin, nshift, nelect
  real(dp), intent( in ) :: e0( nbands, nkpts, nspin, nshift ), fermi
  character(len=5), intent( in ) :: band_style

  integer, intent( out ) :: start_band( nkpts, nspin, nshift )
  !

  integer :: ispin, ik, ielect, ibd
  

  if( band_style .eq. 'bands' ) then
    ielect = nelect / 2 + 1
    do ispin = 1, nspin
      do ik = 1, nkpts
        start_band( ik, ispin, 1 ) = ielect
      enddo
    enddo
    if( nshift .eq. 2 ) then
      do ispin = 1, nspin
        do ik = 1, nkpts
          start_band( ik, ispin, 2 ) = ielect
        enddo
      enddo
    endif
    return
  endif

  do ispin = 1, nspin
    do ik = 1, nkpts
      do ibd = 1, nbands
        if( e0( ibd, ik, ispin, 1 ) .gt. fermi ) then
          start_band( ik, ispin, 1 ) = ibd
          goto 10
        endif
      enddo
10    continue
    enddo
  enddo

  if( nshift .eq. 1 ) return


  do ispin = 1, nspin
    do ik = 1, nkpts
      do ibd = 1, nbands
        if( e0( ibd, ik, ispin, 2 ) .gt. fermi ) then
          start_band( ik, ispin, 2 ) = ibd
          goto 11
        endif
      enddo
11    continue
    enddo
  enddo


  return

end subroutine find_startband
  



subroutine fix_fermi( nbands, nkpts, nspin, nshift, max_val, nelect, ndope, e0, homo, lumo, fermi )
  use kinds, only : dp
  implicit none

  integer, intent( in ) :: nbands, nkpts, nspin, nshift, max_val, nelect, ndope

  real(8), intent( in ) :: e0( nbands, nkpts, nspin, nshift )

  real(8), intent( inout ) :: homo, lumo, fermi



  integer :: tot_electrons
  real(8) :: energy_list( max_val * nkpts * nspin )


  integer :: ispn, ik, ii, jj




  tot_electrons = ( nelect * nkpts * nspin / 2.0 ) + ( ndope )


  do ispn = 1, nspin
    do ik = 1, nkpts
      ii = (ik-1)*max_val + (ispn - 1 ) * nkpts * max_val + 1
      jj = ii + max_val - 1
      energy_list( ii : jj ) = e0( 1 : max_val, ik, ispn, 1 )
    enddo
  enddo

  call sort_energies( max_val * nkpts * nspin, energy_list )


!  write(6,*) homo, 2.0*energy_list( tot_electrons )
!  write(6,*) lumo, 2.0*energy_list( tot_electrons + 1 )
!  write(6,*) fermi
  fermi = ( energy_list( tot_electrons ) + energy_list( tot_electrons + 1 ) ) / 2.0_DP


  homo = 2.0_dp * energy_list( tot_electrons ) 
  lumo = 2.0_dp * energy_list( tot_electrons + 1 )
  write(6,*) fermi, lumo*0.5_dp
  



end subroutine 


subroutine sort_energies( n, energies )
  implicit none
  integer, intent( in ) :: n
  real(8), intent( inout ) :: energies(n)

  integer :: iter, node, node2, top
  real(8) :: temp

      write(6,*) 'sorting'
      top = n
      do iter = n / 2 , 1, -1
        temp = energies( iter )
        node = iter + iter
       node2 = iter
        do
          if ( node .gt. top ) goto 10
          if ( (node .lt. top) .and. (energies( node ) .lt. energies( min(node + 1,top) ) ) ) node = node + 1
          if ( temp .lt. energies( node ) ) then
            energies( node2 ) = energies( node )
            energies( node ) = temp
            node2 = node
            node = node + node
          else
            goto 10
          endif
        enddo
 10 continue
      enddo
      do iter = n, 1, -1
        temp = energies( iter )
        energies( iter ) = energies( 1 )
        node = 2
        node2 = 1
        do
          if ( node .gt. iter - 1) goto 20
          if ( ( node .lt. iter - 1 ) .and. ( energies( node ) .lt. energies( node + 1 ) ) ) node = node + 1
          if ( temp .lt. energies( node ) ) then
            energies( node2 ) = energies( node )
            energies( node ) = temp
            node2 = node
            node = node + node
          else
            goto 20
          endif
        enddo
  20 continue

      enddo

end subroutine
