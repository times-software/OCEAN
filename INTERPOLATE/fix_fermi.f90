! sorts energies



subroutine fix_fermi( nbands, nkpts, nspin, nshift, max_val, nelect, ndope, e0, homo, lumo, fermi )
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


  write(6,*) homo, 2.0*energy_list( tot_electrons )
  write(6,*) lumo, 2.0*energy_list( tot_electrons + 1 )
  write(6,*) fermi
  fermi = ( energy_list( tot_electrons ) + energy_list( tot_electrons + 1 ) ) 
  write(6,*) fermi


  homo = 2.0 * energy_list( tot_electrons ) 
  lumo = 2.0 * energy_list( tot_electrons + 1 )
  



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
          if ( (node .lt. top) .and. (energies( node ) .lt. energies( node + 1 ) ) ) node = node + 1
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
