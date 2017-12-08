program jdftx_density
  implicit none

  real(kind=kind(1.0d0)), allocatable :: inverted_density(:,:,:)
  integer :: nz, ny, nx, iz, iy, ix
  
  open( unit=99,file='nfft',form='formatted',status='old' )
  read( 99, * ) nx, ny, nz
  close( 99 )

  allocate( inverted_density( nz, ny, nx ) )
  open(unit=99,file='scf.n',form='unformatted',access='stream',status='old')
  read(99) inverted_density
  close(99)

  open(unit=99,file='rhoofr',form='formatted',status='unknown')
  write(99,*)"       i1       i2       i3 data"
  do iz = 1, nz
    do iy = 1, ny
      do ix = 1, nx
        write(99,*) ix, iy, iz, inverted_density(iz,iy,ix)
      enddo
    enddo
  enddo
  close(99)

  deallocate( inverted_density )

end program jdftx_density
