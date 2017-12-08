
program jdftx_energy
  implicit none
  real(kind=kind(1.0d0)), allocatable :: energies(:,:), shiftedEnergies(:,:)
  integer :: nkpts, brange(4), ik
  logical :: split_dft

  open(unit=99,file='dft.split',form='formatted',status='old')
  read(99,*) split_dft
  close(99)

  open(unit=99,file='nkpts',form='formatted',status='old')
  read(99,*) nkpts
  close(99)

  open(unit=99,file='brange.ipt',form='formatted',status='old')
  read(99,*) brange(1:4)
  close(99)


  if( split_dft ) then
    allocate( energies( brange(2), nkpts ), shiftedEnergies( brange(4), nkpts ) )
    open(99,file='nscf.eigenvals',form='unformatted',access='stream')
    read(99) energies
    close(99)

    open(99,file='nscf.eigenvals.shift',form='unformatted',access='stream')
    read(99) shiftedEnergies
    close(99)
  else
    allocate( energies( brange(2), nkpts ), shiftedEnergies( brange(4), nkpts ) )
    open(99,file='nscf.eigenvals',form='unformatted',access='stream')
    read(99) shiftedEnergies
    close(99)

    do ik = 1, nkpts
      energies(:,ik) = shiftedEnergies(1:brange(2),ik)
    enddo
  endif

  open(unit=99,file='enkfile', form='formatted',status='unknown')
  do ik = 1, nkpts
    write(99,*) 2.0d0 * energies(brange(1):brange(2), ik )
    write(99,*) 2.0d0 * shiftedEnergies(brange(3):brange(4), ik )
  enddo
  close(99)

  open(unit=99,file='enk_un', form='formatted',status='unknown')
  do ik = 1, nkpts
    write(99,*) 2.0d0 * energies(:, ik )
  enddo
  close(99)

  open(unit=99,file='enk_sh', form='formatted',status='unknown')
  do ik = 1, nkpts
    write(99,*) 2.0d0 * shiftedEnergies(:, ik )
  enddo
  close(99)


  deallocate( energies, shiftedenergies )

end program jdftx_energy
