program xes_spaghetti


  character (LEN=127) :: filename


  integer :: ZNL(3), nalpha, nb, kmesh(3), nkpts, qpt(3)
  real(kind=kind(1.d0)) :: k0(3), tau(3), rr, ri, ir, ii

  integer :: nptot, ntot
  real(kind=kind(1.d0)), allocatable :: pcr(:,:),pci(:,:),mer(:,:),mei(:,:), &
      re_psi(:,:),im_psi(:,:)

  integer :: iter, icml, ialpha, ikpt, iband, ivms

  open(unit=99,file='ZNL',form='formatted',status='old')
  read(99,*) ZNL(:)
  close(99)
  nalpha = 4 * ( 2 * ZNL(3) + 1 )

  open(unit=99,file='nbuse.ipt',form='formatted',status='old')
  read(99,*) nb
  close(99)

  open(unit=99,file='kmesh.ipt',form='formatted',status='old')
  read(99,*) kmesh(:)
  close(99)
  nkpts = kmesh(1)*kmesh(2)*kmesh(3)

  open(unit=99,file='k0.ipt',form='formatted',status='old')
  read(99,*) k0(:)
  close(99)


  read(5,*) filename

  open(unit=99,file=cks_filename,form='unformatted',status='old')
  rewind( 99 )
  read ( 99 ) nptot, ntot
  read ( 99 ) tau( : )
  allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
  read ( 99 ) pcr
  read ( 99 ) pci
  close( unit=99 )

  allocate( mer( nptot, -ZNL(3):ZNL(3) ), mei( nptot, -ZNL(3):ZNL(3) ) )
  read(5,*) filename
  open( unit=99, file=filename, form='formatted', status='old' )
  rewind( 99 )
  do icml = -ZNL(3), ZNL(3)
    do iter = 1, nptot
      read( 99, * ) mer( iter, icml ), mei( iter, icml )
    enddo
  enddo
  close( 99 )


  allocate( re_psi(nb,nkpts), im_psi(nb,nkpts) )

! no spins for now, multiply by 2.0
!  re_psi = 0
!  im_psi = 0

  do icml = -ZNL(3),ZNL(3)
    iter = 0
    do ikpt = 1, nkpts
      do ib = 1, nb
        iter = iter + 1

        rr = dot_product( mer( :, icml ), pcr( :, iter ) )
        ri = dot_product( mer( :, icml ), pci( :, iter ) )
        ir = dot_product( mei( :, icml ), pcr( :, iter ) )
        ii = dot_product( mei( :, icml ), pci( :, iter ) )
        re_psi(ib,ikpt) = 2.0d0 * ( rr - ii )
        im_psi(ib,ikpt) = 2.0d0 * ( -ri -ir )

      enddo
    enddo
  enddo


  allocate(energies(nb,nkpts))

  open(unit=99,file='wvfvainfo',form='unformatted', status='old' )
  rewind(99)
  read(99) energies
  close(99)


  open(unit=99,file='enk.txt',form='formatted',status='unknown')
  rewind(99)

  ikpt = 0
  do xk = 0, kmesh(1)-1
    qpt(1) = k0(1) + dble(xk)/dble(kmesh(1))
    write(99,*) '{'
    do yk = 0, kmesh(2)-1
      qpt(2) = k0(2) + dble(yk)/dble(kmesh(2))
      write(99,*) '{'
      do zk = 0, kmesh(3)-1
        qpt(3) = k0(3) + dble(zk)/dble(kmesh(3))

        ikpt = ikpt + 1
          
        
        
    
    


end program xes_spaghetti
