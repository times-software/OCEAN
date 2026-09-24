! Copyright (C) 2024, 2026 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 10-2024
!
!
program coreExchange
  use ai_kinds, only : DP
  implicit none


  integer :: ZZ, nc, lc, nptot, ntot, nspin, lmin, lmax, l, dumi, npt, maxll, ik, iks, nband, nspin2, ib
  integer :: kgl, kgh, l1, l2, l3, l4, m1, m2, m3, m4, mk, k, ip, jp, istart, brange(4)
  integer :: i, j, nu2, nu3, kmesh(3), nk, isite, iisite, nsite, is, info, i2, kk, nprojmax, mkd, kkmin
  integer, allocatable :: nproj(:)
  real(DP) :: dumr, pi, su, yp( 0 : 1000 ), avecs(3,3), omega, efermi
  real(DP), allocatable :: x(:), w(:), gk(:,:,:), WW(:), rwork(:), gkk(:,:,:,:,:), fkk(:,:,:,:,:), enk(:,:), enkskip(:)
  complex(DP) :: f1, f2, f3, coef
  complex(DP), allocatable :: cks(:,:,:), denMat(:,:), dmat(:,:), VXX(:,:,:), rho(:,:,:), work(:), VF2(:,:)
  logical, parameter :: yes = .true.
  logical, parameter :: no = .false.
  character(len=24) :: str
  character(len=18) :: filnam18
  character(len=10) :: add10
  character(len=2) :: el = 'N_'
  logical :: metal
    !
  include 'sphsetnx.h.f90'
  !
  include 'sphsetx.h.f90'
  ! Currently newgetlym is only programmed for lmax = 5
!  maxll = min( maxll, 5 )
  maxll = 5
  call newgetprefs( yp, maxll, nsphpt, wsph, xsph, ysph, zsph )
  pi = 4.0d0 * atan( 1.0d0 )
  !
  open( unit=99, file='Pquadrature', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) npt
  allocate( x( npt ), w( npt ) )
  su = 0
  do i = 1, npt
     read ( 99, * ) x( i ), w( i )
     su = su + w( i )
  end do
  close( unit=99 )
  w = w * 2 / su
  
  open(unit=99,file='kmesh.ipt',status='old')
  read(99,*) kmesh(:)
  close(99)
  nk = product(kmesh)
!  write(6,*) nk

  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  read(99,*) avecs(:,:)
  close(99)
  call getomega( avecs, omega )

  inquire( file='efermiinrydberg.ipt', exist=metal )
  if( metal ) then
    open(unit=99,file='efermiinrydberg.ipt',form='formatted',status='old')
    read(99,*) efermi
    close(99)

    open(unit=99,file='nspin',form='formatted',status='old')
    read(99,*) nspin2
    close(99)


    open(unit=99,file='brange.ipt', form='formatted',status='old')
    read(99,*) brange(:)
    close(99)
    nband = brange(2) - brange(1) + 1
    allocate( enk(nband,nk*nspin2), enkskip(brange(3):brange(4) ) )
    open(unit=99,file='enkfile',form='formatted',status='old')
    do i = 1, nk*nspin2
      read(99,*) enk(:,i)
      read(99,*) enkskip(:)
    enddo
    close(99)
  endif
    



!  open(unit=99,file='ZNL',form='formatted',status='old')
!  read(99,*) ZZ, nc, lc
!  close(99)
!  write(add10, '(A1,I3.3,A1,I2.2,A1,I2.2)') 'z', ZZ, 'n', nc, 'l', lc

  ! Later maybe add brange, etc, for consistency check
  ! and if we have metals we'll need to figure out occ/unocc
  
!  write(str,'(A8,I3.3)') 'prjfilez', ZZ
!  open(unit=99,file=str,form='formatted',status='old')
!  read(99,*) lmin, lmax, dumi, dumr
!  allocate( nproj(lmin:lmax) )
!  do l=lmin, lmax
!    read(99,*) nproj(l)
!  enddo
!  close(99)


!  write(6,'(A2,A5,X,A14,A14,A14)') '##', 'N L site', 'Real (eV)', 'Imag (eV)', 'Den Trace'
  write(6,'(A2,A,X,A14,A14)') '##', "N L site m_c m'_c spin", 'Real (eV)', 'Imag (eV)'!, 'Den Trace'
  ! 
  !
!  open(unit=99,file='edgelist',form='formatted',status='old')
!  read(99,*) ZZ, nc, lc
!  close(99)

  open(unit=98,file='exx.inp',form='formatted',status='old')
!  open(unit=98,file='sitelist',form='formatted',status='old') 
  read(98,*) nsite
  do isite = 1,nsite
!    read(98,*) el, ZZ, iisite
    read(98,*) el, ZZ, nc, lc, iisite

    write(str,'(A8,I3.3)') 'prjfilez', ZZ
    write(add10, '(A1,I3.3,A1,I2.2,A1,I2.2)') 'z', ZZ, 'n', nc, 'l', lc
    
    open(unit=99,file=str,form='formatted',status='old')
    read(99,*) lmin, lmax, dumi, dumr
    allocate( nproj(lmin:lmax) )
    do l=lmin, lmax
      read(99,*) nproj(l)
    enddo
    close(99)


  write(str,'(A8,A2,I4.4)') 'parcksv.', el, iisite
  open(unit=99,file=str,form='unformatted',access='stream',status='old')
  read(99) nptot, ntot, nspin
!  write(6,*) nptot, ntot, nspin
  allocate(cks(nptot, ntot, nspin) )
  read(99) cks
  close(99)

  allocate( rho( nptot, nptot, nspin ) ) 
  rho(:,:,:) = 0.0_DP

  
  if( metal ) then
    if( nspin .ne. nspin2 ) then
      write(6,*) "parcks and nspin.ipt don't agree"
      stop
    endif
    if( ntot .ne. nk*nband ) then
      write(6,*) "parcks and bands*kpts don't agree"
      stop
    endif

    iks = 0
    do is = 1, nspin
      k = 0
      do ik = 1, nk
        iks = iks + 1
        do ib = 1, nband
          k = k + 1
          if( enk( ib, iks ) .le. efermi ) then
            do j = 1, nptot
              do i = 1, nptot
                rho(i,j,is) = rho(i,j,is) + cks(i,k,is) * conjg(cks(j,k,is))
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  else
    do is = 1, nspin
      do k = 1, ntot
        do j = 1, nptot
          do i = 1, nptot
            rho(i,j,is) = rho(i,j,is) + cks(i,k,is) * conjg(cks(j,k,is))
          enddo
        enddo
      enddo
    enddo
  endif

! Hoist IO
!  i = 0
  kgl = min(lc,lmin)
  kgh = lmax+lc
  nprojmax = 0
  do l = lmin, lmax
    nprojmax = max( nprojmax, nproj(l) )
  enddo
  ! This is way too big, but also still very small
  allocate( gkk( nprojmax, nprojmax, kgl:kgh, lmin:lmax, lmin:lmax ), &
            fkk( nprojmax, nprojmax, kgl:kgh, lmin:lmax, lmin:lmax ) )
  gkk = 0.0_DP
  do l1 = lmin, lmax
    do l2 = lmin, l1
      do kk = 0, min(lc+l1,lc+l2)
        if( ( abs(lc-l1) .gt. kk ) .or. ( abs(lc-l2) .gt. kk )  & !.or. (lc+l1 .lt. kk ) &
            .or. ( mod(lc+l1+kk,2) .ne. 0 ) .or. ( mod(lc+l2+kk,2) .ne. 0 ) ) then
!           .or. (lc+l2 .lt. kk ) .or. ( mod(lc+l1+kk,2) .ne. 0 ) .or. ( mod(lc+l2+kk,2) .ne. 0 ) ) then
          cycle
        endif
        if( l1 .eq. l2 ) then
          write ( filnam18, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'gk', lc, l1, kk, 'z', zz, 'n', nc, 'l', lc
        else
          write ( filnam18, '(1a2,4i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'gk', lc, l1, l2, kk, 'z', zz, 'n', nc, 'l', lc
        endif
        open( unit=99, file=filnam18, form='formatted', status='old' )
        rewind( 99 )
        do i2 = 1, nproj(l2)
          read ( 99, * ) gkk( 1 : nproj(l1), i2, kk, l1, l2 )
!          gkk( 1 : nproj(1), i2, kk, l2, l1 ) = gkk( 1 : nproj(1), i2, kk, l1, l2 )
        end do
        gkk( 1:nproj(l2), 1:nproj(l1), kk, l2, l1 ) = transpose( gkk( 1 : nproj(l1), 1:nproj(l2), kk, l1, l2 ) )
        close( 99 )
      enddo
      do kk = 2, min( 2*lc, l1+l2), 2
        if( ( abs(l1-l2) .gt. kk ) .or. mod(l1+l2,2) .ne. 0 ) then
          cycle
        endif
        write ( filnam18, '(1a2,4i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'fk', lc, l1, l2, kk, 'z', zz, 'n', nc, 'l', lc
!        write(6,*) filnam18
        open( unit=99, file=filnam18, form='formatted', status='old' )
        rewind( 99 )
        do i2 = 1, nproj(l2)
          read ( 99, * ) fkk( 1 : nproj(1), i2, kk, l1, l2 )
!          gkk( 1 : nproj(1), i2, kk, l2, l1 ) = gkk( 1 : nproj(1), i2, kk, l1, l2 )
        end do
        fkk( 1:nproj(l2), 1:nproj(l1), kk, l2, l1 ) = transpose( fkk( 1 : nproj(l1), 1:nproj(l2), kk, l1, l2 ) )
        close( 99 )
      enddo

    enddo
  enddo

    

    
  allocate( VXX( -lc:lc, -lc:lc, nspin ), VF2( -lc:lc, -lc:lc) )
  VXX = 0.0_DP
  VF2 = 0.0_DP
  l1 = lc
  l4 = lc
  do is = 1, nspin
    do m1 = -lc, lc
      do m4 = -lc, lc

        j = 0
        do l2 = lmin, lmax
          do m2 = -l2, l2
            do nu2 = 1, nproj(l2)
              j = j + 1

              i = 0
              do l3 = lmin, lmax
                do m3 = -l3, l3
                  do nu3 = 1, nproj(l3) ! todo nu3
                    i = i + 1

!                   ! This is to match historical limitation in OCEAN
!                   ! actual differences are very small
!                    if( l2 .ne. l3 ) cycle
                
!                    if ( m1 + m2 .ne. m3 + m4 ) cycle


                    if( m1 + m2 .eq. m3 + m4 ) then
                    mk = m1 - m3
                    coef = 0.0_DP
  
                    do kk = 0, min(lc+l3,lc+l2)
                      if( ( abs(lc-l3) .gt. kk ) .or. ( abs(lc-l2) .gt. kk )  &
                          .or. ( mod(lc+l3+kk,2) .ne. 0 ) .or. ( mod(lc+l2+kk,2) .ne. 0 ) ) then
                        cycle
                      endif
                      if ( abs( mk ) .le. kk ) then
                        call threey( l1, m1, kk, mk, l3, m3, no, npt, x, w, yp, f1 )
                        call threey( l2, m2, kk, mk, l4, m4, yes, npt, x, w, yp, f2 )
                        coef = coef + gkk( nu2, nu3, kk, l2, l3 ) * f1 * f2 * ( 4.0_DP * pi / real( 2 * kk + 1, DP ) )
                      endif
                    enddo

                    VXX( m4, m1, is ) = VXX( m4, m1, is ) - coef * rho( i, j, is )
                    endif

                    if( mod( l2+l3,2 ) .eq. 0 ) then
                      coef = 0.0_DP
                      mkd = m4 - m1
                      kkmin = max( 2, abs(l2-l3), abs(mkd) )
                      if( mod( kkmin, 2 ) .ne. 0 ) kkmin = kkmin + 1
                      do kk = kkmin, min( 2*lc, l2+l3), 2
                        call threey( l1, m1, kk, mkd, l4, m4, yes, npt, x, w, yp, f1 )
                        call threey( l2, m2, kk, mkd, l3, m3, yes, npt, x, w, yp, f2 )
                        coef = coef + fkk( nu2, nu3, kk, l2, l3 ) * f1 * f2 * ( 4.0_DP * pi / real( 2 * kk + 1, DP ) )
                      enddo
                      VF2( m4, m1 ) = VF2( m4, m1 ) + coef * rho( i, j, is )
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        write(6,'(A2,1X,I1.1,1X,I1.1,1X,I4.4,1X,I2,1x,I2,1x,I1,1x,F14.6,E14.6)') el, nc, lc, iisite, m4, m1, is, &
                       real(VXX(m4,m1,is),DP)/real(nk,DP)/omega, &
                      aimag(VXX(m4,m1,is))/real(nk,DP)/omega

      enddo
    enddo
  enddo

  
  
  do m1 = -lc, lc
    do m4 = -lc, lc
      write(6,'(A2,1X,I1.1,1X,I1.1,1X,I4.4,1X,I2,1x,I2,1x,I1,1x,F14.6,E14.6)') el, nc, lc, iisite, m4, m1, is, &
                     real(VF2(m4,m1),DP)/real(nk,DP)/omega, &
                    aimag(VF2(m4,m1))/real(nk,DP)/omega
    enddo
  enddo
#if 0
  i = 2*lc + 1
  allocate( work( 2*i), rwork(3*i), WW(i) )
  do is = 1, nspin
    call ZHEEV( 'N', 'U', i, VXX(:,:,is), i, WW, work, 2*i, rwork, info )
    write(6,*) WW(:)/real(nk,DP)/omega
  enddo
  deallocate( work, rwork, WW )
#endif
  

  ! Calculate average VXX, subtract it from the diagonal 
  ! This separates the shift of the edge from the spin/m_l splittings
  su = 0.0_DP
  do is = 1, nspin
    do m1 = -lc, lc
      su = su + VXX(m1,m1,is)
    enddo
  enddo
  su = su / real(nspin*(2*lc+1),DP)
  do is = 1, nspin
    do m1 = -lc, lc
      VXX(m1,m1,is) = VXX(m1,m1,is) - su
    enddo
  enddo
  su = su/(real(nk,DP)*omega)
  write(6,'(A,X,A2,X,I3,X,I1,X,I1.1,X,I4.4,F24.16)') 'AVG', el, ZZ, nc, lc, iisite, su

  write( filnam18 , '(A5,A2,A1,I2.2,A1,I2.2,A1,I4.4)') 'vexx.', el, 'n', nc, 'l', lc, '.', iisite
  open(unit=99,file=filnam18,status='unknown',form='formatted')
  do is = 1, nspin
    do i = -lc, lc
      write(99,*) VXX(:,i,is)/real(nk,DP)/omega
    enddo
  enddo
  close(99)

  write( filnam18 , '(A4,A2,A1,I2.2,A1,I2.2,A1,I4.4)') 'vf2.', el, 'n', nc, 'l', lc, '.', iisite
  open(unit=99,file=filnam18,status='unknown',form='formatted')
  do i = -lc, lc
    write(99,*) VF2(:,i)/real(nk,DP)/omega
  enddo
  close(99)


  deallocate( VXX, rho, gkk, fkk, VF2)

  deallocate( cks, nproj )
  enddo


end program coreExchange
