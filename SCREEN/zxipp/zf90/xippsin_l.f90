program xipostproc
  implicit none
  !
  integer :: npt
  real( kind = kind( 1.0d0 ) ) :: rmax, s1, s2
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vipt( : ), drel( : )
  !
  integer :: nbasis, i, j, k, ii, jj, ibase, jbase, nang, nrad
  real( kind = kind( 1.0d0 ) ) :: q, arg, pi, su, su0, coul, cterm
  real( kind = kind( 1.0d0 ) ), allocatable :: basfcn( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: potfcn( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: qtab( : ), ptab( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: coulmat( :, : ), ximat( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: res( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xibb( :, : ), xifull( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: rhs( : ), pref( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: tmp( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xir( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: nind( : ), nind0( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vind( : ), vind0( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ylm( :, : ), specpnt( :, : ), zeros(:), specwgt( : )
  real( kind = kind( 1.0d0 ) ) :: ang, dumf, output( 10 )
  real(kind=16) :: z, max_z, min_z, eps, diff, real_z
  !
  integer,allocatable :: basis_map( :, : ), nbasisbyl( : )
  integer :: lmax, l, basis_start

  eps = 0.000000000000000001
  write(6,*) eps
  pi = 4.0d0 * atan( 1.0d0 )
  open( unit=99, file='specpnt', form='formatted', status='old' )
  rewind 99
  read( 99, * ) nang
  allocate( specpnt( 3, nang ), specwgt( nang) )
  do i = 1, nang
    read( 99, * ) specpnt( :, i ), specwgt( i )
  enddo
  close( 99 )
  
  !
  open( unit=99, file='rbfile.bin', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) npt, rmax
  allocate( posn( 3, npt ), wpt( npt ), vipt( npt ), drel( npt ) )
  read ( 99 ) posn
  read ( 99 ) wpt
  read ( 99 ) drel
  close ( unit=99 )
  nrad = npt / nang
  call mkvipt( npt, drel, vipt )
  write(6,*) nang, nrad, npt
  !
  lmax = 2
  allocate( basis_map( 2*(lmax+1), 0 : lmax ) )
! For now try something like this
! Reduce the number of basis functions for each l by the m degeneracy 
  allocate( nbasisbyl( 0 : lmax ) )
  nbasisbyl( 0 ) = 30 
  nbasisbyl( 1 ) = 30
  nbasisbyl( 2 ) = 20

  allocate(zeros(nbasisbyl(1)))
    min_z = pi/2.0d0  !pi / 4.0d0
    real_z = pi 
    max_z = 3.0* pi / 2.0d0 
    z = real_z !- dble( i - 1) * pi
  do i = 1, nbasisbyl(1)
    do
      diff = tan(z) - real_z
      if( abs(diff) .lt. eps ) goto 10
      if( diff .lt. 0.0d0 ) then
        min_z = real_z - eps
        real_z = (max_z + min_z) / 2.0d0
      else
        max_z = real_z + eps
        real_z = (max_z + min_z) / 2.0d0
      endif
      z = real_z - dble( i - 1) * pi
    enddo
 10 continue
    zeros( i ) = dble(real_z/rmax)
    min_z =  dble(i+1) * pi 
    max_z = min_z + pi/2.0d0 
    write(6,*) i, dble(real_z), dble(tan(real_z - real( i - 1) * pi)), i*pi
    real_z = (max_z + min_z) / 2.0d0
    z = real_z - dble( i ) * pi
!    write(6,*) min_z,real_z,max_z
!    write(6,*) z, tan(z)
  enddo
11 continue

  nbasis = 0
  do i = 0, lmax
    do j = 1, 2*i+1
      basis_map( j, i ) = nbasis + 1
      nbasis = nbasis + nbasisbyl( i )
    enddo
  enddo
  
!  read ( 5, * ) nbasis
  
  allocate( basfcn( npt, 2*nbasis ), potfcn( npt, 2*nbasis ) )
  allocate( qtab( 2*nbasis ), ptab( 2*nbasis ), pref( 2*nbasis ) )
  allocate( coulmat( 2*nbasis, 2*nbasis ) )
  !
  ! This version of the code uses q = n * pi / rmax, with n = 1, 2, 3, ... 
  ! which is DIFFERENT FROM the q rmax = tan( q rmax ) model.
  !
  ! L = 0
  do i = basis_map( 1, 0 ), basis_map( 1, 0 ) + nbasisbyl( 0 ) !nbasis
     q = pi * dble( i ) / rmax
     pref( i ) = 2.0d0 * pi * rmax / q ** 2 ! this is true for sin( 2 q rmax ) = 0
     pref( i ) = 1.0d0 / sqrt( pref( i ) )
     qtab( i ) = q
     ptab( i ) = pref( i )
     coul = 4.0d0 * pi / q ** 2
     arg = q * rmax
     cterm = cos( arg )
     do j = 1, npt
        arg = q * drel( j )
        basfcn( j, i ) = pref( i ) * sin( arg ) / arg
        potfcn( j, i ) = pref( i ) * coul * ( sin( arg ) / arg - cterm )
     end do
  end do
  ! l = 1
  allocate( ylm( 3, npt ) )
  do i = 1, nang
    ylm( 1, i ) =  sqrt( 3.0d0 ) * specpnt( 1, i )
    ylm( 2, i ) =  sqrt( 3.0d0 ) * specpnt( 2, i )
    ylm( 3, i ) =  sqrt( 3.0d0 ) * specpnt( 3, i )
  enddo
  do i = 1, nrad - 1
    do j = 1, nang
      ylm( :, j + i * nang ) = ylm( :, j )
    enddo
  enddo

  do l = 1,3
    basis_start = basis_map( l, 1 )
    do i = basis_start, basis_start + nbasisbyl( 1 )
       q = pi * dble( i - basis_start + 1 ) / rmax
       q = zeros( i - basis_start + 1 ) 
       pref( i ) = 2.0d0 * pi * rmax / q ** 2 ! this is true for sin( 2 q rmax ) = 0
       pref( i ) = pi * ( 2.0d0 * q * rmax - sin( 2.0d0 * q * rmax ) ) / (q**3)
       pref( i ) = 1.0d0 / sqrt( pref( i ) )
      
       qtab( i ) = q
       ptab( i ) = pref( i )
       coul = 5.0d0 * 4.0d0 * pi / q ** 2
       arg = q * rmax
       do j = 1, npt
          arg = q * drel( j )
!          arg = q * sqrt(posn(1,j)*posn(1,j) + posn(2,j)*posn(2,j)+posn(3,j)*posn(3,j) )
          ! Use real Ylm: for l=1, ie. x/r, y/r, z/r
!          ang = sqrt( 3.0d0 / (4.0d0 * pi ) ) * posn( l, j ) / drel( j )
!          ang = sqrt( 3.0d0 ) * posn( l, j ) / drel( j )
!          ang = sqrt( 3.0d0 ) * posn( l, j ) / sqrt(posn(1,j)*posn(1,j) + posn(2,j)*posn(2,j)+posn(3,j)*posn(3,j) )
          basfcn( j, i ) = pref( i ) * ( sin( arg ) / (arg*arg) - cos(arg)/arg ) * ylm( l, j )
          potfcn( j, i ) = coul * basfcn( j, i )
       end do
    end do
  enddo
  ! l = 2
  deallocate( ylm )
  allocate( ylm( 5, npt ) )
  ylm = 0
!  ylm( :, : ) = sqrt( 1.0d0 / (4.0d0 * pi ) ) 
  do i = 1, nang
    ylm( 1, i ) = sqrt( 5.0d0 ) / 2.0d0 &
                * ( 2.0d0*specpnt(3,i)*specpnt(3,i) - specpnt(1,i)*specpnt(1,i) - specpnt(2,i)*specpnt(2,i) )
    ylm( 2, i ) = sqrt( 15.0d0 ) * specpnt(2,i)*specpnt(3,i)
    ylm( 3, i ) = sqrt( 15.0d0 ) * specpnt(1,i)*specpnt(3,i)
    ylm( 4, i ) = sqrt( 15.0d0 ) * specpnt(2,i)*specpnt(1,i)
    ylm( 5, i ) = sqrt( 15.0d0 ) / 2.0d0 * (specpnt(1,i)*specpnt(1,i) - specpnt(2,i)*specpnt(2,i) )
  enddo
  do i = 1, nrad - 1
    do j = 1, nang
      ylm( :, j + i * nang ) = ylm( :, j )
    enddo
  enddo

!  do i = 1, npt
!    ylm( 1, i ) = sqrt( 5.0d0 / pi ) / 4.0d0 & 
!                * ( 2.0d0 * posn(3,i) * posn(3,i) - posn(1,i)*posn(1,i) - posn(2,i)*posn(2,i) ) &
!                / ( drel( i ) * drel( i ) )
!    ylm( 2, i ) = sqrt( 15.0d0 / pi ) / 2.0d0 &
!                * posn(2,i)*posn(3,i)/(drel(i)*drel(i))
!    ylm( 3, i ) = sqrt( 15.0d0 / pi ) / 2.0d0 &
!                * posn(1,i)*posn(3,i)/(drel(i)*drel(i))
!    ylm( 4, i ) = sqrt( 15.0d0 / pi ) / 2.0d0 &
!                * posn(2,i)*posn(1,i)/(drel(i)*drel(i))
!    ylm( 5, i ) = sqrt( 15.0d0 / pi ) / 4.0d0 &
!                * ( posn(1,i)*posn(1,i) - posn(2,i)*posn(2,i) ) &
!                / ( drel( i ) * drel( i ) )
!  enddo
!  ylm( :, : ) = ylm( :, : ) * 2.0d0 * sqrt( pi )

  do l = 1, 5
    basis_start = basis_map( l, 2 )
    do i = basis_start, basis_start + nbasisbyl( 2 )
       q = pi * dble( i - basis_start + 1 ) / rmax
       q = zeros( i - basis_start + 1 )
       pref( i ) = 2.0d0 * pi * rmax / q ** 2 ! this is true for sin( 2 q rmax ) = 0
       pref( i ) = pi * ( 2.0d0 * q * rmax - sin( 2.0d0 * q * rmax ) ) / (q**3)
       pref( i ) = 1.0d0 / sqrt( pref( i ) )
       qtab( i ) = q
       ptab( i ) = pref( i )
       coul = 5.0d0 * 4.0d0 * pi / q ** 2
       arg = q * rmax
       do j = 1, npt
          arg = q * drel( j )
!          basfcn( j, i ) = pref( i ) * ( ( 3.0d0 / (arg*arg*arg) - 1.0d0/arg ) * sin(arg) & 
!                                        - (3.0d0/(arg*arg))*cos(arg) ) &
          basfcn( j, i ) = pref( i ) * ( ( 3.0d0 / (arg*arg) - 1.0d0 ) * sin(arg)/arg - 3.0d0 * cos(arg)/(arg*arg) ) &
                                     * ylm( l, j)
!          basfcn( j, i ) = pref( i ) * ( sin( arg ) / (arg*arg) - cos(arg)/arg ) * ylm( l, j )
          potfcn( j, i ) = coul * basfcn( j, i )
       end do
    end do
  enddo
!  20 continue


  coulmat = 0
  do i = 1, nbasisbyl( 0 )
     do j = 1, nbasisbyl( 0 )
         coulmat( i, j ) = 8.0d0 * pi * cos( qtab( i ) * rmax ) * cos( qtab( j ) * rmax ) / ( qtab( i ) * qtab( j ) )
     end do
     coulmat( i, i ) = coulmat( i, i ) + 4.0d0 * pi / ( qtab( i ) ** 2 )
  end do
  do i = basis_map( 1, 1 ), nbasis
    coulmat( i, i ) = 4.0d0 * pi / ( qtab( i ) ** 2 )
  enddo
  open ( unit=99, file='chkmat', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbasis
     do j = 1, nbasis
        s1 = 0
        s2 = 0
        do ii = 1, npt
           s1 = s1 + basfcn( ii, i ) * basfcn( ii, j ) * wpt( ii )
           s2 = s2 + basfcn( ii, i ) * potfcn( ii, j ) * wpt( ii )
        end do
        su = 0
        if ( i .eq. j ) su = 1
        write ( 99, '(2i5,4(1x,1e15.8))' ) i, j, s1, s2, su, coulmat( i, j )       
     end do
  end do
  close( unit=99 )
!  goto 111
  !
  allocate( ximat( npt, npt ), res( npt ) )
  open( unit=99, file='ximat', form='unformatted', status='unknown' )
  rewind 99
  do i = 1, npt
     read( 99 ) ximat( :, i )
  end do
  close( unit=99 )
  !
  open( unit=99, file='specpnt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nang
  close( unit=99 )
  nrad = npt / nang
  allocate( xir( nrad, nrad ) )
  xir = 0
  ibase = 1
  do i = 1, nrad
     jbase = 1
     do j = 1, nrad
        su = 0
        do ii = ibase, ibase + nang - 1
           do jj = jbase, jbase + nang - 1
              su = su + ximat( ii, jj )
           end do
        end do
        xir( i, j ) = su
        jbase = jbase + nang
     end do
     ibase =  ibase + nang
  end do
  open( unit=99, file='xirr', form='formatted', status='unknown' )
  rewind 99
  ibase = 1
  do i = 1, nrad
     jbase = 1
     do j = 1, nrad
        write ( 99, '(3(1x,1e15.8))' ) drel( ibase ), drel( jbase ), xir( i, j )
        jbase = jbase + nang
     end do
     ibase =  ibase + nang
     write ( 99, * )
  end do
  close( unit=99 )
  !
  allocate( xibb( nbasis, nbasis ), rhs( nbasis ) )
  rhs( : ) = 0
  open( unit=99, file='rhs', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbasis
     su = 0
     do ii = 1, npt
        su = su + basfcn( ii, i ) * vipt( ii ) * wpt( ii )
     end do
     rhs( i ) = su
     write ( 99, '(1i5,2(1x,1e15.8))' ) i, qtab( i ), rhs( i )
  end do
  close( unit=99 )
  !
  open( unit=99, file='xibb', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbasis
     do jj = 1, npt
        su = 0
        do ii = 1, npt
           su = su + ximat( ii, jj ) * basfcn( ii, i ) * wpt( ii )
        end do
        res( jj ) = su
     end do
     do j = 1, nbasis
        su = 0 
        do jj = 1, npt
           su = su + res( jj ) * basfcn( jj, j ) * wpt( jj )
        end do
        xibb( i, j ) = su
        write ( 99, '(2i5,3(1x,1e15.8))' ) i, j, qtab( i ), qtab( j ), xibb( i, j )
     end do
     write ( 99, * )
  end do
  close( unit=99 )
  !
  allocate( tmp( nbasis, nbasis ), xifull( nbasis, nbasis ) )
  tmp = 0
  do i = 1, nbasis
     do j = 1, nbasis
        do k = 1, nbasis
           tmp( i, j ) = tmp( i, j ) - xibb( i, k ) * coulmat( k, j )
        end do
     end do
     tmp( i, i ) = tmp( i, i ) + 1.0d0
  end do
  call rinvert( nbasis, tmp )
  xifull = 0
  do i = 1, nbasis
     do j = 1, nbasis
        do k = 1, nbasis
           xifull( i, j ) = xifull( i, j ) + tmp( i, k ) * xibb( k, j )
        end do
     end do
  end do
  allocate( nind( npt ), nind0( npt ) )
  allocate( vind( npt ), vind0( npt ) )
  nind = 0
  nind0 = 0
  vind = 0
  vind0 = 0
  open( unit=98, file='basisc', form='formatted', status='unknown' )
  open( unit=99, file='xifull', form='formatted', status='unknown' )
  rewind 98
  rewind 99
  do i = 1, nbasis
     su = 0
     su0 = 0
     do j = 1, nbasis
        su = su + xifull( i, j ) * rhs( j )
        su0 = su0 + xibb( i, j ) * rhs( j )
        write ( 99, '(2i5,3(1x,1e15.8))' ) i, j, qtab( i ), qtab( j ), xifull( i, j )
     end do
     write ( 98, '(1i5,2(1x,1e15.8))' ) i, su, su0
     do ii = 1, npt
        nind0( ii ) = nind0( ii ) + su0 * basfcn( ii, i )
        vind0( ii ) = vind0( ii ) + su0 * potfcn( ii, i )
        nind( ii ) = nind( ii ) + su * basfcn( ii, i )
        vind( ii ) = vind( ii ) + su * potfcn( ii, i )
     end do
     write ( 99, * )
  end do
  close( unit=98 )
  close( unit=99 )
  open( unit=99, file='ninduced', form='formatted', status='unknown' )
  rewind 99
  do ii = 1, npt
     write ( 99, '(5(1x,1f10.5))' ) drel( ii ), nind( ii ), nind0( ii ), vind( ii ), vind0( ii )
  end do
  close( unit=99 )
  open( unit=99, file='vind', form='formatted', status='unknown' )
  rewind 99
!  vind( : ) = vind( : ) * wpt( : )
  specwgt( : ) = specwgt( : ) / (4.0d0 * pi )
  do ii = 1, nrad
    output( : ) = 0.0d0
    do i = 1, nang
      output( 1 ) = output(1) + vind( i + (ii-1)*nang ) * specwgt( i ) 
      output( 2 ) = output(2) + vind( i + (ii-1)*nang ) * sqrt(3.0d0)* specpnt( 1, i ) * specwgt( i )
      output( 3 ) = output( 3 ) + vind( i + (ii-1)*nang ) * sqrt( 3.0d0 ) * specpnt( 2, i )* specwgt( i )
      output( 4 ) = output( 4 ) + vind( i + (ii-1)*nang ) * sqrt( 3.0d0 ) * specpnt( 3, i )* specwgt( i )
      output( 5 ) = output( 5 ) + vind( i + (ii-1)*nang ) * ylm( 1, i )* specwgt( i )
      output( 6 ) = output( 6 ) + vind( i + (ii-1)*nang ) * ylm( 2, i )* specwgt( i )
      output( 7 ) = output( 7 ) + vind( i + (ii-1)*nang ) * ylm( 3, i )* specwgt( i )
      output( 8 ) = output( 8 ) + vind( i + (ii-1)*nang ) * ylm( 4, i )* specwgt( i )
      output( 9 ) = output( 9 ) + vind( i + (ii-1)*nang ) * ylm( 5, i )* specwgt( i )
    enddo
     write ( 99, '(10(1x,1f10.5))' ) drel( 1+(ii-1)*nang ), output( 1 : 9 )
  enddo
  close( 99 )
      
  !
  deallocate( posn, wpt, vipt, drel, qtab, ptab, coulmat, nind )
  !
111 continue
end program xipostproc
