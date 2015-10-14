program xipostproc
  use ISO_FORTRAN_ENV
  implicit none
  !
  integer :: npt
  real( kind = kind( 1.0d0 ) ) :: s1, s2, rmax_old
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vipt( : ), drel( : )
  !
  integer :: nbasis, i, j, k, ii, jj, ibase, jbase, nang, nrad, tot_nbasis
  integer :: lmax, mtot, im, m, l, offset, offsetj, lj, mj, imj
  real( kind = kind( 1.0d0 ) ) :: su, su0, coul, cterm,  tmp_pref1
  real( kind = kind( 1.0d0 ) ), allocatable :: basfcn( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: potfcn( :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: qtab( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: coulmat( :, : ), ximat( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: res( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xibb( :, : ), xifull( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: rhs( : )
  real( kind=16  ), allocatable :: pref( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: tmp( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xir( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: nind( : ), nind0( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vind( : ), vind0( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: Ylm( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: specpnt( :, : )

  real( kind=16 ) :: q, pi, arg, root2piR, rmax
  !
  integer, parameter :: DP = 16
  include 'tanqr.f90'
  !
!  pi = 4.0d0 * atan( 1.0d0 )
  pi = real( 3.1415926535897932384626433832795029, 16 )
  write(6,*) DP, pi, tanqr( 1 )
  !
  open( unit=99, file='rbfile.bin', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) npt, rmax_old
  rmax = rmax_old
  allocate( posn( 3, npt ), wpt( npt ), vipt( npt ), drel( npt ) )
  read ( 99 ) posn
  read ( 99 ) wpt
  read ( 99 ) drel
  close ( unit=99 )
  call mkvipt( npt, drel, vipt )
  !
  read ( 5, * ) nbasis, lmax
  write(6,*)
  write(6,*) nbasis, lmax
!  lmax = 0
  ! mtot contains the total number of m's so we can use sequential, not padded arrays
  mtot = 1
  do i = 1, lmax
    mtot = mtot + ( 2 * i + 1 )
  enddo

  ! the tot_nbasis is the size of the coulmat, xibb
  tot_nbasis = nbasis * mtot

  !
  allocate( basfcn( npt, nbasis, 0 : lmax ), potfcn( npt, nbasis, 0 : lmax ) )
  allocate( qtab( nbasis ), pref( nbasis ) )
  allocate( coulmat( tot_nbasis, tot_nbasis ) )
  basfcn = 0.0_DP
  !
!!!  ! This version of the code uses q = n * pi / rmax, with n = 1, 2, 3, ... 
!!!  ! which is DIFFERENT FROM the q rmax = tan( q rmax ) model.
  ! Changing to q rmax = tan( q rmax ) model!
  !
  root2piR = sqrt( 2.0_DP * pi * rmax )
  do i = 1, nbasis
!    q = tanqr( i ) / rmax 
!    pref( i ) = q / ( root2piR * SIN( tanqr( i ) ) )
!    qtab( i ) = q
!    coul = 4.0_DP * pi / q ** 2
!    cterm = cos( tanqr( i ) )

     q = pi * dble( i ) / rmax
     pref( i ) = 2.0d0 * pi * rmax / q ** 2 ! this is true for sin( 2 q rmax ) = 0
     pref( i ) = 1.0d0 / sqrt( pref( i ) )
     qtab( i ) = q
     coul = 4.0d0 * pi / q ** 2
     arg = q * rmax
     cterm = cos( arg )


    ! l = 0
    do j = 1, npt
      arg = q * drel( j )
!      if( abs( arg ) .lt. 0.01_DP )  then 
!        basfcn( j, i, 0 ) = pref( i ) * ( 1.0_DP - ( arg**2 )/ 6.0_DP + (arg ** 4 ) / 120.0_DP )
!        potfcn( j, i, 0 ) = pref( i ) * coul * ( sin( arg ) / arg - cterm )
!      else
        basfcn( j, i, 0 ) = pref( i ) * sin( arg ) / arg
        potfcn( j, i, 0 ) = pref( i ) * coul * ( sin( arg ) / arg - cterm )
!      endif
    end do
  enddo    

! coulmat is block diagonal
  coulmat = 0.0d0
  do i = 1, nbasis
     do j = 1, nbasis
         coulmat( i, j ) = 8.0d0 * pi * cos( qtab( i ) * rmax ) * cos( qtab( j ) * rmax ) / ( qtab( i ) * qtab( j ) )
     end do
     coulmat( i, i ) = coulmat( i, i ) + 4.0d0 * pi / ( qtab( i ) ** 2 )
  end do
  

  if ( lmax .ge. 1 ) then
  if( .false. ) then
    do i = 1, nbasis
      q = tanqr( i ) / rmax 
      pref( i ) = q / ( root2piR * SIN( tanqr( i ) ) )
      qtab( i ) = q
      coul = 4.0_DP * pi / q ** 2
      cterm = sin( tanqr( i ) ) / ( 3.0_DP * rmax )
      do j = 1, npt
        arg = q * drel( j )
        basfcn( j, i, 1 ) = pref( i ) * ( sin( arg ) / arg ** 2 - cos( arg ) / arg )
        potfcn( j, i, 1 ) = pref( i ) * coul * ( - sin( arg ) / arg ** 2 + cos( arg ) / arg - drel( j ) * cterm )
      enddo
    end do


    offset = 0
    do m = -1, 1
      offset = offset + nbasis
      offsetj = offset
      do i = 1, nbasis
        do j = 1, nbasis
          ! I was missing the 3 here and below?
          coulmat( i + offset, j + offsetj ) = 8.0d0 * pi / ( qtab( i ) * qtab( j ) * 3.0d0)
        enddo
        coulmat( i + offset, i + offset ) = 20.0d0 * pi / (  qtab( i ) ** 2 * 3.0d0 )
      enddo
    enddo
  else
    do i = 1, nbasis
      q = pi * dble( i ) / rmax
      pref( i ) = q / ( root2piR ) 
      qtab( i ) = q
      coul = 4.0_DP * pi / q ** 2
      do j = 1, npt
        arg = q * drel( j )
        basfcn( j, i, 1 ) = pref( i ) * ( sin( arg ) / arg ** 2 - cos( arg ) / arg )
        potfcn( j, i, 1 ) = pref( i ) * coul * ( sin( arg ) / arg ** 2 + cos( arg ) / arg )
      enddo
    end do


    offset = 0
    do m = -1, 1
      offset = offset + nbasis
      offsetj = offset
      do i = 1, nbasis
        do j = 1, nbasis
          ! I was missing the 3 here and below?
!          coulmat( i + offset, j + offsetj ) = 8.0d0 * pi / ( qtab( i ) * qtab( j ) * 3.0d0)
        enddo
        coulmat( i + offset, i + offset ) = 4.0d0 * pi / (  qtab( i ) ** 2 )
      enddo
    enddo

  endif
  endif


  open( unit=99, file='specpnt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nang
  allocate( specpnt( 4, nang ) )
  read( 99, * ) specpnt( :, : )
  close( unit=99 )
  ! Creat Ylm array
  ! one of the mtot is for l=0,m=0
  allocate( Ylm( npt, mtot ) )
  ! I want real Ylm instead of complex
  !
  ! The 4pi normalization is already in the bessel
  tmp_pref1 = sqrt( 3.0_DP )
!    tmp_pref1 = sqrt( 0.75_DP / pi )
!    tmp_pref2 = 0.5_DP * sqrt( 7.5_DP / pi )
!    do i = 1, nang
    Ylm( :, 1 ) = 1.0_DP
    if( lmax .ge. 1 ) then
      do i = 1, npt, nang
        j = i + nang - 1
        Ylm( i:j, 2 ) = tmp_pref1 * specpnt( 1, : )
        Ylm( i:j, 3 ) = tmp_pref1 * specpnt( 2, : )
        Ylm( i:j, 4 ) = tmp_pref1 * specpnt( 3, : )
      enddo
    endif
!    enddo


  open ( unit=99, file='chkmat', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbasis
     do j = 1, nbasis
        s1 = 0.0_DP
        s2 = 0.0_DP
        do ii = 1, npt
           s1 = s1 + basfcn( ii, i, 0 ) * basfcn( ii, j, 0 ) * wpt( ii )
           s2 = s2 + basfcn( ii, i, 0 ) * potfcn( ii, j, 0 ) * wpt( ii )
        end do
        su = 0
        if ( i .eq. j ) su = 1
        write ( 99, '(2i5,4(1x,1e15.8))' ) i, j, s1, s2, su, coulmat( i, j )       
     end do
  end do

  if( lmax .ge. 1 ) then
    offset = nbasis
    do i = 1, nbasis
       do j = 1, nbasis
          s1 = 0
          s2 = 0
          do ii = 1, npt
             s1 = s1 + basfcn( ii, i, 1 ) * basfcn( ii, j, 1 ) * Ylm(ii,2)**2 * wpt( ii )
             s2 = s2 + basfcn( ii, i, 1 ) * potfcn( ii, j, 1 ) * Ylm(ii,2)**2 * wpt( ii )
          end do
          su = 0
          if ( i .eq. j ) su = 1
          write ( 99, '(2i5,4(1x,1e15.8))' ) i, j, s1, s2, su, coulmat( i + offset, j + offset)
       end do
    end do
  endif
  close( unit=99 )
  !

  !
  allocate( ximat( npt, npt ), res( npt ) )
  open( unit=99, file='ximat', form='unformatted', status='unknown' )
  rewind 99
  do i = 1, npt
     read( 99 ) ximat( :, i )
  end do
  close( unit=99 )
  !
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




  !
  allocate( xibb( tot_nbasis, tot_nbasis ), rhs( tot_nbasis ) )
  rhs( : ) = 0
  open( unit=99, file='rhs', form='formatted', status='unknown' )
  rewind 99
  im = 0
  offset = 0
  do l = 0, lmax
    do m = -l, l
      im = im + 1
      do i = 1, nbasis
        su = 0
        do ii = 1, npt
          su = su + basfcn( ii, i, l ) * vipt( ii ) * wpt( ii ) * Ylm( ii, im )
        end do
        rhs( i + offset ) = rhs( i + offset ) + su
!        write ( 99, '(1i5,2(1x,1e15.8))' ) i, qtab( i ), rhs( i )
      end do
      offset = offset + nbasis
    enddo
  enddo
  offset = 0
  do im = 1, mtot
    do i = 1, nbasis
      write ( 99, '(1i5,2(1x,1e15.8))' ) i, qtab( i ), rhs( i + offset)
    enddo
    offset = offset + nbasis
  enddo
  close( unit=99 )
  !


  !
  open( unit=99, file='xibb', form='formatted', status='unknown' )
  rewind 99
  im = 0
  offset = 0
  do l = 0, lmax
    do m = -l, l
      im = im + 1
      do i = 1, nbasis
        do jj = 1, npt
          su = 0
          do ii = 1, npt
            su = su + ximat( ii, jj ) * basfcn( ii, i, l ) * wpt( ii ) * Ylm( ii, im )
          enddo
          res( jj ) = su
        end do
        !
        imj = 0
        offsetj = 0
        do lj = 0, lmax
          do mj = -lj, lj
            imj = imj + 1
            do j = 1, nbasis
              su = 0 
              do jj = 1, npt
                su = su + res( jj ) * basfcn( jj, j, lj ) * wpt( jj )  * Ylm( jj, imj )
              enddo
              xibb( i + offset , j + offsetj) = su
              write ( 99, '(2i5,3(1x,1e15.8))' ) i + offset, j + offsetj, qtab( i ), qtab( j ), & 
                                                 xibb( i + offset, j + offsetj)
            enddo
            offsetj = offsetj + nbasis
          enddo
        enddo
        write ( 99, * )
      enddo
      offset = offset + nbasis
    enddo
  enddo
  close( unit=99 )
  !



  !
  allocate( tmp( tot_nbasis, tot_nbasis ), xifull( tot_nbasis, tot_nbasis ) )
  tmp = 0
  do i = 1, tot_nbasis
     do j = 1, tot_nbasis
        do k = 1, tot_nbasis
           tmp( i, j ) = tmp( i, j ) - xibb( i, k ) * coulmat( k, j )
        end do
     end do
     tmp( i, i ) = tmp( i, i ) + 1.0d0
  end do
  call rinvert( tot_nbasis, tmp )
  xifull = 0
  do i = 1, tot_nbasis
     do j = 1, tot_nbasis
        do k = 1, tot_nbasis
           xifull( i, j ) = xifull( i, j ) + tmp( i, k ) * xibb( k, j )
        end do
     end do
  end do
  !


  !
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
  !
  im = 0
  offset = 0
  do l = 0, lmax
    do m = -l, l
      im = im + 1
      do i = 1, nbasis
        su = 0
        su0 = 0
        offsetj = 0
        do imj = 1, 1 ! Only s-type comes in rhs mtot
          do j = 1, nbasis
            su = su + xifull( i + offset, j + offsetj ) * rhs( j + offsetj )
            su0 = su0 + xibb( i + offset, j + offsetj ) * rhs( j + offsetj )
            write ( 99, '(2i5,3(1x,1e15.8))' ) i + offset, j + offsetj, qtab( i ), qtab( j ), &
                                               xifull( i + offset, j + offsetj )
          enddo
          offsetj = offsetj + nbasis
        enddo
        write ( 98, '(1i5,2(1x,1e15.8))' ) i + offset, su, su0
        do ii = 1, npt
          ! pot functions end up with m -> -m. Because of the Ylm expansion of V(r)
          nind0( ii ) = nind0( ii ) + su0 * basfcn( ii, i, l ) * Ylm( ii, im )
          vind0( ii ) = vind0( ii ) + su0 * potfcn( ii, i, l ) * Ylm( ii, im )
          nind( ii ) = nind( ii ) + su * basfcn( ii, i, l ) * Ylm( ii, im )
          vind( ii ) = vind( ii ) + su * potfcn( ii, i, l ) * Ylm( ii, im )
        enddo
      enddo
      write ( 99, * )
      offset = offset + nbasis
    enddo
  end do

  close( unit=98 )
  close( unit=99 )
  open( unit=99, file='ninduced_full', form='formatted', status='unknown' )
  rewind 99
  do ii = 1, npt
     write ( 99, '(5(1x,1f10.5))' ) drel( ii ), nind( ii ), nind0( ii ), vind( ii ), vind0( ii )
  end do
  close( unit=99 )
  !
  nind = 0
  nind0 = 0
  vind = 0
  vind0 = 0
    !
  im = 0
  offset = 0
  do l = 0, 0
    do m = -l, l
      im = im + 1
      do i = 1, nbasis
        su = 0
        su0 = 0
        offsetj = 0
        do imj = 1, 1
          do j = 1, nbasis
            su = su + xifull( i + offset, j + offsetj ) * rhs( j + offsetj )
            su0 = su0 + xibb( i + offset, j + offsetj ) * rhs( j + offsetj )
          enddo
          offsetj = offsetj + nbasis
        enddo
        do ii = 1, npt
          ! pot functions end up with m -> -m. Because of the Ylm expansion of V(r)
          nind0( ii ) = nind0( ii ) + su0 * basfcn( ii, i, l ) * Ylm( ii, im )
          vind0( ii ) = vind0( ii ) + su0 * potfcn( ii, i, l ) * Ylm( ii, im )
          nind( ii ) = nind( ii ) + su * basfcn( ii, i, l ) * Ylm( ii, im )
          vind( ii ) = vind( ii ) + su * potfcn( ii, i, l ) * Ylm( ii, im )
        enddo
      enddo
      offset = offset + nbasis
    enddo
  end do
  open( unit=99, file='ninduced', form='formatted', status='unknown' )
  rewind 99
  do ii = 1, npt
     write ( 99, '(5(1x,1f10.5))' ) drel( ii ), nind( ii ), nind0( ii ), vind( ii ), vind0( ii )
  end do
  close( unit=99 )


  !
  deallocate( posn, wpt, vipt, drel, qtab, coulmat, nind )
  !
end program xipostproc
