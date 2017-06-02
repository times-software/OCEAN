! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine kcmprs( longr, longi, shortr, shorti, nq, nbd, imark )
  implicit none
  !
  integer nq, nbd, imark
  double precision shortr( nq * nbd ), shorti( nq * nbd )
  double precision longr( nq * nbd )
  double precision longi( nq * nbd )
  !
  integer ii
  !
  do ii=1,nq*nbd
     shortr(ii)=longr(ii)
     shorti(ii)=longi(ii)
  end do
  !
  return
end subroutine kcmprs
!
subroutine kexpnd( longr, longi, shortr, shorti, nq, nbd, imark )
  implicit none
  !
  integer nq, nbd, imark
  double precision shortr( nq * nbd ), shorti( nq * nbd )
  double precision longr( nq * nbd )
  double precision longi( nq * nbd )
  !
  integer ii
  !
  do ii=1,nq*nbd
     longr( ii ) = shortr( ii )
     longi( ii ) = shorti( ii )
  end do
  !
  return
end subroutine kexpnd
!
subroutine nqtoxr( nx, ny, nz, nfft, nv, nq, nbd, indx, vecr, veci, cor, coi, xwrkr, xwrki, tau, & 
                   nn1, zn, wrk, jfft, ur, ui, qpred )
  implicit none
  !
  integer nx, ny, nz, nfft, nv, nq, nbd, nn1, zn( 3 ), jfft
  integer indx( nq )
  !
  real cor( nq * nbd, nv ), coi( nq * nbd, nv )
  real( kind = kind( 1.0d0 ) ) :: qpred( 3, nq )
  !
  double precision vecr( nq * nbd ), veci( nq * nbd ) 
  double precision xwrkr( nfft, nx * ny * nz )
  double precision xwrki( nfft, nx * ny * nz )
  double precision tau( 3 ), wrk( jfft )
  double precision ur( nx * ny * nz, nv ), ui( nx * ny * nz, nv )
  !
  integer i, j, ix
  double precision tmpr, tmpi
  double precision, allocatable :: fr( : ), fi( : )
  !
  double precision, parameter :: pi = 3.14159265358979323846d0
  !
  xwrkr = 0
  xwrki = 0
  !
  allocate( fr( nfft ) , fi( nfft ) )
  do j = 1, nv
     call recoef(nfft,nq,nbd,indx,fr,fi, vecr,veci,cor(1,j),coi(1,j))
     do ix = 1, nx * ny * nz
        tmpr = ur( ix, j )
        tmpi = ui( ix, j )
        do i = 1, nfft
           xwrkr( i, ix ) = xwrkr( i, ix ) + fr( i ) * tmpr - fi( i ) * tmpi
           xwrki( i, ix ) = xwrki( i, ix ) + fr( i ) * tmpi + fi( i ) * tmpr
        end do
     end do
  end do
  deallocate( fr, fi )
  !
  call krphase( nx, ny, nz, nq, tau, +1.d0, indx, qpred, xwrkr, xwrki )
  !
  do i = 1, nx * ny * nz
     call cfft(xwrkr(1,i),xwrki(1,i),nn1, zn(1),zn(2),zn(3),-1,wrk,jfft)
  end do
  !
  return
end subroutine nqtoxr
!
subroutine xrtonq( nx, ny, nz, nn1, zn, jfft, nq, nv, nbd, tau, indx, xwrkr, xwrki, wrk, cor, coi, hvr, hvi, ur, ui, qpred )
  implicit none
  !
  integer nx, ny, nz, nn1, zn( 3 ), jfft, nq, nv, nbd
  integer indx( nq )
  real( kind = kind( 1.0d0 ) ) :: qpred( 3, nq )
  !
  real cor( nq * nbd, nv ), coi( nq * nbd, nv )
  !
  double precision tau( 3 ), wrk( jfft )
  double precision xwrkr( nq, nx * ny * nz )
  double precision xwrki( nq, nx * ny * nz )
  double precision hvr( nq * nbd ), hvi( nq * nbd )
  double precision ur( nx * ny * nz, nv ), ui( nx * ny * nz, nv )
  !
  integer :: i, j, k, ii, ib
  double precision tmpr, tmpi
  ! 
  double precision, allocatable :: ywrkr( : ), ywrki( : )
  double precision, allocatable :: zwrkr( : , : ), zwrki( : , : )
  !
  double precision, parameter :: pi = 3.14159265358979323846d0
  !
  do i = 1, nx * ny * nz
     call cfft(xwrkr(1,i),xwrki(1,i),nn1, zn(1),zn(2),zn(3),+1,wrk,jfft)
  end do
  !
  call krphase( nx, ny, nz, nq, tau, -1.d0, indx, qpred, xwrkr, xwrki )
  !
  allocate( ywrkr( nx * ny * nz ), zwrkr( nv, nq ) )
  allocate( ywrki( nx * ny * nz ), zwrki( nv, nq ) )
  do i = 1, nq
     do j = 1, nx * ny * nz
        ywrkr(j)=xwrkr(indx(i),j)
        ywrki(j)=xwrki(indx(i),j)
     end do
     do j = 1, nv
        tmpr = 0.d0
        tmpi = 0.d0
        do k = 1, nx * ny * nz
           tmpr = tmpr + ur( k, j ) * ywrkr( k ) + ui( k, j ) * ywrki( k )
           tmpi = tmpi + ur( k, j ) * ywrki( k ) - ui( k, j ) * ywrkr( k )
        end do
        zwrkr( j, i ) = tmpr
        zwrki( j, i ) = tmpi
     end do
  end do
  do i=1,nq*nbd
     hvr(i)=0.d0
     hvi(i)=0.d0
  end do
  do j=1,nv
     ii=0
     do i=1,nq
        tmpr = zwrkr( j, i )
        tmpi = zwrki( j, i )
        do ib=1,nbd
           ii=ii+1
           hvr(ii)=hvr(ii)-( cor(ii,j)*tmpr+coi(ii,j)*tmpi )
           hvi(ii)=hvi(ii)-( cor(ii,j)*tmpi-coi(ii,j)*tmpr )
        end do
     end do
  end do
  deallocate(ywrkr,ywrki,zwrkr,zwrki)
  !
  return
end subroutine xrtonq
!
subroutine krphase( nx, ny, nz, nq, tau, way, indx, qpred, xr, xi )
  implicit none
  !
  integer nx, ny, nz, nq
  integer indx( nq )
  real( kind = kind( 1.0d0 ) ) :: qpred( 3, nq )
  !
  double precision tau( 3 ), way
  double precision xr( nq, nx * ny * nz )
  double precision xi( nq, nx * ny * nz )
  !
  integer j, ix, iy, iz, i
  double precision xpt1, xpt2, xpt3, be1, be2, be3
  double precision phase, phsr, phsi, tmpr, tmpi
  double precision, parameter :: pi = 3.14159265358979323846d0
  !
  j = 0
  do ix = 1, nx
     xpt1 = dble( ix - 1 ) / dble( nx )
     do iy = 1, ny
        xpt2 = dble( iy - 1 ) / dble( ny )
        do iz = 1, nz
           xpt3 = dble( iz - 1 ) / dble( nz )
           j = j + 1
           do i = 1, nq
              be1 = qpred( 1, i )
              be2 = qpred( 2, i )
              be3 = qpred( 3, i )
              phase = 2.0d0 * pi * ( be1 * xpt1 + be2 * xpt2 + be3 * xpt3 )
              phsr = cos( phase )
              phsi = sin( phase ) * way
              tmpr = xr( indx(i), j )
              tmpi = xi( indx(i), j )
              xr( indx(i), j ) = tmpr * phsr - tmpi * phsi
              xi( indx(i), j ) = tmpr * phsi + tmpi * phsr
           end do
        end do
     end do
  end do
  !
  return
end subroutine krphase
!
subroutine xrvmult( nx, ny, nz, zn, tau, amet, epsi, xwrkr, xwrki, ptab )
  implicit none
  !
  integer nx, ny, nz, zn( 3 )
  double precision epsi, tau( 3 ), amet( 3, 3 )
  double precision ptab( 100 )
  double precision xwrkr(zn( 1 ) * zn( 2 ) * zn( 3 ), nx * ny * nz )
  double precision xwrki(zn( 1 ) * zn( 2 ) * zn( 3 ), nx * ny * nz )
  !
  integer i, j, ii, jj, ix, iy, iz, k1, k2, k3, kk1, kk2, kk3
  double precision r, alf( 3 ), potn, frac
  double precision, external :: wav
  !
  i = 0
  do ix = 1, nx
     do iy = 1, ny
        do iz = 1, nz
           i = i + 1
           j = 0
           do k3 = 1, zn( 3 )
              kk3 = k3 - 1
              if ( kk3 .ge. zn( 3 ) / 2 ) kk3 = kk3 - zn( 3 )
              do k2 = 1, zn( 2 )
                 kk2 = k2 - 1
                 if ( kk2 .ge. zn( 2 ) / 2 ) kk2 = kk2 - zn( 2 )
                 do k1 = 1, zn( 1 )
                    kk1 = k1 - 1
                    if ( kk1 .ge. zn( 1 ) / 2 ) kk1 = kk1 - zn( 1 )
                    j = j + 1
                    alf(1) = dble( kk1 ) + dble( ix - 1 ) / dble( nx ) - tau(1)
                    alf(2) = dble( kk2 ) + dble( iy - 1 ) / dble( ny ) - tau(2)
                    alf(3) = dble( kk3 ) + dble( iz - 1 ) / dble( nz ) - tau(3)
                    r = 0.d0
                    do ii = 1, 3
                       do jj = 1, 3
                          r = r + alf( ii ) * alf( jj ) * amet( ii, jj ) 
                       end do
                    end do
                    r = dsqrt( r )
                    potn = 0.d0
                    if ( r .ge. 9.9d0 ) then
                       potn = wav( r, epsi )
                    else
                       ii = 1. + 10. * r
                       frac = 10.d0 * ( r - 0.1d0 * dble( ii - 1 ) )
                       potn = ptab( ii ) + frac * ( ptab( ii + 1 ) - ptab( ii ) )
                    end if
                    xwrkr( j, i ) = xwrkr( j, i ) * potn
                    xwrki( j, i ) = xwrki( j, i ) * potn
                 end do
              end do
           end do
        end do
     end do
  end do
  !
  return
end subroutine xrvmult
