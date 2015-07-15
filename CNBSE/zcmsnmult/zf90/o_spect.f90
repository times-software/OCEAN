! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! # Re-calculate spectra from the Haydock a's and b's
! #
! # Useful for changing broadening or shifting energies for chemical shifts
program o_spect
  implicit none

  integer :: ZNL(3), indx, photon
  character(len=2) :: elname, corelevel
  character(len=5) :: calc_type
  character( LEN=21 ) :: lancfile
  character( LEN=21 ) :: abs_filename

  integer :: ie, jdamp, jj, dumi
  real(kind=(kind(1.d0))), external :: gamfcn
  real(kind=(kind(1.d0))) :: e, gam, dr, di, spct( 0 : 1 ), spkk, pi
  complex(kind=(kind(1.d0))) :: rm1, ctmp, disc, delta
  REAL(kind=(kind(1.d0))) :: el, eh, gam0, eps, nval,  ebase
  REAL(kind=(kind(1.d0))) :: ener, kpref, dumf

  integer :: z, n, l, num, pol
  character(len=3) :: calc
  character(len=21) :: lancfil, spectfil

  INTEGER  :: ne, n_recur, i_recur, nruns, iter, run_iter

  real( kind=(kind(1.d0)) ), allocatable :: a(:), b(:)
  
!  real(kind=kind(1.0d0)), parameter :: Ryd2eV = 13.605693009  ! (84) 2014 CODATA
  real(kind=kind(1.0d0)), parameter :: Ha2eV = 27.211386018  ! (168) 2014 CODATA


  open(unit=99,file='spect.in',form='formatted',status='old')
  rewind(99)
  read(99,*) dumi
  read(99,*) dumf
  read(99,*) calc_type
  select case ( calc_type )
    case('hay')
      read(99,*) ne, el, eh, gam0, ebase
    case default
      write(6,*) 'Incompatible calc type'
      goto 112
  end select
  close(99)

  el = el / Ha2eV
  eh = eh / Ha2eV
  gam0 = gam0 / Ha2eV

!  el = (el -ebase)/ 27.2114d0
!  eh = (eh -ebase)/ 27.2114d0


  open(unit=99,file='epsilon',form='formatted',status='old')
  rewind 99
  read(99,*) eps
  close(99)

  open( unit=99, file='nval.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nval
  close( unit=99 )

!  open(unit=99,file='mulfile',form='formatted',status='old')
!  rewind( 99 )
!  read( 99, * ) kpref
!  close( 99 )
! kpref is stored in lanceigs
  

  open(unit=98,file='runlist',form='formatted',status='old')
  read(98,*) nruns

  do run_iter = 1, nruns

! 3  1  0  Li  1s  64  2  XAS
    read(98,*) z, n, l, elname, corelevel, num, pol, calc
    if( calc .eq. 'XAS' ) then
!abslanc_Li.0035_1s_03
      write(lancfil,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)') 'abslanc_', elname, '.', num, '_', corelevel, '_', pol
      write(spectfil,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)') 'absspct_', elname, '.', num, '_', corelevel, '_', pol
    elseif( calc .eq. 'XES' ) then
      write(lancfil,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)') 'xeslanc_', elname, '.', num, '_', corelevel, '_', pol
      write(spectfil,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)') 'xesspct_', elname, '.', num, '_', corelevel, '_', pol
    else
      goto 111
    endif

    write(6,*) lancfil
    write(6,*) spectfil



    open(unit=99,file=lancfil,form='formatted',status='old')
!    open(unit=99,file='lanceigs',form='formatted',status='old')
    rewind(99)
    read(99,*) n_recur, kpref
    allocate(a(0:n_recur),b(n_recur))
    read(99,*) a(0)
    do i_recur = 1, n_recur
      read(99,*) a(i_recur), b(i_recur)
    enddo
    close(99)
    iter = n_recur


    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
    open( unit=99, file=spectfil, form='formatted', status='unknown' )
    rewind 99
    do ie = 1, 2 * ne, 2
       e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne ) 
       do jdamp = 0, 1 
          gam= gam0 + gamfcn( e, nval, eps ) * dble( jdamp )
          ctmp = e - a( iter - 1 ) + rm1 * gam
          disc = sqrt( ctmp ** 2 - 4 * b( iter ) ** 2 )
          di= -rm1 * disc
          if ( di .gt. 0.0d0 ) then
             delta = ( ctmp + disc ) / 2
          else
             delta = ( ctmp - disc ) / 2
          end if
          do jj = iter - 1, 0, -1
             delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
          end do
          dr = delta
          di = -rm1 * delta
          di = abs( di )
          ener = ebase + Ha2eV * e
!          ener = 27.2114d0 * e
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
       end do
       spkk = kpref * dr / ( dr ** 2 + di ** 2 )
       write ( 99, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
    end do
    close(unit=99)
    goto 111

  enddo

111 continue
  deallocate( a, b )
  close(98)
112 continue  


end program
