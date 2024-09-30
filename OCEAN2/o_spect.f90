! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program o_spect
  use AI_kinds
  use OCEAN_constants, only : Hartree2eV, eV2Hartree, alphainv, bohr
  implicit none

  integer :: ZNL(3), indx, photon, rixs_energy, rixs_pol
  character(len=2) :: elname, corelevel
  character(len=5) :: calc_type
  character( LEN=24 ) :: lancfile, photFile
  character( LEN=25 ) :: abs_filename

  integer :: ie, jdamp, jj, dumi, i, j, k, ii
  real(DP), external :: gamfcn
  real(DP), allocatable :: total_spectrum(:,:)
  real(DP) :: e, gam, dr, di, spct( 0 : 1 ), spkk(0:1), pi, avec(3,3), ucvol, ere, reeps, imeps, lossf, fact, mu, reflct
  complex(dp) :: rm1, ctmp, disc, delta, arg, rp, rm, rrr, al, be, eps, refrac, mdelta
  REAL(dp) :: el, eh, gam0, eps0, nval,  ebase, mdr, mdi
  REAL(dp) :: ener, kpref, dumf, indref, indabs, theta, rad, avgWeight, dumv(3)

  INTEGER  :: ne, n_recur, i_recur, nruns, iter, run_iter

  logical :: is_val, ex, do_avg
  real( dp ), allocatable :: a(:), b(:)


  open(unit=99,file='spect.in',form='formatted',status='old')
  rewind(99)
  read(99,*) ne, el, eh, gam0, ebase
  close(99)

  el = el * eV2Hartree !/ 27.2114d0
  eh = eh * eV2Hartree !/ 27.2114d0
  gam0 = gam0 * eV2Hartree !/ 27.2114d0
  ebase = ebase * eV2Hartree

  open(unit=99,file='epsilon',form='formatted',status='old')
  rewind 99
  read(99,*) eps0
  close(99)

  open( unit=99, file='nval.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nval
  close( unit=99 )

  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old')
  rewind 99
  read(99,*) avec(:,:)
  close( 99 )
  ucvol = 0
  do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     ucvol = ucvol + avec( i, 1 ) * avec( j, 2 ) * avec( k, 3 )
     ucvol = ucvol - avec( i, 1 ) * avec( k, 2 ) * avec( j, 3 )
  end do
  ucvol = abs( ucvol )


  open(unit=98,file='runlist',form='formatted',status='old')
  rewind(98)
  read(98,*) nruns

  inquire(file='do_avg.ipt', exist=ex )
  if( ex ) then
    open(unit=99, file='do_avg.ipt', form='formatted', status='old' )
    read(99,* ) do_avg
    close(99 )
  else
    do_avg = .false.
  endif

  if( do_avg ) then
    allocate( total_spectrum( 5, ne ) )
    total_spectrum(:,:) = 0.0_dp
  endif

  do run_iter = 1, nruns
    read(98,*) ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type
    is_val = .false.

    select case ( calc_type)
    case( 'XES' )
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xesspct_', elname, &
          '.', indx, '_', corelevel, '_', photon
      write(lancfile,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xeslanc_', elname, &
          '.', indx, '_', corelevel, '_', photon
    case( 'XAS' )
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', elname, &
          '.', indx, '_', corelevel, '_', photon
      write(lancfile,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'abslanc_', elname, &
          '.', indx, '_', corelevel, '_', photon
    case( 'RXS' )
      backspace( 98 )
      read(98,*)  ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type, &
                  rixs_energy, rixs_pol
      is_val = .true.
      write(abs_filename,'(A8,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxsspct_', elname, &
        '.', corelevel,  '_', photon, '.', rixs_energy, '.', rixs_pol
      write(lancfile, '(A7,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxlanc_', elname, &
        '.', corelevel,  '_', photon, '.', rixs_energy, '.', rixs_pol

    case( 'VAL' )
      write(abs_filename,'(A)' ) 'opcons'
      write(lancfile,'(A)') 'oplanc'
      is_val = .true.

    case default
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', elname, &
          '.', indx, '_', corelevel, '_', photon
      write(lancfile,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'abslanc_', elname, &
          '.', indx, '_', corelevel, '_', photon
    end select

    if( do_avg ) then
      write( photFile, '(A,I0)' ) 'photon', photon
      open( unit=99, file=photFile, form='formatted', status='old' )
      read( 99, * )
      call fancyvector( dumv, dumv, 99 )
      call fancyvector( dumv, dumv, 99 )
      read( 99, * )
      read( 99, * ) avgWeight
      close(99)
    endif

    open(unit=99,file=trim(lancfile),form='formatted',status='old')
    rewind(99)
    read(99,*) n_recur, kpref
    allocate(a(0:n_recur),b(n_recur))
    read(99,*) a(0)
    do i_recur = 1, n_recur
      read(99,*) a(i_recur), b(i_recur)
    enddo
    close(99)

    inquire( file='niter_override', exist=ex)
    if( ex ) then
      open(unit=99, file='niter_override', form='formatted', status='old' )
      read( 99, * ) iter
      close( 99 )
      iter = max( 1, iter )
      iter = min( iter, n_recur )
    else
      iter = n_recur
    endif

    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
    open( unit=99, file=trim(abs_filename), form='formatted', status='unknown' )
    rewind 99

    if( is_val ) then
!    fact = kpref * real( 2 / val_ham_spin, DP ) * ucvol
    fact = kpref * ucvol * 2.0_DP

    write(99,"(a)") "#   omega (eV)      epsilon_1       epsilon_2       n"// &
      "               kappa           mu (cm^(-1))    R"//  &
      "               epsinv"

!p%kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
    do ie = 1, 2 * ne, 2
      ere = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )

      ctmp = cmplx( ere, gam0, DP )

      arg = ( ere - a( iter - 1 ) ) ** 2 - 4.0_dp * b( iter ) ** 2
      arg = sqrt( arg )

      rp = 0.5_dp * ( ere - a( iter - 1 ) + arg )
      rm = 0.5_dp * ( ere - a( iter - 1 ) - arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = rp
      else
        rrr = rm
      endif

      al =  ctmp - a( iter - 1 ) - rrr
      be = -ctmp - a( iter - 1 ) - rrr
      do i = iter-1, 0, -1
        al = ctmp - a( i ) - b( i+1 )**2 / al
        be = -ctmp - a( i ) - b( i+ 1 ) **2 / be
      enddo

      eps = 1.0_dp - fact / al - fact / be

      reeps = dble( eps )
      imeps = aimag( eps )
      lossf = imeps / ( reeps ** 2 + imeps ** 2 )

      rad = sqrt( reeps ** 2 + imeps ** 2 )
      theta = acos( reeps / rad ) / 2
      indref = sqrt( rad ) * cos( theta )
      indabs = sqrt( rad ) * sin( theta )

      refrac = sqrt(eps)
!      reflct = abs((refrac-1.0d0)/(refrac+1.0d0))**2
      reflct = (refrac-1.0d0)/(refrac+1.0d0) * conjg( (refrac-1.0d0)/(refrac+1.0d0) )
      mu = 2.0d0 * ere * Hartree2eV * aimag(refrac) / ( bohr * alphainv * 1000 )

!      write(99,'(8(1E24.16,1X))') ere*Hartree2eV, reeps, imeps, refrac-1.0d0, mu, reflct, lossf
      write(99,'(8(1E24.16,1X))') ere*Hartree2eV, reeps, imeps, indref, indabs, mu, reflct, lossf

    enddo

    else
      write( 99, '(A,1i5,A,1e15.8,A,1e15.8)' ) '#   iter=', iter, '   gam=', gam0, '   kpref=', kpref
      write( 99, '(5(A15,1x))' ) '#   Energy', 'Spect', 'Spect(0)', 'SPKK', 'SPKK(0)'
      ii = 0
      do ie = 0, 2 * ne, 2
         ii = ii + 1
         e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
         do jdamp = 0, 1 
            gam= gam0 + gamfcn( real(e,dp), real(nval,dp),  real(eps0,dp) ) * dble( jdamp )
            ctmp = e - a( iter - 1 ) + rm1 * gam
            disc = sqrt( ctmp ** 2 - 4 * b( iter ) ** 2 )
            di= -rm1 * disc
            if ( di .gt. 0.0d0 ) then
               delta = ( ctmp + disc ) / 2
            else
               delta = ( ctmp - disc ) / 2
            end if
            mdelta = delta - e - rm1*gam -2.0d0*ebase
            do jj = iter - 1, 0, -1
               delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
               mdelta = -e - a( jj ) - rm1 * gam - 2.0d0*ebase - b( jj + 1 ) ** 2 / mdelta
            end do
            dr = delta
            di = -rm1 * delta
            di = abs( di )
  !          ener = ebase + 27.2114d0 * e
!            ener = ebase + Hartree2eV * e

            mdr = mdelta
            mdi = -rm1 * mdelta
            mdi = abs( mdi )
           
!            spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
!            if( ebase .gt. 0.01d0 ) spct( jdamp ) = spct( jdamp ) + kpref * mdi / ( mdr ** 2 + mdi ** 2 )
            eps = - kpref / delta
            if( ebase .gt. 0.01d0 ) eps = eps - kpref / mdelta
            spct( jdamp ) = aimag( eps )
            spkk( jdamp )  = kpref * dr / ( dr ** 2 + di ** 2 )
         end do
         ener = ( e + ebase ) * Hartree2eV
!         if( ebase .gt. 0.01d0 ) spkk = spkk + kpref * mdr / ( mdr ** 2 + mdi ** 2 )
         spkk = - real( eps, DP )
!         write ( 99, '(5(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
         write ( 99, '(5(1e15.8,1x))' ) ener, spct( 1 ), spct( 0 ), spkk( 1 ), spkk( 0 )
         if( do_avg ) then
           total_spectrum(1,ii) = ener
           total_spectrum(2,ii) = total_spectrum(2,ii) + spkk(0)*avgWeight
           total_spectrum(3,ii) = total_spectrum(3,ii) + spct(0)*avgWeight
           total_spectrum(4,ii) = total_spectrum(4,ii) + spkk(1)*avgWeight
           total_spectrum(5,ii) = total_spectrum(5,ii) + spct(1)*avgWeight
         endif

      end do
      close(unit=99)
    endif

    deallocate( a, b )

  enddo
  
  close(98)

  if( do_Avg ) then
    open(unit=99, file='TotAvg', form='formatted' )
    do ie = 1, ne
      write(99, '(5(E24.16,1X))' ) total_spectrum(:,ie)
    enddo
    close(99)
    deallocate( total_spectrum )
  endif

111 continue

end program
