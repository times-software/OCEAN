! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! This can be run in the RIXS directory. It will re-use runlist (from the valence calc)
! and cnbse.gmres.erange. The only new input is spect.in
! 
! spect.in : gives the emission energy range and steps. The energies are with respect to 
! the zero of the abssorption, so something like -50 to 5 would be reasonable for a wide 
! range of materials. 
!
! output: rixs_plane.txt, suitable for plotting with gnuplot
! set pm3d map
! splot 'rixs_plane.txt'
! 
program rixs_plane
  use AI_kinds
  use OCEAN_constants, only : Hartree2eV, eV2Hartree
  implicit none

  integer :: ZNL(3), indx, photon
  character(len=2) :: elname, corelevel
  character(len=5) :: calc_type
  character( LEN=40 ) :: lancfile
  character( LEN=21 ) :: abs_filename

  integer :: ie, jdamp, jj, dumi, i, j, k, nspin, nHERFD
  real(DP), external :: gamfcn
  real(DP) :: e, gam, dr, di, spct( 0 : 1 ), spkk, pi, ein, eloss, omega, avec(3,3)
  complex(DP) :: rm1, ctmp, disc, delta
  REAL(DP) :: el, eh, gam0, ebase, estart, estop, estep
  REAL(DP) :: kpref, dumf, ere, reeps, imeps, lossf, fact, mu, reflct
  complex(DP) :: arg, rp, rm, rrr, al, be, eps, refrac

  INTEGER  :: ne, n_recur, i_recur, nruns, iter, run_iter, rixs_energy, rixs_pol

  real( DP ), allocatable :: a(:), b(:), herfd(:,:)


  open(unit=99,file='spect.in',form='formatted',status='old')
  rewind(99)
  read(99,*) ne, el, eh, gam0
  close(99)

  el = el * eV2Hartree !/ 27.2114d0
  eh = eh * eV2Hartree !/ 27.2114d0
  gam0 = gam0 * eV2Hartree !/ 27.2114d0

  open(unit=99,file='cnbse.gmres.erange',form='formatted',status='old')
  rewind 99
  read(99,*) estart, estop, estep
  close( 99 )

  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
  rewind 99
  read( 99, * ) avec( :, : )
  close( unit=99 )

  open(unit=99, file='nspin', form='formatted',status='old')
  rewind 99
  read(99,*) nspin
  close(99)
  

  omega = 0
  do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     omega = omega + avec( i, 1 ) * avec( j, 2 ) * avec( k, 3 )
     omega = omega - avec( i, 1 ) * avec( k, 2 ) * avec( j, 3 )
  end do
  omega = abs( omega )


  write(6,*) omega


  open(unit=98,file='runlist',form='formatted',status='old')
  rewind(98)
  read(98,*) nruns


  gam = ( ( (eh - el ) * Hartree2eV - ( estop - estart ) ) / estep )
  nHERFD = ceiling( gam )
  allocate( herfd( nruns, nHERFD ) )
  herfd = 0.0_DP

  open( unit=99, file='rixs_plane.txt', form='formatted', status='unknown' )
  rewind 99

  do run_iter = 1, nruns
    read(98,*) ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type, rixs_energy, rixs_pol

    select case ( calc_type)
    case( 'RXS' )
!      write(abs_filename,'(A8,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxsspct_', elname, &
!            '.', corelevel, '_', photon, '.', rixs_energy, '.', rixs_pol
      write(lancfile,'(A7,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxlanc_', elname, &
                '.', corelevel, '_', photon, '.', rixs_energy, '.', rixs_pol

    case( 'C2C' )

      write(lancfile,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'ctclanc_', elname, &
            '.', indx, '_', corelevel, '_', photon, '.', &
            rixs_energy, '.', rixs_pol


    case default
      write(6,*) 'Wrong runtypes in runlist'
      stop
    end select

    write(6,*) lancfile

    open(unit=97,file=lancfile,form='formatted',status='old')
    rewind(97)
    read(97,*) n_recur, kpref
    allocate(a(0:n_recur),b(n_recur))
    read(97,*) a(0)
    do i_recur = 1, n_recur
      read(97,*) a(i_recur), b(i_recur)
    enddo
    close(97)
    iter = n_recur

    fact = kpref * omega * real( 2 / nspin, DP ) 

    ein = estart + estep * dble( run_iter - 1 )
    ein = ein * eV2Hartree

    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )


    select case ( calc_type )
    case( 'RXS' )

      do ie = 1, 2 * ne, 2
         
        e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )

        eloss = e !+ ein 

        ctmp = cmplx( eloss, gam0, DP )

        arg = ( eloss - a( iter - 1 ) ) ** 2 - 4.0_dp * b( iter ) ** 2
        arg = sqrt( arg )

      
        rp = 0.5_dp * ( eloss - a( iter - 1 ) + arg )
        rm = 0.5_dp * ( eloss - a( iter - 1 ) - arg )
        if( aimag( rp ) .lt. 0.0_dp ) then
          rrr = rp
        else
          rrr = rm
        endif

        al =  ctmp - a( iter - 1 ) - rrr
        be = -ctmp - a( iter - 1 ) - rrr

        do i = iter-1, 0, -1
          al =  ctmp - a( i ) - b( i + 1 ) ** 2 / al
          be = -ctmp - a( i ) - b( i + 1 ) ** 2 / be
        enddo

        eps = 1.0_dp - fact / al - fact / be

        if( eloss .le. 0.0_dp ) eps = 0.0d0

        reeps = dble( eps )
        imeps = aimag( eps )


        write ( 99, '(4(1e15.8,1x))' ) e*Hartree2eV, ein*Hartree2eV, imeps, eloss*Hartree2eV
   
      end do

    case ('C2C')

      do ie = 1, 2 * ne, 2
         e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
            ctmp = e - a( iter - 1 ) + rm1 * gam0
            disc = sqrt( ctmp ** 2 - 4 * b( iter ) ** 2 )
            di= -rm1 * disc
            if ( di .gt. 0.0d0 ) then
               delta = ( ctmp + disc ) / 2
            else
               delta = ( ctmp - disc ) / 2
            end if
            do jj = iter - 1, 0, -1
               delta = e - a( jj ) + rm1 * gam0 - b( jj + 1 ) ** 2 / delta
            end do
            dr = delta
            di = -rm1 * delta
            di = abs( di )
            imeps = kpref * di / ( dr ** 2 + di ** 2 )
         spkk = kpref * dr / ( dr ** 2 + di ** 2 )
         write ( 99, '(4(1e15.8,1x))' ) ein*Hartree2eV, e*Hartree2eV, imeps, spkk
      end do


      do ie = 1, nHERFD
        e = el + dble( ie - 1  + run_iter - 1 ) * estep * eV2Hartree

        ctmp = e - a( iter - 1 ) + rm1 * gam0
        disc = sqrt( ctmp ** 2 - 4 * b( iter ) ** 2 )
        di= -rm1 * disc
        if ( di .gt. 0.0d0 ) then
           delta = ( ctmp + disc ) / 2
        else
           delta = ( ctmp - disc ) / 2
        end if
        do jj = iter - 1, 0, -1
          delta = e - a( jj ) + rm1 * gam0 - b( jj + 1 ) ** 2 / delta
        end do
        dr = delta
        di = -rm1 * delta
        di = abs( di )
        imeps = kpref * di / ( dr ** 2 + di ** 2 )
        herfd( run_iter, ie ) = imeps
      enddo
      
    end select

    write(99,*) ''

    deallocate( a, b )

  enddo
  
  close(unit=99)
  close(98)


  do ie = 1, nHERFD
    write(lancfile,'(A5,I5.5)' ) 'herfd', ie
    open( unit=99, file=lancfile, form='formatted', status='unknown' )
    rewind (99 )

    do run_iter = 1, nruns
      ein = estart + estep * dble( run_iter - 1 )

      write(99, '(2(1e15.8,1x))' ) ein, herfd( run_iter, ie )
    enddo
    close( 99 )
  enddo
    

  deallocate( herfd )

111 continue

end program
