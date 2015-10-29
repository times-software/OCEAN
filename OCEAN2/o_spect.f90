program o_spect
  use AI_kinds
  implicit none

  integer :: ZNL(3), indx, photon
  character(len=2) :: elname, corelevel
  character(len=5) :: calc_type
  character( LEN=21 ) :: lancfile
  character( LEN=21 ) :: abs_filename

  integer :: ie, jdamp, jj, dumi
  real(DP), external :: gamfcn
  real(DP) :: e, gam, dr, di, spct( 0 : 1 ), spkk, pi
  complex(DP) :: rm1, ctmp, disc, delta
  REAL(DP) :: el, eh, gam0, eps, nval,  ebase
  REAL(DP) :: ener, kpref, dumf

  INTEGER  :: ne, n_recur, i_recur, nruns, iter, run_iter

  real( DP ), allocatable :: a(:), b(:)


  open(unit=99,file='spect.in',form='formatted',status='old')
  rewind(99)
  read(99,*) dumi
  read(99,*) dumf
  read(99,*) calc_type
  select case ( calc_type )
    case('hay')
      read(99,*) ne, el, eh, gam0, ebase
    case default
      goto 111
  end select
  close(99)

  el = el / 27.2114d0
  eh = eh / 27.2114d0


  open(unit=99,file='epsilon',form='formatted',status='old')
  rewind 99
  read(99,*) eps
  close(99)

  open( unit=99, file='nval.h', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nval
  close( unit=99 )

  open(unit=97,file='mulfile',form='formatted',status='old')
  


  open(unit=98,file='runlist',form='formatted',status='old')
  rewind(98)
  read(98,*) nruns


  do run_iter = 1, nruns
    read(98,*) ZNL(1), ZNL(2), ZNL(3), elname, corelevel, indx, photon, calc_type
    read(97,*) kpref

    write(lancfile,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'lanceig_', elname, &
        '.', indx, '_', '1s', '_', photon
    select case ( calc_type)
    case( 'XES' )
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xesspct_', elname, &
          '.', indx, '_', '1s', '_', photon
    case( 'XAS' )
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', elname, &
          '.', indx, '_', '1s', '_', photon
    case default
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', elname, &
          '.', indx, '_', '1s', '_', photon
    end select


    open(unit=99,file=lancfile,form='formatted',status='old')
    rewind(99)
    read(99,*) n_recur
    allocate(a(0:n_recur),b(n_recur))
    read(99,*) a(0)
    do i_recur = 1, n_recur
      read(99,*) a(i_recur), b(i_recur)
    enddo
    close(99)
    iter = n_recur


    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
    open( unit=99, file=abs_filename, form='formatted', status='unknown' )
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
          ener = ebase + 27.2114d0 * e
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
       end do
       spkk = kpref * dr / ( dr ** 2 + di ** 2 )
       write ( 99, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
    end do
    close(unit=99)

    deallocate( a, b )

  enddo
  
  close(98)
  close(97)

111 continue

end program
