      program whom0
      implicit none
c
      include 'whom.h'
      integer, parameter :: wdat = 30
      integer, parameter :: w0dat = 31
c
      double precision whom( nr, nrs ), d( nr ), rs( nrs )
c
      integer ir, irs, jl, jh, i, irad
c
      double precision radl, radh, num, den, rho, drho, w, brack, r
c
      open( unit=wdat, file='W.els', form='formatted', status='old' )
      rewind wdat
      do irs = 1, nrs
        do ir = 1, nr
          read ( wdat, * ) rs( irs ), d( ir ),  whom( ir, irs )
          whom( ir, irs ) = whom( ir, irs ) - 1.d0
        end do
      end do
      close(unit=wdat)
c
      open( unit=w0dat, file='W0.els', form='formatted',
     &      status='unknown' )
      rewind w0dat
c
      do irs = 1, nrs
c
        open( unit=99, file='prog', form='formatted',
     &        status='unknown' )
        rewind 99
        write ( 99, '(2i5)' ) irs, nrs
        close( unit=99 )
c
        do irad = 1, nrad
          r = rsphl + ( rsphh - rsphl ) *
     &                dble( irad - 1 ) / dble( nrad - 1 )
          num= 0.d0
          den= 0.d0
          drho = 2.d0 * r / dble( nstep )
          jl = 1
          radl = d( 1 )
          jh = 2
          radh = d( 2 )
          do i = 1, nstep
            rho = drho * ( dble( i ) - 0.5d0 )
            do while ( radh .lt. rho )
              jl = jh
              radl = radh
              jh = jh + 1
              radh = d( jh )
            end do
            if ( radl .gt. rho ) then
              w = whom( jl, irs ) * rho / radl
            else
              w = whom( jl, irs ) +
     &            ( whom( jh, irs ) - whom( jl, irs ) ) * 
     &            ( rho - radl ) / ( radh - radl )
            end if
            brack = ( r ** 3 - 0.125d0 * rho ** 3 ) / 3.d0 -
     &              rho * ( r ** 2 - 0.25d0 * rho ** 2 ) / 4.d0
            num = num + drho * brack * rho * w
            den = den + drho * brack * rho ** 2
          end do
          write ( w0dat, '(3f16.10)' )
     &      rs( irs ), r , num / den + 1.2d0 / r
        end do
        write ( w0dat, * )
      end do
c
      close( unit=w0dat )
c
      end
