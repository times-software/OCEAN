      function levlou( q, qf, lam )
      implicit none
      double precision q, qf, lam, levlou
      double precision xp, xm, term1, term2, term3
      double precision tmp, c0, c2, c4
      real( kind = kind( 1.0d0 ) ), parameter :: pi = 3.1415926535897932384d0
!     include 'general.h'
      if ( q .gt. 0.01d0 ) then
        xp = q * ( 2.d0 + q ) / lam
        xm = q * ( 2.d0 - q ) / lam
        term1 =   1.d0 / q ** 2
        term2 = - lam / ( 2.d0 * q ** 3 ) *
     &             ( datan( xp ) + datan( xm ) )
        term3 =   ( lam ** 2 / ( 8.d0 * q ** 5 ) +
     &              1.d0 / ( 2.d0 * q ** 3 ) -
     &              1.d0 / ( 8.d0 * q ) ) *
     &            dlog( ( 1.d0 + xp ** 2 ) / ( 1.d0 + xm ** 2 ) )
        tmp = term1 + term2 + term3
      else
        c0 =   2.6666667d0 / lam ** 2
        c2 = - 6.4d0 / lam ** 4
        c4 =   ( 384.d0 / lam ** 6 - 56.d0 / lam ** 4 ) / 21.d0
        tmp = c0 + c2 * q ** 2 + c4 * q ** 4
      end if
      levlou = 1.d0 / ( 1.d0 + 2.d0 / ( pi * qf ) * tmp )
      return
      end
