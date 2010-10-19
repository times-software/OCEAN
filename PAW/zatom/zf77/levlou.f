      function levlou( q, qf, lam )
      implicit none
      double precision q, qf, lam, levlou, lred, l 
      double precision xp, xm, term1, term2, term3, tmp, c0, c2, c4
      double precision, parameter :: pi = 3.1415926535898d0
      if ( q .gt. 0.01d0 ) then
        lred = lam / q ** 2
        xp = ( 2.d0 / q + 1.d0 ) / lred
        xm = ( 2.d0 / q - 1.d0 ) / lred
        l = dlog( ( 1.d0 + xp ** 2 ) / ( 1.d0 + xm ** 2 ) )
        term1 = 1.d0 / q ** 2
        term2 = - 0.5d0 * lred / q * ( datan( xp ) + datan( xm ) )
        term3 = 0.125d0 * ( lred ** 2 + 4.d0 / q ** 2 - 1.d0 ) * l / q
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
