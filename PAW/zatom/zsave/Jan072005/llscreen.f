      program screen
      implicit none
c
      integer i, idum, nq
c
      double precision q, uq, vq, epsinv
c
      double precision, parameter :: pi = 3.14159265358979323846d0
c
      double precision eps, lam, n, qf, qq, wf, wpsqd
      double precision, external :: levlou
c
      read ( 5, * ) n, eps
      wpsqd = 4.d0 * pi * n
      qf = ( 3.d0 * pi ** 2 * n ) ** 0.33333333d0
      wf = 0.5d0 * qf ** 2
      lam = dsqrt( wpsqd / ( wf ** 2 * ( eps - 1.d0 ) ) )
      open( unit=99, file='ScreenedPot',
     &      form='formatted', status='unknown' )
      rewind 99
      open( unit=98, file='BarePot',
     &      form='formatted', status='unknown' )
      rewind 98
      read ( 98, * ) idum, nq, q, uq, vq
      do i = 1, nq
        qq = q / qf
        epsinv = levlou( qq, qf, lam )
        uq = uq * ( epsinv - 1.d0 )
        vq = vq * ( epsinv - 1.d0 )
        write ( 99, '(2x,2i5,3(2x,1e15.8))' ) i, nq, q, uq, vq
        if ( i .lt. nq ) read ( 98, * ) idum, nq, q, uq, vq
      end do
      close( unit=98 )
      close( unit=99 )
      end
