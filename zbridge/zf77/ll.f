c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      program llft
      implicit none
c
      integer, parameter :: stdin = 5
      integer, parameter :: whomdat = 30
c
      include 'whom.h'
      double precision remd( nq ), e0, rsh, rsl, rsmax, rsmin
c
      integer irs
      double precision rs
      double precision n, qf, wf, wpsqd, lam
c
      integer ir
      double precision r, rang
c
      integer j
      double precision int, q, qq, dq, pi
      double precision cold, sold, cnew, snew, cs, ss
c
      double complex ra, rb, c1, a, b, aa
      double precision eill, eift, beta, gamma, c0, c2, c4, u, w, raw
      double precision, external :: levlou
c
      pi = 4.0d0 * atan( 1.0d0 )
c static dielectric e0
      read ( stdin, * ) e0
      open( unit=99, file='rsval.h', form='formatted',
     &      status='unknown' )
      call rgetval( rsmax, 'rsmax' )
      call rgetval( rsmin, 'rsmin' )
      close( unit=99 )
      rsl = rsmin * 0.95d0
      rsh = rsmax * 1.05d0
c
      open( unit=whomdat, file='W.els', form='formatted',
     &      status='unknown' )
      rewind whomdat 
c
      do irs = 1, nrs
c
        rs = rsl + ( rsh - rsl ) * dble( irs - 1 ) / dble( nrs - 1 )
c move back from r to the density
c The following mimics Hybertsen and Louie Phys Rev B 37. 2733
        n = 1.d0 / ( 4.d0 * pi * rs ** 3 / 3.d0 )
c calculate fermi momentum associated with n
        qf = ( 3.d0 * pi ** 2 * n ) ** ( 1.d0 / 3.d0 )
c calculate the plasma frequency (squared) asspciated with n
        wpsqd = 4.d0 * pi * n
c fermi energy for n
        wf = 0.5d0 * qf ** 2
c lambda(rs) eq. 6 of above
        lam = dsqrt( wpsqd / ( wf ** 2 * ( e0 - 1 ) ) )
c
        c0 =   1.d0 / e0
        c2 =   c0 ** 2 * 64.d0 / ( 5.d0 * pi * lam ** 4 * qf )
        c4 = - 16.d0 / ( 3.d0 * pi * qf )
c
        beta = - c2 * dabs( c4 ) / ( 1.d0 - c0 ) ** 2
        gamma = c4 / ( c0 - 1.d0 )
c
        u = - 0.5d0 * beta
        w = - 0.25d0 * ( beta ** 2 - 4.d0 * gamma )
c
        raw = dsqrt( dabs( w ) )
        if ( w .ge. 0.d0 ) then
          a = dcmplx( u, - raw )
          b = dcmplx( u, + raw )
          aa = dcmplx( 0.d0, - 0.5d0 * c4 / raw )
        else
          a = dcmplx( u + raw, 0.d0 )
          b = dcmplx( u - raw, 0.d0 )
          aa = dcmplx( - 0.5d0 * c4 / raw, 0.d0 )
        end if
        ra = cdsqrt( a )
        if ( dble( ra ) .lt. 0.d0 ) ra = - ra
        rb = cdsqrt( b )
        if ( dble( rb ) .lt. 0.d0 ) rb = - rb
        c1 = dcmplx( 1.d0, 0.d0 )
        dq = dqq * qf
c
        do j = 1, nq
          qq = dqq * ( dble( j ) - 0.5d0 )
          q = qq * qf
          eill = levlou( qq, qf, lam )
          eift = 1.d0 + c4 / ( ( u + qq ** 2 ) ** 2 + w )
          remd( j ) = ( eill - eift ) / q
        end do
c
        do ir = 1, nr
          rang = rl + ( rh - rl ) * dble( ir - 1 ) / dble( nr - 1 )
          r = rang / 0.529177d0 ! get r into bohrs
c
          open( unit=99, file='prog', form='formatted',
     &          status='unknown' )
          rewind 99
          write ( 99, '(2i5,2x,2i5)' ) ir, nr, irs, nrs
          close( unit=99 )
c
          cs = dcos( dq * r )
          ss = dsin( dq * r )
          q = - dq * 0.5d0
          cold = dcos( q * r )
          sold = dsin( q * r )
          int = 0.d0
          do j = 1, nq
            q = q + dq
            cnew = cold * cs - sold * ss
            snew = sold * cs + cold * ss
            int = int + dq * remd( j ) * snew
            sold = snew
            cold = cnew
          end do
          int = int * ( 2.d0 / pi )
          int = int + c1 + aa / a * ( c1 - cdexp( - ra * qf * r ) )
     &                   - aa / b * ( c1 - cdexp( - rb * qf * r ) )
c
          write ( whomdat, '(3e16.8)' ) rs, rang, int
c
        end do
        write ( whomdat, * )
c
      end do
      close( unit=whomdat )
c
      end
