      program backtran
      implicit none
c
      integer i, idum, nq, nr, j, k, ii
      double precision qval, uval, dum, f, g, rval, fact
      double precision rcut, eps, dq, cval, vval
      double precision w( 0 : 4 )
      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, allocatable :: q( : ), uq( : )
c
      read ( 5, * ) rcut, eps
      fact = 1.d0 - 1.d0 / eps
c
      open( unit=99, file='recipfile.ipt',
     &      form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) idum, nq, qval, uval
      allocate( q( nq ), uq( nq ) )
      q( 1 ) = qval
      uq( 1 ) = uval
      do i = 2, nq
        read ( 99, * ) idum, idum, q( i ), uq( i )
      end do
      close( unit=99 )
c
      dq = q( 2 ) - q( 1 )
      w( 0 ) = 14.d0 * dq / 45.d0
      w( 1 ) = 64.d0 * dq / 45.d0
      w( 2 ) = 24.d0 * dq / 45.d0
      w( 3 ) = 64.d0 * dq / 45.d0
      w( 4 ) = 14.d0 * dq / 45.d0
c
      open( unit=98, file='realfile.ipt',
     &      form='formatted', status='unknown' )
      rewind 98
      open( unit=99, file='realfile.opt',
     &      form='formatted', status='unknown' )
      rewind 99
      i = 0
      nr = 1
      do while ( i .lt. nr )
        read ( 98, * ) rval, cval, dum, dum, dum, nr
        if ( rval .le. rcut ) then
          f = 0.d0
          do j = 0, nq - 4, 4
            do k = 0, 4
              ii = j + k
              if ( ii .eq. 0 ) then
                g = rval * fact
              else
                qval = q( ii )
                g = dsin( rval * qval ) / qval * uq( ii )
              end if
              f = f + w( k ) * g
            end do
          end do
          vval = 2.d0 / ( pi * rval ) * f
        else
          vval = fact / rval
        end if
        i = i + 1
        write ( 99, '(2x,4(2x,1e15.8),2x,1i5)' )
     &    rval, vval, cval, - ( cval + vval ), i
      end do
      close( unit=98 )
      close( unit=99 )
c
      end
