      subroutine transp( l, norb, nr, qmax, dq, nq, dl, irc )
      implicit none
c
                            integer :: l, norb, nr, irc
                   double precision :: qmax, dq, dl
c
                            integer :: i, j, m, nq
      double precision, allocatable :: proj(:,:), f(:), r(:), besq(:)
      double precision, allocatable :: q(:), p(:,:)
                      character * 3 :: name
         double precision, external :: spherbes
            
c    
      allocate( r( nr ), proj( nr, norb ), f( nr ), besq( nr ) )
c
      write ( unit=name, fmt='(1a2,1i1)' ) 'pr', l
      open( unit= 97, file=name,
     &      form='formatted', status='unknown' )
      rewind 97
      do i = 1, nr
         read ( 97, * ) r( i ), ( proj( i, j ), j = 1, norb )
      end do
      close( unit=97 )
c
      nq = 1 + qmax / dq
      allocate( q( nq ), p( nq, norb ) )
      q( 1 ) = 0.d0
      do i = 2, nq
         q( i ) = q( i - 1 ) + dq
      end do
c
      write ( unit=name, fmt='(1a2,1i1)' ) 'ft', l
      open( unit=96, file=name, 
     &      form='formatted', status='unknown' )
      rewind 96
      do m = 1, nq
         do i = 1, nr
            besq( i ) = spherbes( q( m ) * r( i ), l )
         end do
         do j = 1, norb
            do i = 1, nr
               f( i ) = proj( i, j ) * besq( i )
            end do
            call bintegrate( nr, r, dl, f, p( m, j ), irc )
         end do
         write ( 96, '(8f22.11)' ) q( m ), ( p( m, j ), j = 1, norb )
      end do
      close ( unit = 96 )
c
      deallocate( r, proj, f, besq, q, p )
c
      return
      end
