      subroutine mkkxfile( oldnel, newnel )
      implicit none
c
                            integer :: oldnel, newnel
c
                            integer :: l, l1, ll, lh, lmax
                            integer :: i, j, nc, idum
                   double precision :: dum
      double precision, allocatable :: lead( : ) 
      double precision, allocatable :: exch( : ) 
      double precision, allocatable :: eint( : ) 
      double precision, allocatable :: rslt( : ) 
c
      nc = oldnel - newnel
      lmax = newnel - 1
      allocate( lead( 0 : lmax ), rslt( 0 : lmax ) )
      allocate( eint( 0 : lmax ) )
c
      open( unit=99, file='psld',
     &      form='formatted', status='unknown' )
      rewind 99
      do l = 0, lmax
         read ( 99, * ) l1
         if ( l1 .ne. l ) stop ' wrong l order !!! '
         read ( 99, * ) dum, lead( l )   
         write ( 6, '(2x,1i5,2x,1e15.8)' ) l, lead( l )
      end do
      close( unit=99 )
      open( unit=99, file='aex',
     &      form='formatted', status='unknown' )
      rewind 99
      do i = 1, oldnel
         do j = 1, oldnel
            read ( 99, * ) idum, idum, ll, lh
            allocate( exch( ll : lh ) )
            do l = ll, lh
               exch( l ) = 0.d0
            end do 
            do l = ll, lh, 2
               read ( 99, '(1e15.8)' ) exch( l )
            end do
            do l = 0, lmax
               if ( ( i .eq. 1 ) .and. ( j .eq. nc + l + 1 ) ) then
                  eint( l ) = exch( l ) / dble( 2 * l + 1 )
                  write ( 6, '(2x,1i5,2x,1e15.8)' ) l, eint( l )
               end if
            end do          
            deallocate( exch )
         end do
      end do
      close( unit=99 )
      do l = 0, lmax
         rslt( l ) = eint( l ) / lead( l ) ** 2
      end do
      write ( 6, '(2x,5(2x,1e15.8))' ) ( rslt( l ), l = 0, lmax )
      open( unit=99, file='rslt',
     &      form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(2x,5(2x,1e15.8))' ) ( rslt( l ), l = 0, lmax )
      close( unit=99 )
      deallocate( rslt, eint, lead )
c
      return
      end
