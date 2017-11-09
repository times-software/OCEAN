  !
  ! load, check hermitian quadrature points
  !
  open( unit=99, file='hqp', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nhqp
  allocate( xhqp( nhqp ), ahqp( nhqp ) )
  do iquad = 1, nhqp
     read ( 99, * ) xhqp( iquad ), ahqp( iquad )
  end do
  close( unit=99 )
  !
  open( unit=99, file='hqpdiag', form='formatted', status='unknown' )
  rewind 99
  tmpquad = 0.5d0 * sqrt( 4.0d0 * atan( 1.0d0 ) )
  do mquad = 0, 15 ! only even powers work for the +ve range
     nquad = 2 * mquad
     suquad = 0
     do iquad = 1, nhqp
        suquad = suquad + xhqp( iquad ) ** nquad * ahqp( iquad )
     end do
     write ( 99, '(2x,1i5,(2(1x,1e15.8)))' ) nquad, suquad, tmpquad
     tmpquad = tmpquad * dble( 2 * mquad + 1 ) / 2.0d0
  end do
  close( unit=99 )
  !
  ! load, check laguerre quadrature points
  !
  open( unit=99, file='lqp', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nlqp
  allocate( xlqp( nlqp ), alqp( nlqp ) )
  do iquad = 1, nlqp
     read ( 99, * ) xlqp( iquad ), alqp( iquad )
  end do
  close( unit=99 )
  !
  open( unit=99, file='lqpdiag', form='formatted', status='unknown' )
  rewind 99
  tmpquad = 1.0d0
  do nquad = 0, 30
     suquad = 0
     do iquad = 1, nlqp
        suquad = suquad + xlqp( iquad ) ** nquad * alqp( iquad )
     end do
     write ( 99, '(2x,1i5,(2(1x,1e15.8)))' ) nquad, suquad, tmpquad
     tmpquad = tmpquad * dble( nquad + 1 )
  end do
  close( unit=99 )
  !
  ! load, check legendre quadrature points
  !
  open( unit=99, file='gauss16', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) npt
  allocate( xpt( npt ), wpt( npt ) )
  suquad = 0
  do iquad = 1, npt
     read ( 99, * ) xpt( iquad ), wpt( iquad )
     xpt( iquad ) = ( 1.0d0 + xpt( iquad ) ) / 2.0d0
     suquad = suquad + wpt( iquad )
  end do
  wpt = wpt / suquad
  close( unit=99 )
  !
  open( unit=99, file='gauss16diag', form='formatted', status='unknown' )
  rewind 99
  do nquad = 0, 50
     suquad = 0
     do iquad = 1, npt
        suquad = suquad + xpt( iquad ) ** nquad * wpt( iquad )
     end do
     write ( 99, '(2x,1i5,2f20.15)' ) nquad, suquad, 1.0d0 / dble( nquad + 1 )
  end do
  close( unit=99 )
