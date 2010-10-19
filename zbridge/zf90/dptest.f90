subroutine dptest
  implicit none
  !
  complex( kind = kind( 1.0d0 ) ) :: lv( 4 ), rv( 4 ), rm1, ctest
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  write ( 6, '(1a6,2f15.12)' ) 'rm1 = ', rm1
  !
  lv( : ) = 1; rv( : ) = 1; ctest = dot_product( lv, rv )
  write ( 6, '(4(2f10.5,5x))' ) lv( : ), rv( : )
  write ( 6, '(1a7,2(1x,1e15.8))' ) 'dp ==> ', ctest
  ctest = ctest - 4.0d0
  if ( abs( ctest ) .gt. 1.0d-12 ) stop 'dot_product fault!'
  !
  lv( : ) = rm1; rv( : ) = 1; ctest = dot_product( lv, rv )
  write ( 6, '(4(2f10.5,5x))' ) lv( : ), rv( : )
  write ( 6, '(1a7,2(1x,1e15.8))' ) 'dp ==> ', ctest
  ctest = ctest - 4.0d0 * ( -rm1 )
  if ( abs( ctest ) .gt. 1.0d-12 ) stop 'dot_product fault!'
  !
  lv( : ) = 1; rv( : ) = rm1; ctest = dot_product( lv, rv )
  write ( 6, '(4(2f10.5,5x))' ) lv( : ), rv( : )
  write ( 6, '(1a7,2(1x,1e15.8))' ) 'dp ==> ', ctest
  ctest = ctest - 4.0d0 * rm1
  if ( abs( ctest ) .gt. 1.0d-12 ) stop 'dot_product fault!'
  !
  lv( : ) = rm1; rv( : ) = rm1; ctest = dot_product( lv, rv )
  write ( 6, '(4(2f10.5,5x))' ) lv( : ), rv( : )
  write ( 6, '(1a7,2(1x,1e15.8))' ) 'dp ==> ', ctest
  ctest = ctest - 4.0d0
  if ( abs( ctest ) .gt. 1.0d-12 ) stop 'dot_product fault!'
  !
  return
end subroutine dptest
