! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program orthog
  implicit none
  !
  integer, parameter :: stdout = 6, u1dat = 27, u2dat = 28
  integer :: ilftl, ilfth, irgtl, irgth, lolim, hilim
  integer :: ix, ik, ib, ibp, nx, nk, nb, nbv, ispin, nspin
  real( kind = kind( 1.0d0 ) ) :: x
  complex( kind = kind( 1.0d0 ) ) :: w
  complex( kind = kind( 1.0d0 ) ), allocatable :: uu(:,:,:), u(:,:), v(:,:)
  !
  call dptest
  !
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nb, 'nb   ' )
  call igetval( nbv, 'nbv  ' )
  call igetval( nx, 'nx   ' )
  call igetval( nk, 'nk   ' )
  call igetval( ilftl, 'ilftl' )
  call igetval( ilfth, 'ilfth' )
  call igetval( irgtl, 'irgtl' )
  call igetval( irgth, 'irgth' )
  call igetval( nspin, 'nspin' )
  close( unit=99 )
  allocate( uu( nx, nb, nk ), u( nx, nb ), v( nx, nb ) )
  ! 
  open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown' )
  rewind u1dat 
  do ik = 1, nk
     do ib = 1, nb
        read ( u1dat ) uu( : , ib, ik )
     end do
  end do
  close( unit=u1dat )
  !
  open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
  rewind u2dat
  do ispin = 1, nspin
    do ik = 1, nk
       open( unit=99, file='prog', form='formatted', status='unknown' )
       rewind 99
       write ( 99, '(2x,2i6,2x,2i2)' ) ik, nk, ispin, nspin
       close( unit=99 )
       u( :, : ) = uu( :, :, ik )
       do ib = 1, nb
          if( ik .eq. 1 ) then
            w = dot_product( u( :, ib ), u( :, ib ) )
            write(6,*) ib, ib, w
          endif
          call normalize( nx, u( :, ib ) )
          v( :, ib ) = u( :, ib )
          if ( ib .le. nbv ) then
             lolim = 1
          else
             lolim = nbv + 1
          end if
          hilim = ib - 1
          if ( lolim .le. hilim ) then
             do ibp = lolim, hilim
                w = dot_product( v( :, ibp ), v( :, ib ) )
                x = dot_product( v( :, ibp ), v( :, ibp ) )
                w = w / x
                v( : , ib ) = v( : , ib ) - w * v( : , ibp )
                if( ik .eq. 1 ) then
                  write(6,*) ib, ibp, w
                endif
             end do
          end if
          call normalize( nx, v( :, ib ) )
          do ix = 1, nx
             write ( u2dat ) ib, ik, ix, v( ix, ib )
          end do
       end do
    end do
  end do
  close( unit=u2dat )
  !
  deallocate( uu, u, v )
  !
end program orthog
