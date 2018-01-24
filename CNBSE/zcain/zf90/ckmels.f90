subroutine ckmels( nr, zz, lc, lmin, lmax, npmax, nproj, phae, phps, r, dl, wc, jl, jlpow )
  implicit none
  !
  integer, intent(in) :: nr, zz, lc, lmin, lmax, npmax
  integer :: nproj( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: dl
  real( kind = kind( 1.0d0 ) ) :: r( nr ), wc( nr ), jl( nr, 0 : lc + lmax ), jlpow( nr, 0 : lc + lmax )
  real( kind = kind( 1.0d0 ) ) :: phae( nr, npmax, lmin : lmax ), phps( nr, npmax, lmin : lmax )
  !
  integer :: l, iener, nener, idum, iphotl, ip, ir
  real( kind = kind( 1.0d0 ) ) :: dum, aemel, rcmel, pomel
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: dr, coeff, ener
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: phiae, phips
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: meltab
  character( len=8 ) :: s8
  character( len=80 ) :: fnam
  !
  allocate( dr( nr ), coeff( npmax ) )
  dr( : ) = 0
  do ir = 1, nr - 4, 4
     dr( ir + 0 ) = dr( ir + 0 ) + dl * r( ir + 0 ) * 14.0d0 / 45.0d0
     dr( ir + 1 ) = dr( ir + 1 ) + dl * r( ir + 1 ) * 64.0d0 / 45.0d0
     dr( ir + 2 ) = dr( ir + 2 ) + dl * r( ir + 2 ) * 24.0d0 / 45.0d0
     dr( ir + 3 ) = dr( ir + 3 ) + dl * r( ir + 3 ) * 64.0d0 / 45.0d0
     dr( ir + 4 ) = dr( ir + 4 ) + dl * r( ir + 4 ) * 14.0d0 / 45.0d0
  end do
  !
  do l = lmin, lmax
     write ( fnam, '(1a4,1i1,1a1,1i3.3)' ) 'phrc', l, 'z', zz
     open( unit=99, file=fnam, form='formatted', status='old' )
     rewind 99
     read ( 99, * ) s8, nener, idum
     allocate( ener( nener ), phiae( nr, nener ), phips( nr, nener ) )
     allocate( meltab( 3, nener ) )
     do iener = 1, nener
        do ir = 1, nr
           read ( 99, * ) dum, ener( iener ), phiae( ir, iener ), dum, phips( ir, iener )
        end do
     end do
     close( unit=99 )
     do iphotl = iabs( lc - l ), lc + l, 2 
        do iener = 1, nener
           coeff( 1 : nproj( l ) ) = 0
           do ip = 1, nproj( l )
              do ir = 1, nr
                 coeff( ip ) = coeff( ip ) + dr( ir ) * ( r( ir ) * phps( ir, ip, l ) ) * phips( ir, iener )
              end do
           end do
           aemel = 0
           rcmel = 0
           pomel = 0
           do ir = 1, nr
              aemel = aemel + dr( ir ) * wc( ir ) * jl( ir, iphotl ) * phiae( ir, iener )
              do ip = 1, nproj( l )
                 rcmel = rcmel + dr( ir ) * wc( ir ) * jl( ir, iphotl ) * coeff( ip ) * ( r( ir ) * phae( ir, ip, l ) )
                 pomel = pomel + dr( ir ) * wc( ir ) * jlpow( ir, iphotl ) * coeff( ip ) * ( r( ir ) * phae( ir, ip, l ) )
              end do
           end do
           meltab( 1, iener ) = aemel
           meltab( 2, iener ) = rcmel
           meltab( 3, iener ) = pomel
        end do
        write ( fnam, '(1a10,1i1,1i2.2)' ) 'pawmeldiag', l, iphotl
        open( unit=99, file=fnam, form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(1a80)' ) '#    energy, mels for allowed l vals, reconstructed mels for allowed l vals'
        do iener = 1, nener
           write ( 99, '(15(1x,1e15.8))' ) ener( iener ), meltab( 1, iener ), meltab( 2, iener ), &
                meltab( 3, iener )
        end do 
        close( unit=99 )
     end do
     deallocate( ener, phiae, phips, meltab )
  end do
  deallocate( dr, coeff )
  !
  return
end subroutine ckmels
