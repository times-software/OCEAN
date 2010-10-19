  integer, allocatable :: iq( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: qphys( : , : )
  real, allocatable :: cor(:,:), coi(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: e0(:)
  real( kind = kind( 1.0d0 ) ), allocatable :: ur(:,:), ui(:,:)
  character * 6 :: str
  integer :: lv, nc, lmin, lmax, nq, nbd, ip, nptot
  integer :: icms, icml, ivms, ivml, nu, ic
  integer, allocatable :: nproj( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: cms( : ), cml( : ), vms( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: hcml( : ), hvml( : ), hcms( : ), hvms( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: mhr( : ), mhi( : )
  integer, allocatable :: hvnu( : )
  integer :: l, m, n
  real( kind = kind( 1.0d0 ) ), allocatable :: pwr( : ), pwi( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: hpwr( : ), hpwi( : )
  real( kind = kind( 1.0d0 ) ) :: pi
  real( kind = kind( 1.0d0 ) ) :: rcut, rzero, grandn, nval, eps, pref, celvol
  real( kind = kind( 1.0d0 ) ) :: ptab( 100 ), amet( 3, 3 ), tau( 3 )
  real( kind = kind( 1.0d0 ) ), allocatable :: v( :, :, : )
  integer :: npmax, itot, jtot
  integer, allocatable :: ibeg( : ), jbeg( : ), mham( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: mpcr( :, :, :, : ), pcr( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: mpci( :, :, :, : ), pci( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ampr( : ), ampi( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: hampr( : ), hampi( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: mpm( :, :, : ), list( : )
  character * 5 :: s5
  real( kind = kind( 1.0d0 ) ), allocatable :: somel( :, :, : )
  !
  integer :: ncore, lc, atno
  character * 4 :: add04
  character * 10 :: add10
  character * 11 :: s11
