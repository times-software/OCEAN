  subroutine process_wf ( xmesh, num_bands_val, num_bands_con   )
  implicit none
  !
  integer, intent(in) :: xmesh(3), num_bands_val, num_bands_con
  complex(kind=kind(1.d0)), intent(out) :: wf_xmesh_kslice( xmesh(3), xmesh(2), xmesh(1), &
      num_bands_val + num_bands_con )
  complex(kind=kind(1.d0)), intent(out) :: tmels_kslice( 
