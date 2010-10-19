subroutine getphase( ng, gvec, q, tau, tauphs )
  implicit none
  !
  integer :: ng
  integer :: gvec( 3, ng )
  double precision :: q( 3 ), tau( 3 )
  double complex :: tauphs( ng )
  !
  integer :: ig
  double precision :: gphase, qphase, phase, pi
  double complex :: rm1
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  pi = 4.0d0 * atan( 1.0d0 )
  !
  !  there is no phase attached to a core state or projector.
  !  dot tau into gvec+q.
  qphase = 2.0d0 * pi * sum( tau( : ) * q( : ) )
  do ig = 1, ng
     gphase = 2.0d0 * pi * sum( tau( : ) * gvec( :, ig ) )
     phase = gphase + qphase
     tauphs( ig ) = cos( phase ) + rm1 * sin( phase )
  end do
  !
  return
end subroutine getphase
