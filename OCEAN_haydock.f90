module OCEAN_action
  use AI_kinds
  implicit none
  private
  save  


  REAL(DP), ALLOCATABLE :: a( : )
  REAL(DP), ALLOCATABLE :: b( : )
 
  REAL(DP) :: inter_scale_threshold = 0.00001
  REAL(DP) :: inter_scale

  REAL(DP) :: el, eh, gam0, eps, nval, eps, ebase
  REAL(DP) :: gres, gprc, fff, ener
  
  INTEGER  :: haydock_niter = 0
  INTEGER  :: ne
  INTEGER  :: nloop
  

  CHARACTER, LEN=3 :: calc_type

  public :: OCEAN_haydock, OCEAN_hayinit

  contains



  



  subroutine OCEAN_haydock( sys, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_long_range
    


    implicit none
    integer, intent( inout ) :: ierr
    type( ocean_system ), intent( in ) :: sys

    integer :: iter
    integer :: num_threads

!$    integer, external :: omp_get_num_threads

    num_threads = 1
!$  num_threads = omp_get_num_threads()
    num_threads = min( 2, num_threads )
    

    do iter = 1, haydock_niter

      if( sys%long_range ) then

        if( sys%obf ) then
          call lr_act_obf( sys, ierr )
        else
          call lr_act( sys, ierr )
        endif

        if( nproc .gt. 1 ) then
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS num_threads

!$OMP SECTIONS

!$OMP SECTION
          call ocean_psi_sum_lr( sys, ierr )

!$OMP SECTION
          if( sys%mult ) call mult_act( sys, ierr )


!$OMP END SECTIONS

!$OMP END PARALLEL
        else
          if( sys%mult ) call mult_act( sys, ierr )
        endif
      else
        if( sys%mult ) call mult_act( sys, ierr )
      endif

      if( sys%e0 ) call energies_act( sys, ierr )

      call ocean_psi_sum( sys, ierr )

      call ocean_psi_ab( sys, a(iter-1), b(iter-1), b(iter), ierr )
!      call ocean_psi_dump( sys, ierr )
      if( myid .eq. 0 ) then
        write ( 6, '(2x,2f10.6,10x,1e11.4)' ) real(a(iter-1)), b(iter), aimag(a(iter-1))
        call haydump( iter, ierr )
      endif

    enddo
    
  end subroutine OCEAN_haydock

  subroutine haydump( iter, ierr )
    use module psi,  only : kpref
    implicit none
    integer, intent( inout ) :: ierr
    integer, intent( in ) :: iter

    integer :: ie, jdamp, jj
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk, pi
    complex(DP) :: rm1, ctmp, disc, delta
    
    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
    open( unit=99, file='absspct', form='formatted', status='unknown' )
    rewind 99
    do ie = 1, 2 * ne, 2
       e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
       do jdamp = 0, 1
          gam= gam0 + gamfcn( e, nval, eps ) * dble( jdamp )
          ctmp = e - a( iter - 1 ) + rm1 * gam
          disc = sqrt( ctmp ** 2 - 4 * b( iter + 1 ) ** 2 )
          di= -rm1 * disc
          if ( di .gt. 0.0d0 ) then
             delta = ( ctmp + disc ) / 2
          else
             delta = ( ctmp - disc ) / 2
          end if
          do jj = iter - 1, 0, -1
             delta = e - a( jj ) + rm1 * gam - b( jj + 1 + 1 ) ** 2 / delta
          end do
          dr = delta
          di = -rm1 * delta
          di = abs( di )
          ener = ebase + 27.2114d0 * e
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
       end do
       spkk = kpref * dr / ( dr ** 2 + di ** 2 )
       write ( 99, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, j, gam, kpref, ne
    end do
    close(unit=99)
    !
    return
  end subroutine haydump

  subroutine OCEAN_hayinit( ierr )
    use OCEAN_mpi
    implicit none

    integer, intent( inout ) :: ierr

    if( myid .eq. root ) then
      open(unit=99,file='mode',form='formatted',status='old')
      rewind(99)
      read(99,*) inter_scale, haydock_niter
      close(99)

      open(unit=99,file='calc_control',form='formatted',status='old')
      rewind(99)
      read(99,*) calc_type
      select case ( calc_type )
        case('hay')
          read(99,*) ne, el, eh, gam0, ebase
        case('inv')
          read(99,*) nloop, gres, gprc, ffff, ener
        case default
          ierr = -1
      end select
      close(99)

    endif

#ifdef MPI
    call MPI_BCAST( inter_scale, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( haydock_niter, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( calc_type, 3, MPI_CHARACTER, root, comm, ierr )

    call MPI_BCAST( nloop, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( gres, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( gprc, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( ffff, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( ener, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

    if( haydock_niter .gt. 0 ) then
      allocate( a( 0 : haydock_niter ) )
      allocate( b( 0 : haydock_niter ) )
      a(:) = 0_DP
      b(:) = 0_DP
    else
      allocate( a(1), b(1) )
    endif

  end subroutine OCEAN_hayinit


end OCEAN_action
