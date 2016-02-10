! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_action
  use AI_kinds
  use OCEAN_timekeeper
  use iso_c_binding
  implicit none
  private
  save  


  REAL(DP), ALLOCATABLE :: a( : )
  REAL(DP), ALLOCATABLE :: b( : )
  REAL(DP), ALLOCATABLE :: real_a( : )
  REAL(DP), ALLOCATABLE :: imag_a( : )
  REAL(DP), ALLOCATABLE :: real_b( : )
  REAL(DP), ALLOCATABLE :: imag_b( : )
 
  REAL(DP) :: inter_scale_threshold = 0.00001
  REAL(DP) :: inter_scale

  REAL(DP) :: el, eh, gam0, eps, nval,  ebase
  REAL(DP) :: gres, gprc, ffff, ener
  REAL(DP) :: e_start, e_stop, e_step
  REAL(DP), ALLOCATABLE :: e_list( : )

  
  INTEGER  :: haydock_niter = 0
  INTEGER  :: ne
  INTEGER  :: nloop
  INTEGER  :: inv_loop
  

  CHARACTER(LEN=3) :: calc_type
  LOGICAL  :: echamp
  LOGICAL  :: is_first = .true.

  public :: OCEAN_hayinit, OCEAN_action_run

  contains


  subroutine OCEAN_hay_dealloc( psi, old_psi, new_psi, mul_psi, lr_psi, ierr )
    use OCEAN_psi
    implicit none
    integer, intent( inout ) :: ierr
    type( ocean_vector ), intent( inout ) :: psi, old_psi, new_psi, mul_psi, lr_psi

    call OCEAN_psi_kill( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( old_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( mul_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( lr_psi, ierr )
    if( ierr .ne. 0 ) return
 
  end subroutine OCEAN_hay_dealloc

  subroutine OCEAN_hay_alloc( sys, hay_vec, psi, old_psi, new_psi, mul_psi, lr_psi, ierr )
    use OCEAN_system
    use OCEAN_psi

    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    type( ocean_vector ), intent( out ) :: psi, old_psi, new_psi, mul_psi, lr_psi


    call OCEAN_psi_new( psi, ierr, hay_vec )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( new_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( mul_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_new( lr_psi, ierr )
    if( ierr .ne. 0 ) return

  end subroutine OCEAN_hay_alloc


  subroutine OCEAN_action_run( sys, hay_vec, lr, ierr )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_long_range
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    type(long_range), intent( inout ) :: lr


    select case ( calc_type )
      case('hay')
        call OCEAN_haydock( sys, hay_vec, lr, ierr )
      case('inv')
        call OCEAN_GMRES( sys, hay_vec, lr, ierr )
    end select
  end subroutine OCEAN_action_run


  subroutine OCEAN_haydock( sys, hay_vec, lr, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    


    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    type(long_range), intent( inout ) :: lr

    real(DP) :: imag_a
    integer :: iter

    character( LEN=21 ) :: lanc_filename

    type( ocean_vector ), target  :: psi1, psi2, psi3
    type( ocean_vector ), pointer :: psi, old_psi, new_psi, tmp_psi
    

    call OCEAN_psi_new( psi1, ierr, hay_vec )
    if( ierr .ne. 0 ) return
    psi => psi1

    call OCEAN_psi_new( psi2, ierr )
    if( ierr .ne. 0 ) return
    old_psi => psi2

    call OCEAN_psi_new( psi3, ierr )
    if( ierr .ne. 0 ) return
    new_psi => psi3

    if( myid .eq. root ) then 
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      write(6,*) inter_scale, haydock_niter
    endif





    do iter = 1, haydock_niter

      ! This can be packaged into act later
      call OCEAN_psi_zero_full( new_psi, ierr )
      if( ierr .ne. 0 ) return

      call OCEAN_psi_ready_buffer( new_psi, ierr )
      if( ierr .ne. 0 ) return
      
      ! assert that psi is good here?

      if( sys%e0 .and. sys%cur_run%have_core) then 
        call OCEAN_tk_start( tk_e0 )
        call ocean_energies_act( sys, psi, new_psi, ierr )
        call OCEAN_tk_stop( tk_e0 )
      endif

      if( sys%mult .and. sys%cur_run%have_core ) then 
        call OCEAN_tk_start( tk_mult )
        call OCEAN_mult_act( sys, inter_scale, psi, new_psi )
        call OCEAN_tk_stop( tk_mult )
      endif

      if( sys%long_range .and. sys%cur_run%have_core ) then
        call OCEAN_tk_start( tk_lr )
        call lr_act( sys, psi, new_psi, lr, ierr )
        call OCEAN_tk_stop( tk_lr )
      endif

      ! end

      


      ! This should be hoisted back up here
      call ocean_hay_ab( sys, psi, hpsi, old_psi, iter, ierr )

!      ! Shuffle around. This round's old will be written to next round as hpsi
!      temp_psi => old_psi
!      old_psi => psi
!      psi => hpsi
!      hpsi => temp_psi

    enddo

    if( myid .eq. 0 ) then
      write(lanc_filename, '(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'lanceig_', sys%cur_run%elname, &
        '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      call haydump( haydock_niter, sys, ierr )
      call redtrid(  haydock_niter, sys, ierr )
    endif

  end subroutine OCEAN_haydock


  subroutine OCEAN_GMRES( sys, hay_vec, lr, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range


    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    type(long_range), intent( inout ) :: lr


    type( ocean_vector ) :: long_range_psi
    type( ocean_vector ) :: multiplet_psi
    type( ocean_vector ) :: psi1
    type( ocean_vector ) :: psi2
    type( ocean_vector ) :: psi3

    type( ocean_vector ), pointer :: hpsi => null()
    type( ocean_vector ), pointer :: psi => null()
    type( ocean_vector ), pointer :: old_psi => null()
    type( ocean_vector ), pointer :: temp_psi => null()

!    type( ocean_vector ) :: prec_psi

    character( LEN=21 ) :: lanc_filename
    character( LEN=3 ) :: technique, req, bs, as
    character( LEN = 9 ) :: ct
    character( LEN=5) :: eval

    character( LEN = 21 ) :: abs_filename
!    character( LEN = 17 ) :: rhs_filename
    character( LEN = 25 ) :: e_filename

    complex( DP ), allocatable, dimension ( : ) :: x, rhs, v1, v2, pcdiv, cwrk

    integer :: i, ntot, iter, iwrk, need, int1, int2
    real( DP ) :: relative_error, f( 2 ), ener
    complex( DP ) :: rm1

    return
  end subroutine OCEAN_GMRES
#ifdef FALSE
    rm1 = -1
    rm1 = sqrt( rm1 )
    

    iwrk = 1
    eval = 'zerox'
    f( 1 ) = ffff

    ! for error checking
    allocate( cwrk( 1 ) )

    call ocean_hay_alloc( sys, hay_vec, psi1, psi2, psi3, multiplet_psi, &
                          long_range_psi, ierr )


    if( myid .eq. root ) then
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      write(6,*) inter_scale, haydock_niter


      select case ( sys%cur_run%calc_type)
      case( 'XES' )
        write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xesspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      case( 'XAS' )
        write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      case default
        write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      end select

      open( unit=76,file=abs_filename,form='formatted',status='unknown' )
      rewind( 76 )
    endif


    ntot = sys%nalpha * sys%nkpts * sys%num_bands
    allocate( rhs( ntot ), v1( ntot ), v2( ntot ), pcdiv( ntot ), x( ntot ) )

    !

    call vtor( sys, hay_vec, rhs )
!    write(rhs_filename,'(A4,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'rhs_', sys%cur_run%elname, &
!            '.', sys%cur_run%indx, '_', '1s', '_', sys%cur_run%photon
!    open(unit=99,file=rhs_filename,form='unformatted',status='unknown')
!    rewind( 99 )
!    write( 99 ) rhs
!    close( 99 )


!    call OCEAN_tk_init()
    do iter = 1, inv_loop
!      ener = ( e_start + ( iter - 1 ) * e_step ) / 27.2114_DP
      ener = e_list( iter ) / 27.2114_DP
      if( myid .eq. root ) write(6,*) ener * 27.2114_DP

!      call OCEAN_action_set_psi( psi )      


      psi%r( :, :, : ) = 1.0_DP
      psi%i( :, :, : ) = 0.0_DP
      call OCEAN_xact( sys, psi, hpsi, multiplet_psi, long_range_psi, lr, ierr )
      ! After OCEAN_xact every proc has the same copy of hpsi (and should still have the same psi)

      call OCEAN_tk_start( tk_inv )


!      call OCEAN_psi_set_prec( sys, ener, gprc, hpsi, prec_psi )
      call vtor( sys, hpsi, v1 )
      do i = 1, ntot
        pcdiv( i ) = ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )
      end do
      ct = 'beginning'
      req = '---'

      do while ( req .ne. 'end' )
        call OCEAN_invdrv( x, rhs, ntot, int1, int2, nloop, need, iwrk, cwrk, v1, v2, bs, as, &
                     req, ct, eval, f )
        select case( req )
        case( 'all ' )
          if( allocated( cwrk ) ) deallocate( cwrk )
          iwrk = need
          allocate( cwrk( need ) )
        case( 'act' ) ! E - S H ... in what follows, v1 must be untouched
          ! v = v1
          call rtov( sys, psi, v1 )
          call OCEAN_tk_stop( tk_inv )
          call OCEAN_xact( sys, psi, hpsi, multiplet_psi, long_range_psi, lr, ierr )
          call OCEAN_tk_start( tk_inv )
          call vtor( sys, hpsi, v2 )
          v2( : ) = ( ener + rm1 * gres ) * v1( : ) - v2( : )
        case( 'prc' )  ! meaning, divide by S(E-H0) ... in what follows, v1 must be untouched
          v2( : ) = v1 ( : ) * pcdiv( : )
          if( myid .eq. root ) then
            write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) int1, int2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( rhs, x )
!             write ( 66, '(1p,2x,3i5,5(1x,1e15.8))' ) int1, int2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( rhs, x )
          endif
        end select
      enddo
      call OCEAN_tk_stop( tk_inv )



      if( myid .eq. 0 ) then
        relative_error = f( 2 ) / ( dimag( - dot_product( rhs, x ) ) ) !* kpref )
        write ( 76, '(1p,1i5,4(1x,1e15.8))' ) int1, ener*27.2114_DP, &
                  ( 1.0d0 - dot_product( rhs, x ) ) * kpref, relative_error
        call flush(76)

        if( echamp ) then
          write(e_filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', iter
          open(unit=99,file=e_filename,form='unformatted',status='unknown')
          rewind( 99 )
          write( 99 ) x
          close( 99 )
        endif
      endif

      write(e_filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'exciton', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', iter
!      call rtov( sys, psi, v1 )
!      call rtov( sys, psi, x )
!      call dump_exciton( sys, psi, e_filename, ierr )

    enddo

    deallocate( rhs, v1, v2, pcdiv, x )

    if( myid .eq. root ) close( 76 )

    call OCEAN_hay_dealloc( psi1, psi2, psi3, multiplet_psi, long_range_psi, ierr )

!    if( myid .eq. root ) call OCEAN_tk_printtimes
#ifdef MPI
    write(6,*) 'OCEAN_GMRES:', myid, root
    call MPI_BARRIER( comm, ierr )
#endif


  end subroutine OCEAN_GMRES

  subroutine vtor( sys, psi, vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(in) :: psi
    complex( DP ), intent ( out ) :: vec( sys%nalpha * sys%nkpts * sys%num_bands )
    !
    integer :: ia, ik, ib, ii
    complex(DP) :: rm1

    rm1 = -1
    rm1 = sqrt(rm1)
    ii = 0 
    do ia = 1, sys%nalpha
      do ik = 1, sys%nkpts
        do ib = 1, sys%num_bands
          ii = ii + 1
!          vec( ii ) = cmplx( psi%r( ib, ik ,ia ), psi%i( ib, ik ,ia ) )
          vec( ii ) = psi%r( ib, ik ,ia ) + rm1 * psi%i( ib, ik ,ia )
        enddo
      enddo
    enddo
  end subroutine vtor
        

  subroutine rtov( sys, psi, vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi
    complex( DP ), intent ( in ) :: vec( sys%nalpha * sys%nkpts * sys%num_bands )
    ! 
    integer :: ia, ik, ib, ii
    complex(DP) :: rm1

    rm1 = -1
    rm1 = sqrt(rm1)

    ii = 0
    do ia = 1, sys%nalpha
      do ik = 1, sys%nkpts
        do ib = 1, sys%num_bands
          ii = ii + 1
          psi%r( ib, ik, ia ) = vec( ii )
          psi%i( ib, ik, ia ) = -rm1 * vec( ii )
!          psi%r( ib, ik ,ia ) = real( vec( ii ) )
!          psi%i( ib, ik ,ia ) = aimag( vec( ii ) )
        enddo
      enddo
    enddo
  end subroutine rtov
#endif



#ifdef FALSE
  subroutine OCEAN_xact( sys, psi, hpsi, multiplet_psi, long_range_psi, lr, ierr )
    use AI_kinds 
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: hpsi,  multiplet_psi, long_range_psi
    type(OCEAN_vector), intent( in ) :: psi
    type(long_range), intent( inout ) :: lr

    !
    if( sys%long_range ) then
      call OCEAN_tk_start( tk_lr )
      call lr_act( sys, psi, long_range_psi, lr, ierr )
!      if( nproc .gt. 1 ) then
!        call ocean_psi_sum_lr( sys, long_range_psi, ierr )
!      endif
      call OCEAN_tk_stop( tk_lr )
    endif

    if( sys%mult ) then
      call OCEAN_tk_start( tk_mult )
      call OCEAN_mult_act( sys, inter_scale, psi, multiplet_psi )
      call OCEAN_tk_stop( tk_mult )
    endif
    if( sys%e0 ) then 
      call OCEAN_tk_start( tk_e0 )
      call ocean_energies_act( sys, psi, hpsi, ierr )
      call OCEAN_tk_stop( tk_e0 )
    endif
    call OCEAN_tk_start( tk_psisum )
    call ocean_psi_sum( hpsi, multiplet_psi, long_range_psi, ierr )
    call OCEAN_tk_stop( tk_psisum )

  end subroutine
#endif

  subroutine OCEAN_hay_ab( sys, psi, hpsi, old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi

    complex(DP) :: ctmp
    real(dp) :: rtmp
    integer :: ialpha, ikpt, crequest

    ! calc ctmp = < hpsi | psi > and begin Iallreduce
    call OCEAN_psi_dot( ctmp, crequest, hpsi, psi )

    ! hpsi -= b(i-1) * psi^{i-1}
    rtmp = -b(iter-1)
    call OCEAN_psi_axpy( rtmp, old_psi, hpsi )

    ! finish allreduce to get ctmp
    call MPI_WAIT( crequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    real_a(iter-1) = dble( ctmp )
    imag_a(iter-1) = aimag( ctmp )
    rtmp = -real_a( iter - 1 )
    call OCEAN_psi_axpy( rtmp, psi, hpsi )

    call OCEAN_psi_nrm( b(iter), hpsi, ierr, crequest )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( crequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    rtmp = 1.0_dp / b(iter)
    call OCEAN_psi_scal( rtmp, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_store2full( psi, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. 0 ) then
!      write ( 6, '(2x,2f10.6,10x,1e11.4,x,f6.3)' ) a(iter-1), b(iter), imag_a, time2-time1
!      write ( 6, '(2x,2f10.6,10x,1e11.4,8x,i6)' ) a(iter-1), b(iter), imag_a, iter
      write ( 6, '(2x,2f20.6,10x,1e11.4,8x,i6)' ) real_a(iter-1), b(iter), imag_a(iter-1), iter
      if( mod( iter, 10 ) .eq. 0 ) call haydump( iter, sys, ierr )
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
      endif

!      call haydump( iter, sys, ierr )
    endif

  end subroutine OCEAN_hay_ab


  subroutine haydump( iter, sys, ierr )
    use OCEAN_psi,  only : kpref
    use OCEAN_system
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter

    integer :: ie, jdamp, jj
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk, pi
    complex(DP) :: rm1, ctmp, disc, delta

    character( LEN=21 ) :: abs_filename
    
    select case ( sys%cur_run%calc_type)
    case( 'XES' )
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xesspct_', sys%cur_run%elname, &
          '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    case( 'XAS' )
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
          '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    case default
      write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
          '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    end select
    
    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
!    open( unit=99, file='absspct', form='formatted', status='unknown' )
    open( unit=99, file=abs_filename, form='formatted', status='unknown' )
    rewind 99
    do ie = 1, 2 * ne, 2
       e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
       do jdamp = 0, 1
          gam= gam0 + gamfcn( e, nval, eps ) * dble( jdamp )
!          ctmp = e - a( iter - 1 ) + rm1 * gam
          ctmp = e - real_a( iter - 1 ) + rm1 * gam 
          disc = sqrt( ctmp ** 2 - 4 * b( iter ) ** 2 )
          di= -rm1 * disc
          if ( di .gt. 0.0d0 ) then
             delta = ( ctmp + disc ) / 2
          else
             delta = ( ctmp - disc ) / 2
          end if
          do jj = iter - 1, 0, -1
!             delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
             delta = e - real_a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
          end do
          dr = delta
          di = -rm1 * delta
          di = abs( di )
          ener = ebase + 27.2114d0 * e
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
       end do
       spkk = kpref * dr / ( dr ** 2 + di ** 2 )
       write ( 99, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
    end do
    close(unit=99)
    !
    return
  end subroutine haydump

  subroutine OCEAN_hayinit( ierr )
    use OCEAN_mpi
    implicit none

    integer, intent( inout ) :: ierr

    integer :: dumi, iter
    character(len=4) :: inv_style
    real :: dumf

    if( .not. is_first ) goto 10
    is_first = .false.

    if( myid .eq. root ) then
      open(unit=99,file='mode',form='formatted',status='old')
      rewind(99)
      read(99,*) inter_scale, haydock_niter
      close(99)

!      open(unit=99,file='calc_control',form='formatted',status='old')
      open(unit=99,file='bse.in',form='formatted',status='old')
      rewind(99)
      read(99,*) dumi
      read(99,*) dumf
      read(99,*) calc_type
      select case ( calc_type )
        case('hay')
          read(99,*) haydock_niter, ne, el, eh, gam0, ebase
          el = el / 27.2114d0
          eh = eh / 27.2114d0
          gam0 = gam0 / 27.2114d0
          inv_loop = 1
          allocate( e_list( inv_loop ) )
        case('inv')
          read(99,*) nloop, gres, gprc, ffff, ener
          read(99,*) inv_style
          select case( inv_style )
            case('list')
              read(99,*) inv_loop
              allocate( e_list( inv_loop ) )
              do iter = 1, inv_loop
                read(99,*) e_list(iter)
              enddo
            case('loop')
              read(99,*) e_start, e_stop, e_step
              inv_loop = floor( ( e_stop - e_start + e_step*.9) / e_step )
              if (inv_loop .lt. 1 ) inv_loop = 1
              allocate( e_list( inv_loop ) )
              do iter = 1, inv_loop
                e_list( iter ) = e_start + ( iter - 1 ) * e_step
              enddo
          end select          

          inquire(file='echamp.inp',exist=echamp)
          if( echamp ) then
            open(98,file='echamp.inp',form='formatted',status='old')
            read(98,*) echamp
            close(98)
          endif
        case default
          ierr = -1
      end select
      close(99)


      open(unit=99,file='epsilon',form='formatted',status='old')
      rewind 99
      read(99,*) eps
      close(99)

      open( unit=99, file='nval.h', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) nval
      close( unit=99 )
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
    call MPI_BCAST( eps, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( nval, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )


    call MPI_BCAST( e_start, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( e_stop, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( e_step, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( inv_loop, 1, MPI_INTEGER, root, comm, ierr )
    if( inv_loop .gt. 0 ) then
      if( myid .ne. root ) allocate( e_list( inv_loop ) )
      call MPI_BCAST( e_list, inv_loop, MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif
#endif

10 continue

    if( allocated( a ) ) deallocate( a )
    if( allocated( b ) ) deallocate( b )
    if( allocated( real_a ) ) deallocate( real_a )
    if( allocated( imag_a ) ) deallocate( imag_a )
    if( allocated( real_b ) ) deallocate( real_b )
    if( allocated( imag_b ) ) deallocate( imag_b )
    if( haydock_niter .gt. 0 ) then
      allocate( a( 0 : haydock_niter ) )
      allocate( b( 0 : haydock_niter ) )
      allocate( real_a( 0 : haydock_niter ) )
      allocate( imag_a( 0 : haydock_niter ) )
      allocate( real_b( 0 : haydock_niter ) )
      allocate( imag_b( 0 : haydock_niter ) )
      a(:) = 0.0_DP
      b(:) = 0.0_DP
      real_a(:) = 0.0_DP
      imag_a(:) = 0.0_DP
      real_b(:) = 0.0_DP
      imag_b(:) = 0.0_DP
    else
      allocate( a(1), b(1) )
    endif

  end subroutine OCEAN_hayinit

  subroutine redtrid(n,sys, ierr)
    use OCEAN_psi,  only : kpref
    use OCEAN_system
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) ::  n
!    real( DP ), intent( in ) :: a( 0 : n ), b( n )
!    character( LEN=21 ), intent( in ) :: lanc_filename
    double precision, allocatable :: ar(:,:),ai(:,:)
    double precision, allocatable :: w(:),zr(:,:),zi(:,:)
    double precision, allocatable :: fv1(:),fv2(:),fm1(:)
    integer :: matz,nm,i,j,nn


    character( LEN=21 ) :: lanc_filename

    select case ( sys%cur_run%calc_type)
    case( 'XES' )
      write(lanc_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xeslanc_', sys%cur_run%elname, &
          '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    case( 'XAS' )
      write(lanc_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'abslanc_', sys%cur_run%elname, &
          '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    case default
      write(lanc_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'abslanc_', sys%cur_run%elname, &
          '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    end select



    nn=n+1
    nm=nn+10
    allocate(ar(nm,nm),ai(nm,nm),w(nm),zr(nm,nm),zi(nm,nm))
    allocate(fv1(nm),fv2(nm),fm1(2*nm))
    do i=1,n+1
      do j=1,n+1
        ar(j,i)=0.d0
        ai(j,i)=0.d0
      end do
    end do
    do i=1,n+1
      ar(i,i)=real_a(i-1)
      if (i.le.n) then
        ar(i+1,i)=b(i)
        ar(i,i+1)=b(i)
      end if
    end do
    matz=0
    call elsch(nm,nn,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
!    open(unit=99,file='lanceigs',form='formatted',status='unknown')
    open(unit=99,file=lanc_filename,form='formatted',status='unknown')
    rewind 99
    write ( 99, '(1i5,1e26.15)' ) n, kpref
    do i = 0, n
      if ( i .eq. 0 ) then
        write ( 99, '(2x,1f20.10)' ) real_a( i )
      else
        write ( 99, '(2x,2f20.10)' ) real_a( i ), b( i )
      end if
    end do
    write (99,'(2x,2i5,1f20.10)') (i,nn,w(i),i=1,nn)
    close(unit=99)
    deallocate(ar,ai,w,zr,zi,fv1,fv2,fm1)
    return
  end subroutine redtrid




end module OCEAN_action
