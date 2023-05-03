! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_haydock
  use AI_kinds
  use OCEAN_timekeeper
  use iso_c_binding
  implicit none
  private
  save  


!  REAL(DP), ALLOCATABLE :: a( : )
!  REAL(DP), ALLOCATABLE :: b( : )
  REAL(DP), ALLOCATABLE :: real_a( : )
  REAL(DP), ALLOCATABLE :: imag_a( : )
  REAL(DP), ALLOCATABLE :: real_b( : )
  REAL(DP), ALLOCATABLE :: imag_b( : )
  REAL(DP), ALLOCATABLE :: real_c( : )
  REAL(DP), ALLOCATABLE :: imag_c( : )
 

  REAL(DP) :: el, eh, gam0, eps, nval,  ebase
  REAL(DP) :: gres, gprc, ffff, ener
  REAL(DP) :: e_start, e_stop, e_step
  REAL(DP), ALLOCATABLE :: e_list( : )

  real(DP) :: eps1Conv( 3 )

  
  INTEGER  :: haydock_niter = 0
  INTEGER  :: ne
  INTEGER  :: nloop
  INTEGER  :: inv_loop
  

  CHARACTER(LEN=3) :: calc_type
  LOGICAL  :: echamp
  LOGICAL  :: project_absspct
  LOGICAL  :: is_first = .true.

  LOGICAL  :: val_loud = .true.
  LOGICAL  :: complex_haydock = .false.

  public :: OCEAN_haydock_setup, OCEAN_haydock_do

  contains

  subroutine OCEAN_haydock_nonHerm_do( sys, hay_vec, ierr )
    use AI_kinds, only : DP
    use OCEAN_energies
    use OCEAN_system, only : o_system
    use OCEAN_psi
    use OCEAN_action, only : OCEAN_xact
    use OCEAN_mpi, only : myid, root, comm

    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec

    real(DP) :: imag_a
    integer :: iter
    type( ocean_vector ) :: psi, old_psi, new_psi
    type( ocean_vector ) :: back_psi, back_old_psi, back_new_psi


    ! Initialization steps
    call OCEAN_psi_new( psi, ierr, hay_vec )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( back_psi, ierr, hay_vec )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( new_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( old_psi, ierr )
    if( ierr .ne. 0 ) return


    call OCEAN_psi_new( back_new_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( back_new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( back_old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( back_old_psi, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) then
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref
      write(6,*) sys%interactionScale, haydock_niter
    endif
    call MPI_BARRIER( comm, ierr )
    !\Initialization



    do iter = 1, haydock_niter
      if( sys%cur_run%have_val ) then
        if( myid .eq. root ) write(6,*)   " iter. no.", iter-1
      endif


      call OCEAN_energies_allow( sys, psi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_allow( sys, back_psi, ierr )
      if( ierr .ne. 0 ) return


      call OCEAN_xact( sys, sys%interactionScale, psi, new_psi, ierr )
      if( ierr .ne. 0 ) return
      ! need the action of the Hermitian conjugate of the Hamiltonian
      !  obviously we are only bothering to do this when H isn't Hermitian
      call OCEAN_xact( sys, sys%interactionScale, back_psi, back_new_psi, ierr, .true. )
      if( ierr .ne. 0 ) return

      ! This should be hoisted back up here
      call haydock_abc_1( sys, psi, new_psi, old_psi, back_psi, back_new_psi, back_old_psi, & 
                          iter, ierr )

    enddo

    call OCEAN_tk_stop( tk_psisum )
    call MPI_BARRIER( comm, ierr )
    if( myid .eq. 0 ) then
      call haydump( haydock_niter, sys, hay_vec%kpref, ierr )
      call redtrid(  haydock_niter, sys, hay_vec%kpref, ierr )
    endif

    call OCEAN_psi_kill( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( old_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( back_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( back_new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( back_old_psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )

  end subroutine OCEAN_haydock_nonHerm_do


  subroutine OCEAN_haydock_do( sys, hay_vec, restartBSE, newEps, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_psi, only : ocean_vector

    integer, intent( inout ) :: ierr
    logical, intent( inout ) :: restartBSE
    real(DP), intent( inout ) :: newEps
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec


    if( complex_haydock ) then
      call OCEAN_haydock_nonHerm_do( sys, hay_vec, ierr )
    else
      call OCEAN_haydock_Herm_do( sys, hay_vec,  restartBSE, newEps, ierr )
    endif

  end subroutine OCEAN_haydock_do


  subroutine OCEAN_haydock_Herm_do( sys, hay_vec, restartBSE, newEps, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_action, only : OCEAN_xact


    implicit none
    integer, intent( inout ) :: ierr
    logical, intent( inout ) :: restartBSE
    real(DP), intent( inout ) :: newEps
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec

    real(DP) :: imag_a, maxDiff, relArea 
    integer :: iter, haydock_niter_actual, iter2

!    character( LEN=21 ) :: lanc_filename

    type( ocean_vector ) :: psi, old_psi, new_psi

    logical :: prevConv
    
    prevConv = .false.

!    if( myid .eq. root ) write(6,*) 'entering haydock'

    call OCEAN_psi_new( psi, ierr, hay_vec )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi'

    call OCEAN_psi_new( new_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'new_psi'

    call OCEAN_psi_new( old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( old_psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'old_psi'

    if( myid .eq. root ) then 
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref
      write(6,*) sys%interactionScale, haydock_niter
    endif
    call MPI_BARRIER( comm, ierr )



!    call OCEAN_tk_start( tk_psisum )

    do iter = 1, haydock_niter
      if( sys%cur_run%have_val ) then
        if( myid .eq. root ) write(6,*)   " iter. no.", iter-1
      endif
!        call OCEAN_energies_allow( sys, psi, ierr )
!        if( ierr .ne. 0 ) return
!      endif

      call OCEAN_xact( sys, sys%interactionScale, psi, new_psi, ierr )
      if( ierr .ne. 0 ) return
!      if( myid .eq. root ) write(6,*) 'Done with ACT'



      ! This should be hoisted back up here
      call ocean_hay_ab_twoterm( sys, psi, new_psi, old_psi, iter, restartBSE, newEps, ierr )
      if( ierr .ne. 0 ) return
      if( restartBSE ) goto 11

      ! test to check convergence
      if( sys%earlyExit .and. iter .gt. sys%haydockConvergeSpacing ) then
        if( myid .eq. root ) then
          iter2 = iter - sys%haydockConvergeSpacing
          call check_convergence( iter, iter2, sys, hay_vec%kpref, maxDiff, relArea )
        endif
        call MPI_BCAST( relArea, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
        
        if( relArea .lt. sys%haydockConvergeThreshold ) then
          if( prevConv ) then
            if( myid .eq. root ) write(6,*) 'Convergence: ', iter, maxDiff, relArea
            haydock_niter_actual = iter
            goto 11
          else
            prevConv = .true.
          endif
        else
          if( myid .eq. root ) write(6,*) 'Not converged', iter, maxDiff, relArea
          prevConv = .false.
        endif
      endif
      haydock_niter_actual = iter
          

    enddo

11  continue

    call OCEAN_tk_stop( tk_psisum )
    call MPI_BARRIER( comm, ierr )

    if( myid .eq. 0 .and. .not. restartBSE ) then
!      write(lanc_filename, '(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'lanceig_', sys%cur_run%elname, &
!        '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      call haydump( haydock_niter_actual, sys, hay_vec%kpref, ierr )
      call redtrid(  haydock_niter_actual, sys, hay_vec%kpref, ierr )
    endif

    call OCEAN_psi_kill( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( new_psi, ierr )
    if( ierr .ne. 0 ) return
    
    call OCEAN_psi_kill( old_psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )

  end subroutine OCEAN_haydock_Herm_do

#if( 0 )
  subroutine OCEAN_GMRES( sys, hay_vec, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_pfy, only : OCEAN_pfy_load, OCEAN_pfy_act
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    use OCEAN_action, only : OCEAN_xact

    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec


    type( ocean_vector ) :: hpsi 
    type( ocean_vector ) :: psi 

!    character( LEN=21 ) :: lanc_filename
    character( LEN=3 ) :: technique, req, bs, as
    character( LEN = 9 ) :: ct
    character( LEN=5) :: eval

    character( LEN = 40 ) :: abs_filename
    character( LEN = 21 ) :: pfy_filename
    character( LEN = 21 ) :: proj_filename
!    character( LEN = 17 ) :: rhs_filename
    character( LEN = 25 ) :: e_filename

    complex( DP ), allocatable, dimension ( : ) :: x, rhs, v1, v2, pcdiv, cwrk

    integer :: i, ntot, iter, iwrk, need, int1, int2
    integer :: project_file_handle, pfy_file_handle
    real( DP ) :: relative_error, f( 2 ), ener, fact
    complex( DP ) :: rm1


    logical :: do_pfy = .false.

!    return
    rm1 = -1
    rm1 = sqrt( rm1 )
    

    iwrk = 1
    eval = 'zerox'
    f( 1 ) = ffff

    ! for error checking
    allocate( cwrk( 1 ) )

    call OCEAN_psi_new( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( hpsi, ierr )
    if( ierr .ne. 0 ) return
  
    if(  sys%cur_run%have_val ) then
      fact = hay_vec%kpref * 2.0_dp * sys%celvol
    else
      fact = hay_vec%kpref
    endif

    if( myid .eq. root ) then
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref
      write(6,*) sys%interactionScale, haydock_niter, inv_loop


      select case ( sys%cur_run%calc_type)
      case( 'XES' )
        write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xesspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      case( 'XAS' )
        write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
        write(pfy_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'pfyspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
        do_pfy = .true.
      case( 'RXS')
        write(abs_filename,'(A8,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxsspct_', sys%cur_run%elname, &
            '.', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol
      case default
        write(abs_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
        write(pfy_filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'pfyspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
        do_pfy = .true.
      end select

      open( unit=76,file=abs_filename,form='formatted',status='unknown' )
      rewind( 76 )

      do_pfy = .false.
      if( do_pfy ) then
        pfy_file_handle = 75
        open(pfy_file_handle,file=pfy_filename,form='formatted',status='unknown' )
        rewind( pfy_file_handle )

        call OCEAN_pfy_load( sys, ierr )
      endif


      ! This all needs to be moved to a module or type when we expand to make it more general
      if( project_absspct ) then
        write(6,*) 'Projections requested!'
!        nproject = 2
!        allocate( project_file_handle )
        write( proj_filename, '(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'project_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
        open( unit=77,file=proj_filename,form='formatted',status='unknown' )
        rewind( 77 )
        project_file_handle = 77

!        write( proj_filename, '(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'spin_dn_', sys%cur_run%elname, &
!            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
!        open( unit=78,file=proj_filename,form='formatted',status='unknown' )
!        rewind( 78 )
!        project_file_handles( 2 ) = 78
      endif
  
    endif


!    ntot = sys%nalpha * sys%nkpts * sys%num_bands
    ntot = OCEAN_psi_size_full( hay_vec )
    allocate( rhs( ntot ), v1( ntot ), v2( ntot ), pcdiv( ntot ), x( ntot ) )

    !

!    call vtor( sys, hay_vec, rhs )
    call OCEAN_psi_vtor( hay_vec, rhs )
!    write(rhs_filename,'(A4,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'rhs_', sys%cur_run%elname, &
!            '.', sys%cur_run%indx, '_', '1s', '_', sys%cur_run%photon
!    open(unit=99,file=rhs_filename,form='unformatted',status='unknown')
!    rewind( 99 )
!    write( 99 ) rhs
!    close( 99 )

!    call OCEAN_tk_init()


    call OCEAN_psi_zero_full( psi, ierr )

    do iter = 1, inv_loop

      ener = e_list( iter ) * eV2Hartree 
      if( myid .eq. root ) write(6,*) ener * Hartree2eV  



      ! After OCEAN_xact every proc has the same copy of hpsi (and should still have the same psi)
      call OCEAN_psi_one_full( psi, ierr )

      if( sys%cur_run%have_val ) then
        call OCEAN_energies_allow( sys, psi, ierr )
        if( ierr .ne. 0 ) return
      endif


      call OCEAN_xact( sys, sys%interactionScale, psi, hpsi, ierr )
      call OCEAN_psi_prep_min2full( hpsi, ierr )
      call OCEAN_psi_start_min2full( hpsi, ierr )
      call OCEAN_psi_finish_min2full( hpsi, ierr )
      ! After OCEAN_xact every proc has the same copy of hpsi (and should still have the same psi)

      call OCEAN_tk_start( tk_inv )


!      call OCEAN_psi_set_prec( sys, ener, gprc, hpsi, prec_psi )
!      call vtor( sys, hpsi, v1 )
      call OCEAN_psi_vtor( hpsi, v1 )
      do i = 1, ntot
        pcdiv( i ) = ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )
      end do
      ct = 'beginning'
      req = '---'

      if( iter .gt. 1 ) then
        if( abs( e_list( iter ) - e_list( iter - 1 ) ) * eV2Hartree .lt. 3.0_dp * gres ) then
          if( myid .eq. 0 ) write( 6,* ) '    Re-use previous x'
          eval = 'havex'
        else
          eval = 'zerox'
        endif
!        v1(:) = x(:)
      endif

      do while ( req .ne. 'end' )

        if( myid .eq. root ) then
          call OCEAN_invdrv( x, rhs, ntot, int1, int2, nloop, need, iwrk, cwrk, v1, v2, bs, as, &
                       req, ct, eval, f )
        endif
        call MPI_BCAST( req, 3, MPI_CHARACTER, root, psi%core_comm, ierr )
        if( ierr .ne. MPI_SUCCESS ) return

        select case( req )
        case( 'all ' )
          if( myid .eq. root ) then
            if( allocated( cwrk ) ) deallocate( cwrk )
            iwrk = need
            allocate( cwrk( need ) )
          endif
        case( 'act' ) ! E - S H ... in what follows, v1 must be untouched
          ! v = v1
!          call rtov( sys, psi, v1 )

          if( myid .eq. root ) then
            call OCEAN_psi_rtov( psi, v1 )
          endif
          call OCEAN_psi_bcast_full( root, psi, ierr )
          !
          call OCEAN_tk_stop( tk_inv )
!          call OCEAN_xact( sys, psi, hpsi, multiplet_psi, long_range_psi, ierr )
          if( sys%cur_run%have_val ) then
            call OCEAN_energies_allow( sys, psi, ierr )
            if( ierr .ne. 0 ) return
          endif

          call OCEAN_xact( sys, sys%interactionScale, psi, hpsi, ierr )
          call OCEAN_psi_prep_min2full( hpsi, ierr )
          call OCEAN_psi_start_min2full( hpsi, ierr )
          call OCEAN_psi_finish_min2full( hpsi, ierr )

          call OCEAN_tk_start( tk_inv )
!          call vtor( sys, hpsi, v2 )
  
          if( myid .eq. root ) then
            call OCEAN_psi_vtor( hpsi, v2 )
            v2( : ) = ( ener + rm1 * gres ) * v1( : ) - v2( : )
          endif
        case( 'prc' )  ! meaning, divide by S(E-H0) ... in what follows, v1 must be untouched

          if( myid .eq. root ) then
            v2( : ) = v1 ( : ) * pcdiv( : )
!            if( myid .eq. root ) then
              write ( 6, '(1p,2x,3i5,5(1x,1e15.8))' ) int1, int2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( rhs, x )
  !             write ( 66, '(1p,2x,3i5,5(1x,1e15.8))' ) int1, int2, nloop, f( 2 ), f( 1 ), ener, 1.0d0 - dot_product( rhs, x )
          endif
        end select
      enddo
      call OCEAN_tk_stop( tk_inv )



      if( myid .eq. root ) then
        relative_error = f( 2 ) / ( aimag( - dot_product( rhs, x ) ) ) !* kpref )
        write ( 76, '(1p,1i5,4(1x,1e15.8))' ) int1, ener * Hartree2eV,  & !*27.2114_DP, &
                  ( 1.0d0 - dot_product( rhs, x ) ) * fact, relative_error
#ifdef __HAVE_F03
        flush(76)
#else
        call flush(76)
#endif

        if( echamp ) then
          write(e_filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', iter
          open(unit=99,file=e_filename,form='unformatted',status='unknown')
          rewind( 99 )
          write( 99 ) x
          close( 99 )
        endif

        if( project_absspct ) then
          call write_projected_absspct( sys, project_file_handle, hay_vec%kpref, ener, ntot, rhs, x, ierr )
          if( ierr .ne. 0 ) return
        endif

      endif

!      call rtov( sys, psi, v1 )
      write(e_filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'exciton', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', iter
!      call dump_exciton( sys, psi, e_filename, ierr )

!      call rtov( sys, hpsi, x )

      if( myid .eq. root ) then
        call OCEAN_psi_rtov( hpsi, x )
      endif

      if( do_pfy ) then
        hpsi%kpref = hay_vec%kpref
        call OCEAN_pfy_act( sys, hpsi, ener, pfy_file_handle, ierr )
      endif

    enddo

    deallocate( rhs, v1, v2, pcdiv, x )

    if( myid .eq. root ) close( 76 )

!    call OCEAN_hay_dealloc( psi1, psi2, psi3, multiplet_psi, long_range_psi, ierr )
    call OCEAN_psi_kill( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( hpsi, ierr )
    if( ierr .ne. 0 ) return

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


  subroutine OCEAN_hay_ab( sys, psi, hpsi, old_psi, iter, restartBSE, newEps, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    logical, intent(inout) :: restartBSE
    real(DP), intent(inout) :: newEps

    real(dp) :: btmp, atmp, aitmp
    integer :: ialpha, ikpt, arequest, airequest, brequest

!    if( sys%cur_run%have_val ) then
!      call OCEAN_energies_allow( sys, hpsi, ierr )
!      if( ierr .ne. 0 ) return
!    endif

    ! calc ctmp = < hpsi | psi > and begin Iallreduce
    call OCEAN_psi_dot( hpsi, psi, arequest, atmp, ierr, airequest, aitmp )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_dot'

    ! hpsi -= b(i-1) * psi^{i-1}
    btmp = -real_b(iter-1)
    ! y:= a*x + y
    call OCEAN_psi_axpy( btmp, old_psi, hpsi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_axpy 1'

    ! finish allreduce to get atmp
    ! we want iatmp too (for output/diagnostics), but that can wait
    call MPI_WAIT( arequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    real_a(iter-1) = atmp
    atmp = -atmp
!    if( myid .eq. root ) write(6,*) 'ab', real_a(iter-1), real_b(iter-1)
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_axpy 2'


!    if( sys%cur_run%have_val ) then
!      call OCEAN_energies_allow( sys, hpsi, ierr )
!      if( ierr .ne. 0 ) return
!    endif

    
    !JTV was checking to see if any different here. 
!    call OCEAN_psi_nrm( btmp, hpsi, ierr, brequest )
    call OCEAN_psi_dot( hpsi, hpsi, brequest, btmp, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_nrm'

    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_copy 1'

    ! Could move prep and copy here for hspi -> psi
    !   just need to include a way to scale full instead of just min

    call MPI_WAIT( brequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    real_b(iter) = sqrt( btmp )
    real_c(iter) = real_b(iter)
    btmp = 1.0_dp / real_b( iter )
    call OCEAN_psi_scal( btmp, hpsi, ierr )
    if( ierr .ne. 0 ) return
    !

    ! copies hpsi onto psi
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( airequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    imag_a(iter-1) = aitmp

    if( myid .eq. 0 ) then
!      write ( 6, '(2x,2f10.6,10x,1e11.4,x,f6.3)' ) a(iter-1), b(iter), imag_a, time2-time1
!      write ( 6, '(2x,2f10.6,10x,1e11.4,8x,i6)' ) a(iter-1), b(iter), imag_a, iter
!      write ( 6, '(2x,2f24.13,10x,1e24.13,8x,i6)' ) real_a(iter-1) * Hartree2eV, real_b(iter) * Hartree2eV, &
!                                                  imag_a(iter-1) * Hartree2eV, iter
      write ( 6, '(1x,6(f20.13,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter

      if( mod( iter, 10 ) .eq. 0 ) then 
        call haydump( iter, sys, psi%kpref, ierr )
      endif
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

      if( sys%convEps .and. sys%cur_run%calc_type .eq. 'VAL' ) then
        call testConvergeEps( iter, sys, psi%kpref, sys%celvol, sys%valence_ham_spin, restartBSE, newEps )
      endif
    endif

    if( sys%convEps ) then
      call MPI_BCAST( restartBSE, 1, MPI_LOGICAL, root, comm, ierr )
      call MPI_BCAST( newEps, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif

    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return
    

  end subroutine OCEAN_hay_ab

  subroutine OCEAN_hay_ab_twoterm( sys, psi, hpsi, old_psi, iter, restartBSE, newEps, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : O_system
    use OCEAN_psi
    use OCEAN_mpi, only : root, myid, comm, &
                          MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_INTEGER, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV

    implicit none
    integer, intent(inout) :: ierr                  
    integer, intent(in) :: iter                     
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi 
    logical, intent(inout) :: restartBSE
    real(DP), intent(inout) :: newEps

    real(dp) :: btmp, atmp, aitmp
    integer :: ialpha, ikpt, arequest, airequest, brequest 
   
    ! hpsi -= b(i-1) * psi^{i-1}
    btmp = -real_b(iter-1)
    ! y:= a*x + y
    call OCEAN_psi_axpy( btmp, old_psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    ! calc ctmp = < hpsi | psi > and begin Iallreduce
    call OCEAN_psi_dot( hpsi, psi, arequest, atmp, ierr, airequest, aitmp )
    if( ierr .ne. 0 ) return

    ! finish allreduce to get atmp
    ! we want iatmp too (for output/diagnostics), but that can wait
    call MPI_WAIT( arequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    real_a(iter-1) = atmp

    ! hpsi -= a(i) * psi^{i}
    atmp = -atmp
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    !
    call OCEAN_psi_dot( hpsi, hpsi, brequest, btmp, ierr )
    if( ierr .ne. 0 ) return

    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( brequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    real_b(iter) = sqrt( btmp )
    real_c(iter) = real_b(iter)
    btmp = 1.0_dp / real_b( iter )
    call OCEAN_psi_scal( btmp, hpsi, ierr )
    if( ierr .ne. 0 ) return

    ! copies hpsi onto psi
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( airequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    imag_a(iter-1) = aitmp

    if( myid .eq. 0 ) then
      write ( 6, '(1x,6(f20.13,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter

      if( mod( iter, 10 ) .eq. 0 ) then 
        call haydump( iter, sys, psi%kpref, ierr )
        if( ierr .ne. 0 ) return
        call write_lanczos( iter, sys, psi%kpref, ierr )
        if( ierr .ne. 0 ) return

        ! need to sync first and maybe need to write out old_psi above where it is (maybe?) 
        ! already distributed?
!        call OCEAN_psi_write( sys, psi, 'psi_', .false., ierr )
      endif
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

      if( sys%convEps .and. sys%cur_run%calc_type .eq. 'VAL' ) then
        call testConvergeEps( iter, sys, psi%kpref, sys%celvol, sys%valence_ham_spin, restartBSE, newEps )
      endif
!      call haydump( iter, sys, ierr )
    endif

    if( sys%convEps ) then
      call MPI_BCAST( restartBSE, 1, MPI_LOGICAL, root, comm, ierr )
      call MPI_BCAST( newEps, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif

    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return


  end subroutine  OCEAN_hay_ab_twoterm

#if 0
! Alternate ordering of orthogonalization as given originially by Chris Paige
  subroutine OCEAN_hay_abc_Paige( sys, psi, hpsi, old_psi, back_psi, back_hpsi, back_old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : o_system
    use OCEAN_psi
    use OCEAN_mpi, only : myid, root, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    type(OCEAN_vector), intent(inout) :: back_psi, back_hpsi, back_old_psi

    complex(dp) :: ctmp
    real(dp) :: rtmp, itmp, atmp, btmp
    integer :: ialpha, ikpt, irequest, rrequest

!    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_allow( sys, back_hpsi, ierr )
      if( ierr .ne. 0 ) return
!    endif

    ! first calculate hpsi = ( hspi - b * old_psi ) 
    !             back_hpsi = ( back_hpsi - c * back_old_psi )

    ! hpsi -= b(i-1) * psi^{i-1}
    ! y:= a*x + y
    atmp = -real_b(iter-1)
    btmp = -imag_b(iter-1)
    call OCEAN_psi_axpy( atmp, old_psi, hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    atmp = -real_c(iter-1)
    btmp = -imag_c(iter-1)
    call OCEAN_psi_axpy( atmp, back_old_psi, back_hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    ! calc ctmp = < hpsi | back_psi > and begin Iallreduce
    call OCEAN_psi_dot( back_psi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
    if( ierr .ne. 0 ) return


    ! finish allreduce to get a
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    atmp = -rtmp
    btmp = -itmp
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    btmp = -itmp
    call OCEAN_psi_axpy( atmp, back_psi, back_hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    real_a(iter-1) = rtmp
    imag_a(iter-1) = itmp

!    if( myid .eq. root ) write(6,*) 'ab', real_a(iter-1), b(iter-1)

!    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_allow( sys, back_hpsi, ierr )
      if( ierr .ne. 0 ) return
!    endif

    call OCEAN_psi_dot( back_hpsi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
    if( ierr .ne. 0 ) return

    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( back_old_psi, back_psi, ierr )
    if( ierr .ne. 0 ) return

    ! Could move prep and copy here for hspi -> psi
    !   just need to include a way to scale full instead of just min

    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    ctmp = sqrt( cmplx( rtmp, itmp, DP ) )
!    ctmp = sqrt( sqrt( rtmp*rtmp + itmp*itmp ) )

    real_c( iter ) = real( ctmp, DP )
    imag_c( iter ) = aimag( ctmp )

    ctmp = cmplx( rtmp, itmp, DP ) / ctmp

    real_b( iter ) = real( ctmp, DP )
    imag_b( iter ) = aimag( ctmp )

   call OCEAN_psi_divide( back_hpsi, ierr, real_b(iter), -imag_b(iter) )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_divide( hpsi, ierr, real_c(iter), imag_c(iter) )
    if( ierr .ne. 0 ) return
    !

    ! copies hpsi onto psi
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( back_psi, back_hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return



    if( myid .eq. 0 ) then
      write ( 6, '(1x,6(f14.8,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter
      if( mod( iter, 10 ) .eq. 0 ) call haydump( iter, sys, psi%kpref, ierr )
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

!      call haydump( iter, sys, ierr )
    endif
    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_finish_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return


  end subroutine OCEAN_hay_abc_Paige
#endif

  subroutine haydock_abc( sys, psi, hpsi, old_psi, back_psi, back_hpsi, back_old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : o_system
    use OCEAN_psi
    use OCEAN_mpi, only : myid, root, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    type(OCEAN_vector), intent(inout) :: back_psi, back_hpsi, back_old_psi

    complex(dp) :: ctmp
    real(dp) :: rtmp, itmp, atmp, btmp
    integer :: irequest, rrequest


    ! Following 
    ! Inner product (x,y) = \sum_{i=1}^{m} x_i \bar{y}_i
    !  This means that (x,y) = \langle y \vert x \rangle

    ! step 0, enforce allow
    call OCEAN_energies_allow( sys, hpsi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_energies_allow( sys, back_hpsi, ierr )
    if( ierr .ne. 0 ) return


    ! step 2:  New vector  ! HERE THERE IS A DIFFERENCE, beta index
    ! step 2A) r = A v_j - \beta_{j-1} v_{j-1}
    atmp = - real_c( iter - 1 )
    btmp = - imag_c( iter - 1 )
    call OCEAN_psi_axpy( atmp, old_psi, hpsi, ierr, btmp )

    ! step 2B) s= A^* u_j - \beta^*_{j-1} u_{j-1}
    atmp = - real_b( iter - 1 )
    btmp =   imag_b( iter - 1 )
    call OCEAN_psi_axpy( atmp, back_old_psi, back_hpsi, ierr, btmp )

    ! $ \alpha_j = u_j^* r
    call OCEAN_psi_dot( back_hpsi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
    ! Now need to make sure alpha is done
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    ! r = r - alpha v_j
    atmp = -rtmp
    btmp = -itmp
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return
    ! s = s - alpha^* u_j
    atmp = -rtmp
    btmp =  itmp
    call OCEAN_psi_axpy( atmp, back_psi, back_hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    real_a(iter-1) = rtmp
    imag_a(iter-1) = itmp
    call OCEAN_psi_dot( hpsi, back_hpsi, rrequest, rtmp, ierr, irequest, itmp )

    ! get ready for next iteration
    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_copy_min( back_old_psi, back_psi, ierr )
    if( ierr .ne. 0 ) return


    ! Wait on 4A
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    !
    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return


    ctmp = cmplx( rtmp, itmp, DP )
    ctmp = sqrt( ctmp )
    real_b(iter) = real(ctmp,DP)
    imag_b(iter) = aimag(ctmp)

    ctmp = cmplx( rtmp, itmp, DP ) / ctmp
    real_c(iter) = real(ctmp,DP)
    imag_c(iter) = -aimag(ctmp)
    
    ! step 5 w_{j+1) = w_{j+1} / \beta^*_{j}
    call OCEAN_psi_divide( hpsi, ierr, real_b(iter), imag_b(iter) )
    if( ierr .ne. 0 ) return

    ! step 6 v_{j+1} = v_{j+1}/ \delta_j
    call OCEAN_psi_divide( back_hpsi, ierr, real_c(iter), -imag_c(iter) )
    if( ierr .ne. 0 ) return

    ! more copies and prep for next iter
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( back_psi, back_hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return


    if( myid .eq. 0 ) then
      write ( 6, '(1x,6(f14.8,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter
      if( mod( iter, 1 ) .eq. 0 ) call haydump( iter, sys, psi%kpref, ierr )
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

!      call haydump( iter, sys, ierr )
    endif
    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_finish_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return

  end subroutine haydock_abc


  subroutine haydock_abc_1( sys, psi, hpsi, old_psi, back_psi, back_hpsi, back_old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : o_system
    use OCEAN_psi
    use OCEAN_mpi, only : myid, root, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    type(OCEAN_vector), intent(inout) :: back_psi, back_hpsi, back_old_psi

    complex(dp) :: ctmp
    real(dp) :: rtmp, itmp, atmp, btmp
    integer :: irequest, rrequest


    ! Following Saad
    ! Inner product (x,y) = \sum_{i=1}^{m} x_i \bar{y}_i
    !  This means that (x,y) = \langle y \vert x \rangle

    ! step 0, enforce allow
    call OCEAN_energies_allow( sys, hpsi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_energies_allow( sys, back_hpsi, ierr )
    if( ierr .ne. 0 ) return


    ! step 1: 
    ! $ \alpha_j = ( A v_j, w_j )
    call OCEAN_psi_dot( back_psi, hpsi, rrequest, rtmp, ierr, irequest, itmp )


    ! step 2:  New vector  ! HERE THERE IS A DIFFERENCE, beta index
    !  v_{j+1} = A v_j - \alpha_j v_j - \beta_{j-1} v_{j-1}
    ! step 2A) v_{j+1} = A v_j - \beta_{j-1} v_{j-1}
    atmp = - real_b( iter - 1 )
    btmp = - imag_b( iter - 1 )
    call OCEAN_psi_axpy( atmp, old_psi, hpsi, ierr, btmp )

    ! step 3: New vector for the back 
    !  w_{j+1} = A^\dagger w_j - \alpha^*_j w_j - \delta^*_{j-1} w_{j-1}
    ! step 2A) w_{j+1} = A^\dagger W_j - \delta^*_{j-1} w_{j-1}
    atmp = - real_c( iter - 1 )
    btmp =   imag_c( iter - 1 )
    call OCEAN_psi_axpy( atmp, back_old_psi, back_hpsi, ierr, btmp )

    ! Now need to make sure alpha is done
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    !
    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    

    ! step 2B:
    ! v_{j+1} -= \alpha_j v_j
    atmp = -rtmp
    btmp = -itmp
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    ! step 3B:
    ! w_{j+1} -= \alpha^*_j w_j
    atmp = -rtmp
    btmp =  itmp
    call OCEAN_psi_axpy( atmp, back_psi, back_hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return


    real_a( iter-1 ) = rtmp
    imag_a( iter-1 ) = itmp

    ! quick allow enforcement
    call OCEAN_energies_allow( sys, hpsi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_energies_allow( sys, back_hpsi, ierr )
    if( ierr .ne. 0 ) return


    ! Step 4A: ( v_{j+1), w_{j+1} ) 
    call OCEAN_psi_dot( back_hpsi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
    
    ! get ready for next iteration
    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return
    !
    call OCEAN_psi_copy_min( back_old_psi, back_psi, ierr )
    if( ierr .ne. 0 ) return    


    ! Wait on 4A
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    !
    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    

    real_c( iter ) = sqrt( (sqrt(rtmp**2 + itmp**2) + rtmp )/2.0_DP) 
    imag_c( iter ) = sign( sqrt( (sqrt(rtmp**2 + itmp**2) - rtmp )/2.0_DP), itmp)  

    ctmp = cmplx( rtmp, itmp, DP ) / cmplx( real_c( iter ), imag_c( iter ), DP )
    real_b( iter ) = real( ctmp, DP )
    imag_b( iter ) = -aimag( ctmp )

!    real_b( iter ) = sqrt( abs( cmplx( rtmp, itmp, DP ) ) )
!    imag_b( iter ) = 0.0_DP

!    ctmp = cmplx( rtmp, itmp, DP ) / real_b(iter)
!    real_c = real( ctmp, DP )
!    imag_c = aimag( ctmp )

#if 0
    ctmp = sqrt( cmplx( rtmp, itmp, DP ) )
    real_c( iter ) = real( ctmp, DP )
    imag_c( iter ) = aimag( ctmp )

    ctmp = cmplx( rtmp, itmp, DP )/cmplx( real_c(iter), imag_c(iter), DP )
    real_b( iter ) = real( ctmp, DP )
    imag_b( iter ) = aimag( ctmp )
#endif

    ! step 5 w_{j+1) = w_{j+1} / \beta^*_{j}
    call OCEAN_psi_divide( back_hpsi, ierr, real_b(iter), -imag_b(iter) )
    if( ierr .ne. 0 ) return    

    ! step 6 v_{j+1} = v_{j+1}/ \delta_j
    call OCEAN_psi_divide( hpsi, ierr, real_c(iter), imag_c(iter) )
    if( ierr .ne. 0 ) return


    ! more copies and prep for next iter
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( back_psi, back_hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return



    if( myid .eq. 0 ) then
      write ( 6, '(1x,6(f14.8,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter
      if( mod( iter, 1 ) .eq. 0 ) call haydump( iter, sys, psi%kpref, ierr )
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

!      call haydump( iter, sys, ierr )
    endif
    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_finish_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return

  end subroutine haydock_abc_1



  subroutine OCEAN_hay_abc( sys, psi, hpsi, old_psi, back_psi, back_hpsi, back_old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : o_system
    use OCEAN_psi
    use OCEAN_mpi, only : myid, root, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    type(OCEAN_vector), intent(inout) :: back_psi, back_hpsi, back_old_psi

    complex(dp) :: ctmp
    real(dp) :: rtmp, itmp, atmp, btmp
    integer :: ialpha, ikpt, irequest, rrequest

!    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_allow( sys, back_hpsi, ierr )
      if( ierr .ne. 0 ) return
!    endif

    ! calc ctmp = < hpsi | back_psi > and begin Iallreduce
!DERP
    call OCEAN_psi_dot( back_psi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
!    call OCEAN_psi_dot( psi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_dot'

    ! hpsi -= b(i-1) * psi^{i-1}
    ! y:= a*x + y
    atmp = -real_b(iter-1)
    btmp = -imag_b(iter-1)
    call OCEAN_psi_axpy( atmp, old_psi, hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

!DERP
!    atmp = -real_c(iter-1)
!    btmp = -imag_c(iter-1)
    atmp = -real_c(iter-1)
    btmp =  imag_c(iter-1)
    call OCEAN_psi_axpy( atmp, back_old_psi, back_hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_axpy 1'

    ! finish allreduce to get a
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    atmp = -rtmp
    btmp = -itmp
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    btmp = itmp
    call OCEAN_psi_axpy( atmp, back_psi, back_hpsi, ierr, btmp )
    if( ierr .ne. 0 ) return

    real_a(iter-1) = rtmp
    imag_a(iter-1) = itmp

!    if( myid .eq. root ) write(6,*) 'ab', real_a(iter-1), b(iter-1)

!    if( sys%cur_run%have_val ) then
      call OCEAN_energies_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_allow( sys, back_hpsi, ierr )
      if( ierr .ne. 0 ) return
!    endif

!DERP
!    call OCEAN_psi_dot( back_hpsi, hpsi, rrequest, rtmp, ierr, irequest, itmp )
    call OCEAN_psi_dot( hpsi, back_hpsi, rrequest, rtmp, ierr, irequest, itmp )
    if( ierr .ne. 0 ) return

    ! copies psi onto old_psi
    call OCEAN_psi_copy_min( old_psi, psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( back_old_psi, back_psi, ierr )
    if( ierr .ne. 0 ) return

    ! Could move prep and copy here for hspi -> psi
    !   just need to include a way to scale full instead of just min

    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    ctmp = sqrt( cmplx( rtmp, itmp, DP ) ) 
    
    real_b( iter ) = real( ctmp, DP )
    imag_b( iter ) = aimag( ctmp )

    ctmp = cmplx( rtmp, itmp, DP ) / ctmp

    real_c( iter ) = real( ctmp, DP )
    imag_c( iter ) = aimag( ctmp )

!    real_b( iter ) = sqrt( rtmp )
!    imag_b( iter ) = 0.0_DP
!    real_c( iter ) = sqrt( rtmp )
!    imag_c( iter ) = 0.0_DP


    call OCEAN_psi_divide( back_hpsi, ierr, real_c(iter), imag_c(iter) )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_divide( hpsi, ierr, real_b(iter), -imag_b(iter) )
    if( ierr .ne. 0 ) return
    !

    ! copies hpsi onto psi
    call OCEAN_psi_copy_min( psi, hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_copy_min( back_psi, back_hpsi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_prep_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_start_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return



    if( myid .eq. 0 ) then
      write ( 6, '(1x,6(f14.8,2x),i6)' ) real_a(iter-1)*Hartree2eV, imag_a(iter-1) * Hartree2eV, &
                                                    real_b(iter) * Hartree2eV, imag_b(iter) * Hartree2eV, &
                                                    real_c(iter) * Hartree2eV, imag_c(iter) * Hartree2eV, iter
      if( mod( iter, 1 ) .eq. 0 ) call haydump( iter, sys, psi%kpref, ierr )
#ifdef __HAVE_F03
      if( ieee_is_nan( real_a(iter-1) ) ) then
#else
      if( real_a(iter-1) .ne. real_a(iter-1) ) then
#endif
        write(6,*) 'NaN detected'
        ierr = -1
        return
      endif

!      call haydump( iter, sys, ierr )
    endif
    ! Might be moved up & out?
    call OCEAN_psi_finish_min2full( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_finish_min2full( back_psi, ierr )
    if( ierr .ne. 0 ) return


  end subroutine OCEAN_hay_abc


  subroutine check_convergence( iter1, iter2, sys, kpref, maxDiff, relArea )
    use OCEAN_system, only : o_system
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter1, iter2
    real(DP), intent( in ) :: kpref
    real(DP), intent( out ) :: maxDiff, relArea

    real(DP), allocatable :: sp1(:,:), sp2(:,:)
    real(DP) :: area1, area2, diff
    integer :: i

    allocate( sp1(3,ne), sp2(3,ne) )

    select case( sys%cur_run%calc_type)
      case( 'XES', 'XAS' )
        call calc_spect_core( sp1, iter1, kpref )
        call calc_spect_core( sp2, iter2, kpref )
      case( 'VAL', 'RXS' )
        call calc_spect_val( sp1, iter1, kpref, sys%celvol, sys%valence_ham_spin )
        call calc_spect_val( sp2, iter2, kpref, sys%celvol, sys%valence_ham_spin )

      case default
        call calc_spect_core( sp1, iter1, kpref )
        call calc_spect_core( sp2, iter2, kpref )

    end select
    
    maxDiff = 0.0_DP
    relArea = 0.0_DP
    area1 = 0.0_DP
    area2 = 0.0_DP

    do i = 1, ne
      diff = abs(sp1(3,i) - sp2(3,i) )
      if( diff .gt. maxDiff ) maxDiff = diff
      relArea = relArea + diff
      area1 = area1 + abs(sp1(3,i))
      area2 = area2 + abs(sp2(3,i))
    enddo

    write(6,*) relArea, area1, area2
    relArea = 2.0_DP * relArea / ( area1 + area2 )


!    open(unit=99,file='check.txt',form='formatted', status='unknown')
!    do i = 1, ne
!      write(99,*) sp1(1,i), sp1(3,i), sp2(3,i)
!    enddo
!    close(99)
  
    deallocate( sp1, sp2 )

  end subroutine check_convergence



  subroutine haydump( iter, sys, kpref, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_filenames, only : OCEAN_filenames_spectrum
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter
    real(DP), intent( in ) :: kpref

    integer :: ie, jdamp, jj
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk, pi
    complex(DP) :: rm1, ctmp, disc, delta

    character( LEN=40 ) :: abs_filename
    
    call OCEAN_filenames_spectrum( sys, abs_filename, ierr )
    if( ierr .ne. 0 ) return
    
!    rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
!    open( unit=99, file='absspct', form='formatted', status='unknown' )
    open( unit=99, file=abs_filename, form='formatted', status='unknown' )
    rewind 99

    select case ( sys%cur_run%calc_type)
      case( 'XES', 'XAS' )
        call write_core( 99, iter, kpref )
      case( 'VAL', 'RXS' )
        call write_val( 99, iter, kpref, sys%celvol, sys%valence_ham_spin )

      case default
        call write_core( 99, iter, kpref )
    
    end select

    close(unit=99)
    !
    return
  end subroutine haydump

  subroutine write_val( fh, iter, kpref , ucvol, val_ham_spin )
    use OCEAN_constants, only : Hartree2eV, bohr, alphainv
    implicit none
    integer, intent( in ) :: fh, iter, val_ham_spin
    real(DP), intent( in ) :: kpref, ucvol
    !
    integer :: ie, i
    real(DP) :: ere, reeps, imeps, lossf, fact, mu, reflct
    complex(DP) :: ctmp, arg, rp, rm, rrr, al, be, eps, refrac

    fact = kpref * real( 2 / val_ham_spin, DP ) * ucvol

    write(fh,"(a)") "#   omega (eV)      epsilon_1       epsilon_2       n"// &
      "               kappa           mu (cm^(-1))    R"//  &
      "               epsinv"

!p%kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
    do ie = 1, 2 * ne, 2
      ere = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
#if(1)
      ctmp = cmplx( ere, gam0, DP )

      arg = ( ere - real_a( iter - 1 ) ) ** 2 - 4.0_dp * real_b( iter ) ** 2
      arg = sqrt( arg )

      rp = 0.5_dp * ( ere - real_a( iter - 1 ) + arg )
      rm = 0.5_dp * ( ere - real_a( iter - 1 ) - arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = rp
      else
        rrr = rm
      endif

      al =  ctmp - real_a( iter - 1 ) - rrr
      be = -ctmp - real_a( iter - 1 ) - rrr
#else
      ctmp = cmplx( ere - real_a( iter - 1 ), gam0 - imag_a( iter -1 ), DP )
      arg = sqrt( ctmp ** 2 - 4.0_dp * cmplx( real_b( iter ), imag_b( iter ), DP ) &
                                     * cmplx( real_c( iter ), -imag_c( iter ), DP ) )
      rp = 0.5_dp * ( ctmp + arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = 0.5_dp * ( ctmp + arg )
      else
        rrr = 0.5_dp * ( ctmp - arg )
      endif

      ctmp = cmplx( ere, gam0, DP )
  
      al = ctmp - cmplx( real_a( iter-1 ), imag_a( iter-1 ), DP ) - rrr
      be = -ctmp - cmplx( real_a( iter-1 ), imag_a( iter-1 ), DP ) - rrr

#endif

      do i = iter-1, 0, -1
!        al =  ctmp - real_a( i ) - real_b( i + 1 ) ** 2 / al
!        be = -ctmp - real_a( i ) - real_b( i + 1 ) ** 2 / be
        al = ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / al
        be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / be
!        al =  ctmp - cmplx( real_a( i ), imag_a( i ), DP ) - real_b( i + 1 ) **2 / al
!        be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) - real_b( i + 1 ) **2 / be
      enddo

      eps = 1.0_dp - fact / al - fact / be

      reeps = dble( eps )
      imeps = aimag( eps )
!      rad = sqrt( reeps ** 2 + imeps ** 2 )
!      theta = acos( reeps / rad ) / 2
!      indref = sqrt( rad ) * cos( theta )
!      indabs = sqrt( rad ) * sin( theta )
!      ref = ( ( indref - 1 ) ** 2 + indabs ** 2 ) /
!   &        ( ( indref + 1 ) ** 2 + indabs ** 2 )
      lossf = imeps / ( reeps ** 2 + imeps ** 2 )

      refrac = sqrt(eps)
      reflct = abs((refrac-1.0d0)/(refrac+1.0d0))**2
      mu = 2.0d0 * ere * Hartree2eV * aimag(refrac) / ( bohr * alphainv * 1000 )

      write(fh,'(8(1E24.16,1X))') ere*Hartree2eV, reeps, imeps, refrac-1.0d0, mu, reflct, lossf

    enddo

  end subroutine write_val

  subroutine write_core( fh, iter, kpref )
    use OCEAN_constants, only : Hartree2eV
    implicit none
    integer, intent( in ) :: fh, iter
    real(DP), intent( in ) :: kpref
    !
    integer :: ie, jdamp, jj
    real(DP), external :: gamfcn
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk( 0 : 1 )
    complex(DP) :: ctmp, disc, delta, rm1
    !
    write( fh, '(A,1i5,A,1e15.8,A,1e15.8)' ) '#   iter=', iter, '   gam=', gam0, '   kpref=', kpref
    write( fh, '(5(A15,1x))' ) '#   Energy', 'Spect', 'Spect(0)', 'SPKK', 'SPKK(0)'
    rm1 = -1; rm1 = sqrt( rm1 )
    do ie = 0, 2 * ne, 2
       e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
       do jdamp = 0, 1
          gam= gam0 + gamfcn( e, nval, eps ) * dble( jdamp )
!          ctmp = e - a( iter - 1 ) + rm1 * gam

!          if( .true. ) then
#if(1)
          ctmp = cmplx( e - real_a( iter - 1 ), gam + imag_a( iter - 1 ), DP )  
          disc = sqrt( ctmp ** 2 - 4 * cmplx( real_b( iter ), imag_b( iter ) ) & 
                                     * cmplx( real_c( iter ), imag_c( iter ) ) )
          if( aimag( disc ) .gt. 0.0d0 ) then
            delta = (ctmp + disc ) / 2.0_dp
          else
            delta = (ctmp - disc ) / 2.0_dp
          endif

#else
            ctmp = e - real_a( iter - 1 ) + rm1 * gam
            disc = sqrt( ctmp ** 2 - 4 * real_b( iter ) ** 2 )
            di= -rm1 * disc
            if ( di .gt. 0.0d0 ) then
               delta = ( ctmp + disc ) / 2
            else
               delta = ( ctmp - disc ) / 2
            end if
#endif

          do jj = iter - 1, 0, -1
!             delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
!           if( .false. ) then
!             delta = e - real_a( jj ) + rm1 * gam - real_b( jj + 1 ) ** 2 / delta
!           else
            ctmp = cmplx( real_b( jj+1 ), imag_b( jj+1 ) ) * cmplx( real_c( jj+1 ), imag_c( jj+1 ) )
            delta = cmplx( e - real_a( jj ), gam + imag_a( jj ) ) - ctmp / delta
!           endif
          end do
          dr = delta
!          di = -rm1 * delta
!          di = abs( di )
          di = abs(aimag( delta ) )
          ener = ebase + Hartree2eV * e
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
         spkk( jdamp ) = kpref * dr / ( dr ** 2 + di ** 2 )
       end do
!       write ( fh, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
       write ( fh, '(5(1e15.8,1x))') ener, spct( 1 ), spct( 0 ), spkk( 1 ), spkk( 0 )
    end do
  end subroutine write_core

  subroutine calc_spect_core( sp, iter, kpref )
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    implicit none
    real(DP), intent( out ) :: sp(:,:)
    integer, intent( in ) :: iter
    real(DP), intent( in ) :: kpref

    integer :: ie, jj
    real(DP) :: e, dr, di
    complex(DP) :: ctmp, disc, delta

    do ie = 1, ne
      e = el + ( eh - el ) * real( 2*(ie-1)+1, DP ) / real( 2 * ne, DP )
      ctmp = cmplx( e - real_a( iter - 1 ), gam0 + imag_a( iter - 1 ), DP )
      disc = sqrt( ctmp ** 2 - 4 * cmplx( real_b( iter ), imag_b( iter ) ) &
                                 * cmplx( real_c( iter ), -imag_c( iter ) ) )
      if( aimag( disc ) .gt. 0.0d0 ) then
        delta = (ctmp + disc ) / 2.0_dp
      else
        delta = (ctmp - disc ) / 2.0_dp
      endif

      do jj = iter -1, 0, -1
        ctmp = cmplx( real_b( jj+1 ), imag_b( jj+1 ) ) * cmplx( real_c( jj+1 ), -imag_c( jj+1 ) )
        delta = cmplx( e - real_a( jj ), gam0 + imag_a( jj ) ) - ctmp / delta
      enddo

      dr = delta
      di = abs( aimag( delta ) )
      sp(1,ie) = ebase + Hartree2eV * e
      sp(2,ie) = kpref * dr / ( dr ** 2 + di ** 2 )
      sp(3,ie) = kpref * di / ( dr ** 2 + di ** 2 )
    enddo

  

  end subroutine calc_spect_core

  subroutine calc_spect_val( sp, iter, kpref, celvol, nspin )
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    implicit none
    real(DP), intent( out ) :: sp(:,:)
    integer, intent( in ) :: iter
    real(DP), intent( in ) :: kpref
    real(DP), intent( in ) :: celvol
    integer, intent( in ) :: nspin

    integer :: ie, jj, i
    real(DP) :: e, dr, di, fact
    complex(DP) :: ctmp, disc, delta, arg, rp, rm , rrr, al, be, eps

    fact = kpref * real( 2 / nspin, DP ) * celvol

    do ie = 1, ne
      e = el + ( eh - el ) * real( 2*(ie-1)+1, DP ) / real( 2 * ne, DP )

      ctmp = cmplx( e, gam0, DP )
      arg =  ( e - real_a( iter - 1 ) )**2 - 4.0_dp * real_b( iter ) ** 2 
      arg = sqrt(arg)

      rp = 0.5_dp * ( e - real_a( iter - 1 ) + arg )
      rm = 0.5_dp * ( e - real_a( iter - 1 ) - arg )
      if( aimag( rp ) .lt. 0.0_dp ) then
        rrr = rp
      else
        rrr = rm
      endif 
      al =  ctmp - real_a( iter - 1 ) - rrr
      be = -ctmp - real_a( iter - 1 ) - rrr

      do i = iter-1, 0, -1
        al = ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / al
        be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
           - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / be
      enddo

      eps = 1.0_dp - fact / al - fact / be
      dr = real( eps, DP )
      di = aimag( eps ) 
      sp(1,ie) = ebase + Hartree2eV * e
      sp(2,ie) = dr
      sp(3,ie) = di
!      sp(2,ie) = fact * dr / ( dr ** 2 + di ** 2 )
!      sp(3,ie) = fact * di / ( dr ** 2 + di ** 2 )
    enddo

    

  end subroutine calc_spect_val

  subroutine OCEAN_haydock_setup( sys, ierr )
    use OCEAN_mpi
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    use OCEAN_system
    implicit none

    type(o_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: dumi, iter, ierr_
    character(len=4) :: inv_style
    real( DP ) :: dumf
    real( DP ), parameter :: default_gam0 = 0.1_DP

    if( .not. is_first ) goto 10
    is_first = .false.

    if( myid .eq. root ) then

      open(unit=99,file='bse.in',form='formatted',status='old')
      rewind(99)
      read(99,*) dumi
      read(99,*) dumf
      read(99,*) calc_type
!      select case ( calc_type )
!        case('hay')
          read(99,*) haydock_niter, ne, el, eh, gam0, ebase
          call checkBroadening( sys, gam0, default_gam0 )

!          el = el / 27.2114d0
!          eh = eh / 27.2114d0
!          gam0 = gam0 / 27.2114d0
          el = el * eV2Hartree
          eh = eh * eV2Hartree
          gam0 = gam0 * eV2Hartree
!          inv_loop = 1
!          allocate( e_list( inv_loop ) )
#if( 0 )
        case('inv')
          read(99,*) nloop, gres, gprc, ffff, ener
          ! if gres is negative fill it with core-hole lifetime broadening
          call checkBroadening( sys, gres, default_gam0 )
          gprc = gprc * eV2Hartree
          gres = gres * eV2Hartree

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
            case default
              write(6,*) 'Error reading bse.in'
              ierr = -1
          end select          

          inquire(file='echamp.inp',exist=echamp)
          if( echamp ) then
            open(98,file='echamp.inp',form='formatted',status='old')
            read(98,*) echamp
            close(98)
          endif

          ! To allow absspct projected by various attributes of either the core or conduction states
          !   For starters only the conduction-band spin will be enabled.
          project_absspct = .false.
          inquire(file='project_absspct.inp',exist=project_absspct)

        case default
          ierr = -1
      end select
#endif
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
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( haydock_niter, 1, MPI_INTEGER, root, comm, ierr )

!    call MPI_BCAST( nloop, 1, MPI_INTEGER, root, comm, ierr )
!    call MPI_BCAST( gres, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( gprc, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( ffff, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( ener, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( eps, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( nval, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )


!    call MPI_BCAST( e_start, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( e_stop, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( e_step, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    call MPI_BCAST( inv_loop, 1, MPI_INTEGER, root, comm, ierr )
!    if( myid .ne. root ) allocate( e_list( inv_loop ) )
!    call MPI_BCAST( e_list, inv_loop, MPI_DOUBLE_PRECISION, root, comm, ierr )


    call MPI_BCAST( echamp, 1, MPI_LOGICAL, root, comm, ierr )
    call MPI_BCAST( project_absspct, 1, MPI_LOGICAL, root, comm, ierr )
#endif

10 continue

!    if( allocated( a ) ) deallocate( a )
!    if( allocated( b ) ) deallocate( b )
    if( allocated( real_a ) ) deallocate( real_a )
    if( allocated( imag_a ) ) deallocate( imag_a )
    if( allocated( real_b ) ) deallocate( real_b )
    if( allocated( imag_b ) ) deallocate( imag_b )
    if( allocated( real_c ) ) deallocate( real_c )
    if( allocated( imag_c ) ) deallocate( imag_c )
    if( haydock_niter .gt. 0 ) then
!      allocate( a( 0 : haydock_niter ) )
!      allocate( b( 0 : haydock_niter ) )
      allocate( real_a( 0 : haydock_niter ) )
      allocate( imag_a( 0 : haydock_niter ) )
      allocate( real_b( 0 : haydock_niter ) )
      allocate( imag_b( 0 : haydock_niter ) )
      allocate( real_c( 0 : haydock_niter ) )
      allocate( imag_c( 0 : haydock_niter ) )
!      a(:) = 0.0_DP
!      b(:) = 0.0_DP
      real_a(:) = 0.0_DP
      imag_a(:) = 0.0_DP
      real_b(:) = 0.0_DP
      imag_b(:) = 0.0_DP
      real_c(:) = 0.0_DP
      imag_c(:) = 0.0_DP
!    else
!      allocate( a(1), b(1) )
    endif

  end subroutine OCEAN_haydock_setup


  subroutine checkBroadening( sys, broaden, default_broaden )
    use OCEAN_corewidths, only : returnLifetime
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    implicit none
    type( o_system ), intent( in ) :: sys
    real(DP), intent( inout ) :: broaden
    real(DP), intent( in ) :: default_broaden
    
    if( broaden .gt. 0.0_dp ) return


    select case ( sys%cur_run%calc_type )
      case( 'VAL' )
        broaden = default_broaden 
      case( 'XAS' , 'XES', 'RXS' ) 
        call returnLifetime( sys%cur_run%ZNL(1), sys%cur_run%ZNL(2), sys%cur_run%ZNL(3), broaden )
        if( broaden .le. 0 ) broaden = default_broaden
        if( myid .eq. root ) write(6,*) 'Default broadening used: ', broaden
      case default
        broaden = default_broaden
    end select
    write(6,*) 'Default requested for broadening: ', broaden

    end subroutine checkBroadening

  subroutine redtrid(n,sys, kpref, ierr)
    use OCEAN_system, only : o_system
    use OCEAN_filenames, only : OCEAN_filenames_lanc
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) ::  n
    real(DP), intent( in ) :: kpref

    double precision, allocatable :: ar(:,:),ai(:,:)
    double precision, allocatable :: w(:),zr(:,:),zi(:,:)
    double precision, allocatable :: fv1(:),fv2(:),fm1(:)
    integer :: matz,nm,i,j,nn


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
      ai(i,i)=imag_a(i-1)
      if (i.le.n) then
        ar(i+1,i)=real_c(i)
        ar(i,i+1)=real_b(i)
        ai(i+1,i)=imag_c(i)
        ai(i,i+1)=imag_b(i)
      end if
    end do
    matz=0
    call elsch(nm,nn,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
    if( n .lt. 0 ) return

    call write_lanczos( n, sys, kpref, ierr, w )

    deallocate(ar,ai,w,zr,zi,fv1,fv2,fm1)
    return
  end subroutine redtrid

  subroutine write_lanczos( n, sys, kpref, ierr, w )
    use OCEAN_system, only : o_system
    use OCEAN_filenames, only : OCEAN_filenames_lanc
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) ::  n
    real(DP), intent( in ) :: kpref
    real(DP), intent( in ), optional :: w(:)

    integer :: i
    character( LEN=40 ) :: lanc_filename

    call OCEAN_filenames_lanc( sys, lanc_filename, ierr )
    if( ierr .ne. 0 ) return

    open(unit=99,file=lanc_filename,form='formatted',status='unknown')
    rewind 99
    write ( 99, '(1i8,1x,1ES24.17)' ) n, kpref
    if( complex_haydock ) then
      write ( 99, '(2(2x,ES24.17))' ) real_a( 0 ), imag_a( 0 )
      do i = 1, n
        write ( 99, '(2x,6ES24.17)' ) real_a( i ), imag_a( i ), real_b( i ), imag_b( i ), &
                                     real_c( i ), imag_c( i )
      enddo
    else
      do i = 0, n
        if ( i .eq. 0 ) then
          write ( 99, '(2x,ES24.17)' ) real_a( i )
        else
          write ( 99, '(2(2x,ES24.17))' ) real_a( i ), real_b( i )
        end if
      end do
    endif

  
    if( present( w ) ) then
      write (99,'(2x,2i5,1f20.10)') (i,n+1,w(i),i=1,n+1)
    endif
    close(unit=99)
  end subroutine write_lanczos


#if( 0 )
  ! When using GMRES we can project out only part of the exciton.
  ! For now this is hard-coded for only doing spin up/down for the conduction band
  ! In the future we should add things like spin orbit, 3d symmetries, etc.
  subroutine write_projected_absspct( sys, project_file_handle, kpref, ener, ntot, rhs, x, ierr )
    use OCEAN_system
    use OCEAN_constants, only : Hartree2eV
    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: project_file_handle, ntot
    real(dp), intent(in) :: kpref, ener
    complex(DP), intent( IN ):: rhs( ntot ), x( ntot )
    !
    complex(DP) :: dot_up, dot_down
    integer :: ialpha, ik, iband, ispn
    integer :: i

    ! Unfortunately rhs and x have only a single index
    ! The structure of the vector is
    ! do core spin
    !   do core angular momentum
    !     do valence spin
    !       do kpts
    !         do bands


!JTV need to doublecheck which is up and which is down
    dot_up = 0.0_DP
    dot_down = 0.0_DP
    i = 0
    
    do ialpha = 1, sys%nalpha / 2
      do ik = 1, sys%nkpts
        do iband = 1, sys%num_bands
          i = i + 1
          dot_up = dot_up + conjg( rhs( i ) ) * x( i )
        enddo
      enddo
      do ik = 1, sys%nkpts
        do iband = 1, sys%num_bands
          i = i + 1
          dot_down = dot_down + conjg( rhs( i ) ) * x( i )
        enddo
      enddo
    enddo
    
    write ( project_file_handle, '(5(1x,1e15.8))' ) ener*Hartree2eV, ( 1.0d0 - dot_up ) * kpref, & 
                                                    ( 1.0d0 - dot_up ) * kpref
    call flush(project_file_handle)

!    write ( project_file_handles( 2 ), '(3(1x,1e15.8))' ) ener*27.2114_DP, ( 1.0d0 - dot_down ) * kpref
!    call flush(project_file_handles(2))
    

  end subroutine write_projected_absspct
#endif

  subroutine testConvergeEps( iter, sys, kpref, ucvol, val_ham_spin, restartBSE, newEps )
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : PI_DP
    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( in ) :: iter, val_ham_spin
    real(DP), intent( in ) :: kpref, ucvol
    logical, intent( inout ) :: restartBSE
    real(DP), intent( out ) :: newEps

    complex(DP) :: ctmp, arg, rrr, rp, rm, al, be, eps
    real(DP) :: tcEps, fact, oldEps, epsErr
    integer :: i

    fact = kpref * real( 2 / val_ham_spin, DP ) * ucvol
    ctmp = cmplx( 0, gam0, DP )

    arg = real_a( iter - 1 )** 2 - 4.0_dp * real_b( iter ) ** 2
    arg = sqrt( arg )
  
    rp = 0.5_dp * ( - real_a( iter - 1 ) + arg )
    rm = 0.5_dp * ( - real_a( iter - 1 ) - arg )
    if( aimag( rp ) .lt. 0.0_dp ) then
      rrr = rp
    else
      rrr = rm
    endif

    al = ctmp - real_a( iter - 1 ) - rrr
    be = -ctmp - real_a( iter - 1 ) - rrr

    do i = iter-1, 0, -1
      al = ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
         - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / al
      be = -ctmp - cmplx( real_a( i ), imag_a( i ), DP ) &
         - cmplx( real_b( i+1 ), imag_b( i+1 ), DP ) * cmplx( real_c( i+1 ), -imag_c( i+1 ), DP ) / be
    enddo

    eps = 1.0_dp - fact / al - fact / be

    tcEps = abs( dble( eps ) )

    if( iter .lt. 3 ) then
      eps1Conv( iter ) = tcEps
!      write(6,*) 'Estimated eps1(0): ', tcEps
    else
      eps1Conv( 1 ) = eps1Conv( 2 )
      eps1Conv( 2 ) = eps1Conv( 3 )
      eps1Conv( 3 ) = tcEps

      tcEps = sum(eps1Conv(:)) / 3.0_DP
      if( max( sys%epsilon0, eps1Conv( 3 ) ) .lt. 100.0d0 ) then
        write(6,'(3(A,F9.4,X))') 'Est. eps1(0): ', eps1Conv( 3 ), ';  Avg: ', tcEps, ';  Current: ', sys%epsilon0
      else
        write(6,'(3(A,E24.12,X))') 'Est. eps1(0): ', eps1Conv( 3 ), ';  Avg: ', tcEps, ';  Current: ', sys%epsilon0
      endif

      ! change to percentage
      if( ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) )/tcEps .gt. 0.05_dp ) return
      if( abs( sys%epsilon0 - tcEps ) / ( sys%epsilon0 + tcEps - 2.0_dp ) &
                  .lt. 10.0_DP * sys%epsConvergeThreshold ) then
        if( ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) )/tcEps .gt. 0.5_dp * sys%epsConvergeThreshold ) then
!            write(6,*) 'C', ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) )/tcEps
            return
        endif
      endif
          

      if( abs( sys%epsilon0 - tcEps ) .gt. ( maxval(eps1Conv(:)) - minval(eps1Conv(:)) ) .and. &
          abs( sys%epsilon0 - tcEps ) / ( sys%epsilon0 + tcEps - 2.0_dp ) & 
                  .gt. 0.5_DP * sys%epsConvergeThreshold ) then  
        newEps = ( 2.0_DP * eps1Conv( 3 ) + eps1Conv(2) ) / 3.0_DP
#if 0
        if( abs( newEps - sys%epsilon0 )/( newEps + sys%epsilon0 ) .gt. 0.15_dp ) then
          newEps = 0.95_dp * newEps + 0.05_DP * sys%epsilon0
        elseif( abs( newEps - sys%epsilon0 )/( newEps + sys%epsilon0 ) .gt. 0.02_dp ) then
          newEps = 0.98_dp * newEps + 0.02_DP * sys%epsilon0
        endif
#else
        epsErr = 20.0_dp * abs( newEps - sys%epsilon0 )/( newEps + sys%epsilon0 - 2.0_dp )
        epsErr = ( 0.2_DP / PI_DP ) * atan( epsErr )
        write(6,*) newEps, sys%epsilon0, epsErr
        newEps = (1.0_DP-epsErr)*newEps + epsErr * sys%epsilon0
#endif
        
        
        restartBSE = .true.
        write(6,*) 'Restart eps1(0): ', newEps, sys%epsilon0
        return
      endif
    endif

  end subroutine testConvergeEps

end module OCEAN_haydock
