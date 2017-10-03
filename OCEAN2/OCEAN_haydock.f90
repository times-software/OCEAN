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


  REAL(DP), ALLOCATABLE :: a( : )
  REAL(DP), ALLOCATABLE :: b( : )
  REAL(DP), ALLOCATABLE :: real_a( : )
  REAL(DP), ALLOCATABLE :: imag_a( : )
  REAL(DP), ALLOCATABLE :: real_b( : )
  REAL(DP), ALLOCATABLE :: imag_b( : )
  REAL(DP), ALLOCATABLE :: real_c( : )
  REAL(DP), ALLOCATABLE :: imag_c( : )
 
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
  LOGICAL  :: project_absspct
  LOGICAL  :: is_first = .true.

  LOGICAL  :: val_loud = .true.

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


    call OCEAN_psi_new( psi, ierr, hay_vec )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( back_psi, ierr, hay_vec, conj=.true. )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( old_psi, ierr )
    if( ierr .ne. 0 ) return


    call OCEAN_psi_new( back_new_psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_new( back_old_psi, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_zero_min( back_old_psi, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) then
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', hay_vec%kpref
      write(6,*) inter_scale, haydock_niter
    endif
    call MPI_BARRIER( comm, ierr )
    do iter = 1, haydock_niter
      if( sys%cur_run%have_val ) then
        if( myid .eq. root ) write(6,*)   " iter. no.", iter-1
        call OCEAN_energies_val_allow( sys, psi, ierr )
        if( ierr .ne. 0 ) return
      endif

      call OCEAN_xact( sys, inter_scale, psi, new_psi, ierr )
      if( ierr .ne. 0 ) return

      ! need the action of the Hermitian conjugate of the Hamiltonian
      !  obviously we are only bothering to do this when H isn't Hermitian
!      call OCEAN_xact( sys, inter_scale, back_psi, back_new_psi, ierr, backwards=.true. )
      if( ierr .ne. 0 ) return

      ! This should be hoisted back up here
      call ocean_hay_abc( sys, psi, new_psi, old_psi, back_psi, back_new_psi, back_old_psi, & 
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






  subroutine OCEAN_haydock_do( sys, hay_vec, ierr )
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
    type( o_system ), intent( in ) :: sys
    !JTV need to figure out a work-around. Right now hay_vec is inout because of
    ! a depndency tracing back to calling copy and possibly copy_min, and
    ! possibly needing to go min->full, copy full, full->min
    type( ocean_vector ), intent( inout ) :: hay_vec

    real(DP) :: imag_a
    integer :: iter

!    character( LEN=21 ) :: lanc_filename

    type( ocean_vector ) :: psi, old_psi, new_psi
    

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
      write(6,*) inter_scale, haydock_niter
    endif
    call MPI_BARRIER( comm, ierr )



!    call OCEAN_tk_start( tk_psisum )

    do iter = 1, haydock_niter
      if( sys%cur_run%have_val ) then
        if( myid .eq. root ) write(6,*)   " iter. no.", iter-1
        call OCEAN_energies_val_allow( sys, psi, ierr )
        if( ierr .ne. 0 ) return
      endif

      call OCEAN_xact( sys, inter_scale, psi, new_psi, ierr )
      if( ierr .ne. 0 ) return
!      if( myid .eq. root ) write(6,*) 'Done with ACT'

!      if( sys%cur_run%have_val ) then
!        call OCEAN_energies_val_allow( sys, new_psi, ierr )
!        if( ierr .ne. 0 ) return
!      endif


      ! This should be hoisted back up here
      call ocean_hay_ab( sys, psi, new_psi, old_psi, iter, ierr )

    enddo

    call OCEAN_tk_stop( tk_psisum )
    call MPI_BARRIER( comm, ierr )

    if( myid .eq. 0 ) then
!      write(lanc_filename, '(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'lanceig_', sys%cur_run%elname, &
!        '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
      call haydump( haydock_niter, sys, hay_vec%kpref, ierr )
      call redtrid(  haydock_niter, sys, hay_vec%kpref, ierr )
    endif

    call OCEAN_psi_kill( psi, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_kill( new_psi, ierr )
    if( ierr .ne. 0 ) return
    
    call OCEAN_psi_kill( old_psi, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BARRIER( comm, ierr )

  end subroutine OCEAN_haydock_do

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
      write(6,*) inter_scale, haydock_niter, inv_loop


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
        call OCEAN_energies_val_allow( sys, psi, ierr )
        if( ierr .ne. 0 ) return
      endif


      call OCEAN_xact( sys, inter_scale, psi, hpsi, ierr )
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
            call OCEAN_energies_val_allow( sys, psi, ierr )
            if( ierr .ne. 0 ) return
          endif

          call OCEAN_xact( sys, inter_scale, psi, hpsi, ierr )
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
        relative_error = f( 2 ) / ( dimag( - dot_product( rhs, x ) ) ) !* kpref )
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


  subroutine OCEAN_hay_ab( sys, psi, hpsi, old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_val_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi

    real(dp) :: btmp, atmp, aitmp
    integer :: ialpha, ikpt, arequest, airequest, brequest

    if( sys%cur_run%have_val ) then
      call OCEAN_energies_val_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
    endif

    ! calc ctmp = < hpsi | psi > and begin Iallreduce
    call OCEAN_psi_dot( hpsi, psi, arequest, atmp, ierr, airequest, aitmp )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_dot'

    ! hpsi -= b(i-1) * psi^{i-1}
    btmp = -b(iter-1)
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
    if( myid .eq. root ) write(6,*) 'ab', real_a(iter-1), b(iter-1)
    call OCEAN_psi_axpy( atmp, psi, hpsi, ierr )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_axpy 2'


    if( sys%cur_run%have_val ) then
      call OCEAN_energies_val_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
    endif

    
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

    b(iter) = sqrt( btmp )
    btmp = 1.0_dp / b( iter )
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
      write ( 6, '(2x,2f24.13,10x,1e24.13,8x,i6)' ) real_a(iter-1) * Hartree2eV, b(iter) * Hartree2eV, &
                                                  imag_a(iter-1) * Hartree2eV, iter
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
    

  end subroutine OCEAN_hay_ab


  subroutine OCEAN_hay_abc( sys, psi, hpsi, old_psi, back_psi, back_hpsi, back_old_psi, iter, ierr )
#ifdef __HAVE_F03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif
    use OCEAN_system, only : o_system
    use OCEAN_psi
    use OCEAN_mpi, only : myid, root, MPI_STATUS_IGNORE
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_energies, only : OCEAN_energies_val_allow
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi
    type(OCEAN_vector), intent(inout) :: back_psi, back_hpsi, back_old_psi

    complex(dp) :: ctmp
    real(dp) :: rtmp, itmp
    integer :: ialpha, ikpt, irequest, rrequest

    if( sys%cur_run%have_val ) then
      call OCEAN_energies_val_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_val_allow( sys, back_hpsi, ierr )
      if( ierr .ne. 0 ) return
    endif

    ! calc ctmp = < hpsi | back_psi > and begin Iallreduce
    call OCEAN_psi_dot( hpsi, back_psi, rrequest, rtmp, ierr, irequest, itmp )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_dot'

    ! hpsi -= b(i-1) * psi^{i-1}
    ! y:= a*x + y
    call OCEAN_psi_axmy( old_psi, hpsi, ierr, real_b(iter-1), imag_b(iter-1) )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_axmy( back_old_psi, back_hpsi, ierr, real_c(iter-1), imag_c(iter-1) )
    if( ierr .ne. 0 ) return
!    if( myid .eq. root ) write(6,*) 'psi_axpy 1'

    ! finish allreduce to get a
    call MPI_WAIT( rrequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( irequest, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_axmy( psi, hpsi, ierr, rtmp, itmp )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_axmy( back_psi, back_hpsi, ierr, rtmp, -itmp )
    if( ierr .ne. 0 ) return

    real_a(iter-1) = rtmp
    imag_a(iter-1) = itmp

!    if( myid .eq. root ) write(6,*) 'ab', real_a(iter-1), b(iter-1)

    if( sys%cur_run%have_val ) then
      call OCEAN_energies_val_allow( sys, hpsi, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_energies_val_allow( sys, back_hpsi, ierr )
      if( ierr .ne. 0 ) return
    endif

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
!      write ( 6, '(2x,2f10.6,10x,1e11.4,x,f6.3)' ) a(iter-1), b(iter), imag_a, time2-time1
!      write ( 6, '(2x,2f10.6,10x,1e11.4,8x,i6)' ) a(iter-1), b(iter), imag_a, iter
      write ( 6, '(2x,2f24.13,10x,1e24.13,8x,i6)' ) real_a(iter-1) * Hartree2eV, b(iter) * Hartree2eV, &
                                                  imag_a(iter-1) * Hartree2eV, iter
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


  end subroutine OCEAN_hay_abc



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
      ctmp = cmplx( ere, gam0, DP )

      arg = ( ere - real_a( iter - 1 ) ) ** 2 - 4.0_dp * b( iter ) ** 2
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

      do i = iter-1, 0, -1
        al =  ctmp - real_a( i ) - b( i + 1 ) ** 2 / al
        be = -ctmp - real_a( i ) - b( i + 1 ) ** 2 / be
      enddo

      eps = 1.0_dp - fact / al - fact / be

      reeps = dble( eps )
      imeps = dimag( eps )
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
    real(DP) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk
    complex(DP) :: rm1, ctmp, disc, delta
    !
    rm1 = -1; rm1 = sqrt( rm1 )
    do ie = 0, 2 * ne, 2
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
          ener = ebase + Hartree2eV * e
          spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
       end do
       spkk = kpref * dr / ( dr ** 2 + di ** 2 )
       write ( fh, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, iter, gam, kpref, ne
    end do
  end subroutine write_core

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

    call MPI_BCAST( inter_scale, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( haydock_niter, 1, MPI_INTEGER, root, comm, ierr )
!    call MPI_BCAST( calc_type, 3, MPI_CHARACTER, root, comm, ierr )

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

  end subroutine OCEAN_haydock_setup


  subroutine checkBroadening( sys, broaden, default_broaden )
    use OCEAN_corewidths, only : returnLifetime
    use OCEAN_system
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
!    real( DP ), intent( in ) :: a( 0 : n ), b( n )
!    character( LEN=21 ), intent( in ) :: lanc_filename
    double precision, allocatable :: ar(:,:),ai(:,:)
    double precision, allocatable :: w(:),zr(:,:),zi(:,:)
    double precision, allocatable :: fv1(:),fv2(:),fm1(:)
    integer :: matz,nm,i,j,nn


    character( LEN=24 ) :: lanc_filename

    call OCEAN_filenames_lanc( sys, lanc_filename, ierr )
    if( ierr .ne. 0 ) return



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



end module OCEAN_haydock
