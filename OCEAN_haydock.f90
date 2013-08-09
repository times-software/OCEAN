module OCEAN_action
  use AI_kinds
  implicit none
  private
  save  


  REAL(DP), ALLOCATABLE :: a( : )
  REAL(DP), ALLOCATABLE :: b( : )
 
  REAL(DP) :: inter_scale_threshold = 0.00001
  REAL(DP) :: inter_scale

  REAL(DP) :: el, eh, gam0, eps, nval,  ebase
  REAL(DP) :: gres, gprc, ffff, ener
  
  INTEGER  :: haydock_niter = 0
  INTEGER  :: ne
  INTEGER  :: nloop
  

  CHARACTER(LEN=3) :: calc_type


  REAL(DP), POINTER, CONTIGUOUS :: mem_psi_r(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_psi_i(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_hpsi_r(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_hpsi_i(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_oldpsi_r(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_oldpsi_i(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_newpsi_r(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_newpsi_i(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_mulpsi_r(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_mulpsi_i(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_lrpsi_r(:,:,:)
  REAL(DP), POINTER, CONTIGUOUS :: mem_lrpsi_i(:,:,:)

  public :: OCEAN_haydock, OCEAN_hayinit

  contains


  subroutine OCEAN_hay_alloc( sys, hay_vec, psi, hpsi, old_psi, new_psi, mul_psi, lr_psi, ierr )
    use OCEAN_system
    use OCEAN_psi
    use iso_c_binding

    implicit none
    include 'fftw3.f03'
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    type( ocean_vector ), intent( out ) :: psi, hpsi, old_psi, new_psi, mul_psi, lr_psi

    type(C_PTR) :: cptr


    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_psi_r, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_psi_i, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
  

    psi%r => mem_psi_r
    psi%i => mem_psi_i
    psi%r(:,:,:) = hay_vec%r(:,:,:)
    psi%i(:,:,:) = hay_vec%i(:,:,:)
    psi%bands_pad = hay_vec%bands_pad
    psi%kpts_pad  = hay_vec%kpts_pad




    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_hpsi_r, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_hpsi_i, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )

    hpsi%r => mem_hpsi_r
    hpsi%i => mem_hpsi_i
    hpsi%r = 0.0_DP
    hpsi%i = 0.0_DP
    hpsi%bands_pad = hay_vec%bands_pad
    hpsi%kpts_pad  = hay_vec%kpts_pad


    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_oldpsi_r, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_oldpsi_i, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
  
    old_psi%r => mem_oldpsi_r
    old_psi%i => mem_oldpsi_i
    old_psi%r = 0.0_DP
    old_psi%i = 0.0_DP
    old_psi%bands_pad = hay_vec%bands_pad
    old_psi%kpts_pad  = hay_vec%kpts_pad


    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_newpsi_r, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_newpsi_i, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )

    new_psi%r => mem_newpsi_r
    new_psi%i => mem_newpsi_i
    new_psi%r = 0.0_DP
    new_psi%i = 0.0_DP
    new_psi%bands_pad = hay_vec%bands_pad
    new_psi%kpts_pad  = hay_vec%kpts_pad


    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_mulpsi_r, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_mulpsi_i, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
 
    mul_psi%r => mem_mulpsi_r
    mul_psi%i => mem_mulpsi_i
    mul_psi%r = 0.0_DP
    mul_psi%i = 0.0_DP
    mul_psi%bands_pad = hay_vec%bands_pad
    mul_psi%kpts_pad  = hay_vec%kpts_pad


    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_lrpsi_r, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
    cptr = fftw_alloc_real( int(hay_vec%bands_pad * hay_vec%kpts_pad * sys%nalpha, C_SIZE_T) )
    call c_f_pointer( cptr, mem_lrpsi_i, [hay_vec%bands_pad, hay_vec%kpts_pad, sys%nalpha ] )
 
    lr_psi%r => mem_lrpsi_r
    lr_psi%i => mem_lrpsi_i
    lr_psi%r = 0.0_DP
    lr_psi%i = 0.0_DP
    lr_psi%bands_pad = hay_vec%bands_pad
    lr_psi%kpts_pad  = hay_vec%kpts_pad


  end subroutine OCEAN_hay_alloc



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
    integer :: num_threads


    type( ocean_vector ) :: old_psi
    type( ocean_vector ) :: long_range_psi
    type( ocean_vector ) :: multiplet_psi
    type( ocean_vector ) :: hpsi
    type( ocean_vector ) :: new_psi
    type( ocean_vector ) :: psi


!$    integer, external :: omp_get_num_threads

    num_threads = 1
!$  num_threads = omp_get_num_threads()
    num_threads = min( 2, num_threads )
    
    call ocean_hay_alloc( sys, hay_vec, psi, hpsi, old_psi, new_psi, multiplet_psi, long_range_psi, ierr )


    if( myid .eq. root ) then 
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      write(6,*) inter_scale, haydock_niter
    endif
    do iter = 1, haydock_niter

      if( sys%long_range ) then

        if( sys%obf ) then
!          call lr_act_obf( sys, ierr )
        else
          call lr_act( sys, psi, long_range_psi, lr, ierr )
        endif

        if( nproc .gt. 1 ) then
! !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS num_threads

! !$OMP SECTIONS

! !$OMP SECTION
          call ocean_psi_sum_lr( sys, long_range_psi, ierr )

! !$OMP SECTION
          if( sys%mult ) call OCEAN_mult_act( sys, inter_scale, psi, multiplet_psi )


! !$OMP END SECTIONS

! !$OMP END PARALLEL
        else
          if( sys%mult ) call OCEAN_mult_act( sys, inter_scale, psi, multiplet_psi )
        endif
      else
        if( sys%mult ) call OCEAN_mult_act( sys, inter_scale, psi, multiplet_psi )
      endif

      if( sys%e0 ) call ocean_energies_act( sys, psi, hpsi, ierr )

      call ocean_psi_sum( sys, hpsi, multiplet_psi, long_range_psi, ierr )

      call ocean_hay_ab( sys, psi, hpsi, old_psi, iter, ierr )
!      call ocean_psi_ab( sys, a(iter-1), b(iter-1), b(iter), imag_a, psi, hpsi, old_psi, &
!                         ierr )
!!      call ocean_psi_dump( sys, ierr )
!      if( myid .eq. 0 ) then
!        write ( 6, '(2x,2f10.6,10x,1e11.4)' ) a(iter-1), b(iter), imag_a
!        call haydump( iter, ierr )
!
!!        call redtrid( iter-1, a(0), b(1) )
!      endif

    enddo

    if( myid .eq. 0 ) then
      call redtrid( haydock_niter-1, a(0), b(1) )
    endif
    
  end subroutine OCEAN_haydock


  subroutine OCEAN_hay_ab( sys, psi, hpsi, old_psi, iter, ierr )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    integer, intent(in) :: iter
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: psi, hpsi, old_psi

    type(OCEAN_vector) :: temp_psi
    real(DP) :: imag_a
    real(DP) :: time1, time2
    integer :: ialpha, ikpt
    
    imag_a = 0.0_DP  
    call cpu_time( time1 )

    do ialpha = 1, sys%nalpha
      do ikpt = 1, sys%nkpts
        a(iter-1) = a(iter-1) + dot_product( hpsi%r(:,ikpt,ialpha), psi%r(:,ikpt,ialpha) )&
                              + dot_product( hpsi%i(:,ikpt,ialpha), psi%i(:,ikpt,ialpha) )
        imag_a = imag_a + dot_product( hpsi%i(:,ikpt,ialpha), psi%r(:,ikpt,ialpha) )&
                        - dot_product( hpsi%r(:,ikpt,ialpha), psi%i(:,ikpt,ialpha) )
      enddo
    enddo

    hpsi%r( :, :, : ) = hpsi%r( :, :, : ) - a(iter-1) * psi%r( :, :, : ) - b(iter-1) * old_psi%r( :, :, : )
    hpsi%i( :, :, : ) = hpsi%i( :, :, : ) - a(iter-1) * psi%i( :, :, : ) - b(iter-1) * old_psi%i( :, :, : )

    do ialpha = 1, sys%nalpha
      do ikpt = 1, sys%nkpts
        b(iter) = b(iter) + sum( hpsi%r( :, ikpt, ialpha )**2 + hpsi%i( :, ikpt, ialpha )**2 )
      enddo
    enddo
    b(iter) = sqrt( b(iter) )

    hpsi%r( :, :, : ) = hpsi%r( :, :, : ) / b(iter)
    hpsi%i( :, :, : ) = hpsi%i( :, :, : ) / b(iter)


!    old_psi%r( :, :, : ) = psi%r( :, :, : )
!    old_psi%i( :, :, : ) = psi%i( :, :, : )
!    psi%r( :, :, : ) = hpsi%r( :, :, : )
!    psi%i( :, :, : ) = hpsi%i( :, :, : )
    temp_psi%r => old_psi%r
    temp_psi%i => old_psi%i
    old_psi%r => psi%r
    old_psi%i => psi%i
    psi%r => hpsi%r
    psi%i => hpsi%i
    hpsi%r => temp_psi%r
    hpsi%i => temp_psi%i

    call cpu_time( time2 )


    if( myid .eq. 0 ) then
      write ( 6, '(2x,2f10.6,10x,1e11.4,x,f6.3)' ) a(iter-1), b(iter), imag_a, time2-time1
      call haydump( iter, ierr )
    endif


  end subroutine OCEAN_hay_ab


  subroutine haydump( iter, ierr )
    use OCEAN_psi,  only : kpref
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
          disc = sqrt( ctmp ** 2 - 4 * b( iter ) ** 2 )
          di= -rm1 * disc
          if ( di .gt. 0.0d0 ) then
             delta = ( ctmp + disc ) / 2
          else
             delta = ( ctmp - disc ) / 2
          end if
          do jj = iter - 1, 0, -1
             delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
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

    integer :: dumi
    real :: dumf

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
          read(99,*) ne, el, eh, gam0, ebase
        case('inv')
          read(99,*) nloop, gres, gprc, ffff, ener
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


end module OCEAN_action
