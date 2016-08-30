! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_energies
  use AI_kinds
  use OCEAN_psi, only : OCEAN_vector
!  use iso_c_binding

  implicit none
  save
  private

  REAL(DP), ALLOCATABLE, public :: energies(:,:,:)
  real(DP), ALLOCATABLE, public :: imag_selfenergy(:,:,:)


  type( OCEAN_vector ) :: p_energy
  type( OCEAN_vector ) :: allow

#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: energies, imag_selfenergy
#endif



  INTEGER :: energy_bands_pad
  INTEGER :: energy_kpts_pad

  LOGICAL :: have_selfenergy
  LOGICAL :: is_init = .false.
  LOGICAL :: is_loaded = .false.
  LOGICAL :: val_init = .false.
  LOGICAL :: val_loaded = .false.

  public :: OCEAN_energies_val_allow, OCEAN_energies_val_sfact, OCEAN_energies_val_act, &
            OCEAN_energies_val_load, OCEAN_energies_act, OCEAN_energies_init, OCEAN_energies_load
  
  contains

  subroutine OCEAN_energies_val_allow( sys, psi, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_mult
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: psi
    !
    call OCEAN_psi_mult( psi, allow, .true. )
  end subroutine

  subroutine OCEAN_energies_val_sfact( sys, psi, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_mult
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: psi
    !
    call OCEAN_psi_mult( psi, allow, .false. )
  end subroutine

  subroutine OCEAN_energies_val_act( sys, psi, hpsi, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_cmult
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: hpsi
    !
    if( psi%val_myid .eq. 0 ) then
      call OCEAN_psi_cmult( psi, hpsi, p_energy, .false. )
    endif

  end subroutine OCEAN_energies_val_act


  subroutine OCEAN_energies_val_load( sys, ierr )
    use OCEAN_system
    use OCEAN_psi, only : OCEAN_vector, OCEAN_psi_new, OCEAN_psi_zero_full
    use OCEAN_val_energy, only : OCEAN_read_energies
    implicit none
    !
    integer, intent( inout ) :: ierr
    type(O_system), intent( in ) :: sys
    !
    if( sys%have_val .eq. .false. ) return
    !

    if( .not. val_init ) then
!      allocate( p_energy, allow )
      call OCEAN_psi_new( p_energy, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_new( allow, ierr )
      if( ierr .ne. 0 ) return

      val_init = .true.
    endif

    if( .not. val_loaded ) then
      call OCEAN_psi_zero_full( p_energy, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_zero_full( allow, ierr )
      if( ierr .ne. 0 ) return
      
      call OCEAN_read_energies( sys, p_energy, allow, ierr )
      if( ierr .ne. 0 ) return

      val_loaded = .true.
    endif

  end subroutine OCEAN_energies_val_load
    

  subroutine OCEAN_energies_init(  sys, ierr )
    use OCEAN_system
!    use OCEAN_mpi, only : myid, root

    implicit none
!    include 'fftw3.f03'
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys

    integer, parameter :: cacheline_by_Z = 1


!    if( associated( energies ) ) then
!      call fftw_free( rcptr )
!      energies => null()
!    endif
!    if( associated( imag_selfenergy ) ) then
!      call fftw_free( icptr )
!      imag_selfenergy => null()
!    endif
!JTV
    if( allocated( energies ) ) deallocate( energies )
    if( allocated( imag_selfenergy ) ) deallocate( imag_selfenergy )
   

    if( mod( sys%num_bands, cacheline_by_Z ) == 0 ) then
      energy_bands_pad = sys%num_bands
    else
      energy_bands_pad =  cacheline_by_Z * ( sys%num_bands / cacheline_by_Z + 1 )
    endif

    if( mod( sys%nkpts, cacheline_by_Z ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      energy_kpts_pad = sys%nkpts
    else
      energy_kpts_pad =  cacheline_by_Z * ( sys%nkpts / cacheline_by_Z + 1 )
    endif


    allocate( energies( energy_bands_pad, energy_kpts_pad, sys%nspn ), STAT=ierr )
!    rcptr = fftw_alloc_real( int( energy_bands_pad * energy_kpts_pad * sys%nspn, C_SIZE_T ) )
!    call c_f_pointer( rcptr, energies, [ energy_bands_pad, energy_kpts_pad, sys%nspn ] )
    energies = 0.0_DP


    allocate( imag_selfenergy( energy_bands_pad, energy_kpts_pad, sys%nspn ), STAT=ierr )
!    icptr = fftw_alloc_real( int( energy_bands_pad * energy_kpts_pad * sys%nspn, C_SIZE_T ) )
!    call c_f_pointer( icptr, imag_selfenergy, [ energy_bands_pad, energy_kpts_pad, sys%nspn ] )
    imag_selfenergy = 0.0_DP
    

  end subroutine OCEAN_energies_init


  subroutine OCEAN_energies_kill( ierr )
    implicit none
!    include 'fftw3.f03'
    integer, intent(inout) :: ierr

!    if( associated( energies ) ) then
!      call fftw_free( rcptr )
!      energies => null()
!    endif
!    if( associated( imag_selfenergy ) ) then
!      call fftw_free( icptr )
!      imag_selfenergy => null()
!    endif

    deallocate( energies )
    deallocate( imag_selfenergy )

  end subroutine OCEAN_energies_kill

  subroutine OCEAN_energies_load( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi
    use mpi

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys

    real(DP), allocatable :: tmp_e0(:,:,:)
    real(DP) :: core_offset
    integer :: nbd, nq, nspn, iter, i, j
    character(len=9) :: infoname
    character(len=4) :: gw_control
    logical :: file_exists, have_gw

    character(len=18) ::clsFile

    if( sys%conduct ) then
      infoname = 'wvfcninfo'
    else
      infoname = 'wvfvainfo'
    endif

    if( myid .eq. root ) then
      write(6,*) 'Loading energies ...'
      open( unit=99, file=infoname, form='unformatted', status='old' )
      rewind(99)
      if( sys%nspn .eq. 2 ) then
         read( 99 ) nbd, nq, nspn
      else
        read ( 99 ) nbd, nq
        nspn = 1
      endif
      if( nbd .ne. sys%num_bands ) then
        write(6,*) 'Band mismatch!', sys%num_bands, nbd
        ierr = 1
        goto 111
      elseif( nq .ne. sys%nkpts ) then
        write(6,*) 'K-point mismatch', sys%nkpts, nq
        ierr = 1
        goto 111
      elseif( nspn .ne. sys%nspn ) then
        write(6,*) 'Spin mismatch', sys%nspn, nspn
        ierr = 1
        goto 111
      endif

      allocate( tmp_e0( sys%num_bands,  sys%nkpts, sys%nspn ), STAT=ierr )
      if( ierr .ne. 0 ) return
      read(99) tmp_e0
      energies( 1 : sys%num_bands, 1 : sys%nkpts, : ) = tmp_e0( :, :, : )
      deallocate( tmp_e0 )
!      read(99) energies( 1 : sys%num_bands, 1 : sys%nkpts, : )
    

      close( 99 )


      inquire(file='gw_control',exist=have_gw)
      if( have_gw ) then
        open(unit=99,file='gw_control',form='formatted',status='old')
        rewind(99)
        read(99,*) gw_control
        close(99)
        select case (gw_control)
        case ('full')
          if( myid .eq. root ) write(6,*) 'full'
          call OCEAN_abinit_fullgw( sys, ierr, .true. )
        case ('real')
          call OCEAN_abinit_fullgw( sys, ierr, .false. )
        case( 'band' )
          call OCEAN_gw_by_band( sys, ierr, .false. )
        case( 'ibnd' )
          call OCEAN_gw_by_band( sys, ierr, .true. )
        case( 'strc' )
          call OCEAN_gw_stretch( sys, ierr )
        case default
          write(6,*) 'Unrecognized gw_control'
        end select
      endif

      write(clsFile,'(1A5,1A2,1I4.4,1A2,1I2.2,1A1,1I2.2)') 'cls.z', sys%cur_run%elname, sys%cur_run%indx, &
          '_n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)
      inquire(file=clsFile,exist=file_exists)
      if( file_exists ) then
        open(unit=99,file=clsFile,form='formatted',status='old')
        read(99,*) core_offset
        close(99)
        write(6,*) 'Core-level shift:', core_offset

!      if( sys%nruns .gt. 1 ) then
!        inquire(file='core_shift.txt',exist=file_exists)
!        core_offset = 0.0_DP
!        if( file_exists ) then
!          open(unit=99,file='core_shift.txt',form='formatted',status='old')
!          do iter = 1, sys%cur_run%indx
!            read(99,*) core_offset
!          enddo
!          close(99)
!          write(6,*) 'Core offset:', core_offset
          core_offset = core_offset / 27.2114d0
          energies(:,:,:) = energies(:,:,:) + core_offset
        endif
!      else
!        inquire(file='core_offset',exist=file_exists)
!        core_offset = 0.0_DP
!        if( file_exists ) then
!          open(unit=99,file='core_offset',form='formatted',status='old')
!          read(99,*) core_offset
!          write(6,*) 'Core offset:', core_offset
!          core_offset = core_offset / 27.2114d0
!          close(99)
!          energies(:,:,:) = energies(:,:,:) + core_offset
!        endif
!      endif

!      open(unit=99,file='energies.txt',form='formatted')
!      rewind(99)
!      do i = 1, sys%nkpts
!        do j = 1, sys%num_bands
!          write(99,*) j, energies(j,i,1)
!        enddo
!      enddo
!      close(99)

    endif
  ! Need to sychronize to test for ierr from above?
#ifdef MPI

    if( myid .eq. root ) write(6,*) energy_bands_pad,energy_kpts_pad,sys%nspn
    call MPI_BARRIER( comm, ierr )
    write(6,*) myid, sys%nspn, energy_bands_pad*energy_kpts_pad*sys%nspn
    call MPI_BARRIER( comm, ierr )

    call MPI_BCAST( energies, energy_bands_pad*energy_kpts_pad*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

  111 continue
  end subroutine OCEAN_energies_load


 !JTV need to move energies into an ocean_vector then:
 ! 1) this can be done w/ min instead of full
 ! 2) we can write an element-wise y(i) = z(i) * x(i) + y(i) in OCEAN_psi
  subroutine OCEAN_energies_act( sys, psi, hpsi, ierr )
    use OCEAN_system
    use OCEAN_psi

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: hpsi

    integer :: ialpha, ikpt, ibd, icms, icml, ivms, val_spin( sys%nalpha )

!    if( sys%nspn .ne. 1 ) ierr = -1 

    ! predefine the valence spins
    if( sys%nspn .eq. 1 ) then
      val_spin( : ) = 1
    else
      ialpha = 0
      do icms = 1, 2
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = 1, 2
            ialpha = ialpha + 1
            val_spin( ialpha ) = ivms
          enddo
        enddo
      enddo
    endif



    do ialpha = 1, sys%nalpha
!      hpsi%r( :, :, ialpha ) = energies( :, :, 1 ) * psi%r( :, :, ialpha ) &
!                             - imag_selfenergy( :, :, 1 ) * psi%i( :, :, ialpha )
!      hpsi%i( :, :, ialpha ) = energies( :, :, 1 ) * psi%i( :, :, ialpha ) &
!                             + imag_selfenergy( :, :, 1 ) * psi%r( :, :, ialpha )
      do ikpt = 1, sys%nkpts
        do ibd = 1, sys%num_bands
          hpsi%r( ibd, ikpt, ialpha ) = energies( ibd, ikpt, val_spin(ialpha) ) * psi%r( ibd, ikpt, ialpha ) &
                               - imag_selfenergy( ibd, ikpt, val_spin(ialpha) ) * psi%i( ibd, ikpt, ialpha )
        enddo
        do ibd = 1, sys%num_bands
          hpsi%i( ibd, ikpt, ialpha ) = energies( ibd, ikpt, val_spin(ialpha) ) * psi%i( ibd, ikpt, ialpha ) &
                               + imag_selfenergy( ibd, ikpt, val_spin(ialpha) ) * psi%r( ibd, ikpt, ialpha )

        enddo
      enddo 
    enddo

  end subroutine OCEAN_energies_act




  complex(DP) function OCEAN_energies_single( ib, ik, ia )
    use OCEAN_system
    implicit none
    integer, intent( in ) :: ib, ik, ia

    OCEAN_energies_single = CMPLX( energies( ib, ik, 1 ), imag_selfenergy( ib, ik, 1 ), DP )
  end function

  subroutine OCEAN_gw_stretch( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys

    real(DP) :: cstr
    logical :: have_gw

    if( myid .ne. root ) return

    write(6,*) 'Attempting GW stretch!'

    inquire(file='gwcstr', exist=have_gw )

    if( .not. have_gw ) then
      write( 6, * ) 'GW corrections requested (stretch style). File gwcstr not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    open( unit=99, file='gwcstr', form='formatted',status='old')
    rewind(99)
    read(99,*) cstr
    close(99)

    if( abs( cstr ) .lt. 0.00000001_DP ) return

    cstr = cstr + 1.0_DP
    energies( :, :, : ) = energies( :, :, : ) * cstr

  end subroutine OCEAN_gw_stretch

  subroutine OCEAN_gw_by_band( sys, ierr, keep_imag )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    logical, intent( in ) :: keep_imag

    integer :: iter, ispn, kiter
    real( DP ), allocatable :: re_se( : ), im_se( : )
    logical :: have_gw

    if( myid .ne. root ) return
    
    inquire(file='GW_band.in',exist=have_gw)
    if( .not. have_gw ) then
      write( 6, * ) 'GW corrections requested (band style). File GW_band.in not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    allocate( re_se( sys%num_bands ), im_se( sys%num_bands ) )
    open(unit=99,file='GW_band.in',form='formatted',status='old')
    rewind(99)

    if( keep_imag ) then
      do iter = 1, sys%num_bands
        read(99,*) re_se( iter ), im_se( iter )
      enddo
    else
      do iter = 1, sys%num_bands
        read(99,*) re_se( iter )
      enddo
      im_se( : ) = 0.0_DP
    endif
    close( 99 )
    re_se( : ) = re_se( : ) / 27.21138506_DP
    im_se( : ) = -im_se( : ) / ( 27.21138506_DP ) 


    do ispn = 1, sys%nspn
      do kiter = 1, sys%nkpts
        energies( :, kiter, ispn ) = energies( :, kiter, ispn ) + re_se( : )
        imag_selfenergy( :, kiter, ispn ) = im_se( : )
      enddo
    enddo

    deallocate( re_se, im_se )


  end subroutine OCEAN_gw_by_band

  subroutine OCEAN_abinit_fullgw( sys, ierr, keep_imag )
! Only run on root
    use OCEAN_system
    use OCEAN_mpi, only : myid, root

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    logical, intent( in ) :: keep_imag

    integer :: iter, biter, kiter, ix, iy, iz, r_kiter, ibd, ispn
    integer :: num_band, gw_nkpt, gw_nspn, start_band, stop_band, skip_band
    integer, allocatable :: kpt_map(:), start_b(:), stop_b(:)
    real(DP), allocatable :: re_se(:), im_se(:)
    real(DP) :: re_min, im_min, re_max, im_max, kpt( 3 ), e0
    logical :: have_kpt_map, have_gw

    if( myid .ne. root ) return

    inquire(file='kpt_map',exist=have_kpt_map)
    if( .not. have_kpt_map ) then
      write( 6, * ) 'GW corrections requested. File kpt_map not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    inquire(file='GWx_GW',exist=have_gw)
    if( .not. have_gw ) then 
      write( 6, * ) 'GW corrections requested. File GWx_GW not found.'
      write( 6, * ) 'No corrections will be done'
      return 
    endif

    allocate(kpt_map( sys%nkpts ) )
    open(unit=99,file='kpt_map',form='formatted',status='old')
    do iter = 1, sys%nkpts
      read(99,*) kpt_map( iter )
    enddo

    im_min = 0.0_DP
    re_min = 0.0_DP
    re_max = 0.0_DP
    im_max = 0.0_DP
    allocate( start_b( sys%nkpts ), stop_b(sys%nkpts ) )

    open(unit=99,file='GWx_GW',form='formatted',status='old')
    read(99,*) gw_nkpt, gw_nspn

!JTV Part of this is dumb because ABINIT won't warn us in the GW file if there is an extra element 
!    at some of the k-points
    do iter = 1, gw_nkpt
      read(99,*) kpt(1)
      read(99,*) num_band
      allocate( re_se( num_band ), im_se( num_band ) )
!JTV
      read(99,*) start_band, e0, re_se( 1 ), im_se( 1 )
! Format is band GW_energy, some_measure_of_conv, imaginary part
!      read(99,*) start_band, re_se( 1 ), e0, im_se( 1 )
      stop_band = start_band + num_band - 1

      do biter = 2, num_band
        read(99,*) ibd, e0, re_se( biter ), im_se( biter )
!        read(99,*) ibd, re_se( biter ), e0, im_se( biter )
      enddo
      re_se( : ) = re_se( : ) / 27.21138506_DP
!JTV
      im_se( : ) = - im_se( : ) / 27.21138506_DP
!      im_se( : ) = -( 2.0_DP * im_se( : ) ) / 27.21138506_DP
!      im_se( : ) = abs( im_se( : ) ) / ( 27.21138506_DP * 2.0_DP )
!      im_se( 1 ) = 0.0
!      im_se( 2 ) = 0.0

!JTV TEST
!      re_se( : ) = 0.0_DP
!      im_se( : ) = 0.1_DP / 27.21138506_DP
      
      skip_band = sys%cur_run%start_band - 1
      kiter = 0
      do ix = 1, sys%kmesh(1)
        do iy = 1, sys%kmesh(2)
          do iz = 1, sys%kmesh(3)
            kiter = kiter + 1
            r_kiter = ix + (iy-1) * sys%kmesh(1) + (iz-1)*sys%kmesh(1)*sys%kmesh(2)

            if( kpt_map( r_kiter ) .eq. iter ) then
!JTV need some spin fixing 
!              im_min = im_min + im_se( 1 ) 
!              re_min = re_min + ( re_se( 1 ) - energies( 1, kiter, 1 ) )
              im_max = im_max + im_se( num_band )
!              re_max = re_max + ( re_se( num_band ) - energies( stop_band - sys%cur_run%start_band + 1 , kiter, 1 ) )
              re_max = re_max + re_se( num_band )

              start_b( kiter ) = start_band
              stop_b( kiter ) = stop_band


! biter is actual counter of bands
! se_re/se_im goes from 1 to num_bands
              do biter = max( start_band, sys%cur_run%start_band ), &
                         min( stop_band, sys%cur_run%start_band + sys%num_bands - 1 )
                do ispn = 1, sys%nspn
!                  energies( biter - sys%cur_run%start_band + 1, kiter, ispn ) = re_se( biter - start_band + 1 )
                  energies( biter - sys%cur_run%start_band + 1, kiter, ispn ) =  &
                      energies( biter - sys%cur_run%start_band + 1, kiter, ispn ) + re_se( biter - start_band + 1 )
                  imag_selfenergy( biter - sys%cur_run%start_band + 1, kiter, ispn ) =  &
                        im_se( biter - start_band + 1 )
                  if( kiter .eq. 1 ) then
                    write(6,*) biter, biter - sys%cur_run%start_band + 1, biter - start_band + 1, im_se(  biter - start_band + 1 ) * 27.2114_DP
                  endif
                enddo
              enddo


!              do biter = 2 + skip_band - start_band, min( num_band, sys%num_bands )
!!!               if( biter + start_band - 1 .lt. 
!                do ispn = 1, sys%nspn
!                  if( biter == 1+ skip_band ) write(6,*) re_se(biter) * 27.21138506_DP
!                  energies( biter + start_band - 1 - skip_band, kiter, ispn ) = re_se( biter ) &
!                      + energies( biter + start_band - 1 - skip_band, kiter, ispn )
!                  imag_selfenergy( biter + start_band - 1 - skip_band, kiter, ispn ) = im_se( biter )
!                enddo
!              enddo
              
            endif

          enddo
        enddo
      enddo
      deallocate( re_se, im_se )
    enddo
!    return

    im_min = im_min / dble(sys%nkpts)
    re_min = re_min / dble(sys%nkpts)
    im_max = im_max / dble(sys%nkpts)
    re_max = re_max / dble(sys%nkpts)
    
    do ispn = 1, sys%nspn
      do kiter = 1, sys%nkpts
!        do biter = 1, start_b( kiter ) - sys%cur_run%start_band - 1
!!          write(6,*) kiter, start_b( kiter ), sys%cur_run%start_band
!          energies( biter, kiter, ispn ) = energies( biter, kiter, ispn ) + re_min
!          imag_selfenergy( biter, kiter, ispn ) = im_min
!        enddo
        do biter = stop_b( kiter ) + 1 - sys%cur_run%start_band, sys%num_bands
          if( kiter == 1 ) then
            write(6,*) kiter, stop_b( kiter ), sys%num_bands, biter, re_max * 27.2114_DP
          endif
          energies( biter, kiter, ispn ) = energies( biter, kiter, ispn ) + re_max
          imag_selfenergy( biter, kiter, ispn ) = im_max
        enddo
      enddo
    enddo

    deallocate( start_b, stop_b, kpt_map )

    if( .not. keep_imag ) then
      imag_selfenergy(:,:,:) = 0.0_DP
    endif 
!    imag_selfenergy(:,:,:) = -0.014_DP

  end subroutine OCEAN_abinit_fullgw
  
  subroutine OCEAN_energies_act_derp( sys, pi, pr, hpi, hpr, ener, se, ierr )
    use OCEAN_system
    use OCEAN_psi 

    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    real(DP), intent(in), dimension(sys%num_bands,sys%nkpts,sys%nalpha) :: pi, pr
    real(DP), intent(in), dimension(energy_bands_pad,energy_kpts_pad,1) :: ener, se
    real(DP), intent(inout), dimension(sys%num_bands,sys%nkpts,sys%nalpha) :: hpi, hpr

    integer :: ialpha, ikpt, ibd

    if( sys%nspn .ne. 1 ) ierr = -1 

    do ialpha = 1, sys%nalpha
      do ikpt = 1, sys%nkpts
        do ibd = 1, sys%num_bands
          hpr( ibd, ikpt, ialpha  ) = ener( ibd, ikpt, 1 ) * pr( ibd, ikpt, ialpha ) &
                                    - se(   ibd, ikpt, 1 ) * pi( ibd, ikpt, ialpha )

        enddo
      enddo 
    enddo

  end subroutine OCEAN_energies_act_derp

end module OCEAN_energies
