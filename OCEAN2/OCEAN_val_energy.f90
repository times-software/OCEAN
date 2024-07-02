! Copyright (C) 2016 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_val_energy
  use AI_kinds
  implicit none

  contains

  subroutine OCEAN_read_energies( sys, p_energy, allow, ierr )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    use OCEAN_constants, only : Hartree2eV, eV2Hartree
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( inout ) :: p_energy, allow
    integer, intent( inout ) :: ierr
    !


    real( DP ), allocatable, dimension(:,:,:,:) :: val_energies, con_energies, im_val_energies, im_con_energies, tmp_e
    real( DP ), allocatable :: val_buffer( : ), con_buffer( : )
    real( DP ) :: efermi, homo, lumo, cliph
    integer :: ik, ibv, ibc, fh, ispn, jspn, ibeta, i, j
    logical :: metal, have_imaginary, did_gw_correction
    integer :: nbv, nbc(2), nk, nbw, ibw
#ifdef MPI
    integer(MPI_OFFSET_KIND) :: offset
#endif

!    real(DP), parameter :: Ha_to_eV = 27.21138386_dp

    nbw = 1
    if( sys%bwflg ) nbw = 2

!    allocate( val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ),  &
    allocate( val_energies( sys%brange(1) : sys%brange(2), sys%nkpts, sys%nspn, nbw ),  &
              con_energies( sys%brange(3) : sys%brange(4), sys%nkpts, sys%nspn, nbw ), STAT=ierr )
    val_energies(:,:,:,:) = 0.0_DP 
    con_energies(:,:,:,:) = 0.0_DP 
    if( ierr .ne. 0 ) return


! Read in energies
       select case( sys%enk_selector )

       case( 0 )
          
          if( myid .eq. root ) then
            allocate( val_buffer( sys%file_brange(1) : sys%file_brange(2) ), &
                      con_buffer( sys%file_brange(3) : sys%file_brange(4) ) )
            open(unit=99,file='enkfile',form='formatted',status='old')
            do ispn=1, sys%nspn
              do ik = 1, sys%nkpts
                read(99,*) val_buffer( : )
                read(99,*) con_buffer( : )

                val_energies( sys%brange(1) : sys%brange(2), ik, ispn, 1 ) = val_buffer( sys%brange(1) : sys%brange(2) ) / 2.0_DP
                con_energies( sys%brange(3) : sys%brange(4), ik, ispn, 1 ) = con_buffer( sys%brange(3) : sys%brange(4) ) / 2.0_DP
                if( sys%bwflg ) then
                  val_energies( sys%brange(1) : sys%brange(2), ik, ispn, 2 ) = con_buffer( sys%brange(1) : sys%brange(2) ) / 2.0_DP
                  con_energies( sys%brange(3) : sys%brange(4), ik, ispn, 2 ) = val_buffer( sys%brange(3) : sys%brange(4) ) / 2.0_DP
                endif
!                   read(99,*) val_energies( :, ik ,ispn)
!                   read(99,*) con_energies( :, ik ,ispn )
!                enddo
!                val_energies( :, : ,ispn ) = val_energies( :, : ,ispn ) / 2.0_dp ! * Hartree2eV !Ha_to_eV
!                con_energies( :, : ,ispn ) = con_energies( :, : ,ispn ) / 2.0_dp ! * Hartree2eV !Ha_to_eV
              enddo
            enddo  !ispn
          endif
#ifdef MPI
          call MPI_BCAST( val_energies, (sys%brange(2)-sys%brange(1)+1)*sys%nkpts*sys%nspn*nbw, MPI_DOUBLE_PRECISION, root, comm, ierr )
!          call MPI_BCAST( val_energies, sys%cur_run%val_bands*sys%nkpts*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
          call MPI_BCAST( con_energies, (sys%brange(4)-sys%brange(3)+1)*sys%nkpts*sys%nspn*nbw, MPI_DOUBLE_PRECISION, root, comm, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
#endif
          
       case( 1 )
#if 0          
#ifdef MPI
          if( myid .eq. root ) then
             open(unit=99,file='tmels.info',form='formatted',status='old')
             read(99,*) nbv, nbc(1), nbc(2), nk
             close(99)
             
             if( nk .ne. sys%nkpts .and. nk .ne. sys%nkpts*2) then
                ierr = -1
                write(6,*) 'tmels.info not compatible with other run information: NKPTS'
                return
             endif
          endif
          
          call MPI_BCAST( nbc, 2, MPI_INTEGER, 0, comm, ierr )
          if( ierr .ne. 0 ) return
          call MPI_BCAST( nbv, 1, MPI_INTEGER, 0, comm, ierr )
          if( ierr .ne. 0 ) return
          
          call MPI_FILE_OPEN( comm, 'val_energies.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
          offset = 0
          call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL,ierr)
          if( ierr .ne. MPI_SUCCESS ) return
          do ispn=1, sys%nspn
             do ik = 1, sys%nkpts
                offset = ( ispn - 1 )* sys%nkpts * nbv + ( ik - 1 ) * nbv
                call MPI_FILE_READ_AT_ALL( fh, offset, val_energies( 1, ik , ispn), sys%cur_run%val_bands, MPI_DOUBLE_PRECISION, &
                     MPI_STATUS_IGNORE, ierr )
                if( ierr .ne. MPI_SUCCESS) return
             enddo
          enddo  !ispn
          call MPI_FILE_CLOSE( fh, ierr )
          if( ierr .ne. MPI_SUCCESS) return
          
          call MPI_FILE_OPEN( comm, 'con_energies.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
          offset = 0
          call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL,ierr)
          if( ierr .ne. MPI_SUCCESS ) return
          do ispn = 1, sys%nspn
             do ik = 1, sys%nkpts
                offset = ( ispn - 1 ) * sys%nkpts * ( nbc( 2 ) - nbc( 1 ) + 1 ) &
                       + ( ik - 1 ) * ( nbc( 2 ) - nbc( 1 ) + 1 ) + ( sys%brange( 3 ) - nbc( 1 ) )
                call MPI_FILE_READ_AT_ALL( fh, offset, con_energies( 1, ik, ispn ), sys%cur_run%num_bands, MPI_DOUBLE_PRECISION, &
                     MPI_STATUS_IGNORE, ierr )
                if( ierr .ne. MPI_SUCCESS) return
             enddo
          enddo  !ispn
          call MPI_FILE_CLOSE( fh, ierr )
          if( ierr .ne. MPI_SUCCESS) return
#else
          ierr = -1
          if( myid .eq. root ) write(6,*) 'MPI required for OBF-style!'
          return
#endif

       case( 2 )

          if( myid .eq. root ) then
            open(unit=99,file='wvfvainfo',form='unformatted',status='old')
            read(99) ibc, ik, ispn
            if( ibc .lt. sys%brange(2) .or. ik .ne. sys%nkpts .or. ispn .ne. sys%nspn ) then
              write( 6, * ) 'Problem with size of wvfvainfo!!!'
              write( 6, * ) 'Expected: ', sys%brange(2), sys%nkpts, sys%nspn
              write( 6, * ) 'Found   : ', ibc, ik, ispn
              ierr = 129412
              close( 99 )
              return
            endif
            ! Allow the number of bands to be too big
            if( ibc .ne. sys%brange(2) ) then
              allocate( tmp_e( ibc, sys%nkpts, sys%nspn ) )
              read(99) tmp_e
              val_energies( 1 : sys%brange(2), :, : ) = tmp_e( 1 : sys%brange(2), :, : )
              deallocate( tmp_e )
            else
              read(99) val_energies
            endif
            close( 99 )

            open(unit=99,file='wvfcninfo',form='unformatted',status='old')
            read(99) ibc, ik, ispn
            if( ibc .lt. sys%cur_run%num_bands .or. ik .ne. sys%nkpts .or. ispn .ne. sys%nspn ) then
              write( 6, * ) 'Problem with size of wvfcninfo!!!'
              write( 6, * ) 'Expected: ', sys%cur_run%num_bands, sys%nkpts, sys%nspn
              write( 6, * ) 'Found   : ', ibc, ik, ispn
              ierr = 129413
              close( 99 )
              return
            endif
            ! Allow the number of bands to be too big
            if( ibc .ne. sys%cur_run%num_bands ) then
              allocate( tmp_e( ibc, sys%nkpts, sys%nspn ) )
              read(99) tmp_e
              con_energies( 1 : sys%cur_run%num_bands, :, : ) = tmp_e( 1 : sys%cur_run%num_bands, :, : ) 
              deallocate( tmp_e ) 
            else
              read(99) con_energies
            endif
            close( 99 )
          endif
#ifdef MPI
          call MPI_BCAST( val_energies, sys%brange(2)*sys%nkpts*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
!          call MPI_BCAST( val_energies, sys%cur_run%val_bands*sys%nkpts*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
          call MPI_BCAST( con_energies, sys%cur_run%num_bands*sys%nkpts*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
#endif


#endif
          
       case default
          ierr = -2
          if( myid .eq. root ) write(6,*) 'Un-supported enk_selector:', sys%enk_selector
          return
          
       end select
    
    if( myid .eq. root ) then
      open( unit=99,file='val_energy_test.txt', form='formatted',status='unknown' )
      write(99,*) val_energies(:,:,:,:)
      close(99)
      write(6,*) val_energies(:,1,1,1)
    endif
!   call GW corrections time

!    allocate( im_val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ), & 
    allocate( im_val_energies( sys%brange(2), sys%nkpts, sys%nspn, nbw ), & 
              im_con_energies( sys%cur_run%num_bands, sys%nkpts, sys%nspn, nbw ), STAT=ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. 0 .and. abs( nint(sys%nelectron) - sys%nelectron ) .gt. 0.001_DP ) then
      write(6,*) 'WARNING: Fractional electron count not well handled by valence code (yet)', &
                  sys%nelectron, nint(sys%nelectron)
    endif
    call find_fermi( sys, val_energies, con_energies, nint(sys%nelectron), efermi, &
                     homo, lumo, cliph, metal, .true., ierr )
    if( ierr .ne. 0 ) return

    
    if( myid .eq. 0 ) write(6,*) val_energies(:,1,1,1)
    call val_gw( sys, homo, lumo, val_energies, con_energies, have_imaginary, im_val_energies, &
                 im_con_energies, did_gw_correction, ierr )
    if( ierr .ne. 0 ) return

    if( did_gw_correction ) then 
      call find_fermi( sys, val_energies, con_energies, nint(sys%nelectron), efermi, &
                       homo, lumo, cliph, metal, .false., ierr )
      if( ierr .ne. 0 ) return
    endif

    if( myid .eq. 0 ) write(6,*) val_energies(:,1,1,1)
! THIS SHOULD BE OK NOW
#if 0 
    ! move this nonsense back
    if( sys%cur_run%val_bands .ne. sys%brange(2) ) then
      allocate(tmp_e( sys%cur_run%val_bands, sys%nkpts, sys%nspn ) )
      tmp_e( 1:sys%cur_run%val_bands, :, : ) = val_energies( sys%brange(1):sys%brange(2), :, : )
      deallocate( val_energies )
      allocate( val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ) )
      val_energies( :, : ,: ) = tmp_e( :, :, : )
      if( have_imaginary ) then
!        write( 6, * ) sys%brange(1), sys%brange(2)
!        write(6,*) im_val_energies( 1:sys%brange(1), 1, 1 ) * Hartree2eV
!        write(6,*) im_val_energies( sys%brange(1):sys%brange(2), 1, 1)*Hartree2eV
        tmp_e( 1:sys%cur_run%val_bands, :, : ) = im_val_energies( sys%brange(1):sys%brange(2), :, : )
        deallocate( im_val_energies )
        allocate( im_val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ) )
        im_val_energies( :, : ,: ) = tmp_e( :, :, : )
      endif
      deallocate( tmp_e )
    endif
#endif
    call energies_allow( sys, val_energies, con_energies, nint(sys%nelectron), efermi, cliph, &
                                allow, metal, ierr )
    if( ierr .ne. 0 ) return


    do ibw = 1, sys%nbw
      ibeta=0
      do i=1,sys%valence_ham_spin
        ispn = min( i, sys%nspn ) 
        do j=1,sys%valence_ham_spin
          jspn = min( j, sys%nspn )
          ibeta=ibeta+1
          do ik = 1, sys%nkpts
            do ibv = 1, sys%cur_run%val_bands
              do ibc = 1, sys%cur_run%num_bands
                if( ibw .eq. 1 ) then
                p_energy%valr( ibc, ibv, ik, ibeta, ibw ) = con_energies( ibc+sys%brange(3)-1, ik, jspn, ibw ) - val_energies( ibv+sys%brange(1)-1, ik, ispn, ibw )
                else
                p_energy%valr( ibc, ibv, ik, ibeta, ibw ) = -con_energies( ibc+sys%brange(3)-1, ik, jspn, ibw ) + val_energies( ibv+sys%brange(1)-1, ik, ispn, ibw )
                endif
              enddo
            enddo
          enddo
        enddo  !jspn
      enddo !ispn
    enddo  
    if( have_imaginary ) then
      if( myid .eq. 0 ) write(6,*) 'HAVE IMAG'
      do ibw = 1, sys%nbw 
        ibeta=0
        do i=1,sys%valence_ham_spin
          ispn = min( i, sys%nspn )
          do j=1,sys%valence_ham_spin
            jspn = min( j, sys%nspn )
            ibeta=ibeta+1
            do ik = 1, sys%nkpts
              do ibv = 1, sys%cur_run%val_bands
                do ibc = 1, sys%cur_run%num_bands
                  p_energy%vali( ibc, ibv, ik, ibeta, ibw ) = -(im_con_energies( ibc, ik, jspn, ibw ) - im_val_energies( ibv, ik, ispn, ibw ) )
!                      if( p_energy%vali( ibc, ibv, ik, ibeta ) .gt. -0.0001_dp ) then
                  if( myid .eq. 0 .and. ik .eq. 1 ) then
                    write(80,'(3(I8,X),3(F24.6))') ik, ibv, ibc, Hartree2eV*p_energy%vali( ibc, ibv, ik, ibeta, ibw ), &
                         Hartree2eV*im_val_energies( ibv, ik, ispn, ibw ), HArtree2eV*im_con_energies( ibc, ik, jspn, ibw)
                  endif
                enddo
              enddo
            enddo
          enddo !jspn
        enddo  !ispn
      enddo ! ibw
    else
      p_energy%vali( :, :, :, :, : ) = 0.0_dp
    endif


    deallocate( val_energies, con_energies, im_val_energies, im_con_energies )
    
    if( myid .eq. 0 ) write(6,*) 'Done with energies'

    !TODO
    ! Nominally, in the previous every process should have gotten the same set of information. 
    ! This appears to not quite be correct
    call OCEAN_psi_bcast_full( root, p_energy, ierr )
    call OCEAN_psi_bcast_full( root, allow, ierr )


  end subroutine OCEAN_read_energies

  ! Master routine for all the gw flavors (for the valence band). 
  !   For now this is just a stub for what I currently need.
  !   In the future need to better merge the valence and X-ray versions of everything
  subroutine val_gw( sys, homo, lumo, val_energies, con_energies, have_imaginary, im_val_energies, &
                     im_con_energies, did_gw_correction, ierr )
    use OCEAN_system
    use OCEAN_mpi
    implicit none
    type( O_system ), intent( in ) :: sys
    real( DP ), intent( in ) :: homo, lumo
!    real( DP ), intent( inout ), dimension( sys%cur_run%val_bands, sys%nkpts,sys%nspn ) :: &
    real( DP ), intent( inout ), dimension( sys%brange(1) : sys%brange(2), sys%nkpts, sys%nspn, sys%nbw ) :: &
        val_energies, im_val_energies
    real( DP ), intent( inout ), dimension( sys%brange(3) : sys%brange(4), sys%nkpts, sys%nspn, sys%nbw ) :: &
        con_energies, im_con_energies

    logical, intent( out ) :: have_imaginary, did_gw_correction
    integer, intent( inout ) :: ierr
    !
    logical :: have_gw
    character( len=4 ) :: gw_control

    did_gw_correction = .false.
    have_imaginary = .false.

    if( myid .eq. root ) then
      inquire( file="gw_control", exist=have_gw )
      if( have_gw ) then
        open(unit=99,file='gw_control',form='formatted',status='old')
        rewind(99)
        read(99,*) gw_control
        close(99)
        select case (gw_control)
!        case ('full')
!          call val_abinit_fullgw( sys, ierr, .true. )
        case ('list')
          write(6,*) 'GW! Will attempt list-style corrections'
          call val_list_gw( sys, val_energies, con_energies, ierr )
        case( 'band' )
          call val_gw_by_band( sys, val_energies, con_energies, ierr, .false. )
        case( 'ibnd' )
          call val_gw_by_band( sys, val_energies, con_energies, ierr, .true., &
                               im_val_energies, im_con_energies )
          have_imaginary = .true.
        case( 'cstr' )
          call val_gw_stretch( sys, homo, lumo, val_energies, con_energies, ierr )
        case default
          write(6,*) 'Unrecognized gw_control:'
          write(6,*) '   ', gw_control
          have_gw = .false.
        end select
      endif
    endif
#ifdef MPI
    call MPI_BCAST( have_gw, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    if( have_gw ) then
!      call MPI_BCAST( val_energies, sys%cur_run%val_bands * sys%nkpts, MPI_DOUBLE_PRECISION, & 
      call MPI_BCAST( val_energies, (sys%brange(2)-sys%brange(1)) * sys%nkpts * sys%nspn * sys%nbw, &
                      MPI_DOUBLE_PRECISION, root, comm, ierr )
      call MPI_BCAST( con_energies, (sys%brange(4)-sys%brange(3)) * sys%nkpts * sys%nspn * sys%nbw, &
                      MPI_DOUBLE_PRECISION, root, comm, ierr )
      call MPI_BCAST( have_imaginary, 1, MPI_LOGICAL, root, comm, ierr )
      call MPI_BCAST( im_val_energies, (sys%brange(2)-sys%brange(1)) * sys%nkpts * sys%nspn * sys%nbw, &
                      MPI_DOUBLE_PRECISION, root, comm, ierr )
      call MPI_BCAST( im_con_energies, (sys%brange(4)-sys%brange(3)) * sys%nkpts * sys%nspn * sys%nbw, &
                      MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif
#endif
    did_gw_correction = have_gw

  end subroutine val_gw

  subroutine val_gw_by_band( sys, val_energies, con_energies, ierr, keep_imag, im_val_energies, im_con_energies )
    use OCEAN_system
    use OCEAN_constants, only : eV2Hartree
    implicit none
    !
    type( O_system ), intent( in ) :: sys
!    real( DP ), intent( inout ), dimension( sys%cur_run%val_bands, sys%nkpts, sys%nspn ) :: val_energies
!    real( DP ), intent( inout ), dimension( sys%cur_run%num_bands, sys%nkpts, sys%nspn ) :: con_energies
!    real( DP ), intent( inout ), dimension( :, :, : ) :: val_energies, con_energies
    real( DP ), intent( inout ), optional, dimension( sys%brange(1) : sys%brange(2), sys%nkpts, sys%nspn, sys%nbw ) :: val_energies
    real( DP ), intent( inout ), optional, dimension( sys%brange(3) : sys%brange(4), sys%nkpts, sys%nspn, sys%nbw ) :: con_energies
    integer, intent( inout ) :: ierr
    logical, intent( in ) :: keep_imag
    real( DP ), intent( inout ), optional, dimension( sys%brange(1) : sys%brange(2), sys%nkpts, sys%nspn, sys%nbw ) :: im_val_energies
    real( DP ), intent( inout ), optional, dimension( sys%brange(3) : sys%brange(4), sys%nkpts, sys%nspn, sys%nbw ) :: im_con_energies
    !
    integer :: nbands, iter, ispn, ikpt, ibw
    real( DP ), allocatable, dimension(:) :: re_se, im_se
    logical :: have_gw
    !
    inquire(file='GW_band.in',exist=have_gw)
    if( .not. have_gw ) then
      write( 6, * ) 'GW corrections requested (band style). File GW_band.in not found.'
      write( 6, * ) 'No corrections will be done'
      return
    endif

    open(unit=99,file='GW_band.in',form='formatted',status='old')
    rewind(99)
    read(99,*) nbands
    allocate( re_se( nbands ), im_se( nbands ) )
    if( keep_imag ) then
      write(6,*) 'Reading IMAG'
      do iter = 1, nbands
        read(99,*) re_se( iter ), im_se( iter )
      enddo
    else
      do iter = 1, nbands
        read(99,*) re_se( iter )
      enddo
      im_se( : ) = 0.0_DP
    endif
    close( 99 )

    re_se( : ) = re_se( : ) * eV2Hartree
    im_se( : ) = -im_se( : ) * eV2Hartree

    do ibw = 1, sys%nbw
      do ispn = 1, sys%nspn
        do ikpt = 1, sys%nkpts
          do iter = sys%brange(1), min( sys%brange(2), nbands )
            val_energies( iter, ikpt, ispn, ibw ) = val_energies( iter, ikpt, ispn, ibw ) + re_se( iter )
          enddo
          do iter = sys%brange(3), min( sys%brange(4), nbands )
            con_energies( iter, ikpt, ispn, ibw ) = con_energies( iter, ikpt, ispn, ibw ) + re_se( iter )
          enddo
        enddo
      enddo
    enddo

    if( keep_imag ) then
      if( present( im_val_energies ) ) then
        write(6,*) 'Store VAL IMAG'
        do ibw = 1, sys%nbw
          do ispn = 1, sys%nspn
            do ikpt = 1, sys%nkpts
              do iter = sys%brange(1), min( nbands, sys%brange(2) )
                im_val_energies ( iter, ikpt, ispn, ibw ) = im_se( iter )
              enddo   
            enddo
          enddo
        enddo
      endif
        if( present( im_con_energies ) ) then
        write(6,*) 'Store CON IMAG'
        do ibw = 1, sys%nbw
          do ispn = 1, sys%nspn
            do ikpt = 1, sys%nkpts
              do iter = sys%brange(3), min( sys%brange(4), nbands )
                im_con_energies( iter, ikpt, ispn, ibw ) = im_se( iter )
              enddo
            enddo
          enddo
        enddo
      endif
    endif

    deallocate( re_se, im_se )

  
  end subroutine val_gw_by_band


  subroutine val_gw_stretch( sys, homo, lumo, val_energies, con_energies, ierr )
    use OCEAN_system
    use OCEAN_constants, only : eV2Hartree
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    real( DP ), intent( in ) :: homo, lumo
    real( DP ), intent( inout ), dimension( sys%brange(1):sys%brange(2), sys%nkpts, sys%nspn, sys%nbw ) :: val_energies
    real( DP ), intent( inout ), dimension( sys%brange(3):sys%brange(4), sys%nkpts, sys%nspn, sys%nbw ) ::  con_energies
    integer, intent( inout ) :: ierr
    !
    real(dp) :: gw_gap_correction, vstr, cstr !stretch
    logical :: abs_gap


#if 0
    open(unit=99,file='gwgap',form='formatted', status='old' )
    read( 99, * ) gw_gap_correction, abs_gap
    close( 99 ) 

    if( abs_gap ) then
      gw_gap_correction = gw_gap_correction * eV2Hartree + homo
    else
      gw_gap_correction = gw_gap_correction * eV2Hartree + lumo
    endif

    open(unit=99,file='gwcstr',form='formatted', status='old' )
    read( 99, * ) stretch
    close( 99 )
    stretch = stretch + 1.0_dp

    con_energies( :, :, : ) = gw_gap_correction + ( con_energies( :, :, : ) - lumo ) * stretch

    open(unit=99,file='gwvstr',form='formatted', status='old' )
    read( 99, * ) stretch
    close( 99 )
    stretch = stretch + 1.0_dp

    val_energies( :, :, : ) = homo + ( val_energies( :, :, : ) - homo ) * stretch
#endif

    open(unit=99,file='gw_val_cstr',form='formatted', status='old' )
    read( 99, * ) gw_gap_correction, abs_gap, vstr, cstr
    close( 99 )

    if( abs_gap ) then
      gw_gap_correction = gw_gap_correction * eV2Hartree + homo
    else
      gw_gap_correction = gw_gap_correction * eV2Hartree + lumo
    endif

    vstr = vstr + 1.0_DP
    val_energies( :, :, :, : ) = homo + ( val_energies( :, :, :, : ) - homo ) * vstr

    cstr = cstr + 1.0_DP
    con_energies( :, :, :, : ) = gw_gap_correction + ( con_energies( :, :, :,: ) - lumo ) * cstr

  end subroutine val_gw_stretch

  subroutine val_list_gw( sys, val_energies, con_energies, ierr )
    use OCEAN_system
    use OCEAN_constants, only : eV2Hartree
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    real( DP ), intent( inout ), dimension( sys%cur_run%val_bands, sys%nkpts ) :: val_energies
    real( DP ), intent( inout ), dimension( sys%cur_run%num_bands, sys%nkpts ) :: con_energies
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: delta_gw
    integer :: band_max, band_loop, ikpt, iband
    logical :: have_gw
    !
    inquire( file="list_val_gw.txt", exist=have_gw )
    if( have_gw ) then
      !
      open( unit=99, file="list_val_gw.txt", form="formatted", status="old" )
      read( 99, * ) band_max
      if( band_max .lt. sys%cur_run%val_bands ) then
        write( 6, * ) 'WARNING: Bands in list_val_gw.txt is less than the run'
        write( 6, * ) band_max, sys%cur_run%val_bands
        write( 6, * ) 'This is almost certainly not what you want, but will continue'
      endif
      
      ! These will be in eV
      band_loop = min( band_max, sys%cur_run%val_bands )
      do iband = 1, band_loop
        do ikpt = 1, sys%nkpts
          read( 99, * ) delta_gw
          val_energies( iband, ikpt ) = val_energies( iband, ikpt ) + delta_gw * eV2Hartree
        enddo
      enddo
      write(6,*) 'last gw value for valence: ', delta_gw
    endif

    delta_gw = 1.0_dp 
    write( 6, * ) 'WARNING WARNING WARNING'
    write( 6, * ) 'Hardwired conduction band GW correction of: ', delta_gw, ' eV'
    delta_gw = delta_gw * eV2Hartree
    do ikpt = 1, sys%nkpts
      do iband = 1, sys%cur_run%num_bands
        con_energies( iband, ikpt ) = con_energies( iband, ikpt ) + delta_gw
      enddo
    enddo
    

  end subroutine val_list_gw

  subroutine energies_allow( sys, val_energies, con_energies, nelectron, efermi, cliph_, &
                                  allow, metal, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_psi
    implicit none
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( inout ) :: allow
    integer, intent( in ) :: nelectron
    real( DP ), intent( in ) :: con_energies( 1:sys%cur_run%num_bands, sys%nkpts, sys%nspn, sys%nbw ), &
                                val_energies( 1:sys%cur_run%val_bands, sys%nkpts, sys%nspn, sys%nbw ),  &
                                efermi, cliph_
    logical, intent( in ) :: metal
    integer, intent( inout ) :: ierr
    !
    integer :: kiter, biter1, biter2, ibeta, i, ispn, j, jspn, ibw
    logical :: have_clipl, have_cliph
    real(dp) :: clipl, cliph
    !
    !
    ! Already zero'd
    allow%valr = 0.0_dp
    allow%vali = 0.0_dp

    cliph = cliph_

    if( myid .eq. root ) then
      inquire(file='clipl',exist=have_clipl )
      clipl = val_energies( 1, 1, 1, 1 ) - 10000.0_dp
      if( have_clipl ) then
        open(unit=99,file='clipl',form='formatted',status='old' )
        read(99,*) clipl
        close(99)
        write(6,*) 'clipl: ', clipl, val_energies( 1, 1, 1, 1 )
      endif
    endif
    call MPI_BCAST( have_clipl, 1, MPI_LOGICAL, root, comm, ierr )
    call MPI_BCAST( clipl, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )


    if( myid .eq. root ) then
      inquire(file='cliph',exist=have_cliph)
      if(have_cliph) then
        open(unit=99,file='cliph',form='formatted',status='old' )
        read(99,*) cliph
        close(99)
      endif
    endif
    call MPI_BCAST( cliph, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )


    if( myid .eq. root ) write(6,*) 'clipl: ', clipl, val_energies( 1, 1, 1,1 )
    if( myid .eq. root ) write(6,*) '  BWFLG', sys%cur_run%bwflg, sys%cur_run%val_bands, sys%cur_run%num_bands
    if( myid .eq. root .and. sys%disable_intraband ) write(6,*) '  NO INTRABAND'


    !
    ! If we have spins from the DFT then we can't rely on band index, need to use Fermi level 
    !  just like for the metals case. 
    !
    ! In storing psi the conduction band index is the fast index
    if( metal .or. sys%nspn .ne. 1 ) then
    
      do ibw = 1, sys%nbw
      ibeta = 0
      do i = 1, sys%valence_ham_spin
        ispn = min( i, sys%nspn )
        do j = 1, sys%valence_ham_spin
          jspn = min( j, sys%nspn )
          ibeta = ibeta + 1
          do kiter = 1, sys%nkpts
            do biter2 = 1, sys%cur_run%val_bands
!              if( kiter .eq. 1 .and. myid .eq. 0 ) write(6,*) 'clipl: ', clipl, val_energies( biter2, kiter, ispn, ibw )
              if ( val_energies( biter2, kiter, ispn, ibw ) .lt. efermi  .and. &
                    ( val_energies( biter2, kiter, ispn, ibw ) .ge. clipl )) then

                do biter1 = 1, sys%cur_run%num_bands
                  if ( ( con_energies( biter1, kiter, jspn, ibw ) .gt. efermi ) .and. &
                       ( con_energies( biter1, kiter, jspn, ibw ) .le. cliph ) ) then
                    if( sys%disable_intraband .and. &
                       ( biter2+sys%brange(1) .eq. biter1+sys%brange(3) ) ) then    
                      ! above: if we are disabling intraband AND the bands are equal
                      !  then we skip this section setting allow
                      if( myid .eq. root ) write(6,*) 'SKIP', kiter, biter2, biter1
                    else
                      allow%valr( biter1, biter2, kiter, ibeta, ibw ) = 1.0_dp
                      if( ibw .eq. 1 ) then
                        allow%vali( biter1, biter2, kiter, ibeta, ibw ) = 1.0_dp
                      else
                        allow%vali( biter1, biter2, kiter, ibeta, ibw ) = -1.0_dp
                      endif
                    endif
                  endif
                enddo 
!              elseif( sys%cur_run%bwflg .and. &
!                    ( val_energies( biter2, kiter, ispn ) .gt. efermi ) .and. &
!                    ( val_energies( biter2, kiter, ispn ) .le. cliph )) then
!
!                do biter1 = 1, sys%cur_run%num_bands
!                  if ( ( con_energies( biter1, kiter, jspn ) .lt. efermi ) .and. &
!                       ( con_energies( biter1, kiter, jspn ) .ge. clipl ) ) then
!                    allow%valr( biter1, biter2, kiter, ibeta ) = 1.0_dp
!                    allow%vali( biter1, biter2, kiter, ibeta ) = -1.0_dp
!                  endif
!                enddo
!              elseif (  kiter .eq. 1 ) then
!                write(6,*) ibeta, kiter, biter2, val_energies( biter2, kiter, ispn )
              endif
            enddo ! biter1
          enddo ! kiter
        enddo ! j
      enddo ! i
      enddo ! ibw
    else ! not metal && sys%nspn==1
      if( sys%backf ) ierr = 413
      do ibw = 1, sys%nbw
        do ibeta = 1, sys%nbeta
          do kiter = 1, sys%nkpts
            do biter2 = 1, ( nelectron / 2 ) - sys%brange( 1 ) + 1
  !            if( kiter .eq. 1 .and. myid .eq. 0 ) write(6,*) 'clipl: ', clipl, val_energies( biter2, kiter, 1 )
  !            if( (.not. have_clipl) .or.  val_energies( biter2, kiter, 1 ) .gt. clipl ) then
              if( val_energies( biter2, kiter, 1, ibw ) .ge. clipl ) then
                do biter1 = 2 - sys%brange( 3 ) + ( nelectron / 2 ), sys%cur_run%num_bands
                  if(  con_energies( biter1, kiter, 1, ibw ) .le. cliph ) then
                    allow%valr( biter1, biter2, kiter, ibeta, ibw ) = 1.0_dp
                    if( ibw .eq. 1 ) then
                      allow%vali( biter1, biter2, kiter, ibeta, ibw ) = 1.0_dp
                    else
                      allow%vali( biter1, biter2, kiter, ibeta, ibw ) = -1.0_dp
                    endif
                  endif
                enddo ! biter2
              elseif ( myid .eq. 0 .and. kiter .eq. 1 ) then
                write(6,*) ibeta, kiter, biter2, val_energies( biter2, kiter, 1, ibw )
              endif
            enddo ! biter1
          enddo ! kiter
        enddo
      enddo ! ibw
    endif
    !
!    open(myid+2000)
!    write(myid+2000,*) allow%valr
!    close(myid+2000)

!    allow%valr = 1.0_dp
!    allow%vali = 1.0_dp

  end subroutine energies_allow

  subroutine find_fermi( sys, val_energies, con_energies, nelectron, efermi, &
                                     homo, lumo, cliph, metal, dft_flag, ierr )
    use OCEAN_mpi!, only : myid, root, comm
    use OCEAN_system
    use OCEAN_constants, only : Hartree2eV
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    integer, intent( in ) :: nelectron 
    real(dp), intent( in ) :: con_energies( sys%brange(3):sys%brange(4), sys%nkpts, sys%nspn, sys%nbw ), &
                              val_energies( sys%brange(1):sys%brange(2), sys%nkpts, sys%nspn, sys%nbw )
!                              val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn )
    real(dp), intent( out ) :: efermi, homo, lumo, cliph
    logical, intent( out ) :: metal
    logical, intent( in ) :: dft_flag
    integer, intent( inout ) :: ierr
    !
    real(dp), allocatable :: simple_energies( : )
    real(dp) :: temp, per_electron_dope
    integer :: i_band, overlap, t_electron, n_electron_dope
    integer :: iter, node, node2, top, kiter, ierr_, ispn, i, ii, ibw
    logical :: doping
    !
    !
    !
    n_electron_dope = 0
    if( myid .eq. root ) then
      ! if the file metal is not present then the value of metal is auto set to false
      inquire( file='metal', exist=metal )
      if( metal ) then
        open( unit=99, file='metal', form='formatted', status='old' )
        read( 99, * ) metal
        close( 99 )
        !
        inquire( file='doping', exist=doping )
        if( doping ) then
          open( unit=99, file='doping', form='formatted', status='old' )
          read( 99, * ) per_electron_dope
          close( 99 )
          n_electron_dope = floor( per_electron_dope * dble( sys%nkpts * sys%nspn * sys%nbw) / 2.d0 )
          write( 6, * ) 'Doping option:'
          write( 6, * ) 'Percent doping: ', per_electron_dope
          write( 6, * ) 'Modifying electron number by ', n_electron_dope
          write( 6, * ) 'Effective doping percent: ', dble( n_electron_dope ) * 2.d0 / ( sys%nkpts * sys%nspn)
        endif
      endif
      !
      !
      if( mod( nelectron, 2 ) .ne. 0 ) then
        if ( metal .eqv. .false. )  then
          write( 6, * ) 'WARNING: for partial occupation we must have a metal!'
          write( 6, * ) 'Setting metal = true and continuing'
          metal = .true.
        endif
        if( mod( sys%nkpts, 2 ) .ne. 0 ) then
          ierr = 80
          write( 6, * ) 'Number of kpts * number of electrons must be even for spinless calc.'
!          goto 111
        endif
      endif
    endif
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr_ .ne. MPI_SUCCESS ) return
    if( ierr .ne. 0 ) return
    call MPI_BCAST( metal, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    call MPI_BCAST( n_electron_dope, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif
    !
    overlap = sys%brange( 2 ) - sys%brange( 3 ) + 1
    if( ( metal .or. sys%valence_ham_spin .gt. 1 ) .and. overlap .gt. 0 ) then
      ii = 0
      allocate( simple_energies( overlap * sys%nkpts * sys%valence_ham_spin * sys%nbw) )
      do ibw = 1, sys%nbw
        do i = 1, sys%valence_ham_spin
          ispn = min( i, sys%nspn )
          do kiter = 1, sys%nkpts
            do iter = 1, overlap
!                 simple_energies( iter + ( kiter - 1 ) * overlap + ( ispn - 1 ) * sys%nkpts * overlap ) = &
                ii = ii + 1
                simple_energies( ii ) = val_energies( iter + sys%brange( 2 ) - overlap, kiter, ispn, ibw )
            enddo
          enddo
        enddo
      enddo
#if DEBUG
      open(99,file='unsorted_simple_energies.txt',form='formatted',status='unknown')
      node = 0
      do i = 1, sys%valence_ham_spin
        do kiter = 1, sys%nkpts
          do iter = 1, overlap
            node = node + 1
            write(99,'(F24.16)') simple_energies( node )
          enddo
        enddo
      enddo
      close(99)
#endif
     ! heap sort
      call MPI_BARRIER( comm, ierr )
      if( myid .eq. root ) write(6,*) 'sorting'
      top = overlap*sys%nkpts*sys%nspn*sys%nbw
      do iter = overlap*sys%nkpts*sys%nspn*sys%nbw / 2 , 1, -1
        temp = simple_energies( iter )
        node = iter + iter
       node2 = iter
        do
          if ( node .ge. top ) goto 10
          if ( (node .lt. top) .and. (simple_energies( node ) .lt. simple_energies( node + 1 ) ) ) node = node + 1
          if ( temp .lt. simple_energies( node ) ) then
            simple_energies( node2 ) = simple_energies( node )
            simple_energies( node ) = temp
            node2 = node
            node = node + node
          else
            goto 10
          endif
        enddo
 10 continue
      enddo
      do iter = overlap*sys%nkpts*sys%nspn*sys%nbw, 3, -1
        temp = simple_energies( iter )
        simple_energies( iter ) = simple_energies( 1 )
        node = 2
        node2 = 1
        do
          if ( node .gt. iter - 1) goto 20
          if ( ( node .lt. iter - 1 ) .and. ( simple_energies( node ) .lt. simple_energies( node + 1 ) ) ) node = node + 1
          if ( temp .lt. simple_energies( node ) ) then
            simple_energies( node2 ) = simple_energies( node )
            simple_energies( node ) = temp
            node2 = node
            node = node + node
          else
            goto 20
          endif
        enddo
  20 continue

      enddo
      if( simple_energies( 2 ) .lt. simple_energies( 1 ) ) then
        temp = simple_energies( 1 )
        simple_energies( 1 ) = simple_energies( 2 )
        simple_energies( 2 ) = temp
      endif
#if DEBUG
      open(99,file='sorted_simple_energies.txt',form='formatted',status='unknown')
      node = 0
      do i = 1, sys%valence_ham_spin
        do kiter = 1, sys%nkpts
          do iter = 1, overlap
            node = node + 1
            write(99,'(F24.16)') simple_energies( node )
          enddo
        enddo
      enddo
      close(99)
#endif
      t_electron = ( ( nelectron * sys%nkpts * sys%nspn * sys%nbw) / 2 ) &
                 - ( ( sys%brange( 3 ) - 1 ) * sys%nkpts * sys%nspn * sys%nbw) + n_electron_dope
      if( doping .and. ( myid .eq. root ) ) then
        write( 6, * ) 'old HOMO = ', simple_energies( t_electron - n_electron_dope )
        write( 6, * ) 'old LUMO = ', simple_energies( t_electron - n_electron_dope  + 1 )
      endif
      homo = simple_energies( t_electron )
      lumo = simple_energies( t_electron + 1 )
      efermi = ( lumo + homo ) / 2.0_dp
    else ! not metal
!      i_band = nelectron / 2 - sys%brange( 1 ) + 1
      i_band = nelectron / 2
      if( myid .eq. root ) write( 6, * ) "i_band = ", i_band
      homo =  val_energies( i_band, 1 , 1, 1)
      do ibw = 1, sys%nbw
        do ispn = 1, sys%nspn
          do kiter = 1, sys%nkpts
            homo = max( val_energies( i_band, kiter , ispn, ibw), homo )
          enddo
        enddo
!        if( sys%nspn .eq. 2 ) homo = max( val_energies( i_band, kiter , 2), homo )
      enddo
!      i_band = nelectron / 2 - sys%brange( 3 ) + 2
      i_band = nelectron / 2 + 1
      if( myid .eq. root ) write( 6, * ) "i_band = ", i_band
      lumo = con_energies( i_band, 1 , 1, 1)
      do ibw = 1, sys%nbw
        do ispn = 1, sys%nspn
          do kiter = 1, sys%nkpts
            lumo = min( con_energies( i_band, kiter , ispn, ibw), lumo )
          enddo
        enddo
!        if( sys%nspn .eq. 2 ) lumo = min( con_energies( i_band, kiter , 1), lumo )
      enddo
      efermi = ( lumo + homo ) / 2.d0
    endif
    !
!    cliph = con_energies( sys%brange( 4 ) - sys%brange( 3 ) + 1, 1 , 1)
    cliph = con_energies( sys%brange( 4 ), 1, 1, 1 )
    do ibw = 1, sys%nbw
      do ispn = 1, sys%nspn
        do kiter = 1, sys%nkpts
          cliph = min( con_energies( sys%brange( 4 ), kiter , ispn, ibw), cliph )
        enddo
      enddo
    enddo
    !
    if( myid .eq. root ) then
      if( dft_flag ) then
        write( 6, * ) ' **** DFT energy summary **** '
      else
        write( 6, * ) ' **** DFT+GW energy summary **** '
      endif
      write( 6, * ) 'HOMO = ', homo * Hartree2eV
      write( 6, * ) 'LUMO = ', lumo * Hartree2eV
      write( 6, * ) 'Fermi Energy = ', efermi * Hartree2eV
      if( dft_flag ) then
        write( 6, * ) 'DFT gap = ', ( lumo - homo ) * Hartree2eV
      else
        write( 6, * ) '    gap = ', ( lumo - homo ) * Hartree2eV
      endif
      write( 6, * ) 'clips = ', efermi* Hartree2eV, cliph* Hartree2eV, (cliph - efermi)* Hartree2eV
      open(unit=99,file='gap',form='formatted',status='unknown')
      write(99,*) ( lumo - homo ) * Hartree2eV
      close(99)
    endif

  111 continue
  !
  end subroutine find_fermi



end module OCEAN_val_energy
