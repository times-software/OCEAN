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
    use OCEAN_constants, only : Hartree2eV
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( inout ) :: p_energy, allow
    integer, intent( inout ) :: ierr
    !


    real( DP ), allocatable, dimension(:,:,:) :: val_energies, con_energies, im_val_energies, im_con_energies
    real( DP ) :: efermi, homo, lumo, cliph
    integer :: ik, ibv, ibc, fh, ispn, jspn, ibeta, i, j
    logical :: metal, have_imaginary
    integer :: nbv, nbc(2), nk
#ifdef MPI
    integer(MPI_OFFSET_KIND) :: offset
#endif

!    real(DP), parameter :: Ha_to_eV = 27.21138386_dp


    allocate( val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ),  &
              con_energies( sys%cur_run%num_bands, sys%nkpts, sys%nspn ), STAT=ierr )
    if( ierr .ne. 0 ) return


! Read in energies
       select case( sys%enk_selector )

       case( 0 )
          
          if( myid .eq. root ) then
             open(unit=99,file='enkfile',form='formatted',status='old')
             do ispn=1, sys%nspn
                do ik = 1, sys%nkpts
                   read(99,*) val_energies( :, ik ,ispn)
                   read(99,*) con_energies( :, ik ,ispn )
                enddo
                val_energies( :, : ,ispn ) = val_energies( :, : ,ispn ) / 2.0_dp ! * Hartree2eV !Ha_to_eV
                con_energies( :, : ,ispn ) = con_energies( :, : ,ispn ) / 2.0_dp ! * Hartree2eV !Ha_to_eV
             enddo  !ispn
          endif
#ifdef MPI
          call MPI_BCAST( val_energies, sys%cur_run%val_bands*sys%nkpts*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
          call MPI_BCAST( con_energies, sys%cur_run%num_bands*sys%nkpts*sys%nspn, MPI_DOUBLE_PRECISION, root, comm, ierr )
          if( ierr .ne. MPI_SUCCESS ) return
#endif
          
       case( 1 )
          
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
          
       case default
          ierr = -2
          if( myid .eq. root ) write(6,*) 'Un-supported enk_selector:', sys%enk_selector
          return
          
       end select
    
    if( myid .eq. root ) then
      open( unit=99,file='val_energy_test.txt', form='formatted',status='unknown' )
      write(99,*) val_energies(:,:,:)
      close(99)
    endif
!   call GW corrections time

    allocate( im_val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ), & 
              im_con_energies( sys%cur_run%num_bands, sys%nkpts, sys%nspn ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call val_gw( sys, val_energies, con_energies, have_imaginary, im_val_energies, &
                 im_con_energies, ierr )
    if( ierr .ne. 0 ) return

    call find_fermi( sys, val_energies, con_energies, sys%nelectron, efermi, &
                     homo, lumo, cliph, metal, ierr )
    if( ierr .ne. 0 ) return

    call energies_allow( sys, val_energies, con_energies, sys%nelectron, efermi, cliph, &
                                allow, metal, ierr )
    if( ierr .ne. 0 ) return



    ibeta=0
    do i=1,sys%valence_ham_spin
       ispn = min( i, sys%nspn ) 
       do j=1,sys%valence_ham_spin
          jspn = min( j, sys%nspn )
          ibeta=ibeta+1
          do ik = 1, sys%nkpts
             do ibv = 1, sys%cur_run%val_bands
                do ibc = 1, sys%cur_run%num_bands
                   p_energy%valr( ibc, ibv, ik, ibeta ) = con_energies( ibc, ik, jspn ) - val_energies( ibv, ik, ispn )
                   !          p_energy%vali( ibc, ibv, ik, 1 ) = 0.0_dp
                enddo
             enddo
          enddo
       enddo  !jspn
    enddo  !ispn
    if( have_imaginary ) then
       ibeta=0
       do i=1,sys%valence_ham_spin
          ispn = min( i, sys%nspn )
          do j=1,sys%valence_ham_spin
             jspn = min( j, sys%nspn )
             ibeta=ibeta+1
             do ik = 1, sys%nkpts
                do ibv = 1, sys%cur_run%val_bands
                   do ibc = 1, sys%cur_run%num_bands
                      p_energy%vali( ibc, ibv, ik, ibeta ) = im_con_energies( ibc, ik, jspn ) - im_val_energies( ibv, ik, ispn )
                   enddo
                enddo
             enddo
          enddo !jspn
       enddo  !ispn
    else
      p_energy%vali( :, :, :, : ) = 0.0_dp
    endif


    deallocate( val_energies, con_energies, im_val_energies, im_con_energies )
    


  end subroutine OCEAN_read_energies

  ! Master routine for all the gw flavors (for the valence band). 
  !   For now this is just a stub for what I currently need.
  !   In the future need to better merge the valence and X-ray versions of everything
  subroutine val_gw( sys, val_energies, con_energies, have_imaginary, im_val_energies, &
                     im_con_energies, ierr )
    use OCEAN_system
    use OCEAN_mpi
    implicit none
    type( O_system ), intent( in ) :: sys
    real( DP ), intent( inout ), dimension( sys%cur_run%val_bands, sys%nkpts ) :: &
        val_energies, im_val_energies
    real( DP ), intent( inout ), dimension( sys%cur_run%num_bands, sys%nkpts ) :: &
        con_energies, im_con_energies

    logical, intent( out ) :: have_imaginary
    integer, intent( inout ) :: ierr
    !
    logical :: have_gw
    character( len=4 ) :: gw_control

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
!        case( 'band' )
!          call val_gw_by_band( sys, ierr, .false. )
!        case( 'ibnd' )
!          call val_gw_by_band( sys, ierr, .true. )
!        case( 'cstr' )
!          call val_gw_stretch( sys, ierr )
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
      call MPI_BCAST( val_energies, sys%cur_run%val_bands * sys%nkpts, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr )
      call MPI_BCAST( con_energies, sys%cur_run%num_bands * sys%nkpts, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr )
    endif
#endif

  end subroutine val_gw

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

  subroutine energies_allow( sys, val_energies, con_energies, nelectron, efermi, cliph, &
                                  allow, metal, ierr )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( inout ) :: allow
    integer, intent( in ) :: nelectron
    real( DP ), intent( in ) :: con_energies( sys%cur_run%num_bands, sys%nkpts, sys%nspn ), &
                                val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn ),  &
                                efermi, cliph
    logical, intent( in ) :: metal
    integer, intent( inout ) :: ierr
    !
    integer :: kiter, biter1, biter2, ibeta, i, ispn, j, jspn
    !
    !
    ! Already zero'd
    allow%valr = 0.0_dp
    allow%vali = 0.0_dp
    !
    ! If we have spins from the DFT then we can't rely on band index, need to use Fermi level 
    !  just like for the metals case. 
    !
    ! In storing psi the conduction band index is the fast index
    if( metal .or. sys%nspn .ne. 1 ) then
      ibeta = 0
      do i = 1, sys%valence_ham_spin
        ispn = min( i, sys%nspn )
        do j = 1, sys%valence_ham_spin
          jspn = min( j, sys%nspn )
          ibeta = ibeta + 1
          do kiter = 1, sys%nkpts
            do biter2 = 1, sys%cur_run%val_bands
              if ( val_energies( biter2, kiter, ispn ) .le. efermi ) then

                do biter1 = 1, sys%cur_run%num_bands
                  if ( ( con_energies( biter1, kiter, jspn ) .ge. efermi ) .and. &
                       ( con_energies( biter1, kiter, jspn ) .le. cliph ) ) then
                    allow%valr( biter1, biter2, kiter, ibeta ) = 1.0_dp
                    allow%vali( biter1, biter2, kiter, ibeta ) = 1.0_dp
                  endif
                enddo 
              elseif ( sys%backf ) then
                ierr = -413
                return
    !           do biter2 = 1, sys%cur_run%val_bands
    !              if ( ( val_energies( biter2, kiter ) .ge. efermi ) .and. &
    !                   ( val_energies( biter2, kiter ) .le. cliph ) ) then
    !                allow%valr( biter2, biter1, kiter, 1 ) = 1.0d0
    !                allow%vali( biter2, biter1, kiter, 1 ) = -1.0d0
    !              endif
!            enddo ! biter2
              endif
            enddo ! biter1
          enddo ! kiter
        enddo ! j
      enddo ! i
    else ! not metal && sys%nspn==1
      if( sys%backf ) ierr = 413
      do ibeta = 1, sys%nbeta
        do kiter = 1, sys%nkpts
          do biter2 = 1, ( nelectron / 2 ) - sys%brange( 1 ) + 1
            do biter1 = 2 - sys%brange( 3 ) + ( nelectron / 2 ), sys%cur_run%num_bands
              if(  con_energies( biter1, kiter, 1 ) .lt. cliph ) then
                    allow%valr( biter1, biter2, kiter, ibeta ) = 1.0_dp
                    allow%vali( biter1, biter2, kiter, ibeta ) = 1.0_dp
              endif
            enddo ! biter2
          enddo ! biter1
        enddo ! kiter
      enddo
    endif
    !
!    open(myid+2000)
!    write(myid+2000,*) allow%valr
!    close(myid+2000)

!    allow%valr = 1.0_dp
!    allow%vali = 1.0_dp

  end subroutine energies_allow

  subroutine find_fermi( sys, val_energies, con_energies, nelectron, efermi, &
                                     homo, lumo, cliph, metal, ierr )
    use OCEAN_mpi!, only : myid, root, comm
    use OCEAN_system
    use OCEAN_constants, only : Hartree2eV
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    integer, intent( in ) :: nelectron 
    real(dp), intent( in ) :: con_energies( sys%cur_run%num_bands, sys%nkpts, sys%nspn ), val_energies( sys%cur_run%val_bands, sys%nkpts, sys%nspn )
    real(dp), intent( out ) :: efermi, homo, lumo, cliph
    logical, intent( out ) :: metal
    integer, intent( inout ) :: ierr
    !
    real(dp), allocatable :: simple_energies( : )
    real(dp) :: temp, per_electron_dope
    integer :: i_band, overlap, t_electron, n_electron_dope
    integer :: iter, node, node2, top, kiter, ierr_, ispn, i
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
          n_electron_dope = floor( per_electron_dope * dble( sys%nkpts * sys%nspn ) / 2.d0 )
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
    if( metal .or. sys%valence_ham_spin .gt. 1 ) then
      overlap = sys%brange( 2 ) - sys%brange( 3 ) + 1
      allocate( simple_energies( overlap * sys%nkpts * sys%nspn ) )
      do i = 1, sys%valence_ham_spin
         ispn = min( i, sys%nspn )
         do kiter = 1, sys%nkpts
            do iter = 1, overlap
               simple_energies( iter + ( kiter - 1 ) * overlap + ( ispn - 1 ) * sys%nkpts * overlap ) = &
                    val_energies( iter + sys%brange( 2 ) - overlap, kiter, ispn )
            enddo
         enddo
      enddo
     ! heap sort
      write(6,*) 'sorting'
      top = overlap*sys%nkpts*sys%nspn
      do iter = overlap*sys%nkpts*sys%nspn / 2 , 1, -1
        temp = simple_energies( iter )
        node = iter + iter
       node2 = iter
        do
          if ( node .gt. top ) goto 10
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
      do iter = overlap*sys%nkpts*sys%nspn, 1, -1
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
      t_electron = ( ( nelectron * sys%nkpts * sys%nspn ) / 2 ) &
                 - ( ( sys%brange( 3 ) - 1 ) * sys%nkpts * sys%nspn ) + n_electron_dope
      if( doping .and. ( myid .eq. root ) ) then
        write( 6, * ) 'old HOMO = ', simple_energies( t_electron - n_electron_dope )
        write( 6, * ) 'old LUMO = ', simple_energies( t_electron - n_electron_dope  + 1 )
      endif
      homo = simple_energies( t_electron )
      lumo = simple_energies( t_electron + 1 )
      efermi = ( lumo + homo ) / 2.d0
    else ! not metal
      i_band = nelectron / 2 - sys%brange( 1 ) + 1
      if( myid .eq. root ) write( 6, * ) "i_band = ", i_band
      homo =  val_energies( i_band, 1 , 1)
      do kiter = 2, sys%nkpts
        homo = max( val_energies( i_band, kiter , 1), homo )
      enddo
      i_band = nelectron / 2 - sys%brange( 3 ) + 2
!      i_band = nelectron / 2 + 1
      if( myid .eq. root ) write( 6, * ) "i_band = ", i_band
      lumo = con_energies( i_band, 1 , 1)
      do kiter = 2, sys%nkpts
        lumo = min( con_energies( i_band, kiter , 1), lumo )
      enddo
      efermi = ( lumo + homo ) / 2.d0
    endif
    !
    cliph = con_energies( sys%brange( 4 ) - sys%brange( 3 ) + 1, 1 , 1)
    do kiter = 2, sys%nkpts
      cliph = min( con_energies( sys%brange( 4 ) - sys%brange( 3 ) + 1, kiter , 1), cliph )
    enddo
    !
    if( myid .eq. root ) then
      write( 6, * ) 'HOMO = ', homo * Hartree2eV
      write( 6, * ) 'LUMO = ', lumo * Hartree2eV
      write( 6, * ) 'Fermi Energy = ', efermi * Hartree2eV
      write( 6, * ) 'LDA gap = ', ( lumo - homo ) * Hartree2eV
      write( 6, * ) 'clips = ', efermi* Hartree2eV, cliph* Hartree2eV, (cliph - efermi)* Hartree2eV
    endif

  111 continue
  !
  end subroutine find_fermi



end module OCEAN_val_energy
