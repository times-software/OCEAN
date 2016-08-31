module OCEAN_val_energy
  use AI_kinds
  implicit none

  contains

  subroutine OCEAN_read_energies( sys, p_energy, allow, ierr )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( inout ) :: p_energy, allow
    integer, intent( inout ) :: ierr
    !


    real( DP ), allocatable :: val_energies( :, : ), con_energies( :, : )
    real( DP ) :: efermi, homo, lumo, cliph
    integer :: ik, ibv, ibc, fh
    logical :: metal
    integer :: nbv, nbc(2), nk
#ifdef MPI
    integer(MPI_OFFSET_KIND) :: offset
#endif

    real(DP), parameter :: Ha_to_eV = 27.21138386_dp


    allocate( val_energies( sys%cur_run%val_bands, sys%nkpts ), con_energies( sys%cur_run%num_bands, sys%nkpts ), STAT=ierr )
    if( ierr .ne. 0 ) return


! Read in energies
    if( .false. ) then
      if( myid .eq. root ) then
        open(unit=99,file='enkfile',form='formatted',status='old')
        do ik = 1, sys%nkpts
          read(99,*) val_energies( :, ik )
          read(99,*) con_energies( :, ik )
        enddo
        val_energies( :, : ) = val_energies( :, : ) * Ha_to_eV
        con_energies( :, : ) = con_energies( :, : ) * Ha_to_eV
      endif
#ifdef MPI
      call MPI_BCAST( val_energies, sys%cur_run%val_bands*sys%nkpts, MPI_DOUBLE_PRECISION, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
      call MPI_BCAST( con_energies, sys%cur_run%num_bands*sys%nkpts, MPI_DOUBLE_PRECISION, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
#endif
    else
#ifdef MPI
      if( myid .eq. root ) then
        open(unit=99,file='tmels.info',form='formatted',status='old')
        read(99,*) nbv, nbc(1), nbc(2), nk
        close(99)

        if( nk .ne. sys%nkpts ) then
          ierr = -1
          write(6,*) 'tmels.info not compatible with other run information: NKPTS'
          return
        endif
      endif
#ifdef MPI
    call MPI_BCAST( nbc, 2, MPI_INTEGER, 0, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( nbv, 1, MPI_INTEGER, 0, comm, ierr )
    if( ierr .ne. 0 ) return
#endif


      call MPI_FILE_OPEN( comm, 'val_energies.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
      offset = 0
      call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL,ierr)
      if( ierr .ne. MPI_SUCCESS ) return
      do ik = 1, sys%nkpts
        offset = ( ik - 1 ) * nbv
        call MPI_FILE_READ_AT( fh, offset, val_energies( 1, ik ), sys%cur_run%val_bands, MPI_DOUBLE_PRECISION, &
                               MPI_STATUS_IGNORE, ierr )
        if( ierr .ne. MPI_SUCCESS) return
      enddo
      call MPI_FILE_CLOSE( fh, ierr )
      if( ierr .ne. MPI_SUCCESS) return

      call MPI_FILE_OPEN( comm, 'con_energies.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
      offset = 0
      call MPI_FILE_SET_VIEW( fh, offset, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL,ierr)
      if( ierr .ne. MPI_SUCCESS ) return
      do ik = 1, sys%nkpts
        offset = ( ik - 1 ) * ( nbc( 2 ) - nbc( 1 ) + 1 ) + ( sys%brange( 3 ) - nbc( 1 ) )
        call MPI_FILE_READ_AT( fh, offset, con_energies( 1, ik ), sys%cur_run%num_bands, MPI_DOUBLE_PRECISION, &
                               MPI_STATUS_IGNORE, ierr )
        if( ierr .ne. MPI_SUCCESS) return
      enddo
      call MPI_FILE_CLOSE( fh, ierr )
      if( ierr .ne. MPI_SUCCESS) return
#else
      ierr = -1
      if( myid .eq. root ) write(6,*) 'MPI required for OBF-style!'
      return
#endif
    endif
    

    call find_fermi( sys, val_energies, con_energies, sys%nelectron, efermi, &
                     homo, lumo, cliph, metal, ierr )
    if( ierr .ne. 0 ) return

    call energies_allow( sys, val_energies, con_energies, sys%nelectron, efermi, cliph, &
                                allow, metal, ierr )
    if( ierr .ne. 0 ) return


!   call GW corrections time


    do ik = 1, sys%nkpts
      do ibv = 1, sys%cur_run%val_bands
        do ibc = 1, sys%cur_run%num_bands
          p_energy%valr( ibc, ibv, ik, 1 ) = con_energies( ibc, ik ) - val_energies( ibv, ik )
          p_energy%vali( ibc, ibv, ik, 1 ) = 0.0_dp
        enddo
      enddo
    enddo

    deallocate( val_energies, con_energies )
    


  end subroutine OCEAN_read_energies



  subroutine energies_allow( sys, val_energies, con_energies, nelectron, efermi, cliph, &
                                allow, metal, ierr )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( inout ) :: allow
    integer, intent( in ) :: nelectron
    real(kind=kind(1.d0)), intent( in ) :: con_energies( sys%cur_run%num_bands, sys%nkpts ), &
                                           val_energies( sys%cur_run%val_bands, sys%nkpts ),  &
                                           efermi, cliph
    logical, intent( in ) :: metal
    integer, intent( inout ) :: ierr
    !
    integer :: kiter, biter1, biter2
    !
    !
    ! Already zero'd
    allow%valr = 0.0_dp
    allow%vali = 0.0_dp
    !
    ! In storing psi the conduction band index is the fast index
    if( metal ) then
      do kiter = 1, sys%nkpts
        do biter2 = 1, sys%cur_run%val_bands
          if ( val_energies( biter2, kiter ) .le. efermi ) then

            do biter1 = 1, sys%cur_run%num_bands
              if ( ( con_energies( biter1, kiter ) .ge. efermi ) .and. &
                   ( con_energies( biter1, kiter ) .le. cliph ) ) then
                allow%valr( biter1, biter2, kiter, 1 ) = 1.0d0
                allow%vali( biter1, biter2, kiter, 1 ) = 1.0d0
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
    else ! not metal 
      if( sys%backf ) ierr = 413
      do kiter = 1, sys%nkpts
        do biter2 = 1, ( nelectron / 2 ) - sys%brange( 1 ) + 1
          do biter1 = 2 - sys%brange( 3 ) + ( nelectron / 2 ), sys%cur_run%num_bands
            if(  con_energies( biter1, kiter ) .le. cliph ) then
                  allow%valr( biter1, biter2, kiter, 1 ) = 1.0d0
                  allow%vali( biter1, biter2, kiter, 1 ) = 1.0d0
            endif
          enddo ! biter2
        enddo ! biter1
      enddo ! kiter
    endif
    !
!    allow%valr = 1.0_dp
!    allow%vali = 1.0_dp

  end subroutine energies_allow

  subroutine find_fermi( sys, val_energies, con_energies, nelectron, efermi, &
                                     homo, lumo, cliph, metal, ierr )
    use OCEAN_mpi, only : myid, root, comm
    use OCEAN_system
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    integer, intent( in ) :: nelectron 
    real(dp), intent( in ) :: con_energies( sys%cur_run%num_bands, sys%nkpts ), val_energies( sys%cur_run%val_bands, sys%nkpts )
    real(dp), intent( out ) :: efermi, homo, lumo, cliph
    logical, intent( out ) :: metal
    integer, intent( inout ) :: ierr
    !
    real(dp), allocatable :: simple_energies( : )
    real(dp) :: temp, per_electron_dope
    integer :: i_band, overlap, t_electron, n_electron_dope
    integer :: iter, node, node2, top, kiter
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
          n_electron_dope = floor( per_electron_dope * dble( sys%nkpts ) / 2.d0 )
          write( 6, * ) 'Doping option:'
          write( 6, * ) 'Percent doping: ', per_electron_dope
          write( 6, * ) 'Modifying electron number by ', n_electron_dope
          write( 6, * ) 'Effective doping percent: ', dble( n_electron_dope ) * 2.d0 / sys%nkpts
        endif
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
        goto 111
      endif
    endif
#ifdef MPI
    call MPI_BCAST( metal, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    call MPI_BCAST( n_electron_dope, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif
    !
    if( metal ) then
      overlap = sys%brange( 2 ) - sys%brange( 3 ) + 1
      allocate( simple_energies( overlap * sys%nkpts ) )
      do kiter = 1, sys%nkpts
        do iter = 1, overlap
          simple_energies( iter + ( kiter - 1 ) * overlap ) = &
                        val_energies( iter + sys%brange( 2 ) - overlap, kiter )
        enddo
      enddo
     ! heap sort
      write(6,*) 'sorting'
      top = overlap*sys%nkpts
      do iter = overlap*sys%nkpts / 2 , 1, -1
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
      do iter = overlap*sys%nkpts, 1, -1
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
      t_electron = ( ( nelectron * sys%nkpts ) / 2 ) &
                 - ( ( sys%brange( 3 ) - 1 ) * sys%nkpts ) + n_electron_dope
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
      homo =  val_energies( i_band, 1 )
      do kiter = 2, sys%nkpts
        homo = max( val_energies( i_band, kiter ), homo )
      enddo
      i_band = nelectron / 2 - sys%brange( 3 ) + 2
      if( myid .eq. root ) write( 6, * ) "i_band = ", i_band
      lumo = con_energies( i_band, 1 )
      do kiter = 2, sys%nkpts
        lumo = min( con_energies( i_band, kiter ), lumo )
      enddo
      efermi = ( lumo + homo ) / 2.d0
    endif
    !
    cliph = con_energies( sys%brange( 4 ) - sys%brange( 3 ) + 1, 1 )
    do kiter = 2, sys%nkpts
      cliph = min( con_energies( sys%brange( 4 ) - sys%brange( 3 ) + 1, kiter ), cliph )
    enddo
    !
    if( myid .eq. root ) then
      write( 6, * ) 'HOMO = ', homo
      write( 6, * ) 'LUMO = ', lumo
      write( 6, * ) 'Fermi Energy = ', efermi
      write( 6, * ) 'LDA gap = ', lumo - homo
      write( 6, * ) 'clips = ', efermi, cliph, cliph - efermi
    endif

  111 continue
  !
  end subroutine find_fermi



end module OCEAN_val_energy
