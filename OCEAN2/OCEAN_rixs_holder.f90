module OCEAN_rixs_holder
  use AI_kinds
  implicit none

  private

  complex(DP), allocatable :: xes_vec(:,:,:,:)

  integer :: local_ZNL(3)
  logical :: is_init

  public :: OCEAN_rixs_holder_load, OCEAN_rixs_holder_clean

  contains

  subroutine OCEAN_rixs_holder_clean
    implicit none

    if( is_init ) deallocate( xes_vec )
    is_init = .false.
  end subroutine OCEAN_rixs_holder_clean

  subroutine OCEAN_rixs_holder_init( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    implicit none
    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    ! run in serial !
    if( myid .ne. root ) return

    ! Check to see if we are still ok
    !   if previously initiated and both Z and L are the same
    if( is_init ) then
      if( local_ZNL(1) .eq. sys%ZNL(1) .and. local_ZNL(2) .eq. sys%ZNL(2) .and. local_ZNL(3) .eq. sys%ZNL(3) ) return
      deallocate( xes_vec )
      is_init = .false.
    endif

    allocate( xes_vec( sys%val_bands, sys%nkpts, 2 * sys%ZNL(2) + 1, sys%nedges ) )
    is_init = .true.

  end subroutine OCEAN_rixs_holder_init

  subroutine OCEAN_rixs_holder_load( sys, p_vec, file_selector, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root, comm

    implicit none
    
    type(O_system), intent( in ) :: sys 
!    complex(DP), intent( inout ) :: p_vec(sys%num_bands, sys%val_bands, sys%nkpts, 1 )
    complex(DP), intent( inout ) :: p_vec(:,:,:,:)
    integer, intent( in ) :: file_selector
    integer, intent( inout ) :: ierr

    integer :: ierr_

    call OCEAN_rixs_holder_init( sys, ierr )
    if( ierr .ne. 0 ) return


    if( myid .eq. root ) then

      call cksv_read( sys, file_selector, ierr )
      if( ierr .ne. 0 ) return

      call rixs_seed( sys, p_vec, file_selector, ierr )
      
    endif

#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. MPI_SUCCESS ) return
    
    call MPI_BCAST( p_vec, size(p_vec), MPI_DOUBLE_COMPLEX, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif


  end subroutine OCEAN_rixs_holder_load


  subroutine rixs_seed( sys, p_vec, file_selector, ierr )
    use OCEAN_system

    implicit none

    type(O_system), intent( in ) :: sys 
    complex(DP), intent( out ) :: p_vec(:,:,:,:)
    integer, intent( in ) :: file_selector
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: rex( :, :, :, : )
    integer :: edge_iter, ic, icms, icml, ivms, ispin, i, j, ik
    character(len=25) :: echamp_file

    allocate( rex( sys%num_bands, sys%nkpts, 4*(2*sys%ZNL(3)+1), sys%nedges ) )

    do edge_iter = 1, sys%nedges
!      write(6,'(A7,I4.4)') 'echamp.', edge_iter
      write(echamp_file,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', edge_iter, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', & 
              sys%cur_run%rixs_energy
!      write(echamp_file, '(A7,I4.4)') 'echamp.', edge_iter
      write(6,*) echamp_file
      open(unit=99,file=echamp_file,form='unformatted',status='old')
      rewind(99)
      read(99) rex(:,:,:,edge_iter)
      close(99)


! JTV block this so psi is in cache
      ic = 0
      do icms = 0, 1
        do icml = 1, sys%ZNL(3)*2 + 1
          do ivms = 1, 3, 2
            ispin = min( icms + ivms, sys%nspn**2 )
            ic = ic + 1
            do ik = 1, sys%nkpts
              do i = 1, sys%val_bands
                do j = 1, sys%num_bands
                  p_vec( j, i, ik, ispin ) = p_vec( j, i, ik, ispin ) + &
                      rex( j, ik, ic, edge_iter ) * xes_vec(i,ik,icml,edge_iter)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo    

    deallocate( rex )

  end subroutine


  subroutine cksv_read( sys, file_selector, ierr )
    use OCEAN_system
    implicit none
    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: file_selector
    integer, intent( inout ) :: ierr

    real(DP), allocatable, dimension(:,:) :: mer, mei, pcr, pci
    real(DP) :: rr, ri, ir, ii, tau(3)
    integer :: nptot, ntot, nptot_check
    integer :: icml, iter, ik, i, edge_iter

    character(len=11) :: cks_filename
    character(len=18) :: mel_filename

    select case ( file_selector )

    case( 1 )
      
      do edge_iter = 1, sys%nedges
        write(6,'(A5,A2,I4.4)' ) 'cksv.', sys%cur_run%elname, edge_iter
        write(cks_filename,'(A5,A2,I4.4)' ) 'cksv.', sys%cur_run%elname, edge_iter
        open(unit=99,file=cks_filename,form='unformatted',status='old')
        rewind( 99 )
        read ( 99 ) nptot, ntot
        read ( 99 ) tau( : )
        write(6,*) tau(:)
        if( edge_iter .eq. 1 ) allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
        read ( 99 ) pcr
        read ( 99 ) pci
        close( unit=99 )

        ! check ntot
        if( ntot .ne. sys%nkpts * sys%val_bands ) then
          write(6,*) 'Mismatch bands*kpts vs ntot'
          ierr = -1
          return
        endif

        if( edge_iter .eq. 1 ) then
          allocate( mer( nptot, -sys%ZNL(3) : sys%ZNL(3) ),  mei( nptot, -sys%ZNL(3) : sys%ZNL(3) ) )

          write(mel_filename,'(A5,A1,I3.3,A1,I2.2,A1,I2.2,A1,I2.2)' ) 'mels.', 'z', sys%ZNL(1), & 
              'n', sys%ZNL(2), 'l', sys%ZNL(3), 'p', sys%cur_run%rixs_pol
           open( unit=99, file=mel_filename, form='formatted', status='old' ) 
          rewind( 99 ) 
          do icml = -sys%ZNL(3), sys%ZNL(3)
            do iter = 1, nptot
              read( 99, * ) mer( iter, icml ), mei( iter, icml ) 
            enddo
          enddo
          close( 99 ) 
          nptot_check = nptot

        else
          if( nptot .ne. nptot_check ) then
            write(6,*) 'nptot inconsistent between cores'
            ierr = -1
            return
          endif
        endif
    

        do icml = -sys%ZNL(3), sys%ZNL(3)
          iter = 0
          do ik = 1, sys%nkpts
            do i = 1, sys%val_bands
              iter = iter + 1
              rr = dot_product( mer( :, icml ), pcr( :, iter ) )
              ri = dot_product( mer( :, icml ), pci( :, iter ) )
              ir = dot_product( mei( :, icml ), pcr( :, iter ) )
              ii = dot_product( mei( :, icml ), pci( :, iter ) )
              xes_vec(i,ik,1 + icml + sys%ZNL(3), edge_iter) = cmplx( rr - ii, ri + ir )
            enddo
          enddo
        enddo
      enddo

      deallocate( mer, mei, pcr, pci )


    case( 0 )
      write(6,*) 'John is lazy'
      ierr = -1
      return
    case default
      ierr = -1
      return
    end select

  end subroutine cksv_read


end module OCEAN_rixs_holder
