module OCEAN_rixs_holder
  use AI_kinds
  implicit none

  private

  complex(DP), allocatable :: xes_vec(:,:,:,:)

  integer :: local_ZNL(3)
  logical :: is_init

  public :: OCEAN_rixs_holder_load, OCEAN_rixs_holder_clean, OCEAN_rixs_holder_ctc

  contains

  subroutine OCEAN_rixs_holder_clean
    implicit none

    if( is_init ) deallocate( xes_vec )
    is_init = .false.
  end subroutine OCEAN_rixs_holder_clean

  subroutine OCEAN_rixs_holder_init( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi!, only : myid, root
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
    use OCEAN_mpi!, only : myid, root, comm

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


!> @brief Loads up the echamp (conduction-band--core-hole exciton intensities) 
!! and does the inner-product with the core-valence overlaps to get the full 
!! RIXS starting point
!
!> @details Loops over all the edges and fills out the full RIXS starting vec. 
!! The various l<sub>m</sub>'s of the core are summed over. If we are running 
!! spin=1 for the valence Hamiltonian (say a K edge of a nonmagnetic system) 
!! then the various spin states are also compressed down (but the up-down and 
!! down-up channels are zero in that case anyway). 
!! \todo Make a naming module so that the naming convention of files like 
!! the echamps are always consistent and only in a single place. (also things 
!! like absspct, etc.).  Consider bringing in the ocean_vector type instead 
!! of a complex array of fixed size. 
  subroutine rixs_seed( sys, p_vec, file_selector, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_filenames, only : OCEAN_filenames_read_ehamp

    implicit none

    type(O_system), intent( in ) :: sys 
    complex(DP), intent( out ) :: p_vec(:,:,:,:)
    integer, intent( in ) :: file_selector
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: rex( :, :, : )
    integer :: edge_iter, ic, icms, icml, ivms, ispin, i, j, ik
    character(len=25) :: echamp_file

    allocate( rex( sys%num_bands, sys%nkpts, 4*(2*sys%ZNL(3)+1) ) ) 

    p_vec(:,:,:,:) = 0.0_DP

    do edge_iter = 1, sys%nedges

      call OCEAN_filenames_read_ehamp( sys, echamp_file, edge_iter, ierr )
      if( ierr .ne. 0 ) return


      write(6,*) echamp_file
      open(unit=99,file=echamp_file,form='unformatted',status='old')
      rewind(99)
      read(99) rex(:,:,:)
      close(99)


! JTV block this so psi is in cache
      ic = 0
      do icms = 1, 3, 2    ! Core-hole spin becomes valence-hole spin
        do icml = 1, sys%ZNL(3)*2 + 1
          do ivms = 0, 1  ! conduction-band spin stays
            ! the order of ispin will be 1, 2, 3, 4 (assuming that sys%nbeta == 4)
            ispin = min( icms + ivms, sys%nbeta )  
            ic = ic + 1
            do ik = 1, sys%nkpts
              do i = 1, sys%val_bands
                do j = 1, sys%num_bands
                  p_vec( j, i, ik, ispin ) = p_vec( j, i, ik, ispin ) + &
                      rex( j, ik, ic ) * xes_vec(i,ik,icml,edge_iter)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo    

    deallocate( rex )

  end subroutine rixs_seed


  ! This is abstracted to allow for better tracking the 'mels' later
  subroutine OCEAN_rixs_holder_ctc( sys, p_vec, ierr )
    use OCEAN_system, only : o_system
    implicit none

    type(O_system), intent( in ) :: sys
    complex(DP), intent( out ) :: p_vec(:,:,:)
    integer, intent( inout ) :: ierr

    call ctc_rixs_seed( sys, p_vec, ierr )

  end subroutine OCEAN_rixs_holder_ctc

  subroutine ctc_rixs_seed( sys, p_vec, ierr )
    use OCEAN_system, only : o_system
    use OCEAN_filenames, only : OCEAN_filenames_read_ehamp

    implicit none

    type(O_system), intent( in ) :: sys
    complex(DP), intent( out ) :: p_vec(:,:,:)
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: rex(:,:,:), mels(:,:)
    integer :: edge_iter, ialpha, icms, icml, ivms, ic, ik, j, ialpha_old, icms_old, icml_old
    character(len=50) :: echamp_file

    p_vec(:,:,:) = 0.0_DP
    write(6,*) 'C2C hack', sys%nedges, sys%cur_run%rixsInputZNL(3), sys%ZNL(3)
    
    ! size of echamp comes from rixsInputZNL
    allocate( rex( sys%num_bands, sys%nkpts, 4*(2*sys%cur_run%rixsInputZNL(3)+1) ) )

    allocate( mels( 2*sys%cur_run%rixsInputZNL(3)+1, sys%ZNL(3)*2 + 1 ) )

    call ctc_mels_hack( sys%cur_run%rixsInputZNL(3), sys%ZNL(3), mels, sys%cur_run%rixs_pol, ierr )
    if( ierr .ne. 0 ) return

    do edge_iter = 1, sys%nedges

!      edge_iter = sys%cur_run%indx
      call OCEAN_filenames_read_ehamp( sys, echamp_file, edge_iter, ierr )
      if( ierr .ne. 0 ) return
      write(6,*) echamp_file

!      write(6,*) echamp_file
      open(unit=99,file=echamp_file,form='unformatted',status='old')
      rewind(99)
      read(99) rex(:,:,:)
      close(99)

#if 0
      !JTV plenty remains to move beyond 1s2p

      ! carefully write out the ordering of rex, mels

      ialpha = 0
      do icms = 1, 3, 2
        do icml = 1, sys%cur_run%ZNL(3)*2 + 1
          do ivms = 0, 1

            ialpha = ialpha + 1
            ic = ivms+icms
            do ik = 1, sys%nkpts
              do j = 1, sys%num_bands
              
                p_vec( j, ik, ialpha ) = p_vec( j, ik, ialpha ) + rex( j, ik, ic ) * mels( icml )
      
              enddo
            enddo

          enddo
        enddo
      enddo
#endif


      ialpha_old = 0
      do icms_old = 1, 2
        do icml_old = 1, sys%cur_run%rixsInputZNL(3)*2 + 1
          do ivms = 0, 1
            ialpha_old = ialpha_old + 1

            do icms = 0, 1
              do icml = 1, sys%cur_run%ZNL(3)*2 + 1
                ialpha = 1 + ivms + ( icml - 1 ) * 2 + icms * ( sys%cur_run%ZNL(3)*2 + 1 ) * 2

                do ik = 1, sys%nkpts
                  do j = 1, sys%num_bands
                    p_vec( j, ik, ialpha ) = p_vec( j, ik, ialpha ) &
                                           + rex( j, ik, ialpha_old ) * mels( icml_old, icml )
                  enddo
                enddo
      
              enddo
            enddo

          enddo
        enddo
      enddo


    enddo

    deallocate( mels )
    deallocate( rex )

  end subroutine ctc_rixs_seed

  subroutine ctc_mels_hack( l_orig, lc, mels, photonNum, ierr )
    integer, intent( in ) :: l_orig
    integer, intent( in ) :: lc
    complex(DP), intent( out ) :: mels( :, : )
    integer, intent( in ) :: photonNum
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable, dimension(:) :: xsph, ysph, zsph, wsph
    real(DP) :: prefs( 0: 1000 ) 
    real(DP) :: su, ehat(3), edot
    complex(DP) :: csu, ylm, ylcmc
    integer :: nsphpt, i, m_orig, mc
    integer, parameter :: lmax = 5
    character(len=10) :: spcttype
    character(len=255) :: spctfile

    open( unit=99, file='sphpts', form='formatted', status='old' )
    rewind 99
    read ( 99, * ) nsphpt
    allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
    su = 0.0_dp
    do i = 1, nsphpt
      read( 99, * ) xsph(i), ysph(i), zsph(i), wsph( i )
      su = su + wsph( i )
    enddo
    close( 99 )
    wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / su )
    write ( 6, * ) nsphpt, ' points with weights summing to four pi '

    call getprefs( prefs, lmax, nsphpt, wsph, xsph, ysph, zsph )

    write( spctfile, '(A6,I0)' ) 'photon', photonNum
    write(6,*) 'ctc_hack: ', spctfile

!    open( unit=99, file='spectfile', form='formatted', status='unknown' )
    open( unit=99, file=spctfile, form='formatted', status='old' )
    rewind 99
    read ( 99, * ) spcttype
    call fancyvector( ehat, su, 99 )
    close( unit=99 )


    do mc = -lc, lc
      do m_orig = - l_orig, l_orig

        csu = 0.0_dp
        
        do i = 1, nsphpt
           call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
           call getylm( l_orig, m_orig, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
           edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
           csu = csu + conjg(ylcmc) * edot * ylm * wsph(i)
        end do
        mels( m_orig + l_orig + 1, mc + lc + 1 ) = csu
      enddo
    enddo

!    mels(:) = 1.0d0

    deallocate( xsph, ysph, zsph, wsph )

  end subroutine ctc_mels_hack

  subroutine cksv_read( sys, file_selector, ierr )
    use OCEAN_system
    implicit none
    type(O_system), intent( in ) :: sys
    integer, intent( in ) :: file_selector
    integer, intent( inout ) :: ierr

    real(DP), allocatable, dimension(:,:) :: mer, mei, pcr, pci
    complex(DP), allocatable :: pcTemp(:,:,:)
    real(DP) :: rr, ri, ir, ii, tau(3)
    integer :: nptot, ntot, nptot_check, nspn
    integer :: icml, iter, ik, i, edge_iter, ivms

    character(len=127) :: cks_filename
    character(len=18) :: mel_filename

    logical :: ex

    select case ( file_selector )

    case( 1 )
      
      do edge_iter = 1, sys%nedges

        write(cks_filename, '(A8,A2,I4.4)' ) 'parcksv.', sys%cur_run%elname, edge_iter
        inquire( file=cks_filename, exist=ex )
        if( ex ) then
          write(6,*) 'Using parallel cks: ', trim(cks_filename)
          open(unit=99, file=cks_filename, form='unformatted', status='old', access='stream' )
          read(99) nptot, ntot, nspn
          if( nspn .ne. 1 ) then
            ierr = -125
            write(6,*) 'Spin RIXS not yet supported'
            close(99)
            return
          endif
          if( edge_iter .eq. 1 ) allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
          allocate( pcTemp( nptot, ntot, nspn ) )
          read( 99 ) pcTemp(:,:,:)
          close(99)
          do ivms = 1, nspn
            do iter = 1, ntot
!              pcr(:,iter,ivms) = real( pcTemp(:,iter,ivms), DP )
!              pci(:,iter,ivms) = aimag( pcTemp(:,iter,ivms) )
              pcr(:,iter) = real( pcTemp(:,iter,ivms), DP )
              pci(:,iter) = aimag( pcTemp(:,iter,ivms) )
            enddo
          enddo
          deallocate( pcTemp )

        else
    
          write(6,'(A5,A2,I4.4)' ) 'cksv.', sys%cur_run%elname, edge_iter
          write(cks_filename,'(A5,A2,I4.4)' ) 'cksv.', sys%cur_run%elname, edge_iter
          open(unit=99,file=cks_filename,form='unformatted',status='old')
          rewind( 99 )
          read ( 99 ) nptot, ntot
          if( edge_iter .eq. 1 ) allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
          read ( 99 ) tau( : )
          write(6,*) tau(:)
          read ( 99 ) pcr
          read ( 99 ) pci
          close( unit=99 )
        endif

        ! check ntot
        if( ntot .ne. sys%nkpts * sys%val_bands ) then
          write(6,*) 'Mismatch bands*kpts vs ntot', ntot, sys%nkpts,  sys%val_bands
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
