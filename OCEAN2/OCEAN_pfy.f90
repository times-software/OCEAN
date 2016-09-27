!
!  To account for angular/geometric dependence in the detection of 
!  partial fluoresence yield (or total?) it is necessary to carry 
!  out a pseudo-rixs-like calculation. We don't need the full RIXS 
!  losses, but just the total cross section that is headed off in 
!  a given outgoing q and polarization.
!
module OCEAN_PFY
  use AI_kinds
  
  implicit none

  private

  real(DP), allocatable :: re_xes(:,:,:,:)
  real(DP), allocatable :: im_xes(:,:,:,:)

  logical :: is_init = .false.
  logical :: is_loaded = .false.


  public :: OCEAN_pfy_init, OCEAN_pfy_load, OCEAN_pfy_act
  

  contains

  subroutine OCEAN_pfy_act( sys, p, energy, fh, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    use OCEAN_psi
    use OCEAN_constants, only : Hartree2eV

    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: p
    real(dp), intent( in ) :: energy
    integer, intent( in ) :: fh
    !
    real(dp), allocatable :: pfy( : )!, re_temp(:,:,:,:), im_temp(:,:,:,:)
    real(dp) :: re_temp, im_temp, scaling
    integer :: i_xes, icms, icml, ii, ialpha, ibeta, ivms, ivvms, ibv, ibc, ik
    character(len=100) :: fmt_string
    !

    if( myid .ne. root ) return
    if( sys%nXES_photon .lt. 1 ) return

    scaling = p%kpref / sqrt( real( sys%nkpts, DP ) * real( sys%cur_run%num_bands, DP ) )

    allocate( pfy( sys%nXES_photon ), STAT=ierr )
    if( ierr .ne. 0 ) return
    pfy(:) = 0.0_dp

    do i_xes = 1, sys%nXES_photon

      ! Need to loop over both conduction and valence spins
      ialpha = 0
      do icms = 0, 1
        do icml = 0, 2 * sys%cur_run%ZNL(3) !-sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          ii = 0
          do ivms = -1, 1, 2
            ialpha = ialpha + 1

            ibeta = 2* ( icms * (2*sys%cur_run%ZNL(3)+1) + icml )

            do ivvms = 0, 1 
              ibeta = ibeta + 1
              ii = ii + 1

              do ik = 1, sys%nkpts
                do ibv = 1, sys%val_bands
                  do ibc = 1, sys%num_bands
                    re_temp = re_xes(ibv,ik,ibeta,i_xes) * p%r(ibc,ik,ialpha) &
                            + im_xes(ibv,ik,ibeta,i_xes) * p%i(ibc,ik,ialpha)
                    !sign doesn't matter because we are squaring this
                    im_temp = im_xes(ibv,ik,ibeta,i_xes) * p%r(ibc,ik,ialpha) &
                            + re_xes(ibv,ik,ibeta,i_xes) * p%i(ibc,ik,ialpha)
                    pfy( i_xes ) = pfy( i_xes ) + re_temp**2 + im_temp**2
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
      enddo

    enddo

    i_xes = min( 99, sys%nXES_photon )
    write(fmt_string,'( "(1e16.9, ", I2, "(e16.9))" )') i_xes
    
    write(fh,fmt_string) energy*Hartree2eV, pfy(:)*scaling

    deallocate( pfy )


  end subroutine OCEAN_pfy_act

  subroutine OCEAN_pfy_load( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root, comm
!    use mpi, only : MPI_BCAST

    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys

    real(DP) :: tau( 3 ), rr, ri, ir, ii
    real(DP), allocatable, dimension(:,:,:) :: pcr, pci
    real(DP), allocatable, dimension(:,:) :: mer, mei
    integer :: nptot, ntot, ialpha, icms, ivms, icml, ikpt, iband, iter, nspn, ierr_, i_xes

    character (LEN=127) :: cks_filename
    character (LEN=18) :: mel_filename


    if( is_init .eqv. .false. ) then
      call OCEAN_pfy_init( sys, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( myid .eq. root ) then
      do i_xes = 1, sys%nXES_photon
        write(cks_filename,'(A5,A2,I4.4)' ) 'cksv.', sys%cur_run%elname, sys%cur_run%indx

        open(unit=99,file=cks_filename,form='unformatted',status='old')
        rewind( 99 )
        read ( 99 ) nptot, ntot, nspn
        read ( 99 ) tau( : )
        allocate( pcr( nptot, ntot, nspn ), pci( nptot, ntot, nspn ) )
        read ( 99 ) pcr
        read ( 99 ) pci
        close( unit=99 )

        if( nspn .ne. sys%nspn ) then
          ierr = -1
          write(6,*) 'Spin mismatch is fatal'
          return
        endif

        allocate( mer( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3) ),  &
                  mei( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3) ) )

        write(mel_filename,'(A5,A1,I3.3,A1,I2.2,A1,I2.2,A1,I2.2)' ) 'mels.', 'z', sys%cur_run%ZNL(1), &
                'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3), 'p', i_xes
        open( unit=99, file=mel_filename, form='formatted', status='old' )
        rewind( 99 )
    !    do is = 1, 2
          do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
            do iter = 1, nptot
              read( 99, * ) mer( iter, icml ), mei( iter, icml )
            enddo
          enddo
    !    enddo
        close( 99 )

        ialpha = 0
        if( sys%nspn == 1 ) then
          do icms = -1, 1, 2
            do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
              do ivms = -1, 1, 2
                ialpha = ialpha + 1
                if( icms .eq. ivms ) then
                  iter = 0
                  do ikpt = 1, sys%nkpts
                    do iband = 1, sys%cur_run%val_bands
                      iter = iter + 1
                      rr = dot_product( mer( :, icml ), pcr( :, iter, 1 ) )
                      ri = dot_product( mer( :, icml ), pci( :, iter, 1 ) )
                      ir = dot_product( mei( :, icml ), pcr( :, iter, 1 ) )
                      ii = dot_product( mei( :, icml ), pci( :, iter, 1 ) )
                      re_xes(iband,ikpt,ialpha,i_xes) = rr - ii
                      im_xes(iband,ikpt,ialpha,i_xes) = -ri - ir
                    enddo
                  enddo
                endif
              enddo
            enddo
          enddo
        else
          do icms = 1, 2
            do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
              do ivms = 1, 2
                ialpha = ialpha + 1
                if( icms .eq. ivms ) then
                  iter = 0
                  do ikpt = 1, sys%nkpts
                    do iband = 1, sys%cur_run%val_bands
                      iter = iter + 1
                      rr = dot_product( mer( :, icml ), pcr( :, iter, ivms ) )
                      ri = dot_product( mer( :, icml ), pci( :, iter, ivms ) )
                      ir = dot_product( mei( :, icml ), pcr( :, iter, ivms ) )
                      ii = dot_product( mei( :, icml ), pci( :, iter, ivms ) )
                      re_xes(iband,ikpt,ialpha,i_xes) = rr - ii
                      im_xes(iband,ikpt,ialpha,i_xes) = -ri - ir
                    enddo
                  enddo
                endif
              enddo
            enddo
          enddo
        endif

        deallocate( pcr, pci, mer, mei )

      enddo
    endif

!    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
!    if( ierr_ .ne. MPI_SUCCESS ) then
!      ierr = ierr_
!    endif

    if( ierr .ne. 0 ) return
    is_loaded = .true.

  end subroutine OCEAN_pfy_load


  subroutine OCEAN_pfy_init( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, root 

    implicit none
    integer, intent( inout ) :: ierr
    type( o_system ), intent( in ) :: sys

    if( is_init ) return


    if( myid .eq. root .and. sys%nXES_photon .gt. 0 ) then
      allocate( re_xes( sys%val_bands, sys%nkpts, sys%nalpha, sys%nXES_photon ), &
                im_xes( sys%val_bands, sys%nkpts, sys%nalpha, sys%nXES_photon ), STAT=ierr )
    else
      allocate( re_xes( 1, 1, 1, 1 ), im_xes( 1, 1, 1, 1 ), STAT=ierr )
    endif
    if( ierr .ne. 0 ) return
  
    is_init = .true.

  end subroutine OCEAN_pfy_init

end module OCEAN_PFY
