module ocean_long_range

  use AI_kinds

  implicit none
 
  private

  type, public :: long_range

    integer( S_INT ) :: my_nxpts
    integer( S_INT ) :: my_nkpts

    integer( S_INT ) :: my_num_bands
    integer( S_INT ) :: my_start_nx
    integer( S_INT ) :: my_nkpts

    real( DP ) :: tau( 3 )


    real( DP ), pointer :: W( :, : )

    complex( DP ), pointer :: bloch_states( :, :, : )

  end type


  public :: create_lr, destroy_lr, lr_populate_W, lr_populate_bloch, lr_act

  contains

  subroutine lr_act_obf( sys, psi, hpsi, lr, obf, ierr )

    type( ocean_system ), intent( in ) :: sys
    type( long_range ), intent( in ) :: lr
    complex( DP ), intent( in ) :: psi( sys%num_bands, sys%nkpts, sys%nalpha )
    complex( DP ), intent( out ) :: hpsi( sys%num_bands, sys%nkpts, sys%nalpha )
    integer, intent( inout ) :: ierr

    !
    complex( DP ), parameter :: one = 1_DP
    complex( DP ), parameter :: zero = 0_DP
    !
    !
    complex( DP ), allocatable :: phi( :, :, : )
    real( DP ), allocatable :: xwrkr( : ), xwrki( : ), wrk( : )
    integer :: jfft, ialpha, ixpt, ikpt

    num_threads = 1

    hpsi( :, :, : ) = zero
    ! prep info for fft
    jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
    !
!    allocate( phi( lr%my_nkpts, sys%nalpha, lr%my_nxpts ) )
    

    allocate( local_beta( local_obf, local_k, sys%nalpha ) )

! Interleave alpha later to get better flops/mem?
! Possibly interleave message passing too 

!    if( ( sys%nkpt .lt. omp_numthreads ) .or. ( mod( sys%nkpt, omp_numthreads ) .ne. 0 ) ) then
! $OMP PARALLEL DO COLLAPSE( 2 )
!      do ikpt = 1, sys%nkpts
!        do iobf = 1, local_obf, 8
!          obf_width = 8
!          if( iobf + 8 .gt. local_obf ) then
!            obf_width = local_obf - iobf + 1 
!          endif
!        enddo
!      enddo
! $OMP END PARALLEL DO
!
!    else

    allocate( beta( sys%num_obf, sys%nkpt, sys%nalpha ) )
! Do obf over mpi
!$OMP PARALLEL DEFAULT(NONE) PRIVATE( ikpt, iobf, ialpha, obf_width, local_beta ) &
!$OMP& SHARED( system%nkpts, system%nobf, obf_cache, sys%nalpha, sys%num_bands, Bink, psi, beta )
    allocate( local_beta( sys%num_obf, sys%nalpha ) )
!$OMP DO
    do ikpt = 1, system%nkpts
      do iobf = 1, system%nobf, obf_cache
        obf_width = min( iobf + obf_cache - 1, system%nobf )
        obf_width = obf_width - iobf + 1
        do ialpha = 1, sys%nalpha, 4
          call ZGEMV( 'T', sys%num_bands, obf_cache, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha ), 1, zero, local_beta(iobf,ialpha), 1 )
          call ZGEMV( 'T', sys%num_bands, obf_cache, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha+1 ), 1, zero, local_beta(iobf,ialpha+1), 1 )
          call ZGEMV( 'T', sys%num_bands, obf_cache, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha+2 ), 1, zero, local_beta(iobf,ialpha+2), 1 )
          call ZGEMV( 'T', sys%num_bands, obf_cache, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha+3 ), 1, zero, local_beta(iobf,ialpha+3), 1 )
        enddo
      enddo
      do ialpha = 1, sys%nalpha, 4
! Communication of beta happens here
        beta( :, ikpt, ialpha ) = local_beta( :, ialpha )
        beta( :, ikpt, ialpha+1 ) = local_beta( :, ialpha+1 )
        beta( :, ikpt, ialpha+2 ) = local_beta( :, ialpha+2)
        beta( :, ikpt, ialpha+3 ) = local_beta( :, ialpha+3 )
      enddo
    enddo
!$OMP END DO
    deallocate( local_beta )
!$OMP END PARALLEL

! Need to sync up full beta


    if( lr%my_nxpts .gt. 8 * num_threads ) then
      xchunk = 8
    elseif(  lr%my_nxpts .gt. 4 * num_threads ) then
      xchunk = 4
    else
      xchunk = 1
    endif          

    allocate( tphi( lr%my_nxpts, lr%my_nkpts, sys%nalpha ) )

!$OMP PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( lr%my_nxpts, sys%num_alpha, sys%num_obf, sys%nkpt, tphi, beta, obf, lr%W, xchunk ) &
!$OMP& PRIVATE( ixpt, ialpha, iobf, xsize, iixpt, phi, xwrkr, xwrki, wrk )
    allocate( xwrkr( sys%nkpts ), xwrki( sys%nkpts ), wrk( jfft ), phi( sys%nkpts, xchunk, sys%nalpha ) )


!    obf_cache_align = 32
!$OMP DO SCHEDULE( STATIC )
    do ixpt = 1, lr%my_nxpts, xchunk
      xsize = xchunk
      if( ixpt - 1 + xchunk .gt. lr%my_nxpts ) then
        xsize = lr%my_nxpts - ixpt + 1
      endif
      do ialpha = 1, sys%num_alpha
!        do iobf = 1, sys%num_obf, obf_cache_align
!          stop_obf = iobf - 1 + obf_cache_align
!          stop_obf = min( stop_obf, sys%num_obf )
!          do ikpt = 1, sys%nkpt
!            phi( ikpt, ialpha, ixpt ) = phi( ikpt, ialpha, ixpt ) &
!               + sum( beta( iobf : stop_obf, ikpt, ialpha ) * obf( iobf : stop_obf , ixpt ) )
!          enddo
!        enddo
        call ZGEMM( 'T', 'N', sys%nkpts, xsize, sys%num_obf, one, beta( 1, 1, ialpha ), sys%num_obf, &
                     obf( 1, ixpt ), sys%num_obf, zero, phi( 1, 1, alpha ), sys%nkpts )


        do iixpt = 1, xsize

          ! For simplicity we are just going with the fft built in
          xwrkr( : ) = real( phi( :, iixpt, ialpha ) )
          xwrki( : ) = aimag( phi( :, iixpt, ialpha ) )
  
          call cfft( xwrkr, xwrki, zn(1), zn(1), zn(2), zn(3), -1, wrk, jfft )
  
          xwrkr( : ) = xwrkr( : ) * lr%W( :, ixpt + iixpt - 1 )
          xwrki( : ) = xwrki( : ) * lr%W( :, ixpt + iixpt - 1 )

          call cfft( xwrkr, xwrki, zn(1), zn(1), zn(2), zn(3), +1, wrk, jfft )

          phi( :, iixpt, ialpha ) = cmplx( xwrkr( : ), xwrki( : ) )

        enddo

        tphi( ixpt : ixpt + xsize, :, ialpha ) = TRANSPOSE( phi( :, 1:xsize, ialpha ) )


      enddo
    enddo
!$OMP END DO

    deallocate( xwrkr, xwrki, wrk, phi )

!$OMP END PARALLEL




! Need to do blocking for nthread > nalpha
!  Also to have better cache hitting?
!   And to distribute communications



  ! If nbands > 1024 then there is no way we fit tBink( obf, nband, ikpt ) in L2 cache
  !    instead just two-step it and hope for good L3?


    if( sys%numbands .gt. 1024 ) then


    elseif( sys%numbands .gt. 512 ) then
      cache_nkpts = 1
      cache_nalpha = min( 4, sys%num_alpha )
      cache_nxpts = min( 16, lr%my_xpts )
      cache_obf = 10
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%numbands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
      cache_obf = 2 * ( remaining_cache /  &
             ( 2 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    elseif( sys%numbands .gt. 128 ) then
      cache_nkpts = 1
      cache_nalpha = min( 12, sys%num_alpha )
      cache_nxpts = min( 16, lr%my_xpts )
      cache_obf = 24
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%numbands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
      cache_obf = 4 * ( remaining_cache /  &
             ( 4 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    elseif( sys%numbands .gt. 64 ) then
      cache_nkpts = min( 2, sys%nkpts )
      cache_nalpha = min( 12, sys%num_alpha )
      cache_nxpts = min( 32, lr%my_xpts )
      cache_obf = 32
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%numbands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
!      cache_obf = 4 * ( remaining_cache /  &
!             ( 4 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    else
      cache_nkpts = min( 4, sys%nkpts )
      cache_nalpha = min( 12, sys%num_alpha )
      cache_nxpts = min( 32, lr%my_nxpts )
      cache_obf = 32
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%numbands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
!      cache_obf = 4 * ( remaining_cache /  &
!             ( 4 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    endif

!$OMP PARALLEL DO COLLAPSE( 2 )
    do ikpt = 1, sys%nkpts, cache_nkpts
!      do block_obf = 1, sys%nobf/cache_obf
      do obf_start = 1, sys%num_obf/cache_obf
        stop_kpt = min( ikpt - 1 + cache_nkpts, sys%nkpts )
        obf_stop = min( obf_start - 1 + cache_obf, sys%num_obf ) 
!        obf_start = (block_obf - 1 ) * cache_obf + 1
!        obf_stop = block_obf  * cache_obf
        do block_alpha = 1, sys%num_alpha/cache_nalpha
          iialpha_start = (block_alpha-1)*cache_nalpha + 1
          iialpha_stop = iialpha_start + cache_nalpha - 1
          do block_xpt = 1, lr%nxpts/cache_nxpts
            x_start = (block_xpt - 1) * cache_nxpts + 1
            x_stop = block_xpt * cache_nxpts
            tphi_subset( :, :, : ) = tphi( x_start : x_stop, ikpt : stop_kpt, iialpha_start : iialpha_stop )
            obf_subset( :, : ) = obf( x_start : x_stop, obf_start : obf_stop )
            do iikpt = 1, cache_nkpts
              tBink( :, :, iikpt ) = TRANSPOSE( conjg( Bink( :, obf_start : obf_stop, ikpt + iikpt - 1 ) ) )
            enddo
            do ialpha = 1, cache_nalpha
              call ZGEMM( 'T', 'N', cache_obf, cache_nkpt, cache_nxpt, one, obf_subset, cache_nxpts, &
                          tphi_subset( 1, 1, ialpha ), zero, beta_subset( obf_start, 1, alpha ), cache_obf )
            enddo
            do ialpha = iialpha_start, iialpha_stop
              do iikpt = 1, cache_nkpts
                do iband = 1, sys%num_bands
                  hpsi( iband, ikpt, ialpha ) = hpsi( iband, ikpt, ialpha ) &
                                              + sum( tBink( :, iband, iikpt ) * beta_subset( :, iikpt, ialpha ) )
                enddo
              enddo
            enddo        

          enddo
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO
          
      

      

      

!    do ialpha = 1, sys%num_alpha
!      do ikpt = 1, 
!        do iobf = 1,
!          do ixpt = 1,
!            local_beta( iobf, ikpt, ialpha ) = local_beta( iobf, ikpt, ialpha ) &
!                                             + tphi( ixpt, ikpt, ialpha ) * obf( ixpt, iobf )
!          enddo
!        enddo
!        do iband = 1,
!          do iobf = 1, 
!            hpsi( iband, ikpt, ialpha ) = hpsi( iband, ikpt, ialpha ) &
!                                        + tBink( iobf, iband, ikpt ) * local_beta( iobf, ikpt, ialpha )
!          enddo
!        enddo
!      enddo
!    enddo
      
     


! Possibly use reduction or something to take care of hpsi

    deallocate( tphi, beta )

  end subroutine lr_act_obf


  subroutine lr_act( sys, psi, hpsi, lr, ierr )
    type( ocean_system ), intent( in ) :: sys
    type( long_range ), intent( in ) :: lr
    complex( DP ), intent( in ) :: psi( sys%num_bands, sys%nkpts, sys%nalpha )
    complex( DP ), intent( out ) :: hpsi( sys%num_bands, sys%nkpts, sys%nalpha )
    integer, intent( inout ) :: ierr

    !
    complex( DP ), parameter :: one = 1_DP
    complex( DP ), parameter :: zero = 0_DP
    !
    !
    complex( DP ), allocatable :: phi( :, :, : )
    real( DP ), allocatable :: xwrkr( : ), xwrki( : ), wrk( : )
    integer :: jfft, ialpha, ixpt, ikpt
    ! For each x-point in the unit cell
    !   Populate \phi(x,k) = \sum_n u(x,k) \psi_n(x,k)
    !   Do FFT for k-points
    !   Calculate W(x,k) x \phi(x,k)
    !   Do FFT back to k-points


    hpsi( :, :, : ) = zero
    ! prep info for fft
    jfft = 2 * max( zn( 1 ) * ( zn( 1 ) + 1 ), zn( 2 ) * ( zn( 2 ) + 1 ), zn( 3 ) * ( zn( 3 ) + 1 ) )
    !
    allocate( phi( lr%my_nkpts, lr%my_nxpts, sys%nalpha ) )

!$OMP PARALLEL DEFAULT( PRIVATE )
    allocate( xwrkr( sys%nkpts ), xwrki( sys%nkpts ), wrk( jfft ) )

!$OMP DO COLLAPSE( 2 )
    do ialpha = 1, sys%nalpha 
      do ixpt = 1, lr%my_nxpts

    ! Populate phi
#ifdef BLAS
        call ZGEMV( 'N', sys%num_bands, lr%my_nkpts, one, psi( 1, 1, ialpha ), sys%num_bands, &
                    lr%bloch_states( 1, 1, ixpt ), 1, zero, phi( 1, ixpt, ialpha ), 1 )
#else
!    phi = zero
        phi( :, ixpt, ialpha ) = zero
        do ikpt = 1, lr%my_nkpts
          phi( ikpt, ixpt, ialpha ) = phi( ikpt, ixpt, ialpha ) &
                            + sum( psi( :, ikpt, ialpha ) * lr%bloch_states( :, ikpt, ixpt ) )
        enddo
#endif          

    ! If there is some k-point division among procs this would be a problem here
        
        ! For simplicity we are just going with the fft built in
        xwrkr( : ) = real( phi( :, ixpt, ialpha ) )
        xwrki( : ) = aimag( phi( :, ixpt, ialpha ) )

        call cfft( xwrkr, xwrki, zn(1), zn(1), zn(2), zn(3), -1, wrk, jfft )

        xwrkr( : ) = xwrkr( : ) * lr%W( :, ixpt )
        xwrki( : ) = xwrki( : ) * lr%W( :, ixpt )

        call cfft( xwrkr, xwrki, zn(1), zn(1), zn(2), zn(3), +1, wrk, jfft )

        phi( :, ixpt, ialpha ) = cmplx( xwrkr( : ), xwrki( : ) )

!  Option 2, interleave phi here
!       tphi( ixpt, :, ialpha ) = cmplx( xwrkr( : ), xwrki( : ) )
      
      enddo
    enddo
!$OMP END DO


    if( .false. ) then
!$OMP DO
      do ikpt = 1, lr%nkpt
        do ialpha = 1, sys%nalpha
          ! do matrix vector now with tphi and u*
          do ixpt = 1, lr%my_nxpt
            hpsi( :, ikpt, ialpha ) = hpsi( :, ikpt, ialpha ) & 
                                  + conjg( lr%bloch_states( :, ikpt, ixpt ) ) * phi( ikpt, ixpt, ialpha ) 
          enddo
        enddo
      enddo
!$OMP END DO
    else
!$OMP DO COLLAPSE( 2 )
      do ikpt = 1, lr%nkpt  ! swap k and and alpha to reduce dependencies
        do ialpha = 1, sys%nalpha
          do ixpt = 1, lr%my_nxpt
            hpsi( :, ikpt, ialpha ) = hpsi( :, ikpt, ialpha ) &
                                  + conjg( lr%bloch_states( :, ikpt, ixpt ) ) * phi( ikpt, ixpt, ialpha )
          enddo
        enddo
      enddo
!$OMP END DO
    endif
      
    deallocate( xwrkr, xwrki, wrk )
    
!$OMP END PARALLEL
! Possibly use reduction or something to take care of hpsi

    deallocate( phi ) 

  end subroutine

  subroutine create_lr( sys, lr, ierr )
    use OCEAN_mpi, only : myid, nproc
    type( ocean_system ), intent( in ) :: sys
    type( long_range ), intent( inout ) :: lr
    integer, intent( inout ) :: ierr 
    !
    integer( S_INT ) :: nx_left, nx_start, nx_tmp

    ! SOP
    lr%my_num_bands = sys%num_bands
    lr%my_nkpts = sys%nkpts
    nx_left = nx
    nx_start = 1
    do i = 0, nproc - 1
      nx_tmp = nx_left / ( nproc - i )
      nx_left = nx_left - nx_tmp
      if( i .eq. myid ) then
        lr%my_nxpts = nx_tmp
        lr%my_start_nx = nx_start
      endif
      nx_start = nx_start + nx_tmp
    enddo


    allocate( lr%W( lr%my_nkpts, lr%my_nxpts ), &
              lr%bloch_states( lr%my_num_bands, lr%my_nkpts, lr%my_nxpts ), &
              STAT=ierr ) 
    if( ierr /= 0 ) then
      write(6,*) 'Failed to allocate lr%W or lr%bloch_states'
      goto 111
    endif


    call lr_populate_W( lr, ierr )
    if( ierr /= 0 ) then
      write(6,*) 'Failed to populate lr%W'
      goto 111
    endif


    call lr_populate_bloch( lr, ierr )
    if( ierr /= 0 ) then
      write(6,*) 'Failed to populate lr%bloch_states'
      goto 111
    endif



111 continue

  end subroutine create_W


  subroutine destroy_W( lr )
    deallocate( lr%W, lr%bloch_states )
  end subroutine destroy_W


  subroutine lr_populate_W( lr, ierr )
    use OCEAN_mpi, only : myid, nproc, comm

    type( long_range ), intent( inout ) :: lr
    integer, intent( inout ) :: ierr


    real( DP ) :: epsi, ptab( 100 ), avec( 3, 3 ), amet( 3, 3 )
    real( DP ) :: fr( 3 ), xk( 3 ), alf( 3 ), r, frac, potn
    integer :: ix, iy, iz, k1, k2, k3, kk1, kk2, kk3, xiter, kiter


    if( myid .eq. 0 ) then
      open(unit=99,file='epsilon',form='formatted', status='old' )
      rewind( 99 )
      read(99,*) epsi
      close( 99 ) 
      epsi = 1_DP / epsi

      open(unit=99,file='rpottrim',form='formatted',status='old' )
      rewind( 99 ) 
      do i = 1, 100
        read(99,*) ptab( i )
      enddo
      close( 99 )

      open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
      rewind( 99 )
      read( 99, * ) avec( :, : )
      close( 99 )
      do i = 1, 3
        do j = 1, 3
         amet( i , j ) = dot_product( avec( :, i ), avec( :, j ) )
        enddo
      enddo
      
    endif

#ifdef MPI    
    call MPI_BCAST( epsi, 1, MPI_DOUBLE, 0, comm, ierr ) 
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( ptab, 100, MPI_DOUBLE, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( amet, 9, MPI_DOUBLE, 0, comm, ierr )
#endif


    ! Slow way to start

    xiter = 0
    do iz = 1, sys%xmesh( 3 )
      fr( 3 ) = dble( iz - 1 ) / dble( sys%xmesh( 3 ) )
      do iy = 1, sys%xmesh( 2 )
        fr( 2 ) = dble( iy - 1 ) / dble( sys%xmesh( 2 ) )
        do ix = 1, sys%xmesh( 1 )
          fr( 1 ) = dble( ix - 1 ) / dble( sys%xmesh( 1 ) )
          xiter = xiter + 1
          kiter = 0
          if( ( xiter .ge. lr%my_start_nx ) .and. ( xiter .lt. lr%my_start_nx + lr%my_nxpts ) ) then
          do k3 = 1, sys%kmesh( 3 )
            kk3 = k3 - 1
            if ( kk3 .ge. sys%kmesh( 3 ) / 2 ) kk3 = kk3 - sys%kmesh( 3 )
            xk( 3 ) = kk3
              do k2 = 1, sys%kmesh( 2 )
                kk2 = k2 - 1
                if ( kk2 .ge. sys%kmesh( 2 ) / 2 ) kk2 = kk2 - sys%kmesh( 2 )
                xk( 2 ) = kk2
                do k1 = 1, sys%kmesh( 1 )
                  kk1 = k1 - 1
                  if ( kk1 .ge. sys%kmesh( 1 ) / 2 ) kk1 = kk1 - sys%kmesh( 1 )
                  xk( 1 ) = kk1
                  kiter = kiter + 1
                  alf( : ) = xk( : ) + fr( : ) - lr%tau( : )
                  r = sqrt( dot_product( alf, matmul( amet, alf ) ) )
                  if ( r .ge. 9.9d0 ) then
                     potn = epsi / r
                  else
                     ii = 1.0d0 + 10.0d0 * r
                     frac = 10.d0 * ( r - 0.1d0 * dble( ii - 1 ) )
                     potn = ptab( ii ) + frac * ( ptab( ii + 1 ) - ptab( ii ) )
                  end if
                  lr%W( kiter, xiter - lr%my_start_nx + 1 ) =  potn
                end do
              end do
            end do
          endif
        enddo
      enddo
    enddo

    

111 continue

  end subroutine lr_populate_W



  subroutine lr_populate_bloch( lr, ierr )
    use OCEAN_mpi, only : myid, nproc, comm

    type( long_range ), intent( inout ) :: lr
    integer, intent( inout ) :: ierr
    !
    integer, parameter :: u2dat = 35
    !
    integer :: nx, ny, nz, nbd, nq, zn( 3 ), nspn
    real( kind = kind(1.0d0)), intent( inout ) :: tau(3) 
    real( kind = kind( 1.0d0 ) ), dimension( nx, ny, nz, nbd ) :: ur, ui
    !
    integer :: iq, ibd, ig, idum( 3 ), ix, iy, iz, ivl, ivh, icl, ich, ispn
    integer :: iq1, iq2, iq3, dumint, icl2, ich2, ivh2, xshift(3)
    integer :: xtarg, ytarg, ztarg, xph, yph, zph
    real( kind = kind( 1.0d0 ) ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
    real( kind = kind( 1.0d0 ) ) :: su, sul, suh
    real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: tmp_ur, tmp_ui
    logical :: metal, normal

      

    ! As per usual, do this the dumbest way first, 
    ! 1) generic copy of original serial method 
    ! 2) MPI

      

    if( myid .eq. 0 ) then
      sul = 1.0d0
      suh = 1.0d0
      open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) ivl, ivh, icl, ich
      close( unit=99 )
      ivh2 = ivh
    !  icl2 = icl
      open( unit=99, file='metal', form='formatted', status='old')
      read( 99, * ) metal
      close( 99 )
      if( metal ) then
        open( unit=36, file='ibeg.h', form='formatted', status='old' )
      endif
      
      if ( nbd .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'
      open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
      rewind u2dat
      write(6,*) 'nspn: ', nspn
      if( nspn .ne. 1 ) then
        write(6,*) 'No spin yet'
        ierr = 1
        goto 111
      endif

      open(unit=99,file='temp_tau',form='fortmatted', status='old')
      rewind( 99 )
      read(99,*) tau(:)
      close(99)

    endif

    

!    do iq = 1, sys%nkpts
    iq = 0
    do iq1 = 1, zn( 1 )
     do iq2 = 1, zn( 2 )
      do iq3 = 1, zn( 3 )
      iq = iq + 1

      if( myid .eq. 0 ) then
        open( unit=99, file='gumatprog', form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(2i8)' ) iq, nq
        close( unit=99 )
        if( metal ) then
          read( 36, * ) dumint, ivh2
          ivh2 = ivh2 - 1 
        endif

  !  Skip all of the occupied bands (and for metals)
        do ibd = ivl, ivh2
          do ig = 1, nx * ny * nz
             read ( u2dat )
          end do 
        end do
        do ibd = 1, nbd
          do ix = 1, nx
             do iy = 1, ny
                do iz = 1, nz
                   read ( u2dat ) idum( 1 : 3 ), ur( ix, iy, iz, ibd ), ui( ix, iy, iz, ibd )
                end do 
             end do
          end do
          su = sum( ur( :, :, :, ibd ) ** 2 + ui( :, :, :, ibd ) ** 2 )
          sul = min( su, sul )
          suh = max( su, suh )
        end do
  ! Adding 22 nov 2010 get rid of the un-used wfns at the top
        do ibd = ivh2 + nbd + 1, ivh - ivl + ich - icl + 2
          do ig = 1, nx * ny * nz
             read ( u2dat )
          end do
        enddo

        !
        pi = 4.0d0 * atan( 1.0d0 )
        !
        ! Reasoning for the below;
        !  If the kmesh is odd, then the FT of the kmesh will give an equal number of
        !  cells on either side of the central one -- the one with the core hole. 
        !  This means that the core hole should be placed as near as (0.5,0.5,0.5)
        !
        !  On the other hand, if the kmesh is even, the bonus cell will be tacked on 
        !  over toward the negative side. This means that we should center the core 
        !  hole as near as possible to (0,0,0)
        !
        !  We test each dimension individually of course. 
        !
        xshift(:) = 0
        if( mod( zn(1), 2 ) .eq. 0 ) then
          xshift( 1 ) = floor( real(nx, kind( 1.0d0 )) * tau(1) )
        else
          xshift( 1 ) = floor( real(nx, kind( 1.0d0 )) * (tau(1)-0.5d0 ) )
        endif
        if( mod( zn(2), 2 ) .eq. 0 ) then
          xshift( 2 ) = floor( real(ny, kind( 1.0d0 )) * tau(2) )
        else
          xshift( 2 ) = floor( real(ny, kind( 1.0d0 )) * (tau(2)-0.5d0 ) )
        endif
        if( mod( zn(3), 2 ) .eq. 0 ) then
          xshift( 3 ) = floor( real(nz, kind( 1.0d0 )) * tau(3) )
        else
          xshift( 3 ) = floor( real(nz, kind( 1.0d0 )) * (tau(3)-0.5d0 ) )
        endif
        ! 
        write(6,*) 'Shifting X-grid by ', xshift(:)
        write(6,*) 'Original tau ', tau(:)
        tau( 1 ) = tau(1) - real(xshift(1), kind( 1.0d0 ))/real(nx, kind( 1.0d0 ))
        tau( 2 ) = tau(2) - real(xshift(2), kind( 1.0d0 ))/real(ny, kind( 1.0d0 ))
        tau( 3 ) = tau(3) - real(xshift(3), kind( 1.0d0 ))/real(nz, kind( 1.0d0 ))
        write(6,*) 'New tau      ', tau(:)

        allocate( tmp_ur( nx, ny, nz, nbd ), tmp_ui( nz, ny, nz, nbd ) )
        do iz = 1, nz
           ztarg = iz - xshift( 3 )
           if( ztarg .gt. nz ) then
             ztarg = ztarg - nz
             zph = -nz
           elseif( ztarg .lt. 1 ) then
             ztarg = ztarg + nz
             zph = nz
           else
             zph = 0
           endif
            do iy = 1, ny
               ytarg = iy - xshift( 2 )
               if( ytarg .gt. ny ) then
                 ytarg = ytarg - ny
                 yph = -ny
               elseif( ytarg .lt. 1 ) then
                 ytarg = ytarg + ny
                 yph = ny
               else
                 yph = 0
               endif
               do ix = 1, nx
                 xtarg = ix - xshift( 1 )
                 if( xtarg .gt. nx ) then
                   xtarg = xtarg - nx
                   xph = -nx
                 elseif( xtarg .lt. 1 ) then
                   xtarg = xtarg + nx
                   xph = nx
                 else
                   xph = 0
                 endif
                  phsx = 2.0d0 * pi * dble( ( xph + ix - 1 ) * ( iq1 - 1 ) ) / dble( nx * zn( 1 ) )
                  phsy = 2.0d0 * pi * dble( ( yph + iy - 1 ) * ( iq2 - 1 ) ) / dble( ny * zn( 2 ) )
                  phsz = 2.0d0 * pi * dble( ( zph + iz - 1 ) * ( iq3 - 1 ) ) / dble( nz * zn( 3 ) )
                  cphs = cos( phsx + phsy + phsz )
                  sphs = sin( phsx + phsy + phsz )
                  do ibd = 1, nbd
                     psir = cphs * ur( ix, iy, iz, ibd ) - sphs * ui( ix, iy, iz, ibd )
                     psii = cphs * ui( ix, iy, iz, ibd ) + sphs * ur( ix, iy, iz, ibd )
                     tmp_ur( xtarg, ytarg, ztarg, ibd ) = psir
                     tmp_ui( xtarg, ytarg, ztarg, ibd ) = psii
                  end do
               end do
            end do
         end do
!         ur( :, :, :, : ) = tmp_ur( :, :, :, : )
!         ui( :, :, :, : ) = tmp_ui( :, :, :, : )
         xiter = 0
         do iz = 1, nz
           do iy = 1, ny
             do ix = 1, nx
               xiter = xiter + 1
               tmp_bloch( :, xiter ) = cmplx( ur( ix, iy, iz, : ), ui( ix, iy, iz : ) )
             enddo
           enddo
         enddo
       endif


       nx_left = sys%nxpts
       nx_start = 1
       do i = 0, nproc - 1
         nx_tmp = nx_left / ( nproc - i )
         nx_left = nx_left - nx_tmp
         if( i .eq. 0 ) then
           lr%bloch_states( iq, :, : ) = tmp_bloch( :, 1 : nx_tmp )
         elseif( myid .eq. 0 ) then
           call MPI_SEND( tmp_bloch(1,nx_start), nbd*nx_tmp, MPI_DOUBLE_COMPLEX, i, i, comm, ierr )
         elseif( myid .eq. i ) then
           call MPI_RECV( tmp_bloch, nbd*nx_tmp, MPI_DOUBLE_COMPLEX, 0, i, comm, my_status, ierr )
           lr%bloch_states( iq, :, : ) = tmp_bloch( :, 1 : nx_tmp )
         endif
!         if( i .eq. myid ) then
!           lr%my_nxpts = nx_tmp
!           lr%my_start_nx = nx_start
!         endif
         nx_start = nx_start + nx_tmp
       enddo



      enddo
     enddo
    enddo




    if( myid .eq. 0 ) write ( 6, '(1a16,2f20.15)' ) 'norm bounds ... ', sul, suh

#ifdef MPI
    call MPI_BCAST(tau, 3, MPI_DOUBLE, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
#endif
    lr%tau( : ) = tau( : )





111 continue



  end subroutine lr_populate_bloch
end module ocean_long_range

