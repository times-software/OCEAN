module OCEAN_multiplet
  use AI_kinds

  private
  save

  real( DP ), pointer, contiguous :: mpcr( :, :, :, :, :, : )
  real( DP ), pointer, contiguous :: mpci( :, :, :, :, :, : )
  real( DP ), pointer, contiguous :: mpm( :, :, : )
  real( DP ), pointer, contiguous :: mhr( : )
  real( DP ), pointer, contiguous :: mhi( : )
  real( DP ), pointer :: cms( : ), cml( : ), vms( : )
  real( DP ), pointer :: hcml( : ), hvml( : ), hcms( : ), hvms( : )
  real( DP ), pointer, contiguous :: somelr( :, : )
  real( DP ), pointer, contiguous :: someli( :, : )


  integer, pointer :: nproj( : ), hvnu(:)
  integer, pointer :: ibeg( : )
  integer, pointer :: jbeg( : )
  integer, pointer :: mham( : )

  real(DP ) :: xi

  integer :: npmax
  integer :: nptot
  integer :: lmin
  integer :: lmax
  integer :: lvl
  integer :: lvh
  integer :: jtot
  integer :: itot

  logical :: is_init = .false.

  public OCEAN_create_central, OCEAN_soprep, OCEAN_mult_act

  contains


  subroutine OCEAN_create_central( sys, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    
    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: lv, ic, icms, icml, ivms, ii, ivml, jj, nu, i, iband, ikpt, ip, ispn
    integer :: l, m, tmp_n
    real( DP ), allocatable :: pcr(:,:,:,:), pci(:,:,:,:), list(:)
    
    character(len=4) :: add04
    character(len=10) :: add10
    character(len=12) :: add12
    character(len=11) :: s11
    character(len=5) :: s5

    character(len=11) :: cks_filename
    character(len=5) :: cks_prefix

    ierr = 0
    if( is_init .and. associated( sys%cur_run%prev_run ) ) then
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_init = .false.
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_init = .false. 
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_init = .false. 
      if( sys%cur_run%indx .ne. sys%cur_run%prev_run%indx ) is_init = .false.
!      if( ( myid .eq. root ) .and. ( .not. is_init ) ) then
      if( .not. is_init )  then
        if( myid .eq. root ) write(6,*) 'Multiplets are reloading'
        deallocate( mpcr, mpci, mpm, mhr, mhi, cms, cml, vms, hcml, hvml, hcms, hvms, somelr, someli )
        deallocate( nproj, hvnu, ibeg, jbeg, mham )
      endif
    endif

    if( is_init ) return

    write(add12 , '(A2,I2.2,A1,A1,I2.2,A1,I2.2)' ) sys%cur_run%elname, sys%cur_run%indx, &
            '_', 'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)

    if( myid .eq. root ) then
!      write(deflinz,'(a6,a1,i2.2,a1,i2.2,a1,i2.2)') 'deflin', 'z', sys%cur_run%ZNL(1), 'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)
!      open(unit=99,file=deflinz,form='formatted',status='old')
!      open(unit=99,file='derp_control',form='formatted', status='old' )
      open(unit=99,file='bse.in',form='formatted', status='old')
      rewind 99
      read( 99, * ) tmp_n, lvl, lvh
      read(99,*) xi
      close( 99 )
      xi = xi / 27.2114d0

      write( add04, '(1a1,1i3.3)' ) 'z', sys%cur_run%ZNL(1)
      write( add10, '(1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'z', sys%cur_run%ZNL(1), &
                  'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3)
      write ( s11, '(1a7,1a4)' ) 'prjfile', add04
      write(6,*) s11
      open( unit=99, file=s11, form='formatted', status='old' )
      rewind 99
      read ( 99, * ) lmin, lmax
      allocate( nproj( lmin : lmax ) )
      read ( 99, * ) nproj
      close( unit=99 )
      if ( ( lvl .lt. lmin ) .or. ( lvh .gt. lmax ) ) stop 'bad lvl lvh'

      allocate( ibeg( lvl : lvh ), jbeg( lvl : lvh ), mham( lvl : lvh ) )
      itot = 0
      jtot = 0
      do lv = lvl, lvh
         ibeg( lv ) = itot + 1
         jbeg( lv ) = jtot + 1
         mham( lv ) = 4 * ( 2 * sys%cur_run%ZNL(3) + 1 ) * ( 2 * lv + 1 ) * nproj( lv )
         itot = itot + mham( lv )
         jtot = jtot + mham( lv ) ** 2
      end do

      allocate( cml( sys%cur_run%nalpha ), cms( sys%cur_run%nalpha ), vms( sys%cur_run%nalpha ) )
      allocate( hcml( itot ), hcms( itot ) )
      allocate( hvml( itot ), hvms( itot ) )
      allocate( hvnu( itot ) )
      allocate( mhr( jtot ), mhi( jtot ) )
!      allocate( pwr( itot ), pwi( itot ) )
!      allocate( hpwr( itot ), hpwi( itot ) )

      open( unit=99, file='channelmap', form='formatted', status='unknown' )
      rewind 99
      ic = 0
      do icms = -1, 1, 2
         do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
            do ivms = -1, 1, 2
               ic = ic + 1
               cms( ic ) = 0.5d0 * icms
               cml( ic ) = icml
               vms( ic ) = 0.5d0 * ivms
               write ( 99, '(3f8.2,2i5,1a11)' ) cms( ic ), cml( ic ), vms( ic ), ic, sys%cur_run%nalpha, 'cms cml vms'
            end do
         end do
      end do
      !
      ii = 0
      do lv = lvl, lvh
         do ic = 1, 4 * ( 2 * sys%cur_run%ZNL(3) + 1 )
            do ivml = -lv, lv
               do nu = 1, nproj( lv )
                  ii = ii + 1
                  hcms( ii ) = cms( ic )
                  hcml( ii ) = cml( ic )
                  hvms( ii ) = vms( ic )
                  hvml( ii ) = ivml
                  hvnu( ii ) = nu
               end do
            end do
         end do
      end do
      !
      call nbsemkcmel( add04, add12 )
      do lv = lvl, lvh
         ii = ibeg( lv )
         jj = jbeg( lv )
!         call nbsemhsetup( sys%cur_run%ZNL(3), lv, nproj( lv ), mham( lv ), hcms( ii:ii+mham( lv )-1 ) , &
!              hcml( ii:ii+mham( lv )-1 ), hvms( ii:ii+mham( lv )-1 ), hvml( ii:ii+mham( lv )-1 ), hvnu( ii:ii+mham( lv )-1 ), &
!              mhr( jj:jj+mham( lv )*mham( lv )-1 ), mhi( jj:jj+mham( lv )*mham( lv )-1 ), add10 )
          call nbsemhsetup( sys%cur_run%ZNL(3), lv, nproj( lv ), mham( lv ), hcms( ii ) , hcml( ii ), hvms( ii ), hvml( ii ), hvnu( ii ), &
          mhr( jj ), mhi( jj ), add10 )
      end do
      mhr = mhr / 27.2114d0
      mhi = mhi / 27.2114d0
      write ( 6, * ) 'multiplet hamiltonian set up'
      write(6,*) 'Number of spins = ', sys%nspn
      write ( 6, * ) 'n, nc, nspn', sys%num_bands*sys%nkpts, sys%cur_run%nalpha, sys%nspn


      npmax = maxval( nproj( lmin : lmax ) )
      allocate( mpcr( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
      allocate( mpci( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
!      allocate( ampr( npmax ), ampi( npmax ), hampr( npmax ), hampi( npmax ) )
      ! 
      ip = 0
      do l = lmin, lmax
         do m = -l, l
            do nu = 1, nproj( l )
               ip = ip + 1
            end do 
         end do
      end do
      nptot = ip
    ! Spin needed here too!
      allocate( pcr( nptot, sys%num_bands, sys%nkpts, sys%nspn ), pci( nptot, sys%num_bands, sys%nkpts, sys%nspn ) )
!JTV
      select case ( sys%cur_run%calc_type)
      case( 'XES' )
        cks_prefix = 'cksv.'
      case( 'XAS' )
        cks_prefix = 'cksc.'
      case default
        cks_prefix = 'cksc.'
      end select
      write(cks_filename,'(A5,A2,I4.4)' ) cks_prefix, sys%cur_run%elname, sys%cur_run%indx
      open(unit=99,file=cks_filename,form='unformatted', status='old' )
!      open( unit=99, file='ufmi', form='unformatted', status='old' )
      rewind 99
      read ( 99 )
      read ( 99 )
      read ( 99 ) pcr 
      read ( 99 ) pci
      do ispn = 1, sys%nspn
        do ikpt = 1, sys%nkpts
        do iband = 1, sys%num_bands 
          ip = 0
          do l = lmin, lmax
            do m = -l, l
               do nu = 1, nproj( l )
                  ip = ip + 1
                  mpcr( iband, ikpt, nu, m, l, ispn ) = pcr( ip, iband, ikpt, ispn )
                  mpci( iband, ikpt, nu, m, l, ispn ) = pci( ip, iband, ikpt, ispn )
               end do 
            end do
          end do
        end do
        enddo
      enddo
      close( unit=99 )
      deallocate( pcr, pci )
      write ( 6, * ) 'projector coefficients have been read in'
      !
      allocate( list( npmax * npmax ), mpm( npmax, npmax, lmin : lmax ) )
      mpm( :, :, : ) = 0.0_DP
      do l = lmin, lmax
         if ( l .gt. 9 ) stop 'l exceeds 9 in multip'
         write ( s5, '(1a4,1i1)' ) 'cmel', l
         open( unit=99, file=s5, form='formatted', status='old' )
         rewind 99
         read( 99, * ) list( 1: nproj( l ) * nproj( l ) )
         i = 0
         do nu = 1, nproj( l ) 
            mpm( 1:nproj(l), nu, l ) = list( i + 1 : i + nproj( l ) )
            i = i + nproj( l )
         end do
         close( unit=99 )
      end do
      mpm( :, :, : ) = mpm( :, :, : ) / 27.2114d0
      deallocate( list )
      write ( 6, * ) 'central projector matrix elements have been read in'
!      write(100,*) mpcr(:,1,1,0,0,1)
!      write(101,*) mpci(1,1,1,0,0,1)
!      write(102,*) mpm(:,1,0)

    endif

#ifdef MPI
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( itot, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( jtot, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( npmax, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lmin, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lmax, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lvl, 1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( lvh, 1, MPI_INTEGER, root, comm, ierr )


    if( myid .ne. root ) then
      allocate( cml( sys%cur_run%nalpha ), cms( sys%cur_run%nalpha ), vms( sys%cur_run%nalpha ) )
      allocate( hcml( itot ), hcms( itot ) )
      allocate( hvml( itot ), hvms( itot ) )
      allocate( hvnu( itot ) )
      allocate( mhr( jtot ), mhi( jtot ) )
      allocate( nproj( lmin : lmax ) )
      allocate( ibeg( lvl : lvh ), jbeg( lvl : lvh ), mham( lvl : lvh ) )
      allocate( mpcr( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
      allocate( mpci( sys%num_bands, sys%nkpts, npmax, -lmax : lmax, lmin : lmax, sys%nspn ) )
      allocate( mpm( npmax, npmax, lmin : lmax ) )
    endif
    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write(6,*) 'HERE'
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( cml, sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( cms, sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( vms, sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hcml, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hcms, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hvml, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hvms, itot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( hvnu, itot, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( mhr, jtot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( mhi, jtot, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( nproj( lmin ), lmax-lmin+1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( ibeg( lvl ), lvh-lvl+1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( jbeg( lvl ), lvh-lvl+1, MPI_INTEGER, root, comm, ierr )
    call MPI_BCAST( mham( lvl ), lvh-lvl+1, MPI_INTEGER, root, comm, ierr )

    call MPI_BCAST( mpcr, sys%num_bands * sys%nkpts * npmax * (2*lmax+1) * (lmax-lmin+1) *sys%nspn, &
                    MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( mpci, sys%num_bands * sys%nkpts * npmax * (2*lmax+1) * (lmax-lmin+1) *sys%nspn, &
                    MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( mpm, npmax*npmax*(lmax-lmin+1), MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'multiplet hamilotian shared'
#endif
  
    call OCEAN_soprep( sys, ierr )

    is_init = .true.

  end subroutine OCEAN_create_central

  subroutine OCEAN_soprep( sys, ierr )
    use OCEAN_mpi
    use OCEAN_system
    implicit none
    
    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: ic, jc, i
    complex( kind = kind( 1.0d0 ) ) :: ctmp, rm1, vrslt( 3 )
    complex( kind = kind( 1.0d0 ) ), external :: jimel
    real( kind = kind( 1.d0 ) )  :: life_time( 2 ), l_alpha, l_beta, delta_so( 2 )
    logical :: broaden_exist
    !
    real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
    integer :: nsphpt, isphpt
    real( kind = kind( 1.0d0 ) ) :: sphsu
    real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: xsph, ysph, zsph, wsph

    ierr = 0
    allocate( somelr( sys%cur_run%nalpha, sys%cur_run%nalpha), someli( sys%cur_run%nalpha, sys%cur_run%nalpha) )

    if( myid .eq. root ) then

      open( unit=99, file='sphpts', form='formatted', status='old' )
      rewind 99
      read ( 99, * ) nsphpt
      allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
      do isphpt = 1, nsphpt
         read ( 99, * ) xsph( isphpt ), ysph( isphpt ), zsph( isphpt ), wsph( isphpt )
      end do
      close( unit=99 )
      sphsu = sum( wsph( : ) )
      wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / sphsu )
      write ( 6, * ) nsphpt, ' points with weights summing to four pi '

      inquire( file='core_broaden.ipt', exist= broaden_exist )
      if( broaden_exist ) then
        write(6,*) 'spin-orbit dependent broadening'
        open( unit=99, file='core_broaden.ipt', form='formatted', status='old' )
        read( 99, * ) life_time( : )
        close( 99 )
        life_time( : ) = life_time( : ) / 27.2114d0

        delta_so( 1 ) = 0.5d0 * ( ( real( sys%cur_run%ZNL(3) ) + 0.5 ) * ( real( sys%cur_run%ZNL(3) ) + 1.5 ) - & 
                                    real( sys%cur_run%ZNL(3) ) * real( sys%cur_run%ZNL(3) + 1 ) - 0.75d0 )
        delta_so( 2 ) = 0.5d0 * ( ( real( sys%cur_run%ZNL(3) ) - 0.5 ) * ( real( sys%cur_run%ZNL(3) ) + 0.5 ) - &
                                    real( sys%cur_run%ZNL(3) ) * real( sys%cur_run%ZNL(3) + 1 ) - 0.75d0 )
        !
        l_alpha = ( life_time( 2 ) - life_time( 1 ) ) / ( - xi * ( delta_so( 2 ) - delta_so( 1 ) ) )
        l_beta = ( life_time( 2 ) * delta_so( 1 ) - life_time( 1 ) * delta_so( 2 ) ) / ( delta_so( 1 ) - delta_so( 2 ) )
        write(6,*) l_alpha, l_beta
      else
        l_alpha = 0.d0
        l_beta  = 0.d0
      endif
      somelr( :, : ) = 0.0_DP
      someli( :, : ) = 0.0_DP
      rm1 = -1; rm1 = sqrt( rm1 )
      do ic = 1, sys%cur_run%nalpha
         do jc = 1, sys%cur_run%nalpha
            if ( vms( ic ) .eq. vms( jc ) ) then
               call limel( sys%cur_run%ZNL(3), nint( cml( jc ) ), nint( cml( ic ) ), vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
               ctmp = 0
               do i = 1, 3
                  ctmp = ctmp + vrslt( i ) * jimel( 0.5d0, cms( jc ), cms( ic ), i )
               end do
               write(20,*) ic, jc, real( ctmp )
               write(21,*)
               ctmp = -xi * ctmp
               somelr( ic, jc ) = real( ctmp ) - aimag( ctmp ) * l_alpha
               someli( ic, jc ) = -aimag( ctmp ) - real( ctmp ) * l_alpha
               if( ic .eq. jc ) then
                 someli( ic, jc ) = someli( ic, jc ) - l_beta !* real( ctmp )
               endif
             end if
         end do
      end do
    endif

#ifdef MPI
    if( nproc .gt. 1 ) then
      call MPI_BCAST( someli, sys%cur_run%nalpha*sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
      call MPI_BCAST( somelr, sys%cur_run%nalpha*sys%cur_run%nalpha, MPI_DOUBLE_PRECISION, root, comm, ierr )
    endif
#endif

  end subroutine OCEAN_soprep
    
  
  subroutine OCEAN_mult_act( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec

    out_vec%r(:,:,:) = 0.0_DP
    out_vec%i(:,:,:) = 0.0_DP
! !$OMP PARALLEL DEFAULT( PRIVATE ) SHARED( sys, inter, in_vec, out_vec )
    call OCEAN_ctact( sys, inter, in_vec, out_vec )
! !$OMP END PARALLEL
    call fgact( sys, inter, in_vec, out_vec )
!    call OCEAN_soact( sys, in_vec, out_vec )


  end subroutine OCEAN_mult_act


  subroutine OCEAN_ctact( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec
    !
!    integer :: sys%cur_run%nalpha, n, nq, nbd, nspn, lmin, lmax, npmax
!    real( DP ) :: celvol, inter
!    integer :: nproj( lmin : lmax )
!    real( DP ) :: mpcr( nq * nbd, npmax, -lmax : lmax, lmin : lmax, nspn )
!    real( DP ) :: mpci( nq * nbd, npmax, -lmax : lmax, lmin : lmax, nspn )
!    real( DP ) :: mpm( npmax, npmax, lmin : lmax )
    !
    integer :: ialpha, l, m, nu, ispn, ikpt, ibnd
    real( DP ) :: mul
    real( DP ), dimension( npmax ) :: ampr, ampi, hampr, hampi
  !


!    write(103,*) mpm(:,1,0)
!    write(104,*) mpcr(:,1,1,0,0,1)

    
    mul = inter / ( dble( sys%nkpts ) * sys%celvol )
!    mul = 1.0_DP / ( dble( sys%nkpts ) * sys%celvol )
!$OMP PARALLEL DO &
!$OMP DEFAULT( NONE ) &
!$OMP PRIVATE( ialpha, l, nu, m, ampr, ampi, hampr, hampi, ispn ) &
!$OMP SHARED( mul, sys, lmin, lmax, nproj ) &
!$OMP SHARED( mpcr, mpci, mpm, in_vec, out_vec )
    do ialpha = 1, sys%cur_run%nalpha
      ispn = 2 - mod( ialpha, 2 )
      if( sys%nspn .eq. 1 ) then
        ispn = 1
      endif
      do l = lmin, lmax
        do m = -l, l
          ampr( : ) = 0
          ampi( : ) = 0
!             do i = 1, n
          do nu = 1, nproj( l )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                ampr( nu ) = ampr( nu ) + in_vec%r( ibnd, ikpt, ialpha ) * mpcr( ibnd, ikpt, nu, m, l, ispn ) & 
                                        - in_vec%i( ibnd, ikpt, ialpha ) * mpci( ibnd, ikpt, nu, m, l, ispn )
                ampi( nu ) = ampi( nu ) + in_vec%r( ibnd, ikpt, ialpha ) * mpci( ibnd, ikpt, nu, m, l, ispn ) & 
                                        + in_vec%i( ibnd, ikpt, ialpha ) * mpcr( ibnd, ikpt, nu, m, l, ispn )
               enddo
            enddo
          enddo


          hampr( : ) = 0
          hampi( : ) = 0
          do nu = 1, nproj( l )
            hampr( : ) = hampr( : ) - mpm( :, nu, l ) * ampr( nu ) 
            hampi( : ) = hampi( : ) - mpm( :, nu, l ) * ampi( nu ) 
          enddo
          hampr( : ) = hampr( : ) * mul
          hampi( : ) = hampi( : ) * mul
          do nu = 1, nproj( l )
            do ikpt = 1, sys%nkpts 
              do ibnd = 1, sys%num_bands
                out_vec%r( ibnd, ikpt, ialpha ) = out_vec%r( ibnd, ikpt, ialpha ) &
                                                + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )  &
                                                + mpci( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu )
                out_vec%i( ibnd, ikpt, ialpha ) = out_vec%i( ibnd, ikpt, ialpha ) &
                                                + mpcr( ibnd, ikpt, nu, m, l, ispn ) * hampi( nu ) &
                                                - mpci( ibnd, ikpt, nu, m, l, ispn ) * hampr( nu )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine OCEAN_ctact



  subroutine fgact( sys, inter, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none
    !
    type(O_system), intent( in ) :: sys
    real( DP ), intent( in ) :: inter
    type( ocean_vector ), intent( in ) :: in_vec
    type( ocean_vector ), intent( inout ) :: out_vec
    !
!    real( DP ), dimension( n, nc, 2 ) :: v, hv
!    real( DP ), dimension( n, npmax, -lmax : lmax, lmin : lmax, nspn ) :: mpcr, mpci
!    real( DP ), dimension( jtot ) :: mhr, mhi
    !
    integer :: lv, ii, ic, ivml, nu, j1, jj, ispn, ikpt, ibnd
    real( DP ) :: mul
    real( DP ), allocatable, dimension( : ) :: pwr, pwi, hpwr, hpwi
    !
    allocate( pwr(itot), pwi( itot), hpwr(itot), hpwi(itot) )
    ! lvl, lvh is not enough for omp
    ! should probably pull thread spawning out of the do loop though

    mul = inter / ( dble( sys%nkpts ) * sys%celvol )
!    write(6,*) mul
    do lv = lvl, lvh
       pwr( : ) = 0
       pwi( : ) = 0
       hpwr( : ) = 0
       hpwi( : ) = 0
!$OMP PARALLEL &
!$OMP PRIVATE( ic, ivml, nu, ii, jj, j1, ispn ) &
!$OMP SHARED( mul, nproj, lv, mham, jbeg, pwr, pwi, mhr, mhi, mpcr, mpci, hpwr, hpwi, in_vec, out_vec, sys ) &
!$OMP DEFAULT( NONE )

!$OMP DO 
      do ic = 1, sys%cur_run%nalpha
        ispn = 2 - mod( ic, 2 )
  !       ispn = 1 + mod( ic, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
            ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                pwr( ii ) = pwr( ii ) &
                          + in_vec%r( ibnd, ikpt, ic ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                          - in_vec%i( ibnd, ikpt, ic ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn ) 
                pwi( ii ) = pwi( ii ) &
                          + in_vec%r( ibnd, ikpt, ic ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn ) &
                          + in_vec%i( ibnd, ikpt, ic ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) 
              enddo
             enddo
           enddo
         enddo
       enddo
!$OMP END DO


!$OMP DO
      do ii = 1, mham( lv )
        do jj = 1, mham( lv )
          j1 = jbeg( lv ) + ( jj -1 ) + ( ii - 1 ) * mham( lv )
          hpwr( ii ) = hpwr( ii ) + mhr( j1 ) * pwr( jj ) - mhi( j1 ) * pwi( jj )
          hpwi( ii ) = hpwi( ii ) + mhr( j1 ) * pwi( jj ) + mhi( j1 ) * pwr( jj )
        end do
        hpwr( ii ) = hpwr( ii ) * mul
        hpwi( ii ) = hpwi( ii ) * mul
      end do
!$OMP END DO


!$OMP DO
      do ic = 1, sys%cur_run%nalpha
        ispn = 2 - mod( ic, 2 )
  !       ispn = 1 + mod( ic, 2 )
        if( sys%nspn .eq. 1 ) then
          ispn = 1
        endif
        do ivml = -lv, lv
          do nu = 1, nproj( lv )
            ii = nu + ( ivml + lv ) * nproj( lv ) + ( ic - 1 ) * ( 2 * lv + 1 ) * nproj( lv )
            do ikpt = 1, sys%nkpts
              do ibnd = 1, sys%num_bands
                out_vec%r( ibnd, ikpt, ic ) = out_vec%r( ibnd, ikpt, ic )  &
                                            + hpwr( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                                            + hpwi( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
                out_vec%i( ibnd, ikpt, ic ) = out_vec%i( ibnd, ikpt, ic )  &
                                            + hpwi( ii ) * mpcr( ibnd, ikpt, nu, ivml, lv, ispn ) &
                                            - hpwr( ii ) * mpci( ibnd, ikpt, nu, ivml, lv, ispn )
              enddo
            enddo
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL 
    end do
    !

  end subroutine fgact

  subroutine OCEAN_soact( sys, in_vec, out_vec )
    use OCEAN_system
    use OCEAN_psi
    implicit none

    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( in ) :: in_vec
    type( OCEAN_vector ), intent( inout ) :: out_vec

    integer :: ic, jc, ikpt
! !$OMP PARALLEL DO &
! !$OMP DEFAULT( NONE ) &
! !$OMP PRIVATE( ic, jc, ikpt ) &
! !$OMP SHARED( sys, in_vec, out_vec, somelr, someli )
    do ic = 1, sys%cur_run%nalpha
      do jc = 1, sys%cur_run%nalpha
        do ikpt = 1, sys%nkpts
          out_vec%r( 1:sys%num_bands, ikpt, ic ) = out_vec%r( 1:sys%num_bands, ikpt, ic )  & 
                                   + somelr( ic, jc ) * in_vec%r( 1:sys%num_bands, ikpt, jc ) &
                                   - someli( ic, jc ) * in_vec%i( 1:sys%num_bands, ikpt, jc )
          out_vec%i( 1:sys%num_bands, ikpt, ic ) = out_vec%i( 1:sys%num_bands, ikpt, ic )  &
                                   + somelr( ic, jc ) * in_vec%i( 1:sys%num_bands, ikpt, jc ) &
                                   + someli( ic, jc ) * in_vec%r( 1:sys%num_bands, ikpt, jc )
        enddo
      enddo
    enddo
! !$OMP END PARALLEL DO

  end subroutine OCEAN_soact


  subroutine nbsemhsetup2( lc, lv, np, mham, cms, cml, vms, vml, vnu, mhr, mhi, add10 )
    implicit none
    !
    integer :: lc, lv, np, mham
    real( kind = kind( 1.0d0 ) ) :: cms( mham ), cml( mham )
    real( kind = kind( 1.0d0 ) ) :: vms( mham ), vml( mham )
    real( kind = kind( 1.0d0 ) ) :: mhr( mham, mham ), mhi( mham, mham )
    integer :: vnu( mham )
    !
    integer :: npt
    real( kind = kind( 1.0d0 ) ) :: pi, su, yp( 0 : 1000 )
    complex( kind = kind( 1.0d0 ) ) :: rm1
    real( kind = kind( 1.0d0 ) ), allocatable :: x( : ), w( : )
    !
    integer :: kfl, kfh, kgl, kgh
    real( kind = kind( 1.0d0 ) ), allocatable :: fk( :, :, : ), scfk( : )
    real( kind = kind( 1.0d0 ) ), allocatable :: gk( :, :, : ), scgk( : )
    !
    integer :: i, i1, i2, nu1, nu2
    integer :: l1, m1, s1, l2, m2, s2, l3, m3, s3, l4, m4, s4, k, mk
    real( kind = kind( 1.0d0 ) ) :: ggk, ffk
    complex( kind = kind( 1.0d0 ) ) :: f1, f2, ctmp
    logical, parameter :: no = .false., yes = .true.
    logical :: tdlda
    !
    character * 10 :: add10
    character * 15 :: filnam
    !
    include 'sphsetnx.h.f90'
    !
    include 'sphsetx.h.f90'
    !
     write(6,*) 'nbsemhsetup', lv
    call newgetprefs( yp, max( lc, lv ), nsphpt, wsph, xsph, ysph, zsph )
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    !
    write ( 6, * ) ' add10 = ', add10 
    open( unit=99, file='Pquadrature', form='formatted', status='old' )
    rewind 99
    read ( 99, * ) npt
    allocate( x( npt ), w( npt ) )
    su = 0
    do i = 1, npt
       read ( 99, * ) x( i ), w( i )
       su = su + w( i )
    end do
    close( unit=99 )
    w = w * 2 / su
    !
    if ( lv .gt. 9 ) stop 'lv must be single digit'
    if ( lc .gt. 9 ) stop 'lc must be single digit'
    !
    ! TDLDA
    inquire( file='tdlda', exist=tdlda )
    if( tdlda ) then
      open(unit=99, file='tdlda', form='formatted', status='old' )
      read( 99, * ) tdlda
      close( 99 )
    endif
    kfh = min( 2 * lc, 2 * lv )
    kfl = 2
    if( tdlda ) then
       kfl = 0
       kfh = 0
       allocate( fk( np, np, kfl : kfh ), scfk( kfl : kfh ) )
       fk = 0; scfk = 0
       do k = kfl, kfh, 2
          write ( filnam, '(1a2,3i1,1a10)' ) 'kk', lc, lv, k, add10
          write( 6, * ) filnam
          open( unit=99, file=filnam, form='formatted', status='old' )
          rewind 99
          read ( 99, * ) fk( :, :, k ), scfk( k )
          close( unit=99 )
       end do
    else
    if ( kfh .ge. kfl ) then
       allocate( fk( np, np, kfl : kfh ), scfk( kfl : kfh ) )
       fk = 0; scfk = 0
       do k = kfl, kfh, 2
          write ( filnam, '(1a2,3i1,1a10)' ) 'fk', lc, lv, k, add10
          open( unit=99, file=filnam, form='formatted', status='old' )
          rewind 99
          read ( 99, * ) fk( :, :, k ), scfk( k )
          close( unit=99 )
       end do
    end if
    end if
    !
    kgh = lc + lv
    kgl = abs( lc - lv )
    if ( kgh .ge. kgl ) then
       allocate( gk( np, np, kgl : kgh ), scgk( kgl : kgh ) )
       gk = 0; scgk = 0
       do k = kgl, kgh, 2
          write ( filnam, '(1a2,3i1,1a10)' ) 'gk', lc, lv, k, add10
          open( unit=99, file=filnam, form='formatted', status='old' )
          rewind 99
          read ( 99, * ) gk( :, :, k ), scgk( k )
        close( unit=99 )
       end do
    end if
    !
    mhr = 0
    mhi = 0
    do i1 = 1, mham
       nu1 = vnu( i1 )
       do i2 = 1, mham
          nu2 = vnu( i2 )
          !
          if( kfh .ge. kfl ) then
             l1 = lc; m1 = nint( cml( i2 ) ); s1 = nint( 2 * cms( i2 ) )
             l2 = lv; m2 = nint( vml( i1 ) ); s2 = nint( 2 * vms( i1 ) )
             l3 = lc; m3 = nint( cml( i1 ) ); s3 = nint( 2 * cms( i1 ) )
             l4 = lv; m4 = nint( vml( i2 ) ); s4 = nint( 2 * vms( i2 ) )
             if ( ( s1 .eq. s3 ) .and. ( s2 .eq. s4 ) ) then
                mk = m1 - m3
                if ( m1 + m2 .eq. m3 + m4 ) then
                   do k = kfl, kfh, 2
                      if ( abs( mk ) .le. k ) then
                         ffk = scfk( k ) * fk( nu1, nu2, k ) 
                         call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                         call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
                         ctmp = - ffk * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                         mhr( i1, i2 ) = mhr( i1, i2 ) + real(ctmp)
                         mhi( i1, i2 ) = mhi( i1, i2 ) + aimag(ctmp) !- ctmp * rm1
                      end if
                   end do
                end if
             end if
          end if
          !
          if( kgh .ge. kgl ) then
             l1 = lc; m1 = nint( cml( i2 ) ); s1 = nint( 2 * cms( i2 ) )
             l2 = lv; m2 = nint( vml( i1 ) ); s2 = nint( 2 * vms( i1 ) )
             l3 = lv; m3 = nint( vml( i2 ) ); s3 = nint( 2 * vms( i2 ) )
             l4 = lc; m4 = nint( cml( i1 ) ); s4 = nint( 2 * cms( i1 ) )
             if ( ( s1 .eq. s3 ) .and. ( s2 .eq. s4 ) ) then
                mk = m1 - m3
                if ( m1 + m2 .eq. m3 + m4 ) then
                   do k = kgl, kgh, 2
                      if ( abs( mk ) .le. k ) then
                         ggk = scgk( k ) * gk( nu1, nu2, k ) 
                         call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                         call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
                         ctmp = ggk * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                         mhr( i1, i2 ) = mhr( i1, i2 ) + real(ctmp)
                         mhi( i1, i2 ) = mhi( i1, i2 ) + aimag(ctmp)! - ctmp * rm1
                      end if
                   end do
                end if
             end if
          end if
          !
       end do
    end do
    !
    return
  end subroutine nbsemhsetup2
end module OCEAN_multiplet
