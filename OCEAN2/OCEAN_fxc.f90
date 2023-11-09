! Copyright (C) 2023 OCEAN collaboration
!             
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!     
!   
! This is the contact attraction Kxc for adiabatic LDA
      
module OCEAN_fxc
    
  use AI_kinds
  implicit none
  private   
  save      
            
  logical :: MATCH_AI2NBSE = .true.
              
  real( DP ), allocatable  :: fxc( :, : )
    
          
  logical :: is_init = .false.
  
  public :: OCEAN_fxc_prep, OCEAN_fxc_act, OCEAN_fxc_clean


  contains

  subroutine OCEAN_fxc_clean( )
    if( allocated( fxc ) ) deallocate( fxc )
    is_init = .false.
  end subroutine OCEAN_fxc_clean
    

  subroutine OCEAN_fxc_prep( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi
    use OCEAN_val_states, only : startx, nxpts
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: rho(:,:,:)

    real(DP) :: n, dn, np2, np1, nz0, nm1, nm2
    real(DP) :: exp2, exp1, exz0, exm1, exm2
    real(DP) :: ecp2, ecp1, ecz0, ecm1, ecm2
    real(DP) :: ux1, ux2, uc1, uc2
    real(DP) :: valp2, valp1, valz0, valm1, valm2

    integer :: iix, ix, iy, iz
    
    if( is_init ) then
      if( myid .eq. root ) write(6,*) '  fxc cached'
      return
    endif


    allocate( rho( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1) ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_get_rho( sys%xmesh, sys%celvol, rho, ierr )
    if( ierr .ne. 0 ) return

    allocate( fxc( nxpts, 1 ), STAT=ierr ) 
    if( ierr .ne. 0 ) return

    iix = 0
    do ix = 1, sys%xmesh( 1 )
      do iy = 1, sys%xmesh( 2 )
        do iz = 1, sys%xmesh( 3 )
          iix = iix + 1
          if( iix .ge. startx .and. ( iix - startx .lt. nxpts ) ) then
             n = rho( iz, iy, ix )
!             su = su + n
             dn = 0.01_dp * n
             np2 = n + 2 * dn
             np1 = n + dn
             nz0 = n
             nm1 = n - dn
             nm2 = n - 2 * dn
             call cacorr( np2, exp2, ecp2, ux1, ux2, uc1, uc2 )
             call cacorr( np1, exp1, ecp1, ux1, ux2, uc1, uc2 )
             call cacorr( nz0, exz0, ecz0, ux1, ux2, uc1, uc2 )
             call cacorr( nm1, exm1, ecm1, ux1, ux2, uc1, uc2 )
             call cacorr( nm2, exm2, ecm2, ux1, ux2, uc1, uc2 )
             valp2 = np2 * ( exp2 + ecp2 )
             valp1 = np1 * ( exp1 + ecp1 )
             valz0 = nz0 * ( exz0 + ecz0 )
             valm1 = nm1 * ( exm1 + ecm1 )
             valm2 = nm2 * ( exm2 + ecm2 )
             fxc(iix-startx+1,1) = ( 16.0_dp * ( valm1 + valp1 ) - ( valm2 + valp2 ) - 30.0_dp * valz0 ) / ( 12.0_dp * dn ** 2 ) &
                                 * dble(sys%nxpts) / ( dble( sys%nkpts ) * sys%celvol )

          endif
        enddo
      enddo
    enddo

    if( .true. .and. myid .eq. root ) then
      open(unit=99,file='rhokxc.txt', form='formatted', status='unknown' )
      iix = 0
      do ix = 1, sys%xmesh( 1 )
        do iy = 1, sys%xmesh( 2 )
          do iz = 1, sys%xmesh( 3 )
            iix = iix + 1
            if( iix .ge. startx .and. ( iix - startx .lt. nxpts ) ) then
               n = rho( iz, iy, ix )
               write( 99, '(1p,2(1x,1e15.8))' ) n, fxc(iix-startx+1,1) * dble( sys%nkpts ) * sys%celvol /dble(sys%nxpts )
            endif
          enddo
        enddo
      enddo
      close( 99 )
    endif

    deallocate( rho )
    
    is_init = .true.

    call MPI_BARRIER( comm, ierr )
  end subroutine OCEAN_fxc_prep



  subroutine OCEAN_fxc_act(sys, psi, psiout, ierr )
    use OCEAN_val_states, only : nkpts, nxpts, nbc, nbv, nxpts_pad, &
                                 re_val, im_val, re_con, im_con, &
                                 cache_double, startx_by_mpiID, nxpts_by_mpiID
    use OCEAN_psi
    use OCEAN_mpi!, only : nproc, myid, root, comm
    use OCEAN_system
    implicit none
    !
    type( O_system ), intent( in ) :: sys
    type( OCEAN_vector ), intent( in ) :: psi
    type( OCEAN_vector ), intent( inout ) :: psiout
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: re_phi(:), im_phi(:), re_amat(:,:,:), im_amat(:,:,:)
    real(DP) :: spin_prefac, minus_spin_prefac
    real(DP), parameter :: minusone = -1.0_DP
    real(DP), parameter :: one = 1.0_DP
    real(DP), parameter :: zero = 0.0_DP
    integer :: ix, xwidth, ibw, ik, ib, dft_spin, dft_spin2, ispin, psi_spin, psi_spin2
    integer :: psi_con_pad

    ! Populate sizes from psi and val_states
    call OCEAN_psi_returnBandPad( psi_con_pad, ierr )
    if( ierr .ne. 0 ) return
    
    !TODO: spins
    if( sys%valence_ham_spin .ne. 1 ) then
      if( myid .eq. root ) write(6,*) 'fxc not implemeneted for spin systems yet'
      call MPI_BARRIER( comm, ierr )
      ierr = -5
      return
    endif
    psi_spin = 1
    psi_spin2 = 1
    dft_spin = 1
    dft_spin2 = 1
    ispin = 1

    spin_prefac = 2.0_DP
    minus_spin_prefac = -2.0_DP
    ! end TODO

    allocate( re_amat( max(1,nxpts_pad), nbv, nkpts ), im_amat( max(nxpts_pad,1), nbv, nkpts ), &
              re_phi( nxpts_pad ), im_phi( nxpts_pad ) )
    

    if( nxpts .gt. 0 ) then
      xwidth = nxpts
      ix = 1
      re_phi(:) = 0.0_DP
      im_phi(:) = 0.0_DP
! $OMP DO COLLAPSE(1)
      do ik = 1, nkpts
        call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik, dft_spin, 1 ), nxpts_pad, &
                    psi%valr( 1, 1, ik, psi_spin, 1 ), psi_con_pad, zero, re_amat( ix, 1, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', xwidth, nbv, nbc, minusone, im_con(ix, 1, ik,dft_spin, 1 ), nxpts_pad, &
                    psi%vali( 1, 1, ik, psi_spin, 1 ), psi_con_pad, one, re_amat( ix, 1, ik ), nxpts_pad )

        call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, im_con(ix, 1, ik,dft_spin, 1 ), nxpts_pad, &
                    psi%valr( 1, 1, ik, psi_spin, 1), psi_con_pad, zero, im_amat( ix, 1, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik,dft_spin, 1 ), nxpts_pad, &
                    psi%vali( 1, 1, ik, psi_spin, 1 ), psi_con_pad, one, im_amat( ix, 1, ik ), nxpts_pad )
      enddo
! $OMP END DO

!TODO chunk over xpoints for OMP
      do ik = 1, nkpts
        do ib = 1, nbv
          re_phi(:) = re_phi(:) &
                    + re_val(:,ib,ik,dft_spin,1) * re_amat(:,ib,ik) &
                    + im_val(:,ib,ik,dft_spin,1) * im_amat(:,ib,ik)

          im_phi(:) = im_phi(:) &
                    + re_val(:,ib,ik,dft_spin,1) * im_amat(:,ib,ik) &
                    - im_val(:,ib,ik,dft_spin,1) * re_amat(:,ib,ik)
        enddo
      enddo

      if( sys%nbw .eq. 2 ) then
! $OMP DO COLLAPSE(1)
        do ik = 1, nkpts
          call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik, dft_spin, 2 ), nxpts_pad, &
                      psi%valr( 1, 1, ik, psi_spin, 2 ), psi_con_pad, zero, re_amat( ix, 1, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, im_con(ix, 1, ik,dft_spin, 2 ), nxpts_pad, &
                      psi%vali( 1, 1, ik, psi_spin, 2 ), psi_con_pad, one, re_amat( ix, 1, ik ), nxpts_pad )

          call DGEMM( 'N', 'N', xwidth, nbv, nbc, minusone, im_con(ix, 1, ik,dft_spin, 2 ), nxpts_pad, &
                      psi%valr( 1, 1, ik, psi_spin, 2), psi_con_pad, zero, im_amat( ix, 1, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', xwidth, nbv, nbc, one, re_con(ix, 1, ik,dft_spin, 2 ), nxpts_pad, &
                      psi%vali( 1, 1, ik, psi_spin, 2 ), psi_con_pad, one, im_amat( ix, 1, ik ), nxpts_pad )
        enddo
! $OMP END DO

!TODO chunk over xpoints for OMP
      do ik = 1, nkpts
        do ib = 1, nbv
          re_phi(:) = re_phi(:) & 
                    + re_val(:,ib,ik,dft_spin,2) * re_amat(:,ib,ik) &
                    - im_val(:,ib,ik,dft_spin,2) * im_amat(:,ib,ik)

          im_phi(:) = im_phi(:) & 
                    + re_val(:,ib,ik,dft_spin,2) * im_amat(:,ib,ik) &
                    + im_val(:,ib,ik,dft_spin,2) * re_amat(:,ib,ik)
        enddo
      enddo
      endif

      !TODO fix spin
      ! probably save to a new array so we can capture both up-up and up-down potentials from an up exciton
      ! and both down-down and down-up from a down exciton
      re_phi(:) = re_phi(:) * fxc(:,ispin)
      im_phi(:) = im_phi(:) * fxc(:,ispin)


!TODO OMP
      do ik = 1, nkpts
        do ib = 1, nbv
          re_amat(:,ib,ik) = re_phi(:) * re_val(:,ib,ik,dft_spin2,1) &
                           - im_phi(:) * im_val(:,ib,ik,dft_spin2,1)
          im_amat(:,ib,ik) = im_phi(:) * re_val(:,ib,ik,dft_spin2,1) &
                           + re_phi(:) * im_val(:,ib,ik,dft_spin2,1)
        enddo
      enddo

      do ik = 1, nkpts
        call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, dft_spin2, 1 ), nxpts_pad, &
                    re_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%valr( 1, 1, ik, psi_spin2, 1 ), psi_con_pad )
        call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, im_con( 1, 1, ik, dft_spin2, 1 ), nxpts_pad, &
                    im_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%valr( 1, 1, ik, psi_spin2, 1 ), psi_con_pad )

        call DGEMM( 'T', 'N', nbc, nbv, nxpts, minus_spin_prefac, im_con( 1, 1, ik, dft_spin2, 1 ), nxpts_pad, &
                    re_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%vali( 1, 1, ik, psi_spin2, 1 ), psi_con_pad )
        call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, dft_spin2, 1 ), nxpts_pad, &
                    im_amat( 1, 1, ik ), nxpts_pad, &
                    one, psiout%vali( 1, 1, ik, psi_spin2, 1 ), psi_con_pad )
      enddo
! $OMP END DO

      if( sys%nbw .eq. 2 ) then
        do ik = 1, nkpts
          do ib = 1, nbv
            re_amat(:,ib,ik) = re_phi(:) * re_val(:,ib,ik,dft_spin2,2) &
                             + im_phi(:) * im_val(:,ib,ik,dft_spin2,2)
            im_amat(:,ib,ik) = im_phi(:) * re_val(:,ib,ik,dft_spin2,2) &
                             - re_phi(:) * im_val(:,ib,ik,dft_spin2,2)
          enddo
        enddo

        do ik = 1, nkpts
          call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, dft_spin2, 2 ), nxpts_pad, & 
                      re_amat( 1, 1, ik ), nxpts_pad, & 
                      one, psiout%valr( 1, 1, ik, psi_spin2, 2 ), psi_con_pad ) 
          call DGEMM( 'T', 'N', nbc, nbv, nxpts, minus_spin_prefac, im_con( 1, 1, ik, dft_spin2, 2 ), nxpts_pad, & 
                      im_amat( 1, 1, ik ), nxpts_pad, & 
                      one, psiout%valr( 1, 1, ik, psi_spin2, 2 ), psi_con_pad ) 

          call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, im_con( 1, 1, ik, dft_spin2, 2 ), nxpts_pad, & 
                      re_amat( 1, 1, ik ), nxpts_pad, & 
                      one, psiout%vali( 1, 1, ik, psi_spin2, 2 ), psi_con_pad ) 
          call DGEMM( 'T', 'N', nbc, nbv, nxpts, spin_prefac, re_con( 1, 1, ik, dft_spin2, 2 ), nxpts_pad, & 
                      im_amat( 1, 1, ik ), nxpts_pad, & 
                      one, psiout%vali( 1, 1, ik, psi_spin2, 2 ), psi_con_pad ) 
        enddo
      endif
    endif
      
    deallocate( re_phi, im_phi, re_amat, im_amat )

  end subroutine OCEAN_fxc_act

      

    


!  exchange correlation routine, by Ceperley-Alder, as parametrized by
!  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
!  between the unpolarized and polarized gas for the correlation part.
!
      subroutine cacorr(xn,ex,ec,ux1,ux2,uc1,uc2)
        use OCEAN_constants, only : PI
        real(DP), intent( in ) :: xn
        real(DP), intent( out ) :: ex, ec, ux1, ux2, uc1, uc2

        real(DP) :: rel, fp, xn1, xn2, rs, zeta, exchfactor, befactor, beta, eta, xl
        real(DP) :: fe1, fu1, fe2, fu2, ex2, b2, rootr, beta1, beta2, denom, ecu, ucu
        real(DP) :: ecp, ucp, xlr, rlr, au, bu, cu, du, ap, bp, cp, ddp, f, dfdz, ex1, gamm
        real(DP), parameter :: trd=1.0_DP/3.0_DP
        real(DP), parameter :: ft=4.0_DP/3.0_DP

      rel = 0
!
!  get the n's, and the rs.
!
!      pi=3.14159265358979d0
      fp=4.d0*pi
      xn1=xn/2.0_DP
      xn2=xn/2.0_DP

!  effect cutoff, to avoid overflow

      if ( xn .lt. 0.00000001d0 ) then

        ex=0.d0
        ec=0.d0
        ux1=0.d0
        ux2=0.d0
        uc1=0.d0      
        uc2=0.d0

      else

        rs=(3.d0/(fp*xn))**trd
        zeta=(xn1-xn2)/xn
!       exchfactor=-0.930525736d0
        exchfactor=-1.5d0*(0.75d0/pi)**trd

        befactor=(9.d0*pi/4.d0)**trd/137.03599976d0
        if (xn1.eq.0.d0) then
          fe1=1.d0
          fu1=1.d0
          ex1=0.d0
          ux1=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=dsqrt(1.d0+b2)
          xl=dlog(beta+eta)
          fe1=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu1=-0.5d0+1.5d0*xl/beta/eta
          ex1=exchfactor*xn1**trd
          ux1=4.d0*ex1/3.d0
        endif
        if (xn2.eq.0.d0) then
          fe2=1.d0
          fu2=1.d0
          ex2=0.d0
          ux2=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=dsqrt(1.d0+b2)
          xl=dlog(beta+eta)
          fe2=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu2=-0.5d0+1.5d0*xl/beta/eta
          ex2=exchfactor*xn2**trd
          ux2=4.d0*ex2/3.d0
        endif
!  these next lines do the Ceperley-Alder correlation
        if (rs.ge.1.d0) then

          rootr=dsqrt(rs)

          gamm=-0.1423d0
          beta1=1.0529d0
          beta2=0.3334d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecu=gamm/denom
          ucu=ecu*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

          gamm=-0.0843d0
          beta1=1.3981d0
          beta2=0.2611d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecp=gamm/denom
          ucp=ecp*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

        else

          xlr=dlog(rs)
          rlr=rs*xlr

          au= 0.0311d0
          bu=-0.048d0
          cu= 0.002d0
          du=-0.0116d0
          ecu=au*xlr+bu+cu*rlr+du*rs
          ucu=au*xlr+(bu-au/3.d0)+2.d0/3.d0*cu*rlr+(2.d0*du-cu)*rs/3.d0

          ap= 0.01555d0
          bp=-0.0269d0
          cp= 0.0007d0
          ddp=-0.0048d0
          ecp=ap*xlr+bp+cp*rlr+ddp*rs
          ucp=ap*xlr+(bp-ap/3.d0)+2.d0/3.d0*cp*rlr+(2.d0*ddp-cp)*rs/3.d0

        endif

!  if we are nonrel, turn off the MacDonald-Vosko correction.

        if (rel.eq.0.d0) then
          fe1=1.d0
          fu1=1.d0
          fe2=1.d0
          fu2=1.d0
        endif

!  interpolate the correlation energies.

        denom=2.d0**ft-2.d0
        f=((1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0)/denom
        dfdz=ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
        ec=ecu+f*(ecp-ecu)
        uc1=ucu+f*(ucp-ucu)+(ecp-ecu)*(1.d0-zeta)*dfdz
        uc2=ucu+f*(ucp-ucu)-(ecp-ecu)*(1.d0+zeta)*dfdz        
!
!  get the final functional and potential.
!
        ex=(xn1*fe1*ex1+xn2*fe2*ex2)/xn
        ux1=fu1*ux1
        ux2=fu2*ux2
        uc1=uc1
        uc2=uc2
      endif
!
      return
      end subroutine cacorr


end module OCEAN_fxc

