program rixsPlane
  use AI_kinds, only : DP
  use OCEAN_constants, only : eV2Hartree
  implicit none

  real(DP) :: dumf, estart, estop, estep, rixs_start, rixs_stop, rixs_broaden, energy, &
              kpref, intensity, eloss, fact, eemission, vb
  complex(DP) :: ctmp, arg, rp, rm ,rrr, al, be, eps
  integer :: dumi, rixs_nstep, npol, iXAS, ie, ipol, iter, i, nstep, nVarBroaden
  character(len=4) :: flg
  character(len=2) :: el, nl
  character( LEN=24 ) :: filnam
  integer, allocatable :: polList(:,:)
  real(DP), allocatable :: xas_energies(:), a(:), b(:), emission_slice(:), varBroaden(:,:)
  logical :: ex, do_var_broaden


  inquire(file='bse.in.xas',exist=ex)
  if( .not. ex ) then
    write(6,*) 'Could not find bse.in.xas'
    goto 111
  endif

  open(unit=99,file='rixs_spect.in',form='formatted',status='old')
  read(99,*) rixs_nstep, rixs_start, rixs_stop, rixs_broaden
  rixs_broaden = rixs_broaden * eV2Hartree
  read(99,*) el, nl
  read(99,*) npol
  allocate( polList( 2, npol ) )
  do i = 1, npol
    read(99,*) polList(:,i)
  enddo
  close(99)

  open(file='bse.in.xas',unit=99, form='formatted', status='old')
  read(99,*) dumi
  read(99,*) dumf
  read(99,*) flg
  if( flg .ne. 'inv' ) then
    write(6,*) 'bse.in.xas was not a GMRES run'
    close(99)
    goto 111
  endif
  read(99,*) dumi
  read(99,*) flg
  select case( flg )
    case ("loop" )
      read(99,*) estart, estop, estep
      nstep = nint( (estop-estart)/estep)
      allocate(xas_energies(nstep))
      energy = estart
      do i=1, nstep
        xas_energies(i) = energy
        energy = energy + estep
      enddo
    case ("list")
      read(99,*) nstep
      allocate(xas_energies(nstep))
      do i = 1, nstep
        read(99,*) xas_energies(i)
      enddo
    case default
      write(6,*) 'Unsupported calculation'
      close(99)
      goto 111
  end select
  close(99)

  inquire(file='rixs_plane_broaden.ipt',exist=do_var_broaden)
  if(do_var_broaden) then
    open(unit=99,file='rixs_plane_broaden.ipt',form='formatted', status='old')
    read(99,*) do_var_broaden 
    close(99)
  endif
  
  if( do_var_broaden ) then
    open(unit=99,file='rixs_var_broaden',form='formatted', status='old')
    read(99,*) nVarBroaden
    allocate(varBroaden(nVarBroaden,2) )
    do i = 1, nVarBroaden
      read(99,*) varBroaden(i,1), varBroaden(i,2)
    enddo
    close(99)
  endif

  if( flg .eq. "loop" ) then
    open(unit=98,file='rixs.txt',status='unknown')
  endif

  allocate( emission_slice(rixs_nstep))

  do iXAS = 1, nstep
    if( flg .eq. "list" ) then
      write(filnam, "(A,I0,A)" ) "rixs.", iXAS, ".txt"
      open(unit=98,file=filnam,status='unknown' )
      rewind(98)
    endif
    emission_slice(:) = 0.0_DP
    do ipol = 1, npol

      write(filnam, '(A7,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxlanc_', el, '.', nl, '_', &
                polList(1,ipol), '.', iXAS, '.', polList(2,ipol)
      open(file=filnam,unit=99,form='formatted',status='old')
      read(99,*) iter, kpref
      fact = kpref ! ignore cell volume/spin normalization
      allocate( a(0:iter), b(iter))
      read(99,*) a(0)
      do i = 1, iter
        read(99,*) a(i), b(i)
      enddo
      close(99)

      do ie = 1, rixs_nstep
        intensity = 0.0_DP
        eemission = rixs_start + (ie-1)*(rixs_stop-rixs_start)/dble(rixs_nstep)  
        eloss = xas_energies(iXAS) - eemission
!        eloss = xas_energies(iXAS) - rixs_start - (ie-1)*(rixs_stop-rixs_start)/dble(rixs_nstep)
        eloss = eloss * eV2Hartree
        if (eloss .gt. 0 ) then
          vb = 0.0_DP
          if( do_var_broaden ) then
            call intval( nVarBroaden, varBroaden(:,1), varBroaden(:,2), eemission, vb, 'cap', 'cap' )
            vb = vb * eV2Hartree
          endif
          ctmp = cmplx(eloss, max(rixs_broaden,vb), DP )
          arg = ( eloss - a( iter - 1 ) ) ** 2 - 4.0_dp * b( iter ) ** 2
          arg = sqrt( arg )

          rp = 0.5_dp * ( eloss - a( iter - 1 ) + arg )
          rm = 0.5_dp * ( eloss - a( iter - 1 ) - arg )
          if( aimag( rp ) .lt. 0.0_dp ) then
            rrr = rp
          else
            rrr = rm
          endif

          al =  ctmp - a( iter - 1 ) - rrr
          be = -ctmp - a( iter - 1 ) - rrr
          do i = iter-1, 0, -1
            al = ctmp - a( i ) - b( i+1 )**2 / al
            be = -ctmp - a( i ) - b( i+ 1 ) **2 / be
          enddo

          eps = 1.0_dp - fact / al - fact / be
          intensity = aimag( eps )
        endif
        emission_slice(ie) = emission_slice(ie) + intensity
      enddo

      deallocate( a, b )
    enddo

    if( flg .eq. "list" ) then
      do ie = 1, rixs_nstep
        eemission = rixs_start + (ie-1)*(rixs_stop-rixs_start)/dble(rixs_nstep)
        write(98,'(2E24.15)') eemission, emission_slice(ie)
      enddo
      close(98)
    elseif( flg .eq. "loop" ) then
      do ie = 1, rixs_nstep
        eemission = rixs_start + (ie-1)*(rixs_stop-rixs_start)/dble(rixs_nstep)
  !      write(98,'(4E24.15)') rixs_start + (ie-1)*(rixs_stop-rixs_start)/dble(rixs_nstep), &
        write(98,'(4E24.15)') eemission,  xas_energies(iXAS), emission_slice(ie), &
                              xas_energies(iXAS) - eemission
  !                            xas_energies(iXAS) - rixs_start - (ie-1)*(rixs_stop-rixs_start)/dble(rixs_nstep)
      enddo
      write(98,*)
    endif
  enddo

  close(98)
  deallocate(emission_slice, polList, xas_energies )
          

  111 continue


end program rixsPlane
