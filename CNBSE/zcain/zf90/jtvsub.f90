! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine jtvsub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, ifcn, stext, ehat, qhat, qmag, &
                   nsphpt, xsph, ysph, zsph, wsph, prefs, nsemi, semifcn, semiReducedEnergy, c2c )
  implicit none
  !
  integer :: nsphpt
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  integer :: lmin, lmax, npmax, lc, powmax, nsemi
  real( kind = kind( 1.0d0 ) ) :: ehat( 3 ), qhat( 3 ), qmag
  character(len=10) :: spcttype, stext
  integer :: nproj( lmin : lmax )
  real( kind = kind(1.0d0)) :: ifcn( npmax, 0 : powmax, lmin : lmax), &
                               semifcn( npmax, 0 : powmax, lmin : lmax, nsemi ), &
                               semiReducedEnergy( nsemi ), c2c( 0 : powmax, nsemi )
  complex( kind = kind( 1.0d0 ) ) :: nbsemel( npmax, -lmax : lmax, lmin : lmax, -lc : lc )
  !
  integer :: i, j, l, m, mc, lpick, lt, lldos, mldos, msign, isemi, li, mi
  real( kind = kind( 1.0d0 ) ) :: edot, qdot, jlmel( npmax, 0 : lmax + lc, lmin : lmax ), pl
  real( kind = kind( 1.0d0 ) ) :: ldosmel( npmax  )
  complex( kind = kind( 1.0d0 )) :: csu( powmax ), ylm, ylcmc, iqmag, rp, rm1, angint, mip, cpr
  logical :: mdos
  !
  logical, parameter :: PCCP2022 = .false. ! The behavior in Vinson, PCCP 2022 is incorrect
  ! This flag reproduces it, but is left as a compile time setting because it is wrong
  ! There are additional flags to test only the Quad contribution, or only the 1st core level
  ! (e.g. 2p when looking at 1s two-photon), or only the 2nd core level (e.g. 3p)
  logical, parameter :: QuadAlone = .false.  ! This is now superceded by use of 'tpq' photon
!  logical, parameter :: test1 = .false.
!  logical, parameter :: test2 = .false.
  logical :: test1, test2
  !
  !
  lpick = -1
  spcttype = stext
  select case( spcttype(1:2) )
    case( 's-' )
      lpick = 0
    case( 'p-' )
      lpick = 1
    case( 'd-' )
      lpick = 2
  end select
  if ( lpick .ge. 0 ) then
     do l = lmin, lmax
        if ( l .ne. lpick ) then
           ifcn( :, :, l ) = 0.0d0
           write ( 6, * ) 'LPICK says zero ... ', l
        else
           write ( 6, * ) 'LPICK says do not zero ... ', l
        end if
     end do
  end if
  !
  if ( spcttype(1:2) .eq. 'qR' ) call jlmatfetch( lc, lmin, lmax, npmax, nproj, qmag, jlmel, powmax )
  if ( spcttype(1:4) .eq. 'ldos' ) then
    spcttype = 'ldos'
    call ldosfetch( stext, lmin, lmax, npmax, nproj, qmag, ldosmel )
    select case (stext(5:5) )
      case( '0' )
        lldos = 0
      case( '1' )
        lldos = 1
      case( '2' )
        lldos = 2
      case( '3' )
        lldos = 3
    end select
    select case (stext(6:6) )
      case ('-')
        i = 7
        msign = -1
      case default 
        i = 6
        msign = 1
    end select
    select case (stext(i:i) )
      case('0')
        mldos = 0
        mdos = .true.
      case('1')
        mldos = 1*msign
        mdos = .true.
      case('2')
        mldos = 2*msign
        mdos = .true.
      case('3')
        mldos = 3*msign
        mdos = .true.
      case default
        mdos = .false.
    end select
  endif
  !
  ! We will need powers of ( -i q )
  rm1 = -1
  rm1 = sqrt( rm1 )
  iqmag = -1
  iqmag = sqrt( iqmag )
  iqmag = - iqmag * qmag
  write ( 6, * ) ' qmag = ', qmag, 'inverse bohrs'
  !      
  do mc = -lc, lc
     do l = lmin, lmax
        do m = -l, l
           select case( spcttype )
              !
           case( 'dipole', 's-dipole', 'p-dipole', 'd-dipole' )
              csu(:) = 0
              do i = 1, nsphpt
                 call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 csu(1) = csu(1) + conjg(ylcmc) * edot * ylm * wsph(i)
              end do
              nbsemel( 1 : nproj( l ), m, l, mc ) = csu(1) * ifcn( 1 : nproj( l ), 1, l )
              !
           case( 'tp', 'tpq', 'tp1', 'tp2', 'tpp' )
              csu(:) = 0
              do i = 1, nsphpt
                 call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 csu(1) = csu(1) + conjg(ylcmc) * edot * edot * ylm * wsph(i)
              end do
              nbsemel( 1 : nproj( l ), m, l, mc ) = csu(1) * ifcn( 1 : nproj( l ), 2, l )
              !
           case( 'quad', 's-quad', 'p-quad', 'd-quad' )
              csu(:) = 0
              do i = 1, nsphpt
                 call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ),  ylcmc, prefs )
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) + zsph( i ) * qhat( 3 )
                 csu(1) = csu(1) + conjg(ylcmc) * edot * ylm * wsph(i)
                 csu(2) = csu(2) + conjg(ylcmc) * edot * ylm * wsph(i) * qdot * iqmag
              enddo
              nbsemel( 1 : nproj( l ), m, l, mc ) = csu( 1 ) * ifcn( 1 : nproj( l ), 1, l ) + &
                   csu( 2 ) * ifcn( 1 : nproj( l ), 2, l ) * 0.5d0 ! Joly factor
              !
           case( 'quadalone', 's-quadalone', 'p-quadalone', 'd-quadalone' )
              csu(:) = 0
              do i = 1, nsphpt
                 call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) + zsph( i ) * qhat( 3 )
                 csu(2) = csu(2) + conjg(ylcmc) * edot * ylm * wsph(i) * qdot * iqmag
              enddo
              nbsemel( 1 : nproj( l ), m, l, mc ) = csu(2) * ifcn( 1 : nproj( l ), 2, l ) * 0.5d0 ! Joly factor
              !
           case( 'NRIXS' )
              csu(:) = 0
              do i = 1, nsphpt
                 call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) + zsph( i ) * qhat( 3 )
                 rp = qdot  
                 do j = 1, powmax
                    rp = rp / dble( j )
                    csu( j ) = csu( j ) + conjg( ylcmc ) * ylm * wsph( i ) * rp
                    rp = rp * ( -rm1 * qmag * qdot )
                 enddo
              enddo
              nbsemel( :, m, l, mc ) = 0
              do j = 1, powmax
                 nbsemel( :, m, l, mc ) = nbsemel( :, m, l, mc ) + csu(j) * ifcn( :, j, l ) 
              enddo
              !
           case( 'qRaman', 'qRs', 'qRp', 'qRd', 'qRf' )
              nbsemel( :, m, l, mc ) = 0

              if( spcttype .ne. 'qRaman' ) then
                if( l .ne. 0 .and. spcttype .eq. 'qRs' ) cycle
                if( l .ne. 1 .and. spcttype .eq. 'qRp' ) cycle
                if( l .ne. 2 .and. spcttype .eq. 'qRd' ) cycle
                if( l .ne. 3 .and. spcttype .eq. 'qRf' ) cycle
              endif
              
              mip = 1.0d0
              do lt = 0, lc + l
                 cpr = mip * dble( 2 * lt + 1 )
                 angint = 0.0d0
                 do i = 1, nsphpt
                    call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                    call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                    qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) + zsph( i ) * qhat( 3 )
                    select case( lt )
                    case( 0 )
                      pl = 1.0d0
                    case( 1 )
                      pl = qdot
                    case( 2 )
                      pl = 0.5d0 * ( 3.0d0 * qdot ** 2 - 1.0d0 )
                    case( 3 )
                      pl = 0.5d0 * ( 5.0d0 * qdot ** 3 - 3.0d0 * qdot )
                    case( 4 )
                      pl = 0.125d0 * ( 35.0d0 * qdot ** 4 - 30.0d0 * qdot ** 2 + 3.0d0 )
                    case( 5 )
                      pl = 0.125d0 * ( 63.0d0 * qdot ** 5 - 70.0d0 * qdot ** 3 + 15.0d0 * qdot )
                    case default
                      pl = 0.0d0
                    end select
                    angint = angint + wsph( i ) * conjg( ylcmc ) * ylm * pl
                 enddo
                 nbsemel( :, m, l, mc ) = nbsemel( :, m, l, mc ) + angint * jlmel( :, lt, l ) * cpr
                 mip = -rm1 * mip
              enddo
              nbsemel( :, m, l, mc ) = nbsemel( :, m, l, mc ) * rm1 / qmag
              !
           case ( 'ldos' )
              if( l .eq. lldos ) then
                if( mdos ) then
                  if( m .eq. 0 .and. m .eq. mldos ) then
                    nbsemel( 1 : nproj( l ), m, l, mc ) = ldosmel( 1 : nproj(l) )
                  elseif( m .eq. -mldos ) then
                    nbsemel( 1 : nproj( l ), m, l, mc ) = 1.0d0/sqrt(2.0d0) * ldosmel( 1 : nproj(l) )
                  elseif( m .eq. mldos ) then
                    if( msign .gt. 0 ) then
                      nbsemel( 1 : nproj( l ), m, l, mc ) = (-1)**l * 1.0d0/sqrt(2.0d0) * ldosmel( 1 : nproj(l) )
                    else
                      nbsemel( 1 : nproj( l ), m, l, mc ) = (-1)**l * (-1.0d0/sqrt(2.0d0) * ldosmel( 1 : nproj(l) ))
                    endif
                  endif
                else
                  nbsemel( 1 : nproj( l ), m, l, mc ) = ldosmel( 1 : nproj(l) ) !/ sqrt( dble(2*l+1) )
                endif
              endif
           case default
              write(6,*) "Photon type was not recognized. Must be one of the following:"
              write(6,*) "      dipole, quad, quadalone, NRIXS, qRaman"
              stop
           end select
        end do
     end do
  end do
  !
  if( spcttype .eq. 'tp' .or. spcttype .eq. 'tp1' .or. spcttype .eq. 'tp2' .or. spcttype .eq. 'tpp') then
    ! Current implementation is hardwired for a K edge
    li = 1

    if( spcttype .eq. 'tp1' ) then
      test1 = .true.
    else
      test1 = .false.
    endif
    if( spcttype .eq. 'tp2' ) then
      test2 = .true.
    else
      test2 = .false.
    endif
    if( spcttype .eq. 'tpp' ) then
      test1 = .true.
      test2 = .true.
    endif

    if( test1 .or. test2 ) then
      nbsemel( :, :, :, : ) = 0.0d0
    endif
    do isemi = 1, nsemi
!      write( filename, '(A8,I3.3,A1,I2.2,A1,I2.2)' ) 'melfilez', z, 'n', i, 'l', li
!      inquire(file=filename, exist=ex )
!      if( .not. ex )
!        write( filename, '(A8,I3.3,A1,I2.2,A1,I2.2)' ) 'melsemiz', z, 'n', i, 'l', li
!        inquire(file=filename, exist=ex )
!      endif
!      if( .not. ex ) cycle
      do mc = -lc, lc
        do mi = -li, li
          csu(1) = 0.0d0
          do i = 1, nsphpt
             call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
             call getylm( li, mi, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
             edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
             csu(1) = csu(1) + conjg(ylcmc) * edot * ylm * wsph(i)
          end do
          do l = lmin, lmax
            do m = -l, l
              csu(2) = 0.0d0
              do i = 1, nsphpt
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                 call getylm( li, mi, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 csu(2) = csu(2) + ylcmc * edot * conjg(ylm) * wsph(i)
              end do
              ! As noted above, the following is incorrect
              if( QuadAlone ) then
              else
              if( test1 ) then
                if( isemi .eq. 1 ) then
                  nbsemel( 1 : nproj( l ), m, l, mc ) = nbsemel( 1 : nproj( l ), m, l, mc ) &
                    + csu(1) * csu(2) * c2c( 1, isemi ) * semifcn( 1 : nproj(l), 1, l, isemi ) &
                    * semiReducedEnergy( isemi )
                endif
              endif
              if( test2 ) then
                if( isemi .eq. 2 ) then
                  nbsemel( 1 : nproj( l ), m, l, mc ) = nbsemel( 1 : nproj( l ), m, l, mc ) &
                    + csu(1) * csu(2) * c2c( 1, isemi ) * semifcn( 1 : nproj(l), 1, l, isemi ) &
                    * semiReducedEnergy( isemi )
                endif
              endif
              if( PCCP2022 ) then
                nbsemel( 1 : nproj( l ), m, l, mc ) = nbsemel( 1 : nproj( l ), m, l, mc ) &
                    - csu(1) * csu(2) * c2c( 1, isemi ) * semifcn( 1 : nproj(l), 1, l, isemi )
              endif
              if( .not. test1 .and. .not. test2 .and. .not. PCCP2022 .and. .not. QuadAlone) then ! The correct behavior
                ! transitions can happen first from any occupied core or semi-core state and 
                ! then from the final target level into the newly unoccupied core. E.g., in 
                ! addition to 1s -> 3d/4s with a quadrupole operator, so can first have 
                ! 2p->3d/4s followed by 1s->2p, both with dipole. 
                ! The array semiReducedEnergy contains the fraction which is a scaling based
                ! on energy of the photon and energy of the core level
                nbsemel( 1 : nproj( l ), m, l, mc ) = nbsemel( 1 : nproj( l ), m, l, mc ) &
                    + ( csu(1) * csu(2) * c2c( 1, isemi ) * semifcn( 1 : nproj(l), 1, l, isemi ) &
                      * semiReducedEnergy( isemi ) )
              endif
              endif
            enddo
          enddo
        enddo
      end do
    end do
  endif
  return
end subroutine jtvsub
