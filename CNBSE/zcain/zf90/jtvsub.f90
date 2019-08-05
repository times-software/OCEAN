! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine jtvsub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, ifcn, stext, ehat, qhat, qmag, nsphpt, xsph, ysph, zsph, wsph, &
     prefs )
  implicit none
  !
  integer :: nsphpt
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 ) 
  integer :: lmin, lmax, npmax, lc, powmax
  real( kind = kind( 1.0d0 ) ) :: ehat( 3 ), qhat( 3 ), qmag
  character(len=10) :: spcttype, stext
  integer :: nproj( lmin : lmax )
  real( kind = kind(1.0d0)) :: ifcn( npmax, 0 : powmax, lmin : lmax)
  complex( kind = kind( 1.0d0 ) ) :: nbsemel( npmax, -lmax : lmax, lmin : lmax, -lc : lc )
  !
  integer :: i, j, l, m, mc, lpick, lt
  real( kind = kind( 1.0d0 ) ) :: edot, qdot, jlmel( npmax, 0 : lmax + lc, lmin : lmax ), pl
  complex( kind = kind( 1.0d0 )) :: csu( powmax ), ylm, ylcmc, iqmag, rp, rm1, angint, mip, cpr
  !
  lpick = -1
  spcttype = stext
  if ( stext .eq. 'ldos0' ) then
     spcttype = 'NRIXS'
     lpick = 0
  end if
  if ( stext .eq. 'ldos1' ) then
     spcttype = 'NRIXS'
     lpick = 1
  end if
  if ( stext .eq. 'ldos2' ) then
     spcttype = 'NRIXS'
     lpick = 2
  end if
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
  write ( 6, * ) 'ck1'
  if ( spcttype .eq. 'qRaman' ) call jlmatfetch( lc, lmin, lmax, npmax, nproj, qmag, jlmel, powmax )
  write ( 6, * ) 'ck2'
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
           case( 'dipole' )
              csu(:) = 0
              do i = 1, nsphpt
                 call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
                 call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 csu(1) = csu(1) + conjg(ylcmc) * edot * ylm * wsph(i)
              end do
              nbsemel( 1 : nproj( l ), m, l, mc ) = csu(1) * ifcn( 1 : nproj( l ), 1, l )
              !
           case( 'quad' )
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
           case( 'quadalone' )
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
           case( 'qRaman' )
              nbsemel( :, m, l, mc ) = 0
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
           case default
              write(6,*) "Photon type was not recognized. Must be one of the following:"
              write(6,*) "      dipole, quad, quadalone, NRIXS, qRaman"
              stop
           end select
        end do
     end do
  end do
  !
  return
end subroutine jtvsub
