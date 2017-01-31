! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine elssub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, ifcn, spcttype, ehat, qhat, q, &
     nsphpt, xsph, ysph, zsph, wsph, prefs )
  implicit none
  !
  integer :: nsphpt
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 ) 
  integer :: lmin, lmax, npmax, lc, powmax
  real( kind = kind( 1.0d0 ) ) :: qhat( 3 ), q
  complex( kind = kind( 1.0d0 ) ) :: ehat( 3 ), edot
  character * 15 :: spcttype
  integer :: nproj( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: ifcn( npmax, 0 : powmax, lmin : lmax )
  complex( kind = kind( 1.0d0 ) ) :: nbsemel( npmax, -lmax : lmax, lmin : lmax, -lc : lc )
  !
  integer :: i, j, l, m, mc
  real( kind = kind( 1.0d0 ) ) :: qdot
  complex( kind = kind( 1.0d0 ) ) :: ylm, ylcmc, rm1, emel, eqmel, qpmel( 0 : powmax ), jp
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  ! calculate full matrix elements
  write ( 6, * ) ehat
  write ( 6, * ) qhat
  nbsemel( :, :, :, : ) = 0.0d0
  do mc = -lc, lc
     do l = lmin, lmax
        do m = -l, l
           select case( spcttype )
           case( 'transverse' )
              emel = 0
              eqmel = 0
           case( 'longitudinal' )
              qpmel( : ) = 0.0d0
           end select
           do i = 1, nsphpt
              call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), ylcmc, prefs )
              call getylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
              select case( spcttype )
              case( 'transverse' )
                 edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) + zsph( i ) * ehat( 3 )
                 qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) + zsph( i ) * qhat( 3 )
                 emel = emel + conjg( ylcmc ) * edot * ylm * wsph( i )
                 eqmel = eqmel - conjg( ylcmc ) * edot * rm1 * q * qdot * ylm * wsph( i )
              case( 'longitudinal' )
                 qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) + zsph( i ) * qhat( 3 )
                 jp = qdot ! we omit the -rm1 to have an analogy to dipole results
                 do j = 1, powmax
                    qpmel( j ) = qpmel( j ) + conjg( ylcmc ) * jp * ylm * wsph( i ) 
                    jp = -rm1 * q * qdot * jp / dble( j + 1 )
                 end do
              end select
           end do
           select case( spcttype )
           case( 'transverse' )
              nbsemel( :, m, l, mc ) = emel * ifcn( :, 1, l )
              if ( powmax .ge. 2 ) then
                 nbsemel( :, m, l, mc ) = nbsemel( :, m, l, mc ) + eqmel * ifcn( :, 2, l )
              end if 
           case( 'longitudinal' )
              nbsemel( :, m, l, mc ) = 0
              do j = 1, powmax
                 nbsemel( :, m, l, mc ) = nbsemel( :, m, l, mc ) + qpmel( j ) * ifcn( :, j, l )
              end do
           end select 
        end do
     end do
  end do 
  !
  return
end subroutine elssub
