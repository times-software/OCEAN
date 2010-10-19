      subroutine jtvsub( lmin, lmax, nproj, npmax, lc, nbsemel, powmax, &
     & ifcn, spcttype, ehat, qhat, qmag, nsphpt, xsph, ysph, zsph,      &
     &  wsph, prefs )
      implicit none
!
      integer :: nsphpt
      real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph,  &
     &     zsph, wsph
      real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 ) 
      integer :: lmin, lmax, npmax, lc, powmax
      real( kind = kind( 1.0d0 ) ) :: ehat( 3 ), qhat( 3 ), qmag
      character * 10 :: spcttype
      integer :: nproj( lmin : lmax )
      real( kind = kind(1.0d0)) :: ifcn( npmax, 0 : powmax, lmin : lmax)
      complex( kind = kind( 1.0d0 ) ) :: nbsemel( npmax, -lmax : lmax,  &
     &      lmin : lmax, -lc : lc )
!
      integer :: i, j, l, m, mc, factorial
      real( kind = kind( 1.0d0 ) ) :: edot, qdot
      complex( kind = kind( 1.0d0 )) :: csu( powmax ), ylm, ylcmc, iqmag
!
! We will need powers of ( -i q )
      iqmag = -1
      iqmag = sqrt( iqmag )
      iqmag = - iqmag * qmag
!      
      do mc = -lc, lc
        do l = lmin, lmax
          do m = -l, l
            select case( spcttype )
!
              case( 'dipole' )
                csu(:) = 0
                do i = 1, nsphpt
                  call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), &
     &                         ylcmc, prefs )
                  call getylm( l, m, xsph( i ), ysph( i ), zsph( i ),   &
     &                        ylm, prefs )
                  edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) +&
     &                   zsph( i ) * ehat( 3 )
                  csu(1) = csu(1) + conjg(ylcmc) * edot * ylm * wsph(i)
                end do
                nbsemel( 1 : nproj( l ), m, l, mc ) =                   &
     &            csu(1) * ifcn( 1 : nproj( l ), 1, l )
!
              case( 'quad' )
                csu(:) = 0
                do i = 1, nsphpt
                  call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), &
     &                         ylcmc, prefs )
                  call getylm( l, m, xsph( i ), ysph( i ), zsph( i ),   &
     &                        ylm, prefs )
                  edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) +&
     &                   zsph( i ) * ehat( 3 )
                  qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) +&
     &                   zsph( i ) * qhat( 3 )
                  csu(1) = csu(1) + conjg(ylcmc) * edot * ylm * wsph(i)
                  csu(2) = csu(2) + conjg(ylcmc) * edot * ylm * wsph(i) &
     &                            * qdot * iqmag * 0.5d0 
! factor of 1/2 above added apr19 2010 jtv
                enddo
                nbsemel( 1 : nproj( l ), m, l, mc ) =                   &
     &            csu( 1 ) * ifcn( 1 : nproj( l ), 1, l ) +             &
     &            csu( 2 ) * ifcn( 1 : nproj( l ), 2, l )
!
              case( 'quadalone' )
                csu(:) = 0
                do i = 1, nsphpt
                  call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), &
     &                         ylcmc, prefs )
                  call getylm( l, m, xsph( i ), ysph( i ), zsph( i ),   &
     &                        ylm, prefs )
                  edot = xsph( i ) * ehat( 1 ) + ysph( i ) * ehat( 2 ) +&
     &                   zsph( i ) * ehat( 3 )
                  qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) +&
     &                   zsph( i ) * qhat( 3 )
                  csu(2) = csu(2) + conjg(ylcmc) * edot * ylm * wsph(i) &
     &                            * qdot * iqmag * 0.5d0 
                enddo
                nbsemel( 1 : nproj( l ), m, l, mc ) =                   &
     &            csu(2) * ifcn( 1 : nproj( l ), 2, l )
!
              case( 'NRIXS' )
                csu(:) = 0
                do i = 1, nsphpt
                  call getylm( lc, mc, xsph( i ), ysph( i ), zsph( i ), &
     &                         ylcmc, prefs )
                  call getylm( l, m, xsph( i ), ysph( i ), zsph( i ),   &
     &                        ylm, prefs )
                  qdot = xsph( i ) * qhat( 1 ) + ysph( i ) * qhat( 2 ) +&
     &                   zsph( i ) * qhat( 3 )
                  factorial = 1
                  do j = 1, powmax
                    factorial = factorial * j
                    csu(j) = csu(j) + conjg(ylcmc) * ylm * wsph(i)      &
     &                     * qdot**(j) * iqmag**(j-1) / dble(factorial) 
                  enddo
                enddo
                nbsemel( 1 : nproj( l ), m, l, mc ) = 0
                do j = 1, powmax
                  nbsemel( 1 : nproj( l ), m, l, mc ) =                 &
     &              nbsemel( 1 : nproj( l ), m, l, mc ) +               &
     &              csu(j) * ifcn( 1 : nproj( l ), j, l ) 
                enddo
!
            end select
          end do
        end do
      end do 
!
      return
      end subroutine jtvsub
