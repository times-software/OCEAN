      subroutine mkkbfile( nel, nl, xnj, ev, dl, nr, r, dr, rsqd, cq,    &
     &                     vi, orb, zorig, rel, njrc, psflag )
      implicit none
!
      integer nel, nr
      double precision dl, zorig
!
      integer nl( nel ), njrc( 4 )
      double precision xnj( nel ), ev( nel ), r( nr ), cq( nr )
      double precision dr( nr ), rsqd( nr )
      double precision vi( nr, 7 ), orb( nr, nel )
!
      integer, allocatable :: twoj( : )
      double precision, allocatable :: vps( : ), xm1( : ), xm2( : )
      double precision, allocatable :: phi( : , : ), psi( : )
!
      integer i, j, k, n, iu, flag, neff, jrt, ief, nn, ie
      double precision e, de, el, eh, rstop, kappa, zeff, rel, plead, x
      double precision ruse, su
      double precision d1, d2, ps0, ps1, norm
!
      logical psflag
!
      allocate( twoj( nel ), vps( nr ), xm1( nr ), xm2( nr ) )
      allocate( phi( nr, -1 : 1 ), psi( nr ) )
!
!
      read (5,*) de,rstop,el,eh,n
      jrt=1.0000001d0+log(rstop/r(1))/dl
      do i=1,nel
        twoj(i)=0.0001d0+2.d0*abs(xnj(i))
      end do
      if ( njrc( 1 ) .ne. 0 ) then
        open(unit=13,file='sepipt',form='formatted',status='unknown')
        rewind 13
        write (13,'(1x,2i5,1f20.10)') nel
        write (13,'(1x,2i5,1f20.10)') (nl(i),twoj(i),ev(i),i=1,nel)
        write (13,'(1x,2i5,1f20.10)') jrt,nr,de
        write (13,'(1x,2(1x,1e22.15))') r(1),r(2)
        do i=1,nr
          write (13,'(1x,8(1x,1e19.12))') (vi(i,k),k=1,7),cq(i)
          write (13,'(1x,4(1x,1e19.12))') (orb(i,j),j=1,nel)
        end do
      end if
      do k=1,nel
        kappa=-1.d0
        if ( njrc( 1 ) .ne. 0 ) then
          iu=40+k
          zeff=0.d0
          ruse=0.d0
        else
          iu = 50 + k
          if (abs(xnj(k)).gt.dble(nl(k))+0.25d0) kappa=-nl(k)-1
          if (abs(xnj(k)).lt.dble(nl(k))-0.25d0) kappa= nl(k)
          zeff = zorig
          ruse = rel 
        end if        
        call gropen(iu)
        flag=1
        neff=1000
        call setqmm(k,orb(1,k),nl(k),xnj(k),flag,vps,zeff,               &
     &              zorig,ruse,nr,r,rsqd,dl,xm1,xm2,njrc,vi,psflag)
        call getplead(nl(k),xnj(k),ruse,kappa,plead,zeff)
        if ( njrc( 1 ) .ne. 0 ) then
          phi( :, : ) = 0.0d0
          do i=-1,1
             e = ev( k ) + de * dble( i )
!            e = ev( k ) * 0.5d0 * dble( 1 - i ) ... this was eager ...
            call integ(e,nl(k),kappa,neff,nn,jrt,ief,x,phi(1,i),
     &                 zeff,vps,xm1,xm2,nr,r,dr,rsqd,dl,ruse,plead)
             su = 0.d0
             do j = 1, jrt - 1
                su = su + r( j ) * phi( j, i ) ** 2 +
     &                    r( j + 1 ) * phi( j + 1, i ) ** 2
             end do
             su = 1.0d0 / sqrt( 0.5d0 * dl * su )
             phi( :, i ) = phi( :, i ) * su
          end do
          do j=1,nr
             write (13,'(1x,1i5,3(2x,1e19.12))') j,(phi(j,i),i=-1,1)
             write ( 23, '(1x,2i5,4(1x,1e22.15))' ) j, nr, 
     &          r( j ), phi( j, 0 ), phi( j, -1 ), phi( j, 1 )
          end do
        end if
        open( unit=99, file='jrtval', 
     &       form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(2i12)' ) jrt, int( dble( n ) / 100.0d0 )
        close( unit=99 )
        e = el
        ie = 0
        do while ( e .lt. eh )
          call integ(e,nl(k),kappa,neff,nn,jrt,ief,x,psi,                &
     &               zeff,vps,xm1,xm2,nr,r,dr,rsqd,dl,ruse,plead)
          ps0 = psi( jrt )
          d1 = ( psi( jrt + 1 ) - psi( jrt - 1 ) ) / ( 2.0d0 * dl )
          d2 = ( psi( jrt + 2 ) - psi( jrt - 2 ) ) / ( 4.0d0 * dl )
          ps1 = ( 4.0d0 * d1 - d2 ) / ( 3.0d0 * r( jrt ) )
          ps1 = ps1 / sqrt( 2.0d0 * e )
          write (iu,'(1x,5(1x,1e15.8))' ) x, e, psi( 1 ), ps0, ps1
          ie = ie + 1
          if ( ie .eq. 100 ) then
             norm = sqrt( ps0 ** 2 + ps1 ** 2 )
             psi( : ) = psi( : ) / norm
             do i = 1, jrt
                write ( 30 + k, '(3(1x,1e15.8))' ) e, r( i ), psi( i )
             end do
             write ( 30 + k, * )
             ie = 0
          end if
          e = e + ( eh - el ) / dble( n )
        end do
      end do
      deallocate( twoj, vps, xm1, xm2, phi, psi )
      return
      end
