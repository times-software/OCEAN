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
      integer i, j, k, n, iu, flag, neff, jrt, ief, nn
      double precision e, de, el, eh, rstop, kappa, zeff, rel, plead, x
      double precision ruse, sum
!
      logical psflag
!
      allocate( twoj( nel ), vps( nr ), xm1( nr ), xm2( nr ) )
      allocate( phi( nr, -1 : 1 ), psi( nr ) )
!
!
      read (5,*) de,rstop,el,eh,n
      jrt=1.0000001d0+dlog(rstop/r(1))/dl
      do i=1,nel
        twoj(i)=0.0001d0+2.d0*dabs(xnj(i))
      end do
      if ( njrc( 1 ) .ne. 0 ) then
        open(unit=13,file='sepipt',form='formatted',status='unknown')
        rewind 13
        write (13,'(1x,2i5,1f20.10)') nel
        write (13,'(1x,2i5,1f20.10)') (nl(i),twoj(i),ev(i),i=1,nel)
        write (13,'(1x,2i5,1f20.10)') jrt,nr,de
        write (13,'(1x,2(1x,1d15.8))') r(1),r(2)
        do i=1,nr
          write (13,'(1x,8(1x,1d19.12))') (vi(i,k),k=1,7),cq(i)
          write (13,'(1x,4(1x,1d19.12))') (orb(i,j),j=1,nel)
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
          if (dabs(xnj(k)).gt.dble(nl(k))+0.25d0) kappa=-nl(k)-1
          if (dabs(xnj(k)).lt.dble(nl(k))-0.25d0) kappa= nl(k)
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
          do i=-1,1
            do j=1,nr
              phi(j,i)=0.d0
            end do
            e=ev(k)+de*dble(i)
            e = ev( k ) * 0.5d0 * dble( 1 - i )
            call integ(e,nl(k),kappa,neff,nn,jrt,ief,x,phi(1,i),
     &                 zeff,vps,xm1,xm2,nr,r,dr,rsqd,dl,ruse,plead)
             sum = 0.d0
             do j = 1, jrt - 1
                sum = sum +
     &                r( j ) * phi( j, i ) ** 2 +
     &                r( j + 1 ) * phi( j + 1, i ) ** 2
             end do
             sum = 0.5d0 * dl * sum
             sum = 1.d0 / dsqrt( sum )
             do j = 1, nr
                phi( j, i ) = phi( j, i ) / sum
             end do
          end do
          do j=1,nr
            write (13,'(1x,1i5,3(2x,1d19.12))') j,(phi(j,i),i=-1,1)
          end do
        end if
        e = el
        do while ( e .lt. eh )
          call integ(e,nl(k),kappa,neff,nn,jrt,ief,x,psi,                &
     &               zeff,vps,xm1,xm2,nr,r,dr,rsqd,dl,ruse,plead)
          write (iu,'(1x,2f20.10)') x, e
          e = e + ( eh - el ) / dble( n )
        end do
      end do
      deallocate( twoj, vps, xm1, xm2, phi, psi )
      return
      end
