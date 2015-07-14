! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
  program excitonplot
  implicit none
!
  integer :: idum, indx, atno, ncore, lc, nc, numbands, totaln, nkpts, dumi, ng, startband, gsize(3)
  character * 2 :: element
  real(kind=kind(1.d0)) :: tau(3), qbase(3), qraw(3), phase, phaser, phasei, pi, avecs(3,3), sino, coso, prefac, coface(3,3), newvecs(3,3), deta
  real(kind=kind(1.d0)), allocatable :: rex(:,:), psiGr(:,:,:), psiGi(:,:,:), zzr(:,:), zzi(:,:), rexbnd(:,:)
  integer, allocatable :: g(:,:)
  integer :: ngfft(2,3), kmesh(3), brange(4), iter, iter2, biter, giter, qiter, ierr, ik1, ik2, ik3
  integer :: xiter, yiter, ziter, nbtot
  logical :: metal
  character *20, allocatable :: listw(:)
  complex(kind=kind(1.d0)) :: rm1
  complex(kind=kind(1.d0)), allocatable ::  psi(:,:,:), COMM(:), phi(:,:,:), ztot(:)
  integer :: sizer, xtics, ytics, ztics
  real :: convcell(3), supercell(3), bvecs(3,3), carcoord(3), redcoord(3)
!
  pi = 4.d0*datan(1.d0)
  rm1 = -1.d0
  rm1 = sqrt(rm1)
  xtics = 10
  ytics = 10
  ztics = 10

! no support for multiple sites per run through
  open( unit=99, file='sitelist', form='formatted', status='old')
  read(99,*) idum
  read(99,*) element, indx
  close(99)
  call snatch( element, indx, tau)
!
  open( unit=99, file='ZNL', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) atno, ncore, lc
  close(99)
!
  nc = 4 * ( 2 * lc + 1 )
  open( unit=99, file='wvfcninfo', form='unformatted', status='old' )
  rewind 99
  read ( 99 ) numbands, nkpts
  totaln = nkpts * numbands
  write ( 6, * ) 'n, nc', totaln, nc, numbands
!
! Rex contains the amplitude info for the exciton
  sizer = totaln*nc
  allocate(rex( sizer, 2), rexbnd( numbands, 2) )
  rex(:,:) = 0.d0
  open( unit=99, file='echamp', form='unformatted', status='old' )
  rewind 99
  read ( 99 ) rex
  close( unit=99 )
!
!
  open( unit=99, file='ngfft', form='formatted', status='old')
  read(99,*) ngfft(:,:)
  close(99)
!
! psiG will store all our info
!  allocate(psiGr(ngfft(1,1):ngfft(2,1), ngfft(1,2):ngfft(2,2), ngfft(1,3):ngfft(2,3) ),  &
  allocate( psi(ngfft(1,1):ngfft(2,1), ngfft(1,2):ngfft(2,2), ngfft(1,3):ngfft(2,3) ),  &
          phi(ngfft(1,1):ngfft(2,1), ngfft(1,2):ngfft(2,2), ngfft(1,3):ngfft(2,3) ) )
  psi(:,:,:) = 0.d0
!
! need to log the band we started from for each of the kpts
! for right now, assume not metal
  open( unit=99, file='metal', form='formatted', status='old')
  rewind 99
  read(99, * ) metal
  close( 99 )
  if (metal) stop 'Metal not supported'
!
  open(unit=99, file='brange.ipt', form='formatted', status='old')
  read(99,*) brange(:)
  close(99)
  nbtot = brange(4)-brange(3)+brange(2)-brange(1)+2
  startband = brange(2)-brange(1)+2
  write(6,*) startband, numbands
!
  allocate(listw(nkpts))
  open(unit=99,file='listwfile',form='formatted',status='old')
  do iter=1,nkpts
    read(99,*) dumi, listw(iter)
  enddo
  close(99)
!
  open( unit=99, file='scaledkzero.ipt', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) qbase( : )
  close( unit=99 )
!
  open( unit=99, file='kmesh.ipt', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) kmesh( : )
  close( unit=99 ) 
!
  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old')
  read(99,*) avecs(:,:)
  close(99)
!
  open( unit=99, file='bvecs', form='formatted', status='old')
  read(99,*) bvecs(:,:)
  close(99)
!
  do xiter=1,3
    do yiter=1,3
      newvecs(xiter,yiter) = avecs(yiter,xiter)
    enddo
  enddo
  deta = avecs(1,1) * (avecs(2,2)*avecs(3,3) - avecs(2,3)*avecs(3,2)) &
       - avecs(1,2) * (avecs(2,1)*avecs(3,3) - avecs(2,3)*avecs(3,1)) &
       + avecs(1,3) * (avecs(2,1)*avecs(3,2) - avecs(2,2)*avecs(3,2))
  coface(1,1) = newvecs(2,2)*newvecs(3,3) - newvecs(2,3)*newvecs(3,2)
  coface(1,2) = -(newvecs(2,1)*newvecs(3,3) - newvecs(2,3)*newvecs(3,1))
  coface(1,3) = newvecs(2,1)*newvecs(3,2) - newvecs(2,2)*newvecs(3,1)
  coface(2,1) = -(newvecs(1,2)*newvecs(3,3) - newvecs(1,3)*newvecs(3,2))
  coface(2,2) = newvecs(1,1)*newvecs(3,3) - newvecs(3,1)*newvecs(1,3)
  coface(2,3) = -(newvecs(1,1)*newvecs(3,2) - newvecs(1,2)*newvecs(3,1))
  coface(3,1) = newvecs(1,2)*newvecs(2,3) - newvecs(1,3)*newvecs(2,2)
  coface(3,2) = -(newvecs(1,1)*newvecs(2,3) - newvecs(1,3)*newvecs(2,1))
  coface(3,3) = newvecs(1,1)*newvecs(2,2) - newvecs(1,2)*newvecs(2,1)
  newvecs(:,:) = coface(:,:) / deta
  write(6,*) avecs(:,:)
  write(6,*)
  write(6,*) newvecs(:,:)
!
  convcell(:) = 7.596698474d0
!
!
  gsize(1) = ngfft(2,1)-ngfft(1,1)+1
  gsize(2) = ngfft(2,2)-ngfft(1,2)+1
  gsize(3) = ngfft(2,3)-ngfft(1,3)+1
!
  allocate(COMM(gsize(1)*gsize(2)*gsize(3)+3*(gsize(1)+gsize(2)+gsize(3)) ) )
!
! Loop over kpoints, reading in wavefunctions and processing
  write ( 6, * ) 'looping over k to load...'
  iter=1
  qiter=1
  do ik1 = 0, kmesh( 1 ) - 1
    do ik2 = 0, kmesh( 2 ) - 1
      do ik3 = 0, kmesh( 3 ) - 1
        qraw( 1 ) = ( qbase( 1 ) + dble( ik1 ) ) / dble( kmesh( 1 ) ) 
        qraw( 2 ) = ( qbase( 2 ) + dble( ik2 ) ) / dble( kmesh( 2 ) ) 
        qraw( 3 ) = ( qbase( 3 ) + dble( ik3 ) ) / dble( kmesh( 3 ) ) 
!
        open(unit=50,file=listw(qiter),form='unformatted',status='old')
        read(50) ng
        allocate( g(ng,3), zzr(ng, nbtot), zzi(ng, nbtot), ztot(ng) )
        ztot(:) = 0.d0
        read(50) g
        read(50) zzr
        read(50) zzi
        close(50)
!
        phi(:,:,:) = 0.d0
!        phase = qraw(1) * tau(1) + qraw(2) * tau(2) + qraw(3) * tau(3)
!        phaser = dcos(2*pi*phase)
!        phasei = -dsin(2*pi*phase)
!
        rexbnd(:,:) = 0.d0
        do biter=1, numbands
          do iter2=1,nc
            rexbnd(biter,1) = rexbnd(biter,1) + rex(iter,1)
            rexbnd(biter,2) = rexbnd(biter,2) + rex(iter,2)
            iter = iter+1
          enddo
        enddo
        do giter=1, ng
          do biter=startband, numbands+startband-1
!            psiGr(g(giter,1), g(giter,2), g(giter,3)) =  psiGr(g(giter,1), g(giter,2), g(giter,2)) &
!                                        + zzr(giter,biter) * phaser * rexbnd(biter-startband+1,1)
!            psiGi(g(giter,1), g(giter,2), g(giter,3)) =  psiGi(g(giter,1), g(giter,2), g(giter,2)) &
!                                        + zzi(giter,biter) * phasei * rexbnd(biter-startband+1,2) 
!!!            phi(g(giter,1), g(giter,2), g(giter,3)) =  zzr(giter,biter) * rexbnd(biter-startband+1,1) &
!!!                                             + rm1 * zzi(giter,biter) * rexbnd(biter-startband+1,2)
            ztot(giter) =  ztot(giter) + zzr(giter,biter) * rexbnd(biter-startband+1,1) &
                                       + rm1 * zzi(giter,biter) * rexbnd(biter-startband+1,2)

          enddo

        enddo
!       
        qiter = qiter+1
!!!        deallocate(g, zzr, zzi )

!!!        call ZFFT3D(1,gsize(1),gsize(2),gsize(3),phi,COMM,ierr)
!!!        phi(:,:,:) = phi(:,:,:) * dsqrt(dble(gsize(1)*gsize(2)*gsize(3)))
!!!        if (ierr .ne. 0) stop 'fft failed'
        do xiter=-xtics,xtics    !ngfft(1,1),ngfft(2,1)
          do yiter=-ytics,ytics    !ngfft(1,2),ngfft(2,2)
            do ziter=-ztics, ztics    !ngfft(1,3),ngfft(2,3)
              do giter=1,ng
                carcoord(1) = xiter / xtics * convcell(1)
                carcoord(2) = yiter / ytics * convcell(2)
                carcoord(3) = ziter / ztics * convcell(3)
                do iter=1,3
                  redcoord(:) = newvecs(:,iter) * carcoord(iter)
                enddo
                
                phase = (g(giter,1) + qraw(1)) * (redcoord(1) - tau(1)) & ! xiter / gsize(1) - tau(1)) &
                      + (g(giter,2) + qraw(2)) * (redcoord(2) - tau(2)) &
                      + (g(giter,3) + qraw(3)) * (redcoord(3) - tau(3))
                phase = phase * 2 * pi
                call FASTSINCOS(phase, sino, coso)
                prefac = coso + rm1*sino
!!!                psi(xiter,yiter,ziter) = psi(xiter,yiter,ziter) + phi(xiter,yiter,ziter)*prefac
                psi(xiter,yiter,ziter) = psi(xiter,yiter,ziter) + ztot(giter)*prefac
!              psi(ngfft(1,1)+xiter-1, ngfft(1,2)+yiter-1, ngfft(1,3)+ziter-1)  &
!                          = psi(ngfft(1,1)+xiter-1, ngfft(1,2)+yiter-1, ngfft(1,3)+ziter-1)  &
!                          + phi(ngfft(1,1)+xiter-1, ngfft(1,2)+yiter-1, ngfft(1,3)+ziter-1)*prefac
              enddo
            enddo
          enddo
        enddo
      
      enddo ! ik3
    enddo ! ik2
    write(6,'(F5.1,A1)') 100.0*dble(1+ik1)/dble(kmesh(1)), '%'
  enddo ! ik1
  write(6,*) 'Bands all loaded' 
!
!
!
!!  allocate(psiG( ngfft(1,1):ngfft(2,1), ngfft(1,2):ngfft(2,2), ngfft(1,3):ngfft(2,3) ) )
!!  rm1 = -1.d0
!!  rm1 = sqrt(rm1)
  
!!  psiG(:,:,:) = psiGr(:,:,:) + rm1*psiGi(:,:,:)
!!  deallocate(psiGr, psiGi)
!
!
!
!!  gsize(1) = ngfft(2,1)-ngfft(1,1)
!!  gsize(2) = ngfft(2,2)-ngfft(1,2)
!!  gsize(3) = ngfft(2,3)-ngfft(1,3)
  write(6,*) gsize(:)
!!  allocate(COMM(gsize(1)*gsize(2)*gsize(3)+3*(gsize(1)+gsize(2)+gsize(3)) ) )
!!  call ZFFT3D(1,gsize(1),gsize(2),gsize(3),psiG,COMM,ierr)
!!  if (ierr .ne. 0) stop 'fft failed'
!!  deallocate(COMM)
!
  write(6,*) 'writing out wf'
  open(unit=99,file='output',form='formatted',status='unknown')
  write(99,*) 'dunno'
  write(99,*) 'bah'
  write(99,'(1i5,3(1x,1e15.8))') 0, -convcell(:)
  write(99,'(1i5,3(1x,1e15.8))') xtics*2+1, convcell(1)/dble(xtics*2+1), 0, 0 !gsize(1), avecs(1,:)/gsize(1)
  write(99,'(1i5,3(1x,1e15.8))') ytics*2+1, 0, convcell(2)/dble(ytics*2+1), 0 !gsize(2), avecs(2,:)/gsize(2)
  write(99,'(1i5,3(1x,1e15.8))') ztics*2+1, 0, 0, convcell(3)/dble(ztics*2+1) !gsize(3), avecs(3,:)/gsize(2)
  do xiter=-xtics, xtics !ngfft(1,1),ngfft(2,1)
    do yiter=-ytics, ytics !ngfft(1,2),ngfft(2,2)
      do ziter=-ztics, ztics !ngfft(1,3),ngfft(2,3)
        write(99,*) real(psi(xiter,yiter,ziter) * conjg(psi(xiter,yiter,ziter)))
      enddo
    enddo 
  enddo
  close(99)
!  deallocate(psiG)

  end program
