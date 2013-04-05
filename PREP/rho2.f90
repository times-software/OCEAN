      program rhotest
      implicit none
!
      include 'fftw3.f'
      complex(kind=kind(1.d0)), allocatable :: rhoofr(:,:,:)
      complex(kind=kind(1.d0)), allocatable ::  rhoofg(:,:,:), trhoofg(:,:,:)
      integer*8 :: plan
      integer :: dims(3), counter1, counter2, counter3, dumint, cter1, cter2, cter3, natom
      character*50 :: lineburn
      real(kind=kind(1.d0)) :: norm, modG, bv1(3), bv2(3), bv3(3), b(3),&
     &      dumf, avecs(3,3)
      logical :: qestyle
!
      inquire( file="system.rho.dat", exist=qestyle )
      
      if( qestyle ) then
        write(6,*) 'Using system.rho.dat'
        open(unit=99,file='system.rho.dat',form='formatted',status='old')
        read(99,*) lineburn
        read(99,*) lineburn
        read(99,*) natom
        read(99,*) dims(1)
        read(99,*) dims(2)
        read(99,*) dims(3)
        do counter1 = 1, natom
          read(99,*) lineburn      
        enddo

        allocate(rhoofr(dims(1),dims(2),dims(3)),rhoofg(dims(1),          &
     &          dims(2),dims(3)))
        call dfftw_plan_dft_3d(plan, dims(1),dims(2),dims(3), rhoofg,     &
     &      rhoofg,FFTW_FORWARD,FFTW_ESTIMATE)
        allocate(trhoofg(dims(3),dims(2),dims(1)))
        read(99,*) trhoofg
        close(99)
        do counter3=1,dims(3)
          do counter2=1,dims(2)
            do counter1=1,dims(1)
              rhoofg(counter1,counter2,counter3) =  &
     &               trhoofg(counter3,counter2,counter1)
            enddo
          enddo
        enddo

      else
        open(unit=99,file='nfft',form='formatted',status='old')
        read(99,*) dims(:)
        close(99)
!
        allocate(rhoofr(dims(1),dims(2),dims(3)),rhoofg(dims(1),          &
     &          dims(2),dims(3)))
        call dfftw_plan_dft_3d(plan, dims(1),dims(2),dims(3), rhoofg,     &
     &      rhoofg,FFTW_FORWARD,FFTW_ESTIMATE)
!

        open(unit=99,file='rhoofr',form='formatted',status='old')
        read(99,*) lineburn
        do counter3=1,dims(3)
          do counter2=1,dims(2)
            do counter1=1,dims(1)
              read(99,*) dumint, dumint, dumint, dumf
              rhoofg(counter1,counter2,counter3) = dumf
            enddo
          enddo
        enddo
        close(99)
!
      endif

      call dfftw_execute_dft(plan, rhoofg, rhoofg)
      call dfftw_destroy_plan(plan)
!      open(unit=99,file='omega.h',form='formatted',status='old')
!      read(99,*) norm
!      close(99)
      open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
      read(99,*) avecs(:,:)
      norm = avecs(1,1)*(avecs(2,2)*avecs(3,3)-avecs(3,2)*avecs(2,3))   &
     &     - avecs(2,1)*(avecs(1,2)*avecs(3,3)-avecs(3,2)*avecs(1,3))   &
     &     + avecs(3,1)*(avecs(1,2)*avecs(2,3)-avecs(2,2)*avecs(1,3))
      close(99)
      norm = norm / dble(dims(1)*dims(2)*dims(3))
!      norm = 1.d0 / norm
      
      open(unit=17,file='bvecs',form='formatted',status='old')
       read(17,*) bv1(1),bv1(2),bv1(3)
       read(17,*) bv2(1),bv2(2),bv2(3)
       read(17,*) bv3(1),bv3(2),bv3(3)
      close(17)

      open(unit=99,file='rhoG2',form='formatted',status='unknown')
      do cter1=1,dims(1)
        do cter2=1,dims(2)
          do cter3=1,dims(3)
            counter1 = cter1
            counter2 = cter2
            counter3 = cter3
            if (cter1 .gt. dims(1)/2 ) counter1 = cter1 - dims(1)
            if (cter2 .gt. dims(2)/2 ) counter2 = cter2 - dims(2)
            if (cter3 .gt. dims(3)/2 ) counter3 = cter3 - dims(3)
            b = (/dble(counter1-1), dble(counter2-1), dble(counter3-1)/)
            modG =  (bv1(1)*b(1) + bv2(1)*b(2) + bv3(1)*b(3))**2 +      &
     &          (bv1(2)*b(1) + bv2(2)*b(2) + bv3(2)*b(3))**2 +          &
     &          (bv1(3)*b(1) + bv2(3)*b(2) + bv3(3)*b(3))**2
            write(99,'(3(I5,X),X,2(E19.12,2X),F25.10)') counter1-1,     &
     &        counter2-1, counter3-1,                                   &
     &        real(rhoofg(cter1,cter2,cter3))* norm,                    &
     &        dimag(rhoofg(cter1,cter2,cter3))* norm,modG
          enddo
        enddo
      enddo
      close(99)

      deallocate(rhoofr,rhoofg)
      end program rhotest
