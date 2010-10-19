  program kgen
  implicit none
!
  real(kind=kind(1.d0)) :: k0(3), qvec(3), ikpt, jkpt, kkpt, qpoint(3)
  integer :: core, kpttotal, coreiter, kptiter, nkpt(3), kptiter2(3), Nfiles, umklapp(3), iter
  logical :: change
  character*9 :: kptfile
!
  open(unit=99,file='k0.ipt',form='formatted',status='old')
  read(99,*) k0(:)
  close(99)
!
  open(unit=99,file='nkpt',form='formatted',status='old')
  read(99,*) nkpt(:)
  close(99)
!
  open(unit=99,file='qinunitsofbvectors.ipt',form='formatted',status='old')
  read(99,*) qvec(:)
  close(99)
!
  open(unit=99,file='core',form='formatted',status='old')
  read(99,*) core
  close(99)
!
  open(unit=99,file='nkpts',form='formatted',status='unknown')
  if ( (abs(qvec(1)) + abs(qvec(2)) + abs(qvec(3)) ) .eq. 0) then
    write(99,*)nkpt(1)*nkpt(2)*nkpt(3)
  else
    write(99,*)nkpt(1)*nkpt(2)*nkpt(3)*2.d0
  endif
  close(99)
!
! Two options, either we have a q or not, for core (or maybe testing) we have
!  no qvector.
  ikpt = k0(1)/dble(nkpt(1))
  jkpt = k0(2)/dble(nkpt(2))
  kkpt = k0(3)/dble(nkpt(3))
!
  kptiter2(:) = 1
  if ( (abs(qvec(1)) + abs(qvec(2)) + abs(qvec(3)) ) .eq. 0 ) then 
    kpttotal = nkpt(1) * nkpt(2) * nkpt(3)
!

!    core = core * ceiling(real(kpttotal)/real(50*core) )
    Nfiles = core
!
    do coreiter=1,core
      write(kptfile,'(A5,I4.4)') 'kpts.', coreiter
      open(unit=99,file=kptfile,form='formatted',status='unknown')
      write(99,*) 'nkpt', int(ceiling(real(kpttotal)/real(core)))
      write(99,*) 'kpt'
!
      do kptiter=1, int(ceiling(real(kpttotal)/real(core)))
        if ( kkpt .gt. 1 ) then
!          jkpt = jkpt + 1.d0/dble(nkpt(2))
          kkpt = kkpt - 1.d0
        endif
        if (jkpt .gt. 1 ) then
!          ikpt = ikpt + 1.d0/dble(nkpt(1))
          jkpt = jkpt - 1.d0
        endif
        if ( ikpt .gt. 1 ) ikpt = ikpt - 1.d0
!
        write(99,*) ikpt, jkpt, kkpt

        if ( kptiter2(3) .lt. nkpt(3) ) then
          kkpt = kkpt + 1.d0/dble(nkpt(3))
          kptiter2(3) = kptiter2(3) + 1
        else
          kkpt = k0(3)/dble(nkpt(3))
          kptiter2(3) = 1
          if (kptiter2(2) .lt. nkpt(2)) then
            jkpt = jkpt + 1.d0/dble(nkpt(2))
            kptiter2(2) = kptiter2(2) + 1
          else
            jkpt = k0(2)/dble(nkpt(2))
            kptiter2(2) = 1
            ikpt = ikpt + 1.d0/dble(nkpt(1))
            kptiter2(1) = kptiter2(1) + 1
          endif
        endif
      enddo
      close(99)
      kpttotal = kpttotal - int(ceiling(real(kpttotal)/real(core)))
      core = core - 1
    enddo
!!!!!!
!
  else ! make both k and k+q
    open(unit=50, file='umklapp', form='formatted', status='unknown')
    kpttotal = nkpt(1) * nkpt(2) * nkpt(3)
!    core = core * ceiling(real(kpttotal)/real(25*core) )
    Nfiles = core
!
    do coreiter=1,core
      write(kptfile,'(A5,I4.4)') 'kpts.', coreiter
      open(unit=99,file=kptfile,form='formatted',status='unknown')
      write(99,*) 'nkpt', 2*int(ceiling(real(kpttotal)/real(core)))
      write(99,*) 'kpt'
!
      do kptiter=1, int(ceiling(real(kpttotal)/real(core)))
        if ( kkpt .gt. 1 ) then
!          jkpt = jkpt + 1.d0/dble(nkpt(2))
          kkpt = kkpt - 1.d0
        endif
        if (jkpt .gt. 1 ) then
!          ikpt = ikpt + 1.d0/dble(nkpt(1))
          jkpt = jkpt - 1.d0
        endif
        if ( ikpt .gt. 1 ) ikpt = ikpt - 1.d0
!
        write(99,'(3(F14.10,X))') ikpt, jkpt, kkpt
        umklapp(:) = 0
        qpoint(1) = ikpt+qvec(1)
        qpoint(2) = jkpt+qvec(2)
        qpoint(3) = kkpt+qvec(3)
        do iter=1,3
!          change = .false.
          do
            change = .false.
            if (qpoint(iter) .gt. 1.d0) then
              change = .true.
              umklapp(iter) = umklapp(iter) + 1
              qpoint(iter) = qpoint(iter) - 1
            elseif  (qpoint(iter) .lt. 0.d0) then
              change = .true.
              umklapp(iter) = umklapp(iter) - 1
              qpoint(iter) = qpoint(iter) + 1
            endif
            if (.not. change) goto 10
          enddo
10        continue
        enddo
        write(99,'(3(F14.10,X))') qpoint(:)
        write(50, * ) umklapp
        if ( kptiter2(3) .lt. nkpt(3) ) then 
          kkpt = kkpt + 1.d0/dble(nkpt(3))
          kptiter2(3) = kptiter2(3) + 1
        else
          kkpt = k0(3)/dble(nkpt(3))
          kptiter2(3) = 1
          if (kptiter2(2) .lt. nkpt(2)) then
            jkpt = jkpt + 1.d0/dble(nkpt(2))
            kptiter2(2) = kptiter2(2) + 1
          else
            jkpt = k0(2)/dble(nkpt(2))
            kptiter2(2) = 1
            ikpt = ikpt + 1.d0/dble(nkpt(1))
            kptiter2(1) = kptiter2(1) + 1
          endif
        endif
      enddo
      close(99)
      kpttotal = kpttotal - int(ceiling(real(kpttotal)/real(core)))
      core = core - 1
    enddo
  
    close(50) 
  endif
  open(unit=99,file='Nfiles',form='formatted',status='unknown')
  write(99,*) Nfiles
  close(99)

  end program kgen
