! Copyright (C) 2014, 2015, 2017, 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program kgen
  implicit none
!
  integer, parameter :: qp = selected_real_kind(33, 4931)
  real(qp) :: k0(3), qvec(3), ikpt, jkpt, kkpt, qpoint(3)
!  real(kind=kind(1.d0)) :: k0(3), qvec(3), ikpt, jkpt, kkpt, qpoint(3)
  integer :: kpttotal, nkpt(3), umklapp(3), iter
  integer :: shift_fh, qe_shift_fh, ikx, iky, ikz
  logical :: change
  character(len=9)  :: kptfile
  character(len=12) :: qekptfile
  real(kind=kind(1.d0)), parameter :: small = 10.0d0 * EPSILON(1.d0)
  
  logical :: have_shift, split_dft, newK0

  newK0 = .true.
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

  if ( (abs(qvec(1)) + abs(qvec(2)) + abs(qvec(3)) ) .gt. small ) then 
    have_shift = .true.
    inquire( file='dft.split', exist=split_dft )
    if( split_dft ) then
      open(unit=99, file='dft.split',form='formatted',status='old')
      read(99,*) split_dft
      close(99)
    endif
    if( split_dft ) write(6,*) 'DFT calc will be split between two files.'
  else
    have_shift = .false.
    split_dft = .false.
  endif

  if( have_shift ) then
    kpttotal = nkpt(1) * nkpt(2) * nkpt(3) * 2
  else
    kpttotal = nkpt(1) * nkpt(2) * nkpt(3)
  endif
  
  if( split_dft ) then
    shift_fh = 98
    qe_shift_fh = 96
    ! override kpt total
    kpttotal = nkpt(1) * nkpt(2) * nkpt(3)
  else
    shift_fh = 99
    qe_shift_fh = 97
  endif
!
  open(unit=99,file='nkpts',form='formatted',status='unknown')
  write(99,*) kpttotal
  close(99)
!




  ! Open files here
  write(kptfile,'(A5,I4.4)') 'kpts.', 1
  open(unit=99,file=kptfile,form='formatted',status='unknown')
  write(99, * ) 'nkpt', kpttotal
  write(99,*) 'kpt'

  write(qekptfile,'(A8,I4.4)') 'kpts4qe.', 1
  open(unit=97,file=qekptfile,form='formatted',status='unknown')

  if( split_dft ) then
    write(kptfile,'(A5,I4.4)') 'kpts.', 2
    open(unit=98,file=kptfile,form='formatted',status='unknown')
    write(98,*) 'nkpt', kpttotal
    write(98,*) 'kpt'

    write(qekptfile,'(A8,I4.4)') 'kpts4qe.', 2
    open(unit=96,file=qekptfile,form='formatted',status='unknown')
  endif


  open(unit=50, file='umklapp', form='formatted', status='unknown')
  ! 

  do ikx = 0, nkpt(1)-1
    if( newK0 ) then
      ikpt = k0(1)/real(nkpt(1),QP) + real(ikx,QP)/real(nkpt(1),QP) - qvec(1)
    else
      ikpt = k0(1)/real(nkpt(1),QP) + real(ikx,QP)/real(nkpt(1),QP)
    endif
    do iky = 0, nkpt(2)-1
      if( newK0 ) then
        jkpt = k0(2)/real(nkpt(2),QP) + real(iky,QP)/real(nkpt(2),QP) - qvec(2)
      else
        jkpt = k0(2)/real(nkpt(2),QP) + real(iky,QP)/real(nkpt(2),QP)
      endif
      do ikz = 0, nkpt(3)-1
        if( newK0 ) then
          kkpt = k0(3)/real(nkpt(3),QP) + real(ikz,QP)/real(nkpt(3),QP) - qvec(3)
        else
          kkpt = k0(3)/real(nkpt(3),QP) + real(ikz,QP)/real(nkpt(3),QP) 
        endif

        write(99,'(e19.10,e19.10,e19.10)') ikpt, jkpt, kkpt
        write(97,'(e19.10,e19.10,e19.10,f19.11)') ikpt, jkpt, kkpt, ( 1.0_QP / real(kpttotal,QP) )


        if( have_shift ) then

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
10          continue
          enddo

          write(shift_fh, '(e19.10,e19.10,e19.10)' ) qpoint(:)
          write(qe_shift_fh, '(e19.10,e19.10,e19.10,f19.11)')  qpoint(:) , ( 1.0_QP / real(kpttotal,QP) )

          write(50, * ) umklapp
        endif
      enddo
    enddo
  enddo


  close(99)
  close(97)
  if( split_dft ) then
    close(96)
    close(98)
  endif
  
  close(50) 

end program kgen
