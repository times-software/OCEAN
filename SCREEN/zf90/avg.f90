! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program avg

  integer ::ng,i,j,numsites,elnum
  integer, allocatable :: gvec(:,:)
  real(kind=kind(1.d0)) :: denr, deni, pi, bvec(3,3), avec(3,3),vol,mag, &
            tau(3),radius,gr,greal(3), magi
  real(kind=kind(1.d0)), allocatable :: rhogr(:),rhogi(:),modrealgvec(:)
  real(kind=kind(1.d0)), allocatable :: modrealgvec2(:)
  integer :: gmax( 3 ), ngmax, iter, rad_int, dumi
  real(kind=kind(1.d0)) :: gmodmax( 3 ), sgmodmax
  character(len=2) elname
  character(len=9) avgname

  gmax( : ) = 0
  gmodmax( : ) = 0.d0
  pi = 4.d0*datan(1.d0)

  open(unit=99,file='rhoofg',form='formatted',status='old')
  rewind 99
  read(99,*) ng
  allocate( gvec( 3, ng ), rhogr( ng ), rhogi( ng ), modrealgvec2( ng ) )
  do i=1,ng
    read(99,*)gvec(:,i),rhogr(i),rhogi(i),modrealgvec2(i)
    do iter = 1, 3
      if( gmax( iter ) .lt. abs( gvec( iter, i ) ) ) then
        gmax( iter ) = abs( gvec( iter, i ) )
        gmodmax( iter ) = modrealgvec2(i)
      else if( gmax( iter ) .eq. abs( gvec( iter, i ) ) ) then
        gmodmax( iter ) = min( gmodmax( iter ), modrealgvec2( i ) )
      endif
    enddo
  enddo
  close(99)
  sgmodmax = maxval( gmodmax )
!  write(6,*) maxval( gmodmax ), sgmodmax
  write(6,*) ng
  ngmax = ng
  do i = 1, ng
    if( modrealgvec2( i ) .gt. sgmodmax ) then
      ngmax = i - 1
!      write(6,*) ngmax
      goto 10
    endif
  enddo
      write(6,*) 'did not work'
10  continue
  ng = ngmax
  write(6,*) ng

  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  rewind 99
  read(99,*)avec(:,:)
  close(99)
  open(unit=99,file='bvecs',form='formatted',status='old')
  rewind 99
  read(99,*)bvec(:,:)
  close(99)
! move to getomega
  call getomega( avec, vol )
!      vol = avec(1,1)* (avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3) ) -   &
!     &      avec(2,1)* (avec(1,2)*avec(3,3) - avec(3,2)*avec(1,3) ) +   &
!     &      avec(3,1)* (avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3) ) 

      write(6,*) "uc vol = ",vol
      rhogr( : ) = rhogr( : ) / vol
      rhogi( : ) = rhogi( : ) / vol

      allocate(modrealgvec(ng))
      do i=1,ng
        mag = 0.d0
        do j =1,3
           greal(j) = gvec(1,i)*bvec(j,1) + gvec(2,i)*bvec(j,2) +       &
     &                gvec(3,i)*bvec(j,3)
        enddo
        mag = greal(1)*greal(1)+greal(2)*greal(2)+greal(3)*greal(3)
        modrealgvec(i) = dsqrt(mag)
        write(70,*) modrealgvec(i), dsqrt(modrealgvec2(i))
      enddo !ng
      modrealgvec( 1 : ng ) = dsqrt( modrealgvec2( 1 : ng ) )
      
      open(unit=98,file='sitelist',form='formatted',status='old')
      read(98,*)numsites

      do i=1,numsites

        read(98,*)elname,dumi,elnum
        write(avgname,"(a3,a2,i4.4)")"avg",elname,elnum
        open(unit=97,file=avgname,form='formatted',status='unknown')
        call snatch(elname,elnum,tau)

!        do radius=0.00001d0,40.0,0.1d0
        do rad_int = 0, 400
          radius = 0.00001d0 + dble( rad_int ) / 10.d0
          denr = 0.d0
          deni = 0.d0
!$OMP PARALLEL DO  &
!$OMP& DEFAULT( NONE ) &
!$OMP& PRIVATE( j, gr, mag, magi ) &
!$OMP& SHARED( radius, ng, rhogr, rhogi, tau, gvec, modrealgvec, pi ) &
!$OMP REDUCTION(+:denr,deni)
          do j=1,ng
            gr = modrealgvec(j) * radius
            if (gr .gt. 0.0001d0 ) then
              mag = (dsin(gr)/gr)*(rhogr(j))
              magi = (dsin(gr)/gr)*(rhogi(j))
            else  
              ! The first two terms of the sin expansion divided by x -> ( 1 - x^2/3! )
              !   the next term is x^4/5! ---->  0.0001^4/120 < 9e-19
              mag  = rhogr(j) * ( 1.0d0 - ( gr*gr / 6.0d0 ) )
              magi = rhogi(j) * ( 1.0d0 - ( gr*gr / 6.0d0 ) )
            endif
 
            denr = denr + dcos( 2.0d0 * pi * &
     &               ( dble( gvec(1,j) ) * tau(1) + &
     &                 dble( gvec(2,j) ) * tau(2) + &
     &                 dble( gvec(3,j) ) * tau(3) ) ) * mag
            denr = denr + dsin( 2.0d0 * pi * &
     &               ( dble( gvec(1,j) ) * tau(1) + &
     &                 dble( gvec(2,j) ) * tau(2) + &
     &                 dble( gvec(3,j) ) * tau(3) ) ) * magi
            deni = deni - dcos( 2.0d0 * pi * &
     &               ( dble( gvec(1,j) ) * tau(1) + &
     &                 dble( gvec(2,j) ) * tau(2) + &
     &                 dble( gvec(3,j) ) * tau(3) ) ) * magi
            deni = deni + dsin( 2.0d0 * pi * &
     &               ( dble( gvec(1,j) ) * tau(1) + &
     &                 dble( gvec(2,j) ) * tau(2) + &
     &                 dble( gvec(3,j) ) * tau(3) ) ) * mag


          enddo !ng
!$OMP END PARALLEL DO
          write(97,"(a2, 1x, i2.2, 1x, e18.11, 1x, e18.11, 1x, e18.11)")elname,elnum,radius,denr,deni
        enddo !radius
        close( 97 )
      enddo !numsites

      close(98)
      end
