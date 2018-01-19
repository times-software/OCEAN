! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program avg
  implicit none

  include 'mpif.h'

  integer ::ng,i,j,numsites
  integer, allocatable :: gvec(:,:),elnum(:)
  real(kind=kind(1.d0)) :: denr, deni, pi, bvec(3,3), avec(3,3),vol,mag, &
            radius,gr,greal(3), magi
  real(kind=kind(1.d0)), allocatable :: rhogr(:),rhogi(:),modrealgvec(:)
  real(kind=kind(1.d0)), allocatable :: modrealgvec2(:), tau(:,:)
  integer :: gmax( 3 ), ngmax, iter, rad_int
  real(kind=kind(1.d0)) :: gmodmax( 3 ), sgmodmax, rad_step

  character(len=2), allocatable :: elname(:)
  character(len=9) :: avgname

  real(kind=kind(1.d0)), allocatable :: Vdenr(:), Vdeni(:)

  integer :: ierr, myrank, pool_size, pool_comm, pool_id, num_pools, my_poolrank, mypool
  integer :: max_rad, rad_remain, start_rad, my_rad, my_start_rad, temp_rad

  logical :: ex


  call MPI_INIT( ierr )
  if( ierr .ne. MPI_SUCCESS ) goto 111

  ! super lazy
  ! Later we can have multiple pools. Divis the sites by pool
  pool_comm = MPI_COMM_WORLD
  pool_id = 0
  num_pools = 1
  mypool = 0

  call MPI_COMM_RANK( pool_comm, myrank, ierr )
  if( ierr .ne. MPI_SUCCESS ) goto 111

  my_poolrank = myrank

  call MPI_COMM_SIZE( pool_comm, pool_size, ierr )
  if( ierr .ne. MPI_SUCCESS ) goto 111

  if( myrank .eq. 0 .and. mypool .eq. 0 ) then
    inquire( file='avg.ipt', exist=ex )
    if( ex ) then
      open( unit=99, file='avg.ipt', form='formatted', status='old' )
      read( 99, * ) max_rad, rad_step
      close( 99 )
    else
      max_rad = 401
      rad_step = 1.0d0/ 10.0d0
    endif
  endif
  

  call MPI_BCAST( max_rad, 1, MPI_INTEGER, 0, pool_comm, ierr )
  call MPI_BCAST( rad_step, 1, MPI_DOUBLE_PRECISION, 0, pool_comm, ierr )


  pi = 4.d0*datan(1.d0)

!  max_rad = 401
  rad_remain = max_rad
  start_rad = 0
  do i = 0, pool_size - 1
    temp_rad = rad_remain / ( pool_size - i )
    if( my_poolrank .eq. i ) then 
      my_rad = temp_rad
      my_start_rad = start_rad
    endif
    start_rad = start_rad + temp_rad
    rad_remain = rad_remain - temp_rad
  enddo

  write(6,*) my_poolrank, my_rad, my_start_rad
    
  if( my_poolrank .eq. 0 ) then
    allocate( Vdenr( max_rad ), Vdeni( max_rad ) )
  else
    allocate( Vdenr( my_rad ), Vdeni( my_rad ) )
  endif


  if( myrank .eq. 0 .and. pool_id .eq. 0 ) then

    gmax( : ) = 0
    gmodmax( : ) = 0.d0

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
!         greal(j) = gvec(1,i)*bvec(1,j) + gvec(2,i)*bvec(2,j) +       &
!   &                gvec(3,i)*bvec(3,j)
         greal(j) = gvec(1,i)*bvec(j,1) + gvec(2,i)*bvec(j,2) +       &
   &                gvec(3,i)*bvec(j,3)
      enddo
      mag = greal(1)*greal(1)+greal(2)*greal(2)+greal(3)*greal(3)
      modrealgvec(i) = dsqrt(mag)
!      write(70,*) modrealgvec(i), dsqrt(modrealgvec2(i))
    enddo !ng
    modrealgvec( 1 : ng ) = dsqrt( modrealgvec2( 1 : ng ) )
    
    open(unit=98,file='sitelist',form='formatted',status='old')
    read(98,*)numsites
    allocate( tau(3,numsites), elname(numsites), elnum(numsites) )
    do i = 1, numsites
      read(98,*)elname(i),elnum(i)
      call snatch( elname(i),elnum(i), tau(:,i) )
    enddo
    close(98)

  endif

  if( pool_size .gt. 1 .or. num_pools .gt. 1 ) then
    call MPI_BCAST( ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    call MPI_BCAST( numsites, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

    if( myrank .ne. 0 .or. mypool .ne. 0 ) then
      allocate( rhogr( ng ), rhogi( ng ), gvec( 3, ng ), modrealgvec( ng ) )
      allocate( tau( 3, numsites ) )
    endif


    call MPI_BCAST( rhogr, ng, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_BCAST( rhogi, ng, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_BCAST( gvec, 3*ng, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    call MPI_BCAST( modrealgvec, ng, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    call MPI_BCAST( tau, 3*numsites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  
  endif


  do i=1,numsites
    Vdenr(:) = 0.0d0
    Vdeni(:) = 0.0d0

!    read(98,*)elname,elnum
!    write(avgname,"(a3,a2,i2.2)")"avg",elname,elnum
!    open(unit=97,file=avgname,form='formatted',status='unknown')
!    call snatch(elname,elnum,tau)

!!        do radius=0.00001d0,40.0,0.1d0
!    do rad_int = 0, 400
    do rad_int = 1, my_rad 
      radius = 0.00001d0 + dble( rad_int - 1 + my_start_rad) * rad_step
      denr = 0.d0
      deni = 0.d0
!$OMP PARALLEL DO  &
!$OMP& DEFAULT( NONE ) &
!$OMP& PRIVATE( j, gr, mag, magi ) &
!$OMP& SHARED( radius, ng, rhogr, rhogi, tau, gvec, modrealgvec, pi ) &
!$OMP REDUCTION(+:denr,deni)
      do j=1,ng
        gr = modrealgvec(j) * radius
        if (gr .gt. 0.00001 ) then
          mag = (dsin(gr)/gr)*(rhogr(j))
          magi = (dsin(gr)/gr)*(rhogi(j))
        else
          mag = rhogr(j) * (1.0d0 - gr*gr/4.0 )
          magi = rhogi(j) *(1.0d0 - gr*gr/4.0 )
        endif

        denr = denr + dcos( 2.0d0 * pi * &
 &               ( dble( gvec(1,j) ) * tau(1,i) + &
 &                 dble( gvec(2,j) ) * tau(2,i) + &
 &                 dble( gvec(3,j) ) * tau(3,i) ) ) * mag
        denr = denr + dsin( 2.0d0 * pi * &
 &               ( dble( gvec(1,j) ) * tau(1,i) + &
 &                 dble( gvec(2,j) ) * tau(2,i) + &
 &                 dble( gvec(3,j) ) * tau(3,i) ) ) * magi
        deni = deni - dcos( 2.0d0 * pi * &
 &               ( dble( gvec(1,j) ) * tau(1,i) + &
 &                 dble( gvec(2,j) ) * tau(2,i) + &
 &                 dble( gvec(3,j) ) * tau(3,i) ) ) * magi
        deni = deni + dsin( 2.0d0 * pi * &
 &               ( dble( gvec(1,j) ) * tau(1,i) + &
 &                 dble( gvec(2,j) ) * tau(2,i) + &
 &                 dble( gvec(3,j) ) * tau(3,i) ) ) * mag


      enddo !ng
!$OMP END PARALLEL DO
      Vdenr( rad_int ) = denr
      Vdeni( rad_int ) = deni
!      write(97,"(a2, 1x, i2.2, 1x, e17.11, 1x, e17.11, 1x, e17.11)")elname,elnum,radius,denr,deni
      
    enddo !radius
!    close( 97 )
  rad_remain = max_rad
  start_rad = 0
  do j = 0, pool_size - 1
    temp_rad = rad_remain / ( pool_size - j )
    if( my_poolrank .eq. 0 ) then
      if( j .ne. 0 ) then
        call MPI_RECV( Vdenr( start_rad+1 ), temp_rad, MPI_DOUBLE_PRECISION, j, j, pool_comm, MPI_STATUS_IGNORE, ierr )
        call MPI_RECV( Vdeni( start_rad+1 ), temp_rad, MPI_DOUBLE_PRECISION, j, j, pool_comm, MPI_STATUS_IGNORE, ierr )
      endif
    elseif( my_poolrank .eq. j ) then
      call MPI_SEND(Vdenr, temp_rad, MPI_DOUBLE_PRECISION, 0, j, pool_comm, ierr )
      call MPI_SEND(Vdeni, temp_rad, MPI_DOUBLE_PRECISION, 0, j, pool_comm, ierr )
    endif
    start_rad = start_rad + temp_rad
    rad_remain = rad_remain - temp_rad
  enddo

  if( my_poolrank .eq. 0 ) then
    write(avgname,"(a3,a2,i4.4)")"avg",elname(i),elnum(i)
    open(unit=97,file=avgname,form='formatted',status='unknown')
    do rad_int = 0, max_rad-1
      radius = 0.00001d0 + dble( rad_int ) * rad_step
      write(97,"(a2, 1x, i2.2, 1x, e17.11, 1x, e17.11, 1x, e17.11)")elname(i),elnum(i),radius,Vdenr(rad_int+1),Vdeni(rad_int+1)
    enddo
    close(97)
  endif



  enddo !numsites

  close(98)

111 continue

  call MPI_FINALIZE( ierr )

  end
