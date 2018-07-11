! Copyright (C) 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2018
!
!
module screen_opf
  use ai_kinds, only : DP

  implicit none
  private

  type opf_holder
    integer :: Zee
    integer :: lMin
    integer :: lMax
    integer :: nRad
    integer, allocatable :: nprojPerChannel( : )

    real( DP ) :: rMax
    real( DP ), allocatable :: rad( : )
    real( DP ), allocatable :: aeProj( :, :, : )
    real( DP ), allocatable :: psProj( :, :, : )
!!    real( DP ), allocatable :: aMat( : , :, : )

  end type opf_holder

  type( opf_holder ), allocatable, save :: FullTable( : )

  public :: screen_opf_init, screen_opf_clean, screen_opf_makeNew, screen_opf_loadAll, &
            screen_opf_lbounds, screen_opf_getNCutoff, screen_opf_nprojForChannel, &
            screen_opf_interpProjs, screen_opf_makeAMat, screen_opf_maxNproj, screen_opf_AltInterpProjs

  contains

  subroutine screen_opf_makeAMat( np, nr, rad, drad, proj, amat, ierr )
    integer, intent( in ) :: np, nr
    real(DP), intent( in ) :: rad(:), drad(:), proj(:,:)
    real(DP), intent( out ) :: amat(:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tmat(:,:), rad2drad(:)
    real(DP) :: su
    integer :: i, j, k
    

    allocate( rad2drad( nr ), tmat( np, np ) )

    do i = 1, nr
      rad2drad( i ) = rad( i ) * rad( i ) * drad( i )
    enddo

    do j = 1, np
      do i = 1, j
        
        su = 0.0_DP
        do k = 1, nr
          su = su + rad2drad( k ) * proj( k, i ) * proj( k, j )
        enddo
        tmat( i, j ) = su
        tmat( j, i ) = su
      enddo
    enddo

!    write(6,*) np
!    write(6,*) tmat(:,:)
    call rinvert( np, tmat, ierr )
    amat( 1:np, 1:np ) = tmat( :, : )
    deallocate( tmat, rad2drad )

  end subroutine screen_opf_makeAMat
    

  pure function isRightTarg( zee, targ ) result( isRight )
    integer, intent( in ) :: zee, targ
    logical :: isright

    if( targ .lt. 0 .or. targ .gt. size( FullTable ) ) then
      isRight = .false.
    else
      isRight = ( zee .eq. FullTable( targ )%Zee )
    endif
  end function
    
  pure function getRightTarg( zee ) result( targ )
    integer, intent( in ) :: zee
    integer :: targ
    integer :: i

    do i = 1, size( FullTable )
      if( FullTable( i )%Zee .eq. zee ) then
        targ = i
        return
      endif
    enddo
    targ = -1
  end function

  subroutine screen_opf_maxNproj( zee, maxNproj, ierr, itarg )
    integer, intent( in ) :: zee
    integer, intent( out ) :: maxNproj
    integer, intent( inout ) :: ierr
    integer, intent( inout ), optional :: itarg
    !
    integer :: targ, lmin, lmax, temp, l

    if( present( itarg ) ) then
      targ = itarg
    else
      targ = -1
    endif
  
    call screen_opf_lbounds( zee, lmin, lmax, ierr, targ )
    if( ierr .ne. 0 ) return

    maxNproj = 0
    do l = lmin, lmax
      call screen_opf_nprojForChannel( zee, l, temp, ierr, targ )
      if( ierr .ne. 0 ) return
      maxNproj = max( temp, maxNproj )
    enddo

    ! zero means something has gone horribly wrong
    if( maxNproj .eq. 0 ) ierr = 1
  end subroutine screen_opf_maxNproj
    

  subroutine screen_opf_lbounds( zee, lmin, lmax, ierr, itarg )
    integer, intent( in ) :: zee
    integer, intent( out ) :: lmin, lmax
    integer, intent( inout ) :: ierr
    integer, intent( inout ), optional :: itarg
    !
    integer :: targ

#if 0
    targ = 0
    if( present( itarg ) ) then
      if( itarg .gt. 0 .and. itarg .lt. size( FullTable ) ) then
        if( FullTable( itarg )%Zee .eq. zee ) then 
          targ = itarg
          goto 111
        endif
      endif
    endif

    do i = 1, size( FullTable )
      if( FullTable( itarg )%Zee .eq. zee ) then
        targ = i
        goto 111
      endif
    enddo
    if( targ .eq. 0 ) then
      ierr = 1
    else
#endif

    if( present( itarg ) ) then
      if( isRightTarg( zee, itarg ) ) then
        targ = itarg
      else
        targ = getRightTarg( zee )
      endif
    else
      targ = getRightTarg( zee )
    endif
    
    if( targ .lt. 1 ) then
      ierr = 1
      return
    else
      lmin = FullTable( targ )%lMin
      lmax = FullTable( targ )%lMax
      if( present( itarg ) ) itarg = targ
    endif

  end subroutine screen_opf_lbounds

  subroutine screen_opf_nprojForChannel( zee, l, nproj, ierr, itarg )
    integer, intent( in ) :: zee, l
    integer, intent( out ) :: nproj
    integer, intent( inout ) :: ierr
    integer, intent( inout ), optional :: itarg
    !
    integer :: targ

    if( present( itarg ) ) then
      if( isRightTarg( zee, itarg ) ) then
        targ = itarg
      else
        targ = getRightTarg( zee )
      endif
    else
      targ = getRightTarg( zee )
    endif

    if( targ .lt. 1 ) then
      ierr = 1
      return
    else
      if( l .lt. FullTable( targ )%lMin .or. l .gt. FullTable( targ )%lMax ) then
        ierr = 2
        return
      endif
      nproj = FullTable( targ )%nprojPerChannel( l )
      if( present( itarg ) ) itarg = targ
    endif

  end subroutine screen_opf_nprojForChannel

  subroutine screen_opf_getNCutoff( zee, nr, rad, ierr, itarg )
    integer, intent( in ) :: zee
    integer, intent( out ) :: nr
    real(DP), intent( in ) :: rad( : ) 
    integer, intent( inout ) :: ierr
    integer, intent( inout ), optional :: itarg
    !
    integer :: targ, i

    if( present( itarg ) ) then
      if( isRightTarg( zee, itarg ) ) then
        targ = itarg
      else
        targ = getRightTarg( zee )
      endif
    else
      targ = getRightTarg( zee )
    endif

    if( targ .lt. 1 ) then
      ierr = 1
      return
    else
      if( present( itarg ) ) itarg = targ

      nr = -1
      do i = 1, size( rad )
        if( rad( i ) .gt. FullTable( targ )%rMax ) then
          nr = i - 1
          exit
        endif
      enddo

      if( nr .eq. -1 ) ierr = 1
    endif
  end subroutine screen_opf_getNCutoff

  subroutine screen_opf_altinterpProjs(  zee, l, rad, psproj, aeproj, ierr, itarg )
    integer, intent( in ) :: zee, l
    real(DP), intent( in ) :: rad( : )
    real(DP), intent( out ) :: psproj( :, : ), aeproj( :, : )
    integer, intent( inout ) :: ierr
    integer, intent( inout ), optional :: itarg
    !
    integer :: targ, i, p

    if( present( itarg ) ) then
      if( isRightTarg( zee, itarg ) ) then
        targ = itarg
      else
        targ = getRightTarg( zee )
      endif
    else
      targ = getRightTarg( zee )
    endif

    if( targ .lt. 1 ) then
      ierr = 1
      return
    endif

    if( present( itarg ) ) itarg = targ

    if( size( psProj, 2 ) .lt. FullTable( targ )%nprojPerChannel( l ) ) then
      ierr = 2
      return
    endif

    do p = 1, FullTable( targ )%nprojPerChannel( l )
      do i = 1, size( psProj, 1 )
        call intval( FullTable( targ )%nrad, FullTable( targ )%rad, FullTable( targ )%psProj( :, p, l ), &
                     rad( i ), psProj( i, p ), 'err', 'err' )
        call intval( FullTable( targ )%nrad, FullTable( targ )%rad, FullTable( targ )%aeProj( :, p, l ), &
                     rad( i ), aeProj( i, p ), 'err', 'err' )
      enddo
    enddo

  end subroutine screen_opf_altinterpProjs


  subroutine screen_opf_interpProjs(  zee, l, rad, psproj, diffproj, ierr, itarg )
    use OCEAN_mpi, only : myid, root
    integer, intent( in ) :: zee, l
    real(DP), intent( in ) :: rad( : )
    real(DP), intent( out ) :: psproj( :, : ), diffproj( :, : )
    integer, intent( inout ) :: ierr
    integer, intent( inout ), optional :: itarg
    !
    integer :: targ, i, p
    character(len=20 ) :: filnam
    character(len=100) :: formatting

    if( present( itarg ) ) then
      if( isRightTarg( zee, itarg ) ) then
        targ = itarg
      else
        targ = getRightTarg( zee )
      endif
    else
      targ = getRightTarg( zee )
    endif

    if( targ .lt. 1 ) then
      ierr = 1
      return
    endif

    if( present( itarg ) ) itarg = targ

    if( size( psProj, 2 ) .lt. FullTable( targ )%nprojPerChannel( l ) ) then
      ierr = 2
      return
    endif

    do p = 1, FullTable( targ )%nprojPerChannel( l )
      do i = 1, size( psProj, 1 )
        call intval( FullTable( targ )%nrad, FullTable( targ )%rad, FullTable( targ )%psProj( :, p, l ), &
                     rad( i ), psProj( i, p ), 'err', 'err' )
        call intval( FullTable( targ )%nrad, FullTable( targ )%rad, FullTable( targ )%aeProj( :, p, l ), &
                     rad( i ), diffProj( i, p ), 'err', 'err' )
        diffProj( i, p ) = diffProj( i, p ) - psProj( i, p )
      enddo
    enddo

    if( myid .eq. root ) then
      write(filnam, '(A,I2.2,I1.1)' ) 'test1.', zee, l
      write(formatting, '("("I0"(F20.10))")' ) FullTable( targ )%nprojPerChannel( l )+1
      open(unit=99,file=filnam)
      do i = 1, size( psProj, 1 )
        write( 99, formatting ) rad(i), psProj( i, : )
      enddo
      close( 99 )

      write(filnam, '(A,I2.2,I1.1)' ) 'test2.', zee, l
      write(formatting, '("("I0"(F20.10))")' ) FullTable( targ )%nprojPerChannel( l )+1
      open(unit=99,file=filnam)
      do i = 1, size( psProj, 1 )
        write( 99, formatting ) rad(i), diffProj( i, : )
      enddo
      close( 99 )

    endif

  end subroutine screen_opf_interpProjs

  subroutine screen_opf_loadAll( ierr )
    use OCEAN_mpi, only : myid, root, comm, MPI_INTEGER
    integer, intent( inout ) :: ierr

    integer :: maxUnique, Zee, i

    if( myid .eq. root ) then
      open(unit=99,file='zeelist',form='formatted',status='old')
      read(99,*) maxUnique
    endif
#ifdef MPI
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( maxUnique, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    call screen_opf_init( maxUnique, ierr )
    if( ierr .ne. 0 ) return

    do i = 1, maxUnique
      if( myid .eq. root ) then
        read(99,*) Zee
      endif
#ifdef MPI
      call MPI_BARRIER( comm, ierr )
      call MPI_BCAST( Zee, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return
#endif

      call screen_opf_makeNew( Zee, ierr )
      if( ierr .ne. 0 ) return
    enddo

    if( myid .eq. root ) close( 99 )

  end subroutine screen_opf_loadAll



  subroutine screen_opf_init( maxUnique, ierr )
    integer, intent( in ) :: maxUnique
    integer, intent( inout ) :: ierr

    integer :: i

    allocate( FullTable( maxUnique ), stat=ierr )
    do i = 1, maxUnique
      FullTable( i )%Zee = 0
    enddo

  end subroutine screen_opf_init

  subroutine screen_opf_clean( )
    integer :: i

    if( allocated( FullTable ) ) then
      do i = 1, size( FullTable )
        if( allocated( FullTable( i )%rad ) ) deallocate( FullTable( i )%rad ) 
        if( allocated( FullTable( i )%aeProj ) ) deallocate( FullTable( i )%aeProj ) 
        if( allocated( FullTable( i )%psProj ) ) deallocate( FullTable( i )%psProj ) 
      enddo
    
      deallocate( FullTable )
    endif
  end subroutine screen_opf_clean


  subroutine screen_opf_makeNew( Zee, ierr )
    integer, intent( in ) :: Zee
    integer, intent( inout ) :: ierr

    integer :: i, targ

    targ = 0
    if( allocated( FullTable ) ) then
      do i = 1, size( FullTable )
        ! Skip if loaded
        if( FullTable( i )%Zee .eq. Zee ) return
        if( FullTable( i )%Zee .eq. 0 ) then
          targ = i
          exit
        endif
      enddo
    else
      ierr = 1
      return
    endif

    ! if targ is zero we've run out of space!
    if( targ .eq. 0 ) then
      write(6,*) "Logic error in screen_opf. Too few unique Z's were passed in for initialization"
      ierr = 2
      return
    endif
      
    call screen_opf_loadNew( Zee, FullTable( targ ), ierr )

!    call screen_opf_makeAmat( FullTable( targ ), ierr )

  end subroutine screen_opf_makeNew

#if 0
  subroutine screen_opf_makeAmat(  oh, ierr )
    use OCEAN_mpi, only : myid, root
    type( opf_holder ), intent( inout ) :: oh
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tMat( :, : ), r(:), dr(:), r2(:)
    real(DP) :: su, rmin, rmax, rat, dl, xrat, xr1

    integer :: l, i, j, n, ir, nr
  
    real(DP), parameter :: rmifac = 0.00000001d0
    real(DP), parameter :: rmafac = 800.d0

    nr = 2000
    allocate( r( nr ), dr( nr ), r2( nr ) )
    rmin = rmifac / dble( oh%Zee )
    rmax = rmafac / sqrt(dble( oh%Zee ))
    rat = rmax/rmin
    dl=dlog(rat)/dble(nr)
    xrat=dexp(dl)
    xr1=dsqrt(xrat)-dsqrt(1.d0/xrat)
    do i=1,nr
      r(i)=rmin*xrat**dble(i)
      dr(i)=r(i)*xr1
      ! r2 is the measure (r^2 dr) not r^2
      r2(i)=r(i)*r(i)*r(i)*xr1
    end do

    do i = 1, oh%nrad
      write(2000,*) r(i), oh%rad( i )
    enddo

    do l = oh%lMin, oh%lMax
      n = oh%nprojPerChannel( l )
      allocate( tMat( n, n ) )

      do i = 1, n
        do j = 1, i

          su = 0.0_DP
          do ir = 1, oh%nRad
            su = su + r2( ir ) * oh%psProj( ir, i, l ) * oh%psProj( ir, j, l )
          enddo
          tmat( i, j ) = su
          tmat( j, i ) = su

        enddo
      enddo

      if( myid .eq. root ) then
        write( 6, * ) "*** l = ", l
        do i = 1, n
          write(6,*) tmat( :, i )
        enddo
      endif

      call rinvert( n, tmat, ierr )
      if( ierr .ne. 0 ) return
      oh%aMat( 1 : n, 1 : n, l ) = tmat( :, : )
      deallocate( tmat )
    end do

    deallocate( r, dr, r2 )

  end subroutine  screen_opf_makeAmat
#endif


  subroutine screen_opf_loadNew( Zee, oh, ierr )
    use OCEAN_mpi, only : myid, comm, root, MPI_INTEGER, MPI_DOUBLE_PRECISION
    integer, intent( in ) :: Zee
    type( opf_holder ), intent( inout ) :: oh
    integer, intent( inout ) :: ierr

    integer :: ierr_, temp( 3 )

    oh%Zee = Zee

    if( myid .eq. root ) then
      call screen_opf_readPrjFile( oh, ierr )
      if( ierr .ne. 0 ) goto 111
      call screen_opf_readRad( oh, ierr )
      if( ierr .ne. 0 ) goto 111
      temp(1) = oh%lMin
      temp(2) = oh%lMax
      temp(3) = oh%nRad
    endif

111   continue

#ifdef MPI
    call MPI_BARRIER( comm, ierr_ )
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to sync ierr in screen_sites_count_sitelist'
      ierr = ierr_
      return
    endif
    call MPI_BCAST( temp, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( oh%rMax, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return

    if( myid .ne. root ) then 
      oh%lMin = temp(1)
      oh%lMax = temp(2)
      oh%nRad = temp(3)
      allocate( oh%nprojPerChannel( oh%lMin : oh%lMax ) )
    endif

    temp(1) = oh%lMax - oh%lMin + 1
    call MPI_BCAST( oh%nprojPerChannel, temp(1), MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif


    call screen_opf_allocate( oh, ierr )

    if( myid .eq. root ) then
      call screen_opf_readProjectors( oh, ierr )
    endif

#ifdef MPI
    call MPI_BCAST( oh%rad, oh%nrad, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return

    temp(1) = maxval( oh%nprojPerChannel )
    temp(1) = temp(1)*oh%nrad*(oh%lMax - oh%lMin + 1)
    call MPI_BCAST( oh%aeProj, temp(1), MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( oh%psProj, temp(1), MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine screen_opf_loadNew


  subroutine screen_opf_readProjectors( oh, ierr )
    type( opf_holder ), intent( inout ) :: oh
    integer, intent( inout ) :: ierr

    integer :: l, i
    character( len=16 ) :: filnam
    logical :: ex
    real(DP), allocatable :: readbuf(:)

    do l = oh%lMin, oh%lMax

      write( filnam, '(A9,A2,I1.1,A1,I3.3)' ) 'zpawinfo/', 'ae', l, 'z', oh%Zee
      inquire( file=filnam, exist=ex )
      if( ex .eqv. .false. ) then
        write( 6, * ) 'Could not find file: ', filnam
        ierr = 2
        return
      endif

      allocate( readbuf( oh%nprojPerChannel( l ) ) )
      open( unit=99, file=filnam, form='formatted', status='old' )
      do i = 1, oh%nRad
        read( 99, * ) oh%rad( i ), readbuf( : )
        oh%aeProj( i, 1 : oh%nprojPerChannel( l ), l ) = readbuf( : )
!        read( 99, * ) oh%rad( i ), oh%aeProj( i, 1 : oh%nprojPerChannel( l ), l )
      enddo
      close( 99 )

      write( filnam, '(A9,A2,I1.1,A1,I3.3)' ) 'zpawinfo/', 'ps', l, 'z', oh%Zee
      inquire( file=filnam, exist=ex )
      if( ex .eqv. .false. ) then
        write( 6, * ) 'Could not find file: ', filnam
        ierr = 2
        return
      endif

      open( unit=99, file=filnam, form='formatted', status='old' )
      do i = 1, oh%nRad
        read( 99, * ) oh%rad( i ), readbuf( : )
        oh%psProj( i, 1 : oh%nprojPerChannel( l ), l ) = readbuf( : )
!        read( 99, * ) oh%rad( i ), oh%psProj( i, 1 : oh%nprojPerChannel( l ), l )
      enddo
      close( 99 )

      deallocate( readbuf ) 

    enddo


  end subroutine screen_opf_readProjectors


  subroutine screen_opf_allocate( oh, ierr )
    type( opf_holder ), intent( inout ) :: oh
    integer, intent( inout ) :: ierr

    integer :: i

    i = maxval( oh%nprojPerChannel( : ) )
    allocate( oh%aeProj( oh%nRad, i, oh%lMin : oh%lMax ), &
              oh%psProj( oh%nRad, i, oh%lMin : oh%lMax ), &
              oh%rad( oh%nRad ), STAT=ierr )
!              oh%aMat( i, i, oh%lMin : oh%lMax ), oh%rad( oh%nRad ), STAT=ierr )
    if( ierr .ne. 0 ) return

  end subroutine screen_opf_allocate


  subroutine screen_opf_readRad( oh, ierr )
    type( opf_holder ), intent( inout ) :: oh
    integer, intent( inout ) :: ierr

    integer :: dumi
    character( len=20 ) :: filnam
    logical :: ex


    if( oh%Zee .lt. 1 .or. oh%Zee .gt. 109 ) then
      ierr = 1
      return
    endif

    write(filnam, '(A17,I3.3)' ) 'zpawinfo/radfilez', oh%Zee
    inquire( file=filnam, exist=ex )
    if( ex .eqv. .false. ) then
      write( 6, * ) 'Could not find file: ', filnam
      ierr = 2
      return
    endif

    open( unit=99, file=filnam, form='formatted', status='old' )
    read( 99, *  ) oh%rMax, dumi, oh%nRad
    close( 99 )

  end subroutine screen_opf_readRad


  subroutine screen_opf_readPrjFile( oh, ierr )
    type( opf_holder ), intent( inout ) :: oh
    integer, intent( inout ) :: ierr

    integer :: i
    character( len=20 ) :: filnam
    logical :: ex

    if( oh%Zee .lt. 1 .or. oh%Zee .gt. 109 ) then
      ierr = 1
      return
    endif

    write(filnam, '(A17,I3.3)' ) 'zpawinfo/prjfilez', oh%Zee
    inquire( file=filnam, exist=ex )
    if( ex .eqv. .false. ) then
      write( 6, * ) 'Could not find file: ', filnam
      ierr = 2
      return
    endif

    open( unit=99, file=filnam, form='formatted', status='old' )
    read( 99, *  ) oh%lMin, oh%lMax
    if( oh%lMax .lt. oh%lMin .or. oh%lMin .lt. 0 ) then
      write( 6, * ) 'Bad first line of ', filnam
      ierr = 3
      close( 99 ) 
      return
    endif

    allocate( oh%nprojPerChannel( oh%lMin : oh%lMax ) )
    do i = oh%lMin, oh%lMax
      read( 99, * ) oh%nprojPerChannel( i ) 
    enddo
    close( 99 )

  end subroutine
    
  subroutine rinvert( n, smat, ierr )
    integer, intent( in ) :: n
    real( DP ), intent( inout) :: smat( n, n )
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, k, ii
    real( DP ), allocatable :: suse( :, : )
    real( DP ) :: ratio, swap, rc1, rc2
    !
    ! copy smat, reset smax to identity to become inverse
    !
    allocate( suse( n, n ) )
    suse = smat
    smat = 0
    do i = 1, n
       smat( i, i ) = 1
    end do
    !
    ! do inversion by pivoted Gaussian elimination
    !
    do i = 1, n
       ii = i
       do j = i + 1, n
          rc1 = dabs( suse( j, i ) )
          rc2 = dabs( suse( ii, i ) )
          if ( rc1 .gt. rc2 ) ii = j
       end do
       if ( ii .gt. i ) then
          do j = i, n
             swap=suse( i, j )
             suse( i, j ) = suse( ii, j )
             suse( ii, j ) = swap
          end do
          do j = 1, n
             swap = smat( i, j )
             smat( i, j ) = smat( ii, j )
             smat( ii, j ) = swap
          end do
       end if
       if ( suse( i, i ) .eq. 0.0d0 ) then
          write ( 6, * ) 'ZERO DETERMINANT...'
          ierr = -i
          return
       end if
       do j = 1, n
          if ( j .ne. i ) then
             ratio = - suse( j, i ) / suse( i, i )
          else
             ratio = 1.d0 / suse( i, i ) - 1.0d0
          endif
          do k = i, n
             suse( j, k ) = suse( j, k ) + ratio * suse( i, k )
          end do
          do k = 1, n
             smat( j, k ) = smat( j, k ) + ratio * smat( i, k )
          end do
       end do
    end do
    !
    return
  end subroutine rinvert

  subroutine intval( n, xtab, ytab, x, y, lopt, hopt )
    implicit none
    !
    integer, intent( in ) :: n
    real(DP), intent( in) :: xtab( n ), ytab( n ), x
    real(DP), intent( out ) :: y
    character(len=3), intent( in ) :: lopt, hopt
    !
    integer :: ii, il, ih
    real(DP) :: rat
    logical :: below, above, interp
    !
    below = ( x .lt. xtab( 1 ) )
    above = ( x .gt. xtab( n ) )
    if ( below .or. above ) then
       interp = .false.
       if ( below ) then
          select case( lopt )
          case( 'ext' )
             ii = 1
             interp = .true.
          case( 'cap' )
             y = ytab( 1 )
          case( 'err' )
             stop 'error ... we are below!'
          end select
       else
          select case( hopt )
          case( 'ext' )
             ii = n - 1
             interp = .true.
          case( 'cap' )
             y = ytab( n )
          case( 'err' )
             stop 'error ... we are above!'
          end select
       end if
    else
       interp = .true.
       il = 1
       ih = n - 1
       do while ( il + 3 .lt. ih )
          ii = ( il + ih ) / 2
          if ( xtab( ii ) .gt. x ) then
             ih = ii - 1
          else
             il = ii
          end if
       end do
       ii = il
       do while ( xtab( ii + 1 ) .lt. x )
          ii = ii + 1
       end do
    end if
    if ( interp ) then
       rat = ( x - xtab( ii ) ) / ( xtab( ii + 1 ) - xtab( ii ) )
       y = ytab( ii ) + rat * ( ytab( ii + 1 ) - ytab( ii ) )
    end if
    !
    return
  end subroutine intval


end module screen_opf
