module screen_centralPotential
  use ai_kinds, only : DP

  implicit none
  private
  save

  type potential
    real(DP), allocatable :: pot(:)
    real(DP), allocatable :: rad(:)
    integer :: z, n, l
  end type potential

  public :: potential
  public :: screen_centralPotential_newScreenShell, screen_centralPotential_loadAll
  public :: screen_centralPotential_prepAll

  contains

  subroutine screen_centralPotential_loadAll( znlPot, znl, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    type( potential ), intent( out ) :: znlPot( : )
    integer, intent( in ) :: znl(:,:)
    integer, intent( inout ) :: ierr

    integer :: i, Nznl, ierr_

    Nznl = size( znlPot, 1 )
    ierr_ = 0

    do i = 1, Nznl
      if( myid .eq. root ) then
        call screen_centralPotential_load( znl(1,i), znl(2,i), znl(3,i), znlPot(i), ierr )
      endif
#ifdef MPI
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
#endif
      if( ierr .ne. 0 ) return
      if( ierr_ .ne. 0 ) then
        ierr = ierr_
        return
      endif
      call screen_centralPotential_share( znlPot(i), myid, root, comm, ierr )
      if( ierr .ne. 0 ) return
    enddo
  end subroutine screen_centralPotential_loadAll

  subroutine screen_centralPotential_share( pot, myid_, root_, comm_, ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_DOUBLE_PRECISION
    type( potential ), intent( inout ) :: pot
    integer, intent( in ) :: myid_, root_, comm_
    integer, intent( inout ) :: ierr
    !
    integer :: intArray( 4 )

    if( myid_ .eq. root_ ) then
      intArray(4) = size( pot%pot, 1 )
      intArray(1) = pot%z
      intArray(2) = pot%n
      intArray(3) = pot%l
    endif

    call MPI_BCAST( intArray, 4, MPI_INTEGER, root_, comm_, ierr )
    if( ierr .ne. 0 ) return
    if( myid_ .ne. root_ ) then
      allocate( pot%pot( intArray(4) ), pot%rad( intArray(4) ) )
    endif
    call MPI_BCAST( pot%pot, intArray(4), MPI_DOUBLE_PRECISION, root_, comm_, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( pot%rad, intArray(4), MPI_DOUBLE_PRECISION, root_, comm_, ierr )
    if( ierr .ne. 0 ) return

  end subroutine screen_centralPotential_share

      

  subroutine screen_centralPotential_prepAll( znl, n, ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER
    integer, intent( out ) :: znl(:,:)
    integer, intent( out ) :: n
    integer, intent( inout ) :: ierr

    integer :: i, maxSize, j

    if( myid .ne. root ) then
#ifdef MPI
      call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif
      return
    endif

    maxSize = size( znl, 2 )
    open( unit=99, file='edgelist', form='formatted', status='old' )

    do i = 1, maxSize
      read(99,*,iostat=j) znl(1:3, i )
      if( j .lt. 0 ) then
        n = i - 1
        close( 99 )
#ifdef MPI
        call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif
        return
      elseif( j .gt. 0 ) then
        ierr = j
        return
      endif
    enddo
    close( 99 )

    ! if this isn't here the other mpi will hang above ...
#ifdef MPI
        call MPI_BCAST( n, 1, MPI_INTEGER, root, comm, ierr )
#endif

    write(6,'(A,I8)') 'Too many unique edges! Greater than: ', maxSize
    ierr = 1
  end subroutine 

  subroutine screen_centralPotential_newScreenShell( pot, newPot, rad, ierr )
    type( potential ), intent( in ) :: pot
    type( potential ), intent( out ) :: newPot
    real(DP), intent( in ) :: rad
    integer, intent( inout ) :: ierr

    integer :: i, restart
    real(DP) :: invRad

    if( ( .not. allocated( pot%pot ) ) .or. ( .not. allocated( pot%rad ) ) ) then
      ierr = 5
      return
    endif

    allocate( newPot%pot( size( pot%pot ) ), newPot%rad( size( pot%rad ) ) )
    newPot%rad(:) = pot%rad(:)
    
    invRad = 1.0_DP / rad
    do i = 1, size( pot%pot ) 
      if( newPot%rad( i ) .lt. rad ) then
        newPot%pot( i ) = pot%pot( i ) + invRad
      else
        restart = i
        exit
      endif
    enddo

    newPot%pot( restart : size( newPot%pot ) ) = 0.0_DP
    newPot%z = pot%z
    newPot%n = pot%n
    newPot%l = pot%l

  end subroutine screen_centralPotential_newScreenShell

  subroutine screen_centralPotential_load( z, n, l, pot, ierr )
    integer, intent( in ) :: z, n, l
    type( potential ), intent( out ) :: pot
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tmpPot(:), tmpRad(:), tmp(:)
    integer :: maxLength, curLength, fh

    character(len=26) :: fileName

    pot%z = z
    pot%n = n
    pot%l = l
    
    write(fileName,'(A17,I3.3,A1,I2.2,A1,I2.2)') 'zpawinfo/vc_barez', z, 'n', n, 'l'
    fh = 99
    open( unit=fh, file=fileName, form='formatted', status='unknown' )
    rewind( 99 )

    maxLength = 1000
    curLength = 1

    allocate( tmpPot( maxLength ), tmpRad( maxLength ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do while( maxLength .le. 512000 )

      call doRead( fh, curLength, MaxLength, tmpRad, tmpPot, ierr )

      if( curLength .lt. 0 ) goto 11

      allocate( tmp( maxLength ), stat=ierr )
      if( ierr .ne. 0 ) return
      tmp(:) = tmpPot(:)
      deallocate( tmpPot )
      allocate( tmpPot( maxLength*2 ), stat=ierr )
      if( ierr .ne. 0 ) return
      tmpPot(1:maxLength) = tmp(:)
      tmp(:) = tmpRad(:)
      deallocate( tmpRad )
      allocate( tmpRad( maxLength*2 ), stat=ierr )
      tmpRad(1:maxLength) = tmp(:)

      curLength = maxLength + 1
      maxLength = 2*maxLength
      deallocate( tmp )
    enddo
    !    If we exit the loop then we've given up not reached the end of the file
    ierr = -1
    return

11  continue
    curLength = abs( curLength )
    allocate( pot%rad( curLength ), pot%pot( curLength ), stat=ierr )
    if( ierr .ne. 0 ) return
    pot%rad(:) = tmpRad(1:curLength)
    pot%pot(:) = tmpPot(1:curLength)
    deallocate( tmpPot, tmpRad )

    close( fh )

  end subroutine screen_centralPotential_load

  subroutine doRead( fh, curLength, maxLength, rad, pot, ierr )
    integer, intent( in ) :: fh, maxLength
    integer, intent( inout ) :: curLength
    real(DP), intent( inout ) :: rad(:), pot(:)
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, start

    start = curLength

    do i = start, maxLength
      read(fh,*,iostat=j) rad(i), pot(i)
      if( j .lt. 0 ) then
        curLength = -( i - 1 )
        return
      elseif( j .gt. 0 ) then
        ierr = j
        return
      endif
      ! if j .eq. 0 then keep going!
    enddo

  end subroutine doRead

end module screen_centralPotential
