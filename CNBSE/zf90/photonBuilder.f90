program photonBuilder
  implicit none

  integer, parameter :: DP = selected_real_kind(15, 307)
  integer, parameter :: OH4(4,6) = reshape( (/ 1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 0, 1, 0, -1, 0, 1, &
        0, 0, 1, 1, 0, 0, -1, 1 /), shape(OH4) )
  integer, parameter :: OH6(4,18  )  = reshape( (/ 1, 0, 0, 1, -1, 0, 0, 1, 0, 1, 0, 1, 0, -1, 0, 1, &
        0, 0, 1, 1, 0, 0, -1, 1, 1, 1, 0, 2, 1, -1, 0, 2, -1, 1, 0, 2, -1, -1, 0, 2, 1, 0, 1, 2, &
        1, 0, -1, 2, -1, 0, 1, 2, -1, 0, -1, 2, 0, 1, 1, 2, 0, 1, -1, 2, 0, -1, 1, 2, 0, -1, -1, 2 /), shape(OH6 ) )
  integer, parameter :: OH8(4, 26 ) = reshape( (/ 1, 1, 1, 27, 1, 1, -1, 27, 1, -1, 1, 27, 1, -1, -1, 27 , &
        -1, 1, 1, 27, -1, 1, -1, 27, -1, -1, 1, 27, -1, -1, -1, 27, 1, 1, 0, 32, 1, -1, 0, 32, &
        -1, 1, 0, 32, -1, -1, 0, 32, 1, 0, 1, 32, 1, 0, -1, 32, -1, 0, 1, 32, -1, 0, -1, 32, &
        0, 1, 1, 32, 0, 1, -1, 32, 0, -1, 1, 32, 0, -1, -1, 32, 1, 0, 0, 40, -1, 0, 0, 40, &
        0, 1, 0, 40, 0, -1, 0, 40, 0, 0, 1, 40, 0, 0, -1, 40 /), shape( OH8 ) )

  integer, allocatable :: OH(:,:)
  integer, allocatable :: photon1( :, : ), sym(:,:,:), photon2(:,:), iphoton(:,:)
  real(DP), allocatable :: rphoton(:,:), rphoton2(:,:)
  real(DP) :: rphot(3), su, rres(3)
  integer :: i, norm, ifound, j, k, res(3), nsym, iphot, res2(3), jfound, l
  logical :: isUnique

  integer :: oorder
  character(len=16) :: photonType, outname
  real(DP) :: energy, qlen
  logical :: twoStep, PolFirst, printQ

  open(unit=99,file='photon_control.ipt', form='formatted', status='old' )
  read(99,*) photonType
  read(99,*) energy
  read(99,*) qlen
  read(99,*) oorder  ! order can be 4, 6, or 8 (only 8 right now/is ignored)
  close(99)

  select case( photonType )
  case( 'dipole' )
    twoStep = .false.
    PolFirst = .true.
    printQ = .false.
  case( 'quad' )
    twoStep = .true.
    PolFirst = .true.
    printQ = .false.
  case( 'quadalone' )
    twoStep = .false.
    PolFirst = .true.
    printQ = .false.
  case( 'NRIXS', 'qRaman' )
    twoStep = .false.
    PolFirst = .true.
    printQ = .true.
  case default
    write( 6, * ) 'Un-supported photon type: ', photonType
  end select

  


  select case( oorder )
    case( 8 )
      ! Here we could swap OH8, OH6, OH4, ...
      allocate( OH( 4, size( OH8, 2 ) ), photon1( 4, size( OH8, 2 ) ) )
      OH(:,:) = OH8(:,:)
    case( 6 ) 
      allocate( OH( 4, size( OH6, 2 ) ), photon1( 4, size( OH6, 2 ) ) )
      OH(:,:) = OH6(:,:)
    case( 4 )
      allocate( OH( 4, size( OH4, 2 ) ), photon1( 4, size( OH4, 2 ) ) )
      OH(:,:) = OH4(:,:)
    case default
      stop
  end select

  norm = 0
  do i = 1, size( OH, 2 ) 
    norm = norm + OH( 4, i )
  enddo

  write( 6, * ) norm, size( OH8, 2 ) 

  if( .false. ) then
  allocate( sym( 3, 3, 3 ) )
  nsym = 3
  sym(:,:,:) = 0
  sym(1,1,1) = -1
  sym(2,2,1) = -1
  sym(3,3,1) = -1

  sym(1,2,2) = 1
  sym(2,1,2) = 1
  sym(3,3,2) = 1

  sym(1,2,3) = -1
  sym(2,1,3) = -1
  sym(3,3,3) = -1
  else
  open(unit=99,file='sym.txt',form='formatted',status='old')
  read(99,*) nsym
  allocate( sym( 3, 3, nsym ) )
  do i = 1, nsym
    do j = 1, 3
      read( 99, * ) sym(:,j,i)
    enddo 
  enddo
  endif

  ifound = 0
  do i = 1, size( OH, 2 )
    isUnique = .true.

    do j = 1, ifound
      do k = 1, nsym
        res(:) = matmul( sym(:,:,k), photon1(1:3,j) )
        if( res(1) .eq. OH(1,i) .and. res(2) .eq. OH(2,i) .and. res(3) .eq. OH(3,i) ) then
          photon1(4,j) = photon1(4,j) + OH(4,i )
          isUnique = .false.
          goto 100
        endif
      enddo
    enddo
100 continue
    if( isUnique ) then
      ifound = ifound + 1
      photon1( :, ifound ) = OH( :, i )
    endif
  enddo
  
  do i = 1, ifound
    write( 6, '(3I3,I5)' ) photon1(:,i)
  enddo

  iphot = 0

  allocate( photon2( 4, 4 ), rphoton( 3, 4 ), iphoton(3,4), rphoton2(3,4) )
  if( twoStep ) norm = norm*4

  do i = 1, ifound 
    if( twoStep ) then
      res2 = (/ 1, 0, 0 /)
      call cross( photon1(1:3,i), res2, res )
      if( res(1) .eq. 0 .and. res(2) .eq. 0 .and. res(3) .eq. 0 ) then
        res2 = (/ 0, 1, 0 /)
        call cross( photon1(1:3,i), res2, res )
      endif

      call cross(  photon1(1:3,i), res, res2 )

      iphoton(:,1) = res(:)
      rphoton(:,1) = real(res(:),DP)
      iphoton(:,2) = res2(:)
      rphoton(:,2) = real(res2(:),DP)
      iphoton(:,3) = res(:) + res2(:)
      rphoton(:,3) = rphoton(:,1) + rphoton(:,2)
      iphoton(:,4) = res(:) - res2(:)
      rphoton(:,4) = rphoton(:,1) - rphoton(:,2)

      do j = 1, 4
        su = dot_product( rphoton(:,j), rphoton(:,j) )
        rphoton(:,j) = rphoton(:,j) / sqrt( su )
      enddo
        

      jfound = 0
      do j = 1, 4
        isUnique = .true.
        do k = 1, jfound
          do l = 1, nsym
            rres(:) = matmul( sym(:,:,l), rphoton2(:,k) )
            if( abs( rres(1) - rphoton(1,j) ) .lt. 0.001 .and. &
                abs( rres(2) - rphoton(2,j) ) .lt. 0.001 .and. &
                abs( rres(3) - rphoton(3,j) ) .lt. 0.001 ) then
              isUnique = .false.
              photon2( 4, k ) = photon2( 4, k ) + photon1(4,i)
              goto 12
            endif
          enddo
        enddo
12      continue
        if( isUnique ) then
          jfound = jfound + 1
          photon2( 1:3, jfound ) = iphoton( :, j )
          photon2( 4, jfound ) = photon1(4,i)
          rphoton2(:,jfound) = rphoton(:,j)
        endif
      enddo

      do j = 1, jfound
        iphot = iphot + 1
        write( outname, '(A,I0)' ) 'photon', iphot
        open( unit=99, file=outname, form='formatted' )
        write( 99, '(A)' ) photonType
        write( 99, '(A,3(I5))' ) 'cartesian', photon1(1:3,i)
        write( 99, '(A)' ) 'end'
        write( 99, '(A,3(I5))' ) 'cartesian', photon2(1:3,j)
        write( 99, '(A)' ) 'end'
        write( 99, * ) energy
        write(99, * ) real( photon2(4, j ), DP )/ real(norm,dp), photon1(4,i), norm
        close( 99 )
      enddo
    
    else
      iphot = iphot + 1
      write( outname, '(A,I0)' ) 'photon', iphot
      open( unit=99, file=outname, form='formatted' )
      write( 99, '(A)' ) photonType
      write( 99, '(A,3(I5))' ) 'cartesian', photon1(1:3,i)
      write( 99, '(A)' ) 'end'
      write( 99, '(A,3(I5))' ) 'cartesian', photon1(1:3,i)
      write( 99, '(A)' ) 'end'
      write( 99, * ) energy
      write(99, * ) real( photon1(4, i ), DP )/ real(norm,dp), photon1(4,i), norm
      close( 99 )

    endif
  enddo

111   continue
  contains

  subroutine cross( xx, yy, zz )
    integer, intent( in ) :: xx(3), yy(3)
    integer, intent( out ) :: zz(3)

    zz(1) = xx(2)*yy(3) - xx(3)*yy(2)
    zz(2) = yy(1)*xx(3) - xx(1)*yy(3)
    zz(3) = xx(1)*yy(2) - xx(2)*yy(1)
  end subroutine cross

end program photonBuilder
