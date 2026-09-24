program coreScreen
  use ai_kinds, only : DP
  implicit none
  


  integer :: ZZ, nc, lc, nptot, ntot, nspin, lmin, lmax, l, dumi, npt, maxll
  integer :: kgl, kgh, l1, l2, l3, l4, m1, m2, m3, m4, mk, k, ip, jp, istart
  integer :: i, j, nu2, nu4, kmesh(3), nk, isite, iisite, nsite
  integer, allocatable :: nproj(:)
  real(DP) :: dumr, pi, su, yp( 0 : 1000 ), avecs(3,3), omega
  real(DP), allocatable :: x(:), w(:), gk(:,:,:)
  complex(DP) :: f1, f2, f3
  complex(DP), allocatable :: cks(:,:,:), denMat(:,:), dmat(:,:)
  logical, parameter :: yes = .true.
  logical, parameter :: no = .false.
  character(len=24) :: str
  character(len=10) :: add10
  character(len=2) :: el = 'N_'

  real(DP) :: rc
  integer :: nrmax, nr, ir, iptot, m, nu, nang, nspecpnt, il, nspcpnt, iang
  real(DP), allocatable :: den( :,: ), proj( :, : ), aeproj( :, : ), ang( :, : )
  real(DP), allocatable :: psproj( :, : ), prefs(:)
  complex(DP), allocatable :: ylm( :, : ), psi( :, :, :, : ), psiae( :, :, :, : )


  nspcpnt = 5
  write(str,'(A8,I0)') 'specpnt.', nspcpnt
  open( unit=99,file=str, form='formatted',status='old')
  read(99,*) nang
  allocate( ang(4,nang) )
  do i = 1, nang
    read(99,*) ang(:,i)
  enddo
  close(99)

  allocate( prefs(0:1000) )
  call getprefs( prefs )

  open(unit=98,file='exx.inp',form='formatted',status='old')

  read(98,*) nsite
  do isite = 1,nsite

    read(98,*) el, ZZ, nc, lc, iisite

    write(str,'(A8,I3.3)') 'prjfilez', ZZ
    write(add10, '(A1,I3.3,A1,I2.2,A1,I2.2)') 'z', ZZ, 'n', nc, 'l', lc

    open(unit=99,file=str,form='formatted',status='old')
    read(99,*) lmin, lmax, dumi, dumr
    allocate( nproj(lmin:lmax) )
    npt = 0
    do l=lmin, lmax
      read(99,*) nproj(l)
      npt = npt + nproj(l)
    enddo
    close(99)

    write(str,'(A8,I3.3)') 'radfilez', ZZ
    open(unit=99,file=str,form='formatted',status='old')
    read(99,*) rc, nrmax, nr
    close( 99 )
    allocate( den( 3, nr ) )
    den(:,:) = 0.0_DP

    allocate( psproj( nr, npt ), aeproj( nr, npt ) )
    allocate( ylm( nang, (lmax+1)**2 ) )
    il = 0
    do l = lmin, lmax
      do m = -l, l
        il = il + 1
        do iang = 1, nang
          call ylmeval( l, m, ang(1,iang), ang(2,iang), ang(3,iang), ylm(iang,il), prefs )
        enddo
      enddo
    enddo

    npt = 1
    do l = lmin, lmax
      write(str, '(A2,I1,A1,I3.3)') 'ps', l, 'z', ZZ
      open(unit=99,file=str,form='formatted',status='old')
      do ir = 1, nr
        read(99,*) den(1,ir), psproj(ir, npt:+npt+nproj(l)-1)
      enddo
      close( 99 )
      write(str, '(A2,I1,A1,I3.3)') 'ae', l, 'z', ZZ
      open(unit=99,file=str,form='formatted',status='old')
      do ir = 1, nr
        read(99,*) den(1,ir), aeproj(ir,npt:+npt+nproj(l)-1)
      enddo
      close( 99 )
      npt = npt + nproj(l)
    enddo

    write(str,'(A8,A2,I4.4)') 'parcksv.', el, iisite
    open(unit=99,file=str,form='unformatted',access='stream',status='old')
    read(99) nptot, ntot, nspin
    allocate(cks(nptot, ntot, nspin) )
    read(99) cks
    close(99)

    allocate( psi(nang,nr,ntot,nspin), psiae(nang,nr,ntot,nspin) )
    psi(:,:,:,:) = 0.0_DP
    psiae(:,:,:,:) = 0.0_DP

    iptot = 0
    npt = 0
    il = 0
    do l = lmin, lmax
!      allocate( proj( nproj(l), nr ), aeproj( nproj(l), nr ) )
!      write(str, '(A2,I1,A1,I3.3)') 'ps', l, 'z', ZZ
!      open(unit=99,file=str,form='formatted',status='old')
!      do ir = 1, nr
!        read(99,*) den(1,ir), proj(:,ir)
!      enddo
!      close( 99 )

!      write(str, '(A2,I1,A1,I3.3)') 'ae', l, 'z', ZZ
!      open(unit=99,file=str,form='formatted',status='old')
!      do ir = 1, nr
!        read(99,*) den(1,ir), aeproj(:,ir)
!      enddo
!      close( 99 )
      

      do m = -l, l
        il = il + 1
        do nu = 1, nproj(l)
          iptot = iptot + 1
          su = 0.0_DP
          do j = 1, nspin
            do i = 1, ntot
!              su = su + real( cks(iptot, i, j ) * conjg( cks( iptot, i, j ) ), DP )
              do ir = 1, nr
                do iang = 1, nang
                  psi(iang,ir,i,j) = psi(iang,ir,i,j) &
                                   + psproj(ir,npt+nu) * cks( iptot, i, j ) * ylm(iang,il)
                  psiae(iang,ir,i,j) = psiae(iang,ir,i,j) &
                                   + aeproj(ir,npt+nu) * cks( iptot, i, j ) * ylm(iang,il)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      npt = npt + nproj(l)
    enddo

    do j = 1, nspin
      do i = 1, ntot
        do ir = 1, nr
          su = 0.0_DP
          do iang = 1, nang
            su = su + psi(iang,ir,i,j) * conjg( psi(iang,ir,i,j) ) * ang(4,iang)
          enddo
          den(2,ir) = den(2,ir) + su
          su = 0.0_DP
          do iang = 1, nang
            su = su + psiae(iang,ir,i,j) * conjg( psiae(iang,ir,i,j) ) * ang(4,iang)
          enddo
          den(3,ir) = den(3,ir) + su
        enddo
      enddo
    enddo

    write(str,'(A4,I4.4)' ) 'avg.', iisite
    open(unit=99,file=str,form='formatted',status='unknown')
    do ir = 1, nr
      write(99,*) den(:,ir)
    enddo
    close(99)

    deallocate( nproj, den, cks, psproj, aeproj )
    deallocate( psi, psiae, ylm )
  enddo

  deallocate( prefs, ang )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains


  subroutine ylmeval( l, m, x, y, z, ylm, prefs )
    implicit none                                    
    !
    integer, intent( in )  :: l, m                   
    !
    real( DP ), intent( in ) :: x, y, z, prefs( 0 : 1000 )
    complex( DP ), intent( out ) :: ylm
    !  
    integer :: lam, j, mm
    real( DP ) :: r, rinv, xred, yred, zred, f       
    real( DP ) :: u, u2, u3, u4, u5
    complex( DP ) :: rm1
    !
    if ( l .gt. 5 ) stop 'l .gt. 5 not yet allowed'
    !
    r = sqrt( x ** 2 + y ** 2 + z ** 2 )
    if ( r .eq. 0.d0 ) r = 1
    rinv = 1 / r
    xred = x * rinv
    yred = y * rinv
    zred = z * rinv
    !
    u = zred
    u2 = u * u
    u3 = u * u2
    u4 = u * u3
    u5 = u * u4
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    !
    mm = abs( m ) + 0.1
    lam = 10 * mm + l
    !
    select case( lam )
       !
    case( 00 )
       f =   1                                       !00
       !
    case( 11 )
       f = - 1                                       !11
    case( 01 )
       f =   u                                       !10
       !
    case( 22 )
       f =   3                                       !22
    case( 12 )
       f = - 3 * u                                   !21
    case( 02 )
       f =   ( 3 * u2 - 1 ) / 2                      !20
       !
    case( 33 )
       f = - 15                                      !33
    case( 23 )
       f =   15 * u                                  !32
    case( 13 )
       f = - ( 15 * u2 - 3 ) / 2                     !31
    case( 03 )
       f =   ( 5 * u3 - 3 * u ) / 2                  !30
       !
    case( 44 )
       f =   105                                     !44
    case( 34 )
       f = - 105 * u                                 !43
    case( 24 )
       f =   ( 105 * u2 - 15 ) / 2                   !42
    case( 14 )
       f = - ( 35 * u3 - 15 * u ) / 2                !41
    case( 04 )
       f =   ( 35 * u4 - 30 * u2 + 3 ) / 8           !40
       !
    case( 55 )
       f = - 945                                     !55
    case( 45 )
       f =   945 * u                                 !54
    case( 35 )
       f = - ( 945 * u2 - 105 ) / 2                  !53
    case( 25 )
       f =   ( 315 * u3 - 105 * u ) / 2              !52
    case( 15 )
       f = - ( 315 * u4 - 210 * u2 + 15 ) / 8        !51
    case( 05 )
       f =   ( 63 * u5 - 70 * u3 + 15 * u ) / 8      !50
       !
    end select
    !
    ylm = prefs( lam ) * f
    if ( m .gt. 0 ) then
       do j = 1, m
          ylm = ylm * ( xred + rm1 * yred )
       end do
    end if
    if ( m .lt. 0 ) then
       do j = 1, mm
          ylm = - ylm * ( xred - rm1 * yred )
       end do
    end if
    !
    return
  end subroutine ylmeval
  subroutine getprefs( prefs )
    implicit none
    !           
    real( DP ), intent(out) :: prefs( 0 : 1000 )
    !           
    integer l, m, lam, lamold
    real( DP ) :: pi 
    !
    pi = 4.0d0 * atan( 1.0d0 )
    !                                      
    do l = 0, 5
       prefs( l ) = dble( 2 * l + 1 ) / ( 4.0d0 * pi )
       lamold = l 
       do m = 1, l
          lam = 10 * m + l
          prefs( lam ) = prefs( lamold ) / dble( ( l - m + 1 ) * ( l + m ) )
          lamold = lam
       end do
    end do
    !
    do l = 0, 5
       do m = 0, l
          lam = 10 * m + l
          prefs( lam ) = sqrt( prefs( lam ) )
       end do
    end do
    !
    return
  end subroutine getprefs 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program
