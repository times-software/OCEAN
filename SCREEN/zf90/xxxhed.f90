! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program xxxhed
  implicit none
  !
  ! these variables are for non-x-ray calcs.
  integer :: nv, ng, npass, nvnew, i, j, k, nvtape, i1, i2, i3
  real( kind = kind( 1.0d0 ) ) :: bmet( 3, 3 ),bvec( 3, 3 ),avec( 3, 3 ),uu( 3 ),temp
  character * 20 :: com
  integer, allocatable :: kvc( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vr( : , : ),vi( : , : ),t2( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: tr( : , : , : ),ti( : , : , : ),t3( : , : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: kr( : , : , : , : , : ), ki( : , : , : , : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: r2( : , : )
  real( kind = kind( 1.0d0 ) ), allocatable :: cfr( : , : ),cfi( : , : )
  !
  ! these below should be only for x-ray calcs
  !
  integer :: zn( 3 ), step
  !
  read ( 5, * ) nv
  open( unit = 99,file = 'kdotp',form = 'unformatted',status = 'unknown' )
  rewind 99
  read ( 99 ) bmet, bvec,avec,uu,nvtape,ng
  if ( nv .gt. nvtape ) stop 'bad nv'
  if ( nv .le. 0 ) nv = nvtape
  allocate( t3( 0 : 3, nvtape, nvtape ), tr( 0 : 3, nv, nv ), ti( 0 : 3, nv, nv ) )
  read ( 99 ) ( ( ( t3( k, i, j ), k = 0, 3 ), i = 1, nvtape ), j = 1, nvtape )
  do j = 1, nv
     do i = 1, nv
        do k = 0, 3
           tr( k, i, j ) = t3( k, i, j )
        end do
     end do
  end do
  read ( 99 ) ( ( ( t3( k, i, j ), k = 0, 3 ), i = 1, nvtape ), j = 1, nvtape )
  do j = 1, nv
     do i = 1, nv
        do k = 0, 3
           ti( k, i, j ) = t3( k, i, j )
        end do
     end do
  end do
  deallocate( t3 )
  allocate( t2( nvtape, nvtape ), vr( nv, nv ), vi( nv, nv ) )
  read ( 99 ) ( ( t2( i, j ), i = 1, nvtape ), j = 1, nvtape )
  do j = 1, nv
     do i = 1, nv
        vr( i, j ) = t2( i, j )
     end do
  end do
  read ( 99 ) ( ( t2( i, j ), i = 1, nvtape ), j = 1, nvtape )
  do j = 1, nv
     do i = 1, nv
        vi( i, j ) = t2( i, j )
     end do
  end do
  deallocate( t2 )
  allocate( r2( nvtape, nvtape ) )
  allocate( kr( 0 : 4, 0 : 4, 0 : 4, nv, nv ), ki( 0 : 4, 0 : 4, 0 : 4, nv, nv ) )
  do i1 = 0, 4
     do i2 = 0, 4
        do i3 = 0, 4
           read ( 99 ) ( ( r2( i, j ), i = 1, nvtape ), j = 1, nvtape )
           do j = 1, nv
              do i = 1, nv
                 kr( i3, i2, i1, i, j ) = r2( i, j )
              end do
           end do
           read ( 99 ) ( ( r2( i, j ), i = 1, nvtape ), j = 1, nvtape )
           do j = 1, nv
              do i = 1, nv
                 ki( i3, i2, i1, i, j ) = r2( i, j )
              end do
           end do
        end do
     end do
  end do
  deallocate( r2 )
  allocate( kvc( 3, ng ) )
  read ( 99 ) ( ( kvc( j, i ), j = 1, 3 ), i = 1, ng )
  read ( 99 ) npass
  if ( npass .ne. ng ) stop 'npass not ng'
  allocate( cfr( nv, ng ), cfi( nv, ng ) )
  do j = 1, nv
     read ( 99 ) ( cfr( j, i ), cfi( j, i ), i = 1, npass )
  end do
  close ( unit = 99 )
  write ( 6, * ) 'The file has been read from disk.'
  if ( npass.gt.ng ) then
     write ( 6, * ) 'This is a two-component system'
  else
     write ( 6, * ) 'This is a one-component system'
  end if
  !
  do
     read ( 5, * ) com
     write ( 6, * ) ' command = ', com
     if ( com .eq. 'polargrid' ) then
        call realbsutil( nv, ng, cfr, cfi, kvc, bvec, avec )
     end if
     if ( com .eq. 'P' ) then
        call gcdmat( nv, ng, avec, bvec, bmet, uu, kvc, tr, ti, vr, vi, kr, ki, cfr, cfi )
     end if
     if ( com .eq. 'A' ) then
        call gconv2( nv, ng, avec, bvec, bmet, uu, kvc, tr, ti, vr, vi, kr, ki, cfr, cfi )
     end if
     if ( com .eq. 'R' ) then
        call gcunocc( nv, ng, avec, bvec, bmet, uu, kvc, tr, ti, vr, vi, kr, ki, cfr, cfi )
     end if
     if ( com .eq. 'gbuilder' ) then
        call gbuilder( nv, ng, avec, bvec, bmet, uu, kvc, tr, ti, vr, vi, kr, ki, cfr, cfi )
     end if
     if ( com .eq. 'a' ) call amdec( nv, ng, avec, bvec, bmet, uu, npass, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki )
     if ( com .eq. 'b' ) call bands( nv, bmet, bvec, avec, uu, vr, vi, tr, ti, kr, ki )
     if ( com .eq. 'crit' ) call critan( nv, bvec, bmet, avec, vr, vi, tr, ti, kr, ki, uu )
     if ( com .eq. 'B' ) call bbands( nv, bmet, bvec, avec, uu, vr, vi, tr, ti, kr, ki )
     if ( com .eq. 'c' ) then
        do i = 1, ng
           temp = 0.d0
           do j = 1, 3
              do k = 1, 3
                 temp = temp + dble( kvc( j, i ) * kvc( k, i ) ) * bmet( j, k )
              end do
           end do
           write ( 6, '( 2x, 5i6, 1f10.5 )' ) i, ng, ( kvc( j, i ), j = 1, 3 ), temp
        end do
     end if
     if ( com .eq. 'E' ) call evu( ng, kvc, avec, bvec, bmet )
     if ( com .eq. 'g' ) call lgmain( ng, kvc, nv, npass, cfr, cfi, bvec, bmet )
     if ( com .eq. 'f' ) then
        read ( 5, * ) step
        select case( step )
        case( 1 )
           call vget( bmet, bvec, avec, uu, vr, vi, tr, ti, kr, ki, cfr, cfi, kvc, nv, ng )
           ! case( 2 )
           ! call fullvget( bmet, bvec, avec )
        end select
     end if
     if ( com .eq. 'F' ) then
        write ( 6, * ) 'goodfluo disabled pending zone 1 compliance'
        ! call goodfluo( nv, bvec )
     end if
     if ( com .eq. 'Y' ) call recon( nv, ng, avec, bvec, bmet, uu, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki )
     if ( com .eq. 'L' ) call eld
     if ( com .eq. 'M' ) call padanal
     if ( com .eq. 'm' ) call redpad( cfr, cfi, kvc, nv, npass, ng, bvec, avec, bmet, vr, vi, tr, ti, kr, ki, uu )
     if ( com .eq. 'pchek' ) call pchek( avec, bvec, bmet, nv, ng, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki, uu )
     select case( com )
     case( 'D' )
        call gridbsutil( nv, ng, cfr, cfi, kvc, bvec, avec, bmet, tr, ti, vr, vi, uu, kr, ki )
     case( 'Grid' )
        call getkgrid( zn )
     case( 'k' )
        call kripes( nv, npass, ng, kvc, bmet, avec, bvec, uu, vr, vi, tr, ti, kr, ki, cfr, cfi )
     case( 'K' )
         call ksetup( nv, ng, avec, bvec, bmet, uu, zn, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki )
     case( 'Multiplet' )
        deallocate ( kr, ki, tr, ti, vr, vi )
        call multip( nv, avec, bvec, bmet, zn, cfr, cfi, ng, kvc )
        stop
     case( 'Newmult' )
        deallocate ( kr, ki, tr, ti, vr, vi )
        call newmultip( nv, avec, bvec, bmet, zn, cfr, cfi, ng, kvc )
        stop
     case( 'fsassem' )
        call fsassem( nv, avec, bvec, bmet, zn, cfr, cfi, ng, kvc )
     case( 'quit' )
        write ( 6, * ) 'reader terminus ordered'
     case( 'reducebasis' )
        read ( 5, * ) nvnew
        if ( nvnew.gt.nv ) then
           write ( 6, * ) 'nv is only ', nv
           stop
        end if
        nv = nvnew
     case( 's' )
        call stackplot( nv, ng, avec, bvec, bmet, uu, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki )
!    case( 'T' )
!       call test( nv, bmet, bvec, avec, vr, vi, tr, ti, kr, ki, uu, ng, kvc )
     case( 'u' )
        call uget( bmet, bvec, avec, uu, vr, vi, tr, ti, kr, ki, cfr, cfi, kvc, nv, ng )
     case( 'U' )
        call newuget( bmet, bvec, avec, uu, vr, vi, tr, ti, kr, ki, cfr, cfi, kvc, nv, ng )
     case( 'v' )
        call wvvu( nv, ng, bmet, avec, bvec, kvc, vr, vi, tr, ti, kr, ki, cfr, cfi, uu )
     case( 'w' )
        call psivuf( nv, ng, avec, bvec, bmet, uu, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki )
     case( 'x' )
        call ximat( bmet, avec, bvec, kvc, ng, nv, npass, vr, vi, tr, ti, kr, ki, cfr, cfi, uu )
     case( 'y' )
        call epsopt( avec, npass, kvc, bmet )
     case( 'z' )
        call chkonp( bmet, bvec, avec, uu, nv, ng, npass, tr, ti, vr, vi, kr, ki, cfr, cfi, kvc )
     case( 'Z' )
        call pdatst( bmet, bvec, avec, ng, nv, npass, vr, vi, tr, ti, kr, ki, cfr, cfi, uu )
     case( 'zeromoment' )
        call gridmomz( nv, ng, cfr, cfi, kvc, bvec, avec, bmet, tr, ti, vr, vi, uu, kr, ki )
     end select
     if ( com .eq. 'quit' ) exit  
  end do
  !
  write ( 6, * ) 'reader terminus executed'
end program xxxhed
!
!
! outmoded
!    if ( com .eq. 'e' ) call nexafs( nv, ng, avec, bvec, bmet, uu, cfr, cfi, kvc, tr, ti, vr, vi, kr, ki )
!    if ( com .eq. 'H' ) call gmels( nv, ng, kvc, cfr, cfi, bmet )
!    if ( com .eq. 'N' ) call edcint( bvec )
!    if ( com .eq. 'n' ) call rededc( cfr, cfi, kvc, nv, npass, ng, bvec, avec, bmet, vr, vi, tr, ti, kr, ki, uu )
!    if ( com .eq. 'R' ) then
!       write ( 6, * ) 'excitvu is disabled pending zone1 compliance'
!       ! call excitvu( ng, kvc, avec, bvec, nv, cfr, cfi, zn )
!    end if
!    case( 'X' )
!       deallocate ( kr, ki, tr, ti, vr, vi )
!       call xxxhrm( nv, avec, bvec, bmet, zn, cfr, cfi, ng, ng1, kvc, sigr, sigi )
!       stop
