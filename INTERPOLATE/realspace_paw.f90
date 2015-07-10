subroutine realspace_paw( npt, lmin, lmax, nptot, nproj, posn, drel, ae_func, ps_func, radfunc, & 
                          ae_psi, ps_psi, delta_psi, prefs, loud, wpt )
  use kinds, only : dp
  implicit none
  integer, intent( in ) :: npt
  real(dp), intent( in ) :: posn( 3, npt ), drel( npt )
  integer, intent( in ) :: lmin, lmax, nproj(lmin:lmax), nptot

  complex(dp), intent( out ) :: delta_psi( npt, nptot )
  complex(dp), intent( out ) :: ae_psi( npt, nptot )
  complex(dp), intent( out ) :: ps_psi( npt, nptot )
  real(dp), intent( in ) :: prefs( 0 : 1000 ), wpt( npt )
  real(dp), intent( in ) :: ae_func( 2000, nptot ), ps_func( 2000, nptot ), radfunc( 2000 )
  logical, intent( in ) :: loud

  complex(dp), allocatable :: ylm( : )
  integer :: l, m, iproj, iter, iter2, iptot, iprev, drel_map( npt ), ii
  real(dp) :: prev, frac, ae_su, ps_su

  allocate( ylm( npt ) )

  iptot = 0
  do l = lmin, lmax
!    allocate( aefunc( 2000, nproj(l) ), psfunc( 2000, nproj(l) ), radfunc( 2000 ) )

!    write(filnam,'(a2,i1.1,a1,i3.3)') 'ae', l, 'z', z
!    open(unit=99,file=filnam,form='formatted', status='old' )
!    rewind(99)
!    do i = 1, 2000
!      read(99,*) radfunc(i), aefunc(i,:)
!    enddo
!    close(99)
!
!    write(filnam,'(a2,i1.1,a1,i3.3)') 'ps', l, 'z', z
!    open(unit=99,file=filnam,form='formatted', status='old' )
!    rewind(99)
!    do i = 1, 2000
!      read(99,*) radfunc(i), psfunc(i,:)
!    enddo
!    close(99)

!    call find_drel_align( 2000, npt, radfunc, drel, drel_map )
    prev = -1.0
    iprev = 1
    do iter = 1, npt
      if( drel(iter) .eq. prev )  then
        drel_map(iter) = iprev
        cycle
      endif
      do iter2 = 2, 2000
        if( radfunc(iter2) .ge. drel(iter ) ) then
          drel_map(iter) = iter2 - 1
          goto 111
        endif
      enddo
111   continue
      prev = drel(iter)
      iprev = drel_map(iter)
    enddo

    do m = -l, l
      do iter = 1, npt
        call getylm( l, m, posn(1,iter),posn(2,iter),posn(3,iter), ylm(iter), prefs )
      enddo
      do iproj = 1, nproj(l)
        iptot = 1
        do ii = lmin, l-1
          iptot = iptot + nproj(1)
        enddo
        ae_su = 0
        ps_su = 0
        do iter = 1, npt
          frac = ( drel( iter ) - radfunc( drel_map( iter ) ) ) &
               / ( radfunc( drel_map( iter ) + 1 ) - radfunc( drel_map( iter ) ) )
!          if( frac .lt. 0 ) stop
          if( loud .and. ( mod(iter-1,36) .eq. 0 ) .and. drel( iter ) .lt. 3.d0 .and. l .eq. 0 .and. iproj .eq. 1 ) then
            write(6,*) drel(iter), radfunc( drel_map( iter ) + 1 ), radfunc( drel_map( iter ) )
            write(6,*) ae_func( drel_map(iter) + 1, iptot ), ae_func( drel_map(iter), iptot )
            write(6,*) frac, &
             frac * ( ae_func( drel_map(iter) + 1, iptot ) - ae_func( drel_map(iter), iptot ) ) &
                  + ae_func( drel_map(iter), iptot )
          endif

!          ae_psi( iter, iptot ) = ae_func( drel_map(iter), iptot )
!          ps_psi( iter, iptot ) = ps_func( drel_map(iter), iptot )
          
          ae_psi( iter, iptot ) = frac * ( ae_func( drel_map(iter) + 1, iptot ) &
                                         - ae_func( drel_map(iter), iptot ) ) &
                                + ae_func( drel_map(iter), iptot ) * ylm( iter )
          ps_psi( iter, iptot ) = frac * ( ps_func( drel_map(iter) + 1, iptot ) &
                                         - ps_func( drel_map(iter), iptot ) ) &
                                +  ps_func( drel_map(iter), iptot ) * ylm( iter )
          delta_psi( iter, iptot ) = ae_psi( iter, iptot ) - ps_psi( iter, iptot )*0
          ae_su = ae_su + ae_psi( iter, iptot ) * conjg( ae_psi( iter, iptot ) ) * wpt( iter )
          ps_su = ps_su + ps_psi( iter, iptot ) * conjg( ps_psi( iter, iptot ) ) * wpt( iter )
        enddo

        if( loud ) write(6,'(3i4,4f12.5)') l, m, iproj, ae_su, ps_su, ylm(1)
      enddo
    enddo

!    deallocate( aefunc, psfunc, radfunc )
  enddo
  

end subroutine realspace_paw

