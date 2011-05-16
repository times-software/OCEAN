program wan2u1
implicit none

  integer :: kpt, ngx, ngy, ngz, nband, iter, ik, nx, biter
  integer :: ix, inx, iny, inz, nkpt( 3 ), ikx, iky, ikz, torecp, idwrk
  integer :: outx, outy, outz, giter, xmesh( 3 ), x_fac( 3 )
  integer, allocatable :: gvecs( :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: cres( :, :, : )
  double precision, allocatable :: fofr( :, :, : ), fofi(:, :, : ), wrk ( : ), &
          fofr_big( :, :, :, : ), fofi_big( :, :, :, : )
  real( kind=kind(1.d0) ) :: normrecp
  character *10 :: filename
  character *12 ::  psiname
  logical :: test, loud

  loud = .false.
!  read(5,*) kpt
!  read(5,*) test

  open(unit=98,file='u1.dat', form='unformatted' )
  rewind(98)

  open(unit=99,file='nkpt',form='formatted',status='old')
  read(99,*) nkpt(:)
  close(99 )
  kpt = nkpt(1)*nkpt(2)*nkpt(3)
  open(unit=99,file='masterwfile', form='formatted')
  rewind(99)
  write(99,*) kpt
  close( 99 )

  open(unit=99,file='xmesh.ipt',form='formatted',status='old')
  read(99,*) xmesh( : )
  close(99)

  open(unit=97,file='listwfile',form='formatted' )
  rewind( 97 )


  call chkfftrecp( torecp, normrecp, loud )


!  do ikx = 1, nkpt( 1 )
!  do iky = 1, nkpt( 2 )
!  do ikz = 1, nkpt( 3 )
!  iter = ( ikz - 1 ) * nkpt( 1 )  * nkpt( 2 ) + (iky - 1 ) * nkpt( 1 ) + ikx
  do iter = 1, kpt
    write(filename, '(A3,I5.5,A2)' ) 'UNK', iter, '.1'
    write(psiname, '(A4,I8.8)') '.Psi', iter
    write(6,*) filename, psiname
    write(97,*) iter, psiname
    open(96,file=psiname, form='unformatted' )
    rewind(96)

    open(99, file=filename, form='unformatted', status='old' )
    read(99) ngx, ngy, ngz, ik, nband
    write(6,*) ngx, ngy, ngz, ik, nband
    idwrk = 2 * max( ngx, ngy, ngz ) * ( max( ngx, ngy, ngz ) + 1 )
    allocate( wrk( idwrk ) )
!    if( test .ne. .true. ) then
    if( ( mod( ngx, xmesh( 1 ) ) .ne. 0 ) .or. &
        ( mod( ngy, xmesh( 2 ) ) .ne. 0 ) .or. &
        ( mod( ngz, xmesh( 3 ) ) .ne. 0 ) )  &
      goto 111

    x_fac( 1 ) = ngx / xmesh( 1 )
    x_fac( 2 ) = ngy / xmesh( 2 )
    x_fac( 3 ) = ngz / xmesh( 3 )

    allocate( cres( xmesh(3), xmesh( 2 ), xmesh( 1 ) ), fofr( ngx, ngy, ngz ), fofi( ngx, ngy, ngz ), &
              fofr_big( ngx, ngy, ngz, nband), fofi_big( ngx, ngy, ngz, nband) )
      do biter = 1, nband
        read(99) (((fofr( inx, iny, inz ), fofi( inx, iny, inz ), &
                 inx = 1, ngx), iny = 1, ngy), inz = 1, ngz )
      
        ix = 1
        do inx = 1, xmesh( 1 )
          do iny = 1, xmesh( 2 )
            do inz = 1, xmesh( 3 )
              cres( inz, iny, inx ) = &
       cmplx( fofr( 1 + ( inx - 1 ) * x_fac( 1 ), 1 + ( iny - 1 ) * x_fac( 2 ), 1 + ( inz -1 ) * x_fac( 3 )), &
              fofi( 1 + ( inx - 1 ) * x_fac( 1 ), 1 + ( iny - 1 ) * x_fac( 2 ), 1 + ( inz -1 ) * x_fac( 3 )) )
              ix = ix + 1
            enddo
          enddo
        enddo
        write(98) cres

        call cfft( fofr, fofi, ngx, ngy, ngz, ngz, torecp, wrk, idwrk )
        fofr = fofr / dble( ngx * ngy * ngz ) !** normrecp
        fofi = fofi / dble( ngx * ngy * ngz ) !** normrecp
        fofr_big( :, :, :, biter ) = fofr( :, :, : )
        fofi_big( :, :, :, biter ) = fofi( :, :, : )
 
      enddo
      giter = 1
      allocate( gvecs( ngx * ngy * ngz, 3 ) )
      do inz = 1, ngz
        outz= inz - 1
        if( outz .gt. floor( dble( ngz ) / 2.d0 ) ) outz = outz - ngz
        do iny = 1, ngy
          outy= iny - 1
          if( outy .gt. floor( dble( ngy ) / 2.d0 ) ) outy = outy - ngy
          do inx = 1, ngx
            outx = inx - 1
            if( outx .gt. floor( dble( ngx ) / 2.d0 ) ) outx = outx - ngx
            gvecs( giter, 1 ) = outx
            gvecs( giter, 2 ) = outy 
            gvecs( giter, 3 ) = outz
            giter = giter + 1
          enddo
        enddo
      enddo
      write(96) ngx*ngy*ngz
      write(96) gvecs
      write(96) fofr_big
      write(96) fofi_big
      deallocate( cres, fofi, fofr, gvecs, wrk, fofr_big, fofi_big )
!    endif
!  enddo
!  enddo
  enddo

  close(98)
  close(97)
!  write(98 ) cres

111 continue
end program wan2u1


