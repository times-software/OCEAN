! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
 program shifted_conugtoux
   implicit none
   !
   integer, parameter :: u1dat = 22, iwf = 23
   integer :: nk, nb, nx( 3 ), ng, i, nwfile, iwfile, j, k, nspin, ispin, nbv, nbc, nbands, ishift, brange(4), nshift
   integer :: ig, jg, sh_ng, ib, jb
   integer, allocatable :: g( :, : ), ikl( :, : ), flip_g( :, : )
   integer, allocatable :: sh_flip_g( :, : ), map_g( : ), sh_g( :, : )
   real( kind = kind( 1.0d0 ) ), allocatable :: zr( :, : ), zi( :, : )
   real( kind = kind( 1.0d0 ) ), allocatable :: sh_zr( :, : ), sh_zi( :, : )
   real(kind=kind(1.0d0)) :: su, qval, q1, q2, q3, bv1(3), bv2(3), bv3(3), rr, ii, ri, ir, orthcr, orthci
   character(len=20), allocatable :: wfnam( : )
   character(len=20) :: filePrefix
   character(len=40) :: fileName
   character(len=3) :: DFT
   logical :: is_jdftx, split_dft, shift_dft
   integer, parameter :: tmels = 41
   !
   open( unit=99, file='nkpts', form='formatted', status='old')
   read(99,*) nk
   close(99)
   open( unit=99, file='xmesh.ipt', form='formatted', status='old')
   read(99,*) nx(1), nx(2), nx(3)
   close(99)
   open( unit=99, file='nspin', form='formatted', status='old')
   read(99,*) nspin
   close(99)
   open( unit=99, file='brange.ipt', form='formatted', status='old')
   read(99,*) brange(1:4)
   close(99)

   open(unit=99,file='dft.split', form='formatted', status='old')
   read(99,*) split_dft
   close(99)

   open(unit=99,file='qinunitsofbvectors.ipt',form='formatted', status='old')
   read(99,*) q1,q2,q3
   close(99)

   open(unit=99,file='bvecs', form='formatted', status='old')
     read(99,*) bv1(1),bv1(2),bv1(3)
     read(99,*) bv2(1),bv2(2),bv2(3)
     read(99,*) bv3(1),bv3(2),bv3(3)
   close(99)


   ! at the moment split and shift are identical
   shift_dft = split_dft
   if( shift_dft ) then
     nshift = 2
     open( unit=tmels,file='tmels',form='formatted',status='unknown')
      qval =dsqrt( (q1*bv1(1)+q2*bv2(1)+q3*bv3(1))**2                   &
     &            +(q1*bv1(2)+q2*bv2(2)+q3*bv3(2))**2                   &
     &            +(q1*bv1(3)+q2*bv2(3)+q3*bv3(3))**2)

   else
     nshift = 1 
     qval = 1.0d0
   endif


   !
   inquire( file='dft', exist=is_jdftx )
   if( is_jdftx ) then
     open( unit=99, file='dft', form='formatted', status='old' )
     read(99,*) DFT
     close(99)
     if( DFT .eq. 'jdf' ) then
       is_jdftx = .true. 
     else 
       is_jdftx = .false.
     endif
   endif
   if( is_jdftx .eqv. .false. ) stop 1
   !
   open( unit=99, file='masterwfile', form='formatted', status='unknown' )
   rewind 99
   read ( 99, * ) nwfile
   close( unit=99 )
   allocate( wfnam( nwfile ), ikl( 2, nwfile ) )
   open( unit=99, file='listwfile', form='formatted', status='unknown' )
   rewind 99
   do i = 1, nwfile
      read ( 99, * ) ikl( 1, i ), wfnam( i )
      if ( i .gt. 1 ) ikl( 2, i - 1 ) = ikl( 1, i ) - 1
   end do
   close( unit=99 )
   ikl( 2, nwfile ) = nk
   !
   open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown' )
   rewind u1dat
   iwfile = 1
   k = 0
   do ispin = 1, nspin
     do i = 1, nk
        k = k + 1
        open( unit=99, file='prog', form='formatted', status='unknown' )
        rewind 99
        write ( 99, '(2x,1a6,3x,4i5)' ) 'conv03g', i,nk, ispin, nspin
        close( unit=99 )

        ! here we double-up, first read shifted, then non-shifted
!        filePrefix = 'unshifted/'

        if( split_dft ) then
          nbands = brange( 2 )
        else
          nbands = brange( 4 )
        endif

         write(fileName,'(A,A)') 'unshifted/', wfnam( iwfile )
!         write(6,*) fileName
         if( is_jdftx ) then
!           open ( unit=iwf, file=fileName, form='unformatted',status='old', access='stream' )
           open ( unit=iwf, file=fileName, form='unformatted', status='old', access='stream' )
           read( iwf ) ng
         else
           open ( unit=iwf, file=fileName, form='unformatted',status='old' )
           rewind iwf
           read ( iwf ) ng
         endif
         write(6,*) i, ng, nbands
         allocate( g( ng, 3 ), zr( ng, nbands ), zi( ng, nbands ) )
         if( is_jdftx ) then
           allocate( flip_g( 3, ng ) )
!           do ig = 1, ng
!             read( iwf ) flip_g(3,ig), flip_g(2,ig), flip_g(1,ig)
!           enddo
           read( iwf ) flip_g
           g = transpose( flip_g )
!           deallocate( flip_g )
         else
           read ( iwf ) g
         endif

        read ( iwf ) zr
        read ( iwf ) zi
        close( iwf )

!          if( i .eq. 1 .and. is_jdftx ) then
!            do j = 1, ng
!              write(31,'(3I6,4X,E23.16,4X,E23.16)') g(j,1), g(j,2), g(j,3), zr(j,1), zi(j,1)
!            enddo
!          endif
!        if( i .eq. 1 ) then
        do j = 1, nbands
          su = 0.0d0
          do ig = 1, ng
            su = su + zr(ig,j)**2 + zi(ig,j)**2
          enddo
          write(6,'(I8,I8,A,E22.16)') i, j, '    Norm:  ', su
        enddo

        call gentoreal( nx, nbands, zr, zi, ng, g, u1dat, ( ( i .eq. 1) .and. ( ispin .eq. 1 ) ) )
    
        if( split_dft ) then
!          filePrefix = 'shifted/'
          nbands = brange(4)
          allocate( map_g( ng ) )
          map_g(:) = 0
          write(fileName,'(A,A)') 'shifted/', wfnam( iwfile ) 
!          write(6,*) fileName
          if( is_jdftx ) then
            open ( unit=iwf, file=fileName, form='unformatted',status='old', access='stream' ) 
            read( iwf ) sh_ng
          else
            open ( unit=iwf, file=fileName, form='unformatted',status='old' ) 
            rewind iwf
            read ( iwf ) sh_ng
          endif
          write(6,*) i, sh_ng, nbands
          allocate( sh_g( sh_ng, 3 ), sh_zr( sh_ng, nbands ), sh_zi( sh_ng, nbands ) )
          if( is_jdftx ) then
            allocate( sh_flip_g( 3, sh_ng ) ) 
!            do ig = 1, sh_ng
!              read( iwf ) sh_flip_g(3,ig), sh_flip_g(2,ig), sh_flip_g(1,ig)
!            enddo
            read( iwf ) sh_flip_g
            sh_g = transpose( sh_flip_g ) 
            do ig = 1, min(ng,sh_ng)
              if( flip_g( 1, ig ) .eq. sh_flip_g( 1, ig ) .and. &
                  flip_g( 2, ig ) .eq. sh_flip_g( 2, ig ) .and. &
                  flip_g( 3, ig ) .eq. sh_flip_g( 3, ig ) ) then
                 map_g( ig ) = ig
              else
                do jg = 1, sh_ng !max(1,i-10), sh_ng
                  if( flip_g( 1, ig ) .eq. sh_flip_g( 1, jg ) .and. &
                      flip_g( 2, ig ) .eq. sh_flip_g( 2, jg ) .and. &
                      flip_g( 3, ig ) .eq. sh_flip_g( 3, jg ) ) then
                     map_g( ig ) = jg
                     goto 11
                  endif
                enddo
11              continue
              endif
            enddo
            if( i .eq. 1 .and. ispin .eq. 1 ) then
              do ig = 1, ng
                write(100,*) flip_g(:,ig)
              enddo
              do ig = 1, sh_ng
                write(101,*) sh_flip_g(:,ig)
              enddo
!              do ig = 1, ng
!                if( map_g( ig ) .gt. 0 ) then
!                  write(6,'(8I8)') ig, map_g( ig ), flip_g( :, ig ), sh_flip_g( :, map_g( ig ) )
!                else
!                  write(6,'(5I8)') ig, map_g( ig ), flip_g( :, ig )
!                endif
!              enddo
            endif
            deallocate( sh_flip_g )
          else
            stop
            read ( iwf ) g 
          endif

          read ( iwf ) sh_zr
          read ( iwf ) sh_zi
          close( iwf )

          if( i .le. 2 .and. is_jdftx ) then
            do j = 1, sh_ng
              write(31,'(3I6,4X,E23.16,4X,E23.16)') sh_g(j,1), sh_g(j,2), sh_g(j,3), sh_zr(j,1), sh_zi(j,1)
            enddo
          endif

      
          do ib=brange(3),brange(4)
            do jb=brange(1),brange(2)
              rr = 0.0d0
              ii = 0.0d0
              ri = 0.0d0
              ir = 0.0d0
              do ig = 1, ng
                if( map_g(ig) .gt. 0 ) then
                  rr = rr + zr( ig, jb ) * sh_zr( map_g( ig ), ib )
                  ii = ii + zi( ig, jb ) * sh_zi( map_g( ig ), ib )
                  ri = ri + zr( ig, jb ) * sh_zi( map_g( ig ), ib )
                  ir = ir + zi( ig, jb ) * sh_zr( map_g( ig ), ib )
                endif
              enddo
              orthcr = rr + ii
              orthci = ri - ir
              write(tmels,'(8(1x,1e22.15))')orthcr/qval,orthci/qval,0.0,0.0, &
                                  0.0,0.0,0.0,0.0
            enddo
          enddo

          do j = brange(3), brange(4)
            su = 0.0d0
            do ig = 1, sh_ng
              su = su + sh_zr(ig,j)**2 + sh_zi(ig,j)**2
            enddo
            write(6,'(I8,I8,A,E22.16)') i, j, '    Norm:  ', su
          enddo

          nbands = brange(4)-brange(3)+1
          call gentoreal( nx, nbands, sh_zr(:,brange(3):), sh_zi(:,brange(3):), sh_ng, sh_g, u1dat, & 
                          ( ( i .eq. 1) .and. ( ispin .eq. 1 ) ) )


          deallocate( sh_g, sh_zr, sh_zi )
          deallocate( map_g )
        endif

        deallocate( g, zr, zi )
        if( is_jdftx ) deallocate( flip_g )
        iwfile = iwfile + 1
      end do
   end do
   close( unit=u1dat )
   !
   deallocate( wfnam, ikl )
  !
end program shifted_conugtoux
