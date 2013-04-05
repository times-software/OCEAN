subroutine loaduxzerotau( nx, ny, nz, nbd, nq, nspn, zn, tau, ur, ui )
  implicit none
  !
  integer, parameter :: u2dat = 35
  !
  integer :: nx, ny, nz, nbd, nq, zn( 3 ), nspn
  real( kind = kind(1.0d0)), intent( inout ) :: tau(3)
  real( kind = kind( 1.0d0 ) ), dimension( nx, ny, nz, nbd, nq, nspn ) :: ur, ui
  !
  integer :: iq, ibd, ig, idum( 3 ), ix, iy, iz, ivl, ivh, icl, ich, ispn
  integer :: iq1, iq2, iq3, dumint, icl2, ich2, ivh2, xshift(3)
  integer :: xtarg, ytarg, ztarg, xph, yph, zph
  real( kind = kind( 1.0d0 ) ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
  real( kind = kind( 1.0d0 ) ) :: su, sul, suh
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: tmp_ur, tmp_ui
  logical :: metal, normal
  !
  sul = 1.0d0
  suh = 1.0d0
  open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ivl, ivh, icl, ich
  close( unit=99 )
  ivh2 = ivh
!  icl2 = icl
  open( unit=99, file='metal', form='formatted', status='old')
  read( 99, * ) metal
  close( 99 )
  if( metal ) then
    open( unit=36, file='ibeg.h', form='formatted', status='old' )
  endif
  !
  open(unit=99,file='cks.normal',form='formatted',status='old' )
  read(99,*) normal
  close( 99 )

  if (.not. normal ) then
   write(6,*) 'XES'
   ur = 0.0d0
   ui = 0.0d0
   return
  endif

  if ( nbd .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'
  open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
  rewind u2dat
  write(6,*) 'nspn: ', nspn
  do ispn = 1, nspn
    do iq = 1, nq
      open( unit=99, file='gumatprog', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(2i8)' ) iq, nq
      close( unit=99 )
      if( metal ) then
        read( 36, * ) dumint, ivh2
        ivh2 = ivh2 - 1
      endif

!  Skip all of the occupied bands (and for metals)
      do ibd = ivl, ivh2
        do ig = 1, nx * ny * nz
           read ( u2dat )
        end do
      end do

!!  Skip bands below the fermi level (for metals)
!     do ibd = icl, icl2 - 1
!        do ig = 1, nx * ny * nz
!           read ( u2dat )
!        end do
!     enddo
!
      do ibd = 1, nbd 
        do ix = 1, nx
           do iy = 1, ny
              do iz = 1, nz
                 read ( u2dat ) idum( 1 : 3 ), ur( ix, iy, iz, ibd, iq, ispn ), ui( ix, iy, iz, ibd, iq, ispn ) 
              end do
           end do
        end do
        su = sum( ur( :, :, :, ibd, iq, ispn ) ** 2 + ui( :, :, :, ibd, iq, ispn ) ** 2 )
        sul = min( su, sul )
        suh = max( su, suh )
      end do
! Adding 22 nov 2010 get rid of the un-used wfns at the top
      do ibd = ivh2 + nbd + 1, ivh - ivl + ich - icl + 2
        do ig = 1, nx * ny * nz
           read ( u2dat )
        end do
      enddo
    enddo
  end do
  if( metal ) then
    close( 36 )
  endif
  close( unit=u2dat )
  write ( 6, '(1a16,2f20.15)' ) 'norm bounds ... ', sul, suh
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  ! Reasoning for the below;
  !  If the kmesh is odd, then the FT of the kmesh will give an equal number of
  !  cells on either side of the central one -- the one with the core hole. 
  !  This means that the core hole should be placed as near as (0.5,0.5,0.5)
  !
  !  On the other hand, if the kmesh is even, the bonus cell will be tacked on 
  !  over toward the negative side. This means that we should center the core 
  !  hole as near as possible to (0,0,0)
  !
  !  We test each dimension individually of course. 
  !
  xshift(:) = 0
  if( mod( zn(1), 2 ) .eq. 0 ) then
    xshift( 1 ) = floor( real(nx, kind( 1.0d0 )) * tau(1) )
  else
    xshift( 1 ) = floor( real(nx, kind( 1.0d0 )) * (tau(1)-0.5d0 ) )
  endif
  if( mod( zn(2), 2 ) .eq. 0 ) then
    xshift( 2 ) = floor( real(ny, kind( 1.0d0 )) * tau(2) )
  else
    xshift( 2 ) = floor( real(ny, kind( 1.0d0 )) * (tau(2)-0.5d0 ) )
  endif
  if( mod( zn(3), 2 ) .eq. 0 ) then
    xshift( 3 ) = floor( real(nz, kind( 1.0d0 )) * tau(3) )
  else
    xshift( 3 ) = floor( real(nz, kind( 1.0d0 )) * (tau(3)-0.5d0 ) )
  endif
  ! 
  write(6,*) 'Shifting X-grid by ', xshift(:)
  write(6,*) 'Original tau ', tau(:)
  tau( 1 ) = tau(1) - real(xshift(1), kind( 1.0d0 ))/real(nx, kind( 1.0d0 ))
  tau( 2 ) = tau(2) - real(xshift(2), kind( 1.0d0 ))/real(ny, kind( 1.0d0 ))
  tau( 3 ) = tau(3) - real(xshift(3), kind( 1.0d0 ))/real(nz, kind( 1.0d0 ))
  write(6,*) 'New tau      ', tau(:)

  allocate( tmp_ur( nx, ny, nz, nbd ), tmp_ui( nz, ny, nz, nbd ) )
  

  do ispn = 1, nspn
  iq = 0
  do iq1 = 1, zn( 1 )
     do iq2 = 1, zn( 2 )
        do iq3 = 1, zn( 3 )
           iq = iq + 1
           do iz = 1, nz
             ztarg = iz - xshift( 3 )
             if( ztarg .gt. nz ) then
               ztarg = ztarg - nz
               zph = -nz
             elseif( ztarg .lt. 1 ) then
               ztarg = ztarg + nz
               zph = nz
             else
               zph = 0
             endif
              do iy = 1, ny
                 ytarg = iy - xshift( 2 )
                 if( ytarg .gt. ny ) then
                   ytarg = ytarg - ny
                   yph = -ny
                 elseif( ytarg .lt. 1 ) then
                   ytarg = ytarg + ny
                   yph = ny           
                 else
                   yph = 0
                 endif
                 do ix = 1, nx
                   xtarg = ix - xshift( 1 )
                   if( xtarg .gt. nx ) then
                     xtarg = xtarg - nx
                     xph = -nx
                   elseif( xtarg .lt. 1 ) then
                     xtarg = xtarg + nx
                     xph = nx
                   else
                     xph = 0
                   endif
                    phsx = 2.0d0 * pi * dble( ( xph + ix - 1 ) * ( iq1 - 1 ) ) / dble( nx * zn( 1 ) )
                    phsy = 2.0d0 * pi * dble( ( yph + iy - 1 ) * ( iq2 - 1 ) ) / dble( ny * zn( 2 ) )
                    phsz = 2.0d0 * pi * dble( ( zph + iz - 1 ) * ( iq3 - 1 ) ) / dble( nz * zn( 3 ) )
                    cphs = cos( phsx + phsy + phsz )
                    sphs = sin( phsx + phsy + phsz )
                    do ibd = 1, nbd
                       psir = cphs * ur( ix, iy, iz, ibd, iq, ispn ) - sphs * ui( ix, iy, iz, ibd, iq, ispn )
                       psii = cphs * ui( ix, iy, iz, ibd, iq, ispn ) + sphs * ur( ix, iy, iz, ibd, iq, ispn )
!                       ur( ix, iy, iz, ibd, iq, ispn ) = psir
!                       ui( ix, iy, iz, ibd, iq, ispn ) = psii
                       tmp_ur( xtarg, ytarg, ztarg, ibd ) = psir
                       tmp_ui( xtarg, ytarg, ztarg, ibd ) = psii
                    end do
                 end do
              end do
           end do
           ur( :, :, :, :, iq, ispn ) = tmp_ur( :, :, :, : )
           ui( :, :, :, :, iq, ispn ) = tmp_ui( :, :, :, : )
        end do
     end do
  end do
  end do
  deallocate( tmp_ur, tmp_ui )
  !
  return
end subroutine loaduxzerotau
