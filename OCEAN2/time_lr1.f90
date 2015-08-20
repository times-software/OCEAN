! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program time

  implicit none


  complex( kind = kind(1.d0)), allocatable :: Bink( :, :, : ), beta( :, :, : ), psi( :, :, : )
  complex(kind=kind(1.d0)) :: maa, mab, mac, mad, aaa, aab, aac, aad, derpy(8)

  integer, parameter :: num_bands = 1024
  integer, parameter :: num_obf = 16
  integer, parameter :: nkpts = 8
  integer, parameter :: nalpha = 4
  integer, parameter :: niter = 50000

  complex( kind = kind(1.d0) ), parameter :: one = 1.0d0
  complex( kind = kind(1.d0) ), parameter :: zero = 0.0d0

  integer :: ikpt, iobf, ialpha, iband, iter

  character*10 :: td(3)
  integer :: start_time(8), old_time(8), new_time(8), total_times(8)


  allocate( beta( num_obf, nkpts, nalpha ), psi( num_bands, nkpts, nalpha ), &
            Bink( num_bands, num_obf, nkpts ) )


  Bink = 1.0d0
  psi = 1.5d0

  call date_and_time( td(1), td(2), td(3), old_time )
  total_times( : ) = old_time( : ) 
  write(6,'(A10,X,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') 'Trial 1: ', total_times( 5 ), &
            ':', total_times( 6 ), ':', total_times( 7 ), '.', total_times( 8 )

  write(6,*) niter*dble( nkpts * nalpha ) * dble( num_bands * num_obf ) / (2493.856d0*1000*1000)
  
  call date_and_time( td(1), td(2), td(3), old_time )

  do iter =1,niter 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( ialpha, ikpt ) SHARED( old_time, new_time )
  do ialpha = 1, nalpha
    do ikpt = 1, nkpts
      call ZGEMV( 'T', num_bands, num_obf, one, Bink( 1, 1, ikpt ), num_bands, &
                  psi( 1, ikpt, ialpha ), 1, zero, beta( 1, ikpt, ialpha ), 1 )
    enddo
  enddo
!$OMP END PARALLEL DO
  enddo

  call date_and_time( td(1), td(2), td(3), new_time )

  total_times(:) = new_time(:) - old_time(:)
  if( total_times(8) .lt. 0 ) then
      total_times(7) = total_times(7) - 1
      total_times(8) = total_times(8) + 1000
  else if ( total_times(8) .gt. 999 ) then
      total_times(7) = total_times(7) + 1
      total_times(8) = total_times(8) - 1000
  endif
  if( total_times(7) .lt. 0 ) then
      total_times(6) = total_times(6) - 1
      total_times(7) = total_times(7) + 60
  else if ( total_times(7) .gt. 59 ) then
      total_times(6) = total_times(6) + 1
      total_times(7) = total_times(7) - 60
  endif
  if( total_times(6) .lt. 0 ) then
      total_times(5) = total_times(5) - 1
      total_times(6) = total_times(6) + 60
  else if ( total_times(6) .gt. 59 ) then
      total_times(5) = total_times(5) + 1
      total_times(6) = total_times(6) - 60
  endif

  write(6,'(A10,X,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') 'Trial 1: ', total_times( 5 ), &
            ':', total_times( 6 ), ':', total_times( 7 ), '.', total_times( 8 )


  Bink = 1.0d0
  psi = 1.5d0

  call date_and_time( td(1), td(2), td(3), old_time )
 
  do iter =1 , niter
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( ialpha, ikpt )
  do ikpt = 1, nkpts
    do ialpha = 1, nalpha 
      call ZGEMV( 'T', num_bands, num_obf, one, Bink( 1, 1, ikpt ), num_bands, &
                  psi( 1, ikpt, ialpha ), 1, zero, beta( 1, ikpt, ialpha ), 1 )
    enddo 
  enddo 
!$OMP END PARALLEL DO
  enddo
 
  call date_and_time( td(1), td(2), td(3), new_time )

  total_times(:) = new_time(:) - old_time(:)
  if( total_times(8) .lt. 0 ) then
      total_times(7) = total_times(7) - 1
      total_times(8) = total_times(8) + 1000
  else if ( total_times(8) .gt. 999 ) then
      total_times(7) = total_times(7) + 1
      total_times(8) = total_times(8) - 1000
  endif
  if( total_times(7) .lt. 0 ) then
      total_times(6) = total_times(6) - 1
      total_times(7) = total_times(7) + 60
  else if ( total_times(7) .gt. 59 ) then
      total_times(6) = total_times(6) + 1
      total_times(7) = total_times(7) - 60
  endif
  if( total_times(6) .lt. 0 ) then
      total_times(5) = total_times(5) - 1
      total_times(6) = total_times(6) + 60
  else if ( total_times(6) .gt. 59 ) then
      total_times(5) = total_times(5) + 1
      total_times(6) = total_times(6) - 60
  endif

  write(6,'(A10,X,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') 'Trial 2: ', total_times( 5 ), &
            ':', total_times( 6 ), ':', total_times( 7 ), '.', total_times( 8 )

  Bink = 1.0d0
  psi = 1.5d0

  call date_and_time( td(1), td(2), td(3), old_time )

  do iter = 1, niter
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( ialpha, ikpt )
  do ialpha = 1, nalpha
    do ikpt = 1, nkpts
      do iobf = 1, num_obf
!        aaa = 0
!        aab = 0
!        aac = 0
!        aad = 0
!        maa = 0
!        mab = 0
!        mac = 0
!        mad = 0
        derpy = 0d0
        do iband = 1, num_bands, 8
!          maa = maa + Bink( iband, iobf, ikpt ) * psi( iband, ikpt, ialpha )
!          mab = mab + Bink( iband+1, iobf, ikpt ) * psi( iband+1, ikpt, ialpha )
!          mac = mac + Bink( iband+2, iobf, ikpt ) * psi( iband+2, ikpt, ialpha )
!          mad = mad + Bink( iband+3, iobf, ikpt ) * psi( iband+3, ikpt, ialpha )
!          aaa = aaa + Bink( iband+4, iobf, ikpt ) * psi( iband+4, ikpt, ialpha )
!          aab = aab + Bink( iband+5, iobf, ikpt ) * psi( iband+5, ikpt, ialpha )
!          aac = aac + Bink( iband+6, iobf, ikpt ) * psi( iband+6, ikpt, ialpha )
!          aad = aad + Bink( iband+7, iobf, ikpt ) * psi( iband+7, ikpt, ialpha )
          derpy( : ) = derpy( : ) + Bink( iband:iband+7, iobf, ikpt ) * psi( iband:iband+7, ikpt, ialpha )
!          beta( iobf, ikpt, ialpha ) = beta( iobf, ikpt, ialpha ) &
!            + sum( Bink( iband:iband+4, iobf, ikpt ) * psi( iband:iband+4, ikpt, ialpha ) )
        enddo
        beta( iobf, ikpt, ialpha ) = sum(derpy(:)) !aaa + aab + aac + aad + maa + mab + mac + mad
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
  enddo

  call date_and_time( td(1), td(2), td(3), new_time )

  total_times(:) = new_time(:) - old_time(:)
  if( total_times(8) .lt. 0 ) then
      total_times(7) = total_times(7) - 1
      total_times(8) = total_times(8) + 1000
  else if ( total_times(8) .gt. 999 ) then
      total_times(7) = total_times(7) + 1
      total_times(8) = total_times(8) - 1000
  endif
  if( total_times(7) .lt. 0 ) then
      total_times(6) = total_times(6) - 1
      total_times(7) = total_times(7) + 60
  else if ( total_times(7) .gt. 59 ) then
      total_times(6) = total_times(6) + 1
      total_times(7) = total_times(7) - 60
  endif
  if( total_times(6) .lt. 0 ) then
      total_times(5) = total_times(5) - 1
      total_times(6) = total_times(6) + 60
  else if ( total_times(6) .gt. 59 ) then
      total_times(5) = total_times(5) + 1
      total_times(6) = total_times(6) - 60
  endif

  write(6,'(A10,X,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') 'Trial 3: ', total_times( 5 ), &
            ':', total_times( 6 ), ':', total_times( 7 ), '.', total_times( 8 )




  deallocate( Bink, beta, psi )

end program time
