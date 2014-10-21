!
! Copyright (C) 2010 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================
  PROGRAM qeband
  !========================
  !
  ! Adapted from QuantumESPRESSO QEXML example
  ! KG; July 2012
  !
  ! General comments:
  !
  !   - first init the library
  !
  !   - for each data required, first get the
  !     dimensions, then allocate the target array
  !     and finally read the data with a second call
  !     to the proper qexml read routine
  !     (shown below)
  ! 
  !   - data that don't need any dynamical allocation
  !     (scalar or small arrays) can be read directly
  !
  !   - explicit error handling through the use of IERR arguments
  !     is required
  !
  USE qexml_module
  IMPLICIT NONE
  !
  ! parameters
  !
  INTEGER, PARAMETER  :: iunit = 10
  INTEGER, PARAMETER  :: stdin = 5
  INTEGER, PARAMETER  :: stdout = 6
  INTEGER, PARAMETER  :: DP=KIND(1.0d0)
  !
  REAL(DP), PARAMETER :: ha2ev  = 27.211396132
  REAL(DP), PARAMETER :: ha2ryd = 2.000000000
  REAL(DP), PARAMETER :: ev2ryd = 0.073498618
  REAL(DP), PARAMETER :: ryd2ev = 13.605698066
  !
  ! local vars
  !
  CHARACTER(7)    :: subname='example'
  CHARACTER(9)    :: units
  CHARACTER(14)   :: wfname
  CHARACTER(256)  :: dirname, filename, str_units
  !
  INTEGER      :: i, ivbm, icbm
  INTEGER      :: nbnd, ib, nkpts, ik, nval
  INTEGER      :: ierr
  !
  REAL(DP)     :: efermi, egap
  REAL(DP), ALLOCATABLE    :: sorted_e(:), sorted_o(:)
  REAL(DP), ALLOCATABLE    :: enk(:,:), occ(:,:)
  !
!
!----------------------------------------
! main Body
!----------------------------------------
!

  WRITE(stdout, "(/,'< Wavefunction Parser > ')" )

  !
  !==========================
  ! init the qexml library
  !==========================
  !
!  WRITE(stdout, "(/, 'Init QEXML library...')" )
  !
  dirname = TRIM("./Out") // '/' // TRIM("system") // '.save/'
  CALL qexml_init( iunit, DIR=dirname )

  filename = TRIM(dirname) // "data-file.xml"
  !
!  write(6,*) " dirname = ", dirname
!  write(6,*) " filename = ", filename
  !
  CALL qexml_openfile( filename, "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'opening dftdata file',ABS(ierr))



  !
  !==========================
  ! read band info
  !==========================
  !
  !  SUBROUTINE qexml_read_bands_info( nbnd, num_k_points, nspin, noncolin,
  !                                    natomwfc, nelec, ef, energy_units, k_units, ierr )
  !
!  write(6,*) ":Reading the band info"
  !
  CALL qexml_read_bands_info( NBND=nbnd, NUM_K_POINTS=nkpts, EF=efermi, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
!  write(6,*) ' nbnd  = ', nbnd
!  write(6,*) ' nkpts = ', nkpts
!  write(6,*) ' Efermi = ', efermi

  open( unit=91, file='efermi',form='formatted',status='unknown')
!  write(91,*) efermi * ha2ryd     ! in Ry
  write(91,*) efermi * ha2ev      ! in eV
  close( 91 )

  allocate( enk(nbnd,nkpts), occ(nbnd,nkpts) )
  nval = nbnd * nkpts


  !
  !==========================
  ! read band info
  !==========================
  !
  ! SUBROUTINE qexml_read_bands( ik, ispin, nbnd, eig, energy_units, occ, ef, ierr )
  !
  ! eig, occ
  !
!  write(6,*) ":Reading the band info"
  !
#ifdef __QE51 
  write(6,*) 'QE 5.1!!'
  CALL qexml_read_bands_pw( nkpts, nbnd, nkpts, .false., .true., filename, ET=enk, WG=occ, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
#else
  do ik = 1, nkpts
     CALL qexml_read_bands( IK=ik, EIG=enk(:,ik), OCC=occ(:,ik), ENERGY_UNITS=units, IERR=ierr )
     IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  end do
#endif
  !
!  write(6,*) " BAND ENERGIES ARE IN ", trim( units ), " !!!"
  !
  allocate( sorted_e(nval), sorted_o(nval) )
  !
  call sort( enk, occ, sorted_e, sorted_o, nbnd, nkpts, nval )
!  open( unit=91, file='test',form='formatted',status='unknown')
!  open( unit=92, file='test_o',form='formatted',status='unknown')
!  do i = 1, nval
!    write(91,*) sorted_e(i) * ha2ryd   ! in Ry
!!    write(91,*) sorted_e(i) * ha2ev    ! in eV
!    write(92,*) sorted_o(i)
!  end do
!  close( 92 )
!  close( 91 )
  !
  call getgap( sorted_o, ivbm, icbm, nval )
!  egap = ( sorted_e(icbm) - sorted_e(ivbm) ) * ev2ryd
!  efermi = ( sorted_e(ivbm) * ev2ryd ) + 0.5 * egap
  egap = ( sorted_e(icbm) - sorted_e(ivbm) ) * ha2ev   ! in eV
  efermi = ( sorted_e(ivbm)*ha2ev + 0.5*egap ) * ev2ryd    ! in Ry
  open( unit=91, file='efermiinrydberg.ipt',form='formatted',status='unknown')
  write(91,*) efermi  ! in Ry
  close( 91 )
  !
  open( unit=91, file='ldagap',form='formatted',status='unknown')
  write(91,*) egap    ! in eV
  close( 91 )
  !
  open( unit=91, file='eshift.ipt',form='formatted',status='unknown')
  write(91,*) -sorted_e(icbm) * ha2ev   ! in eV
  close( 91 )
  !
!  call getclips( sorted_o, ivbm, icbm )
  open( unit=91, file='ldaclips',form='formatted',status='unknown')
  write(91,*) efermi * ryd2ev,  sorted_e(nval) * ha2ev   ! in eV
  close( 91 )
  !


  open( unit=91, file='bands.out',form='formatted',status='unknown')
  if ( nkpts .gt. 1 ) then
    do ib = 1, nbnd
      write(91,*) occ(ib,2)
    end do
  else
    do ib = 1, nbnd
      write(91,*) occ(ib,1)
    end do
  end if
  close( 91 )


  !
  ! local cleanup
  !
  deallocate( enk, occ )
  deallocate( sorted_e, sorted_o )


  !
  !==========================
  ! finalize the qexml read
  !==========================
  !
  WRITE(stdout, "(/,'Finalize QEXML...')" )
  !
  CALL qexml_closefile ( "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'closing dftdata file',ABS(ierr))
  

CONTAINS


SUBROUTINE sort( enk, occ, sorted_e, sorted_o, nbnd, nkpts, nval )
  !
  INTEGER, INTENT(in) :: nbnd, nkpts, nval
  REAL(DP), INTENT(in) :: enk(nbnd,nkpts), occ(nbnd,nkpts)
  !
  REAL(DP), INTENT(out) :: sorted_e(nval), sorted_o(nval)
  !
  !
  INTEGER  :: i, ib, ik, iv1, iv2
  REAL(DP) :: temp_e, temp_o

  ! initialize the list
  i = 1
  do ib = 1, nbnd
    do ik = 1, nkpts
      sorted_e(i) = enk(ib,ik)
      sorted_o(i) = occ(ib,ik)
      i = i + 1
    end do
  end do  


  ! bubble sort
  !
  do iv1 = nval, 2, -1
    do iv2 = 1, iv1
      if ( sorted_e(iv1) .lt. sorted_e(iv2) ) then
        temp_e = sorted_e(iv1)
        temp_o = sorted_o(iv1)
        sorted_e(iv1) = sorted_e(iv2)
        sorted_o(iv1) = sorted_o(iv2)
        sorted_e(iv2) = temp_e
        sorted_o(iv2) = temp_o
      end if
    end do
  end do

END SUBROUTINE sort


SUBROUTINE getgap( sorted_o, ivbm, icbm, nval )
  !
  INTEGER, INTENT(in)  :: nval
  REAL(DP), INTENT(in) :: sorted_o(nbnd)
  !
  INTEGER, INTENT(out) :: ivbm, icbm
  !
  
  icbm = 1
  do while ( sorted_o(icbm) .eq. 1.0d0 )
    icbm = icbm + 1
  end do
  ivbm = icbm - 1

END SUBROUTINE getgap


! SUBROUTINE getclips()

! END SUBROUTINE getclips


!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE errore( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output: 
  ! ... if ierr <= 0 it does nothing, 
  ! ... if ierr  > 0 it stops.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routinee
    ! the output messagee
  INTEGER,          INTENT(IN) :: ierr
    ! the error flag
    !
  !
  !
  IF ( ierr <= 0 ) RETURN
  !
  ! ... the error message is written un the "*" unit
  !
  WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
  WRITE( UNIT = *, &
         FMT = '(5X,"from ",A," : error #",I10)' ) calling_routine, ierr
  WRITE( UNIT = *, FMT = '(5X,A)' ) message
  WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
  !
  WRITE( *, '("     stopping ...")' )
  !
  STOP 2
  !
  RETURN
  !
END SUBROUTINE errore

END PROGRAM qeband

