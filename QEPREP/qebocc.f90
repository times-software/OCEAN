!
! Copyright (C) 2010 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================
  PROGRAM qebocc
  !========================
  !
  ! Simple example to show how to use the QEXML library
  ! to read data from the .save directory written by QE
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
  REAL(DP), PARAMETER :: ev2ryd = 0.073498618

  !
  ! local vars
  !
  CHARACTER(7)    :: subname='example'
  CHARACTER(14)   :: wfname
  CHARACTER(256)  :: dirname, filename, str_units
  !
  INTEGER      :: i, ivbm, icbm
  INTEGER      :: nbnd, ib, nkpts, ik, nval
  INTEGER      :: ierr
  !
  REAL(DP), ALLOCATABLE    :: occ(:,:)
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
  WRITE(stdout, "(/, 'Init QEXML library...')" )
  !
  ! we should probably read in work_dir and prefix here
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
  ! SUBROUTINE qexml_read_bands( ik, ispin, nbnd, eig, energy_units, occ, ef, ierr )
  !
  ! occ
  !
!  write(6,*) ":Reading the band info"
  !
#ifdef __QE51
  CALL qexml_read_bands_info( NUM_K_POINTS=nkpts, NBND=nbnd, IERR=ierr )
#else
  CALL qexml_read_bands_info( NBND=nbnd, IERR=ierr )
#endif
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  !
#ifdef __QE51
  allocate( occ(nbnd,nkpts) )
! num_k_points and nkstot actually need to be the same
  CALL qexml_read_bands_pw( nkpts, nbnd, nkpts, .false., .true., filename, WG=occ, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
#else
  allocate( occ(nbnd,1) )
  CALL qexml_read_bands( IK=1, OCC=occ(:,1), IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
#endif
  !
  !
  open( unit=91, file='bands.out',form='formatted',status='unknown')
  do ib = 1, nbnd
    write(91,*) occ(ib,1)
  end do
  close( 91 )

  !
  ! local cleanup
  !
  deallocate( occ )


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

END PROGRAM qebocc

