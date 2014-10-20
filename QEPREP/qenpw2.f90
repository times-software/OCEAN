!
! Copyright (C) 2010 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!========================
  SUBROUTINE qenpw2(prefix,work_dir,ik,maxnpw)
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
  INTEGER, PARAMETER :: iunit = 10
  INTEGER, PARAMETER :: stdin = 5
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: DP=KIND(1.0d0)

  !
  ! input variables
  !
  CHARACTER(256), intent(in)     :: prefix       ! used to locate the .save directory
  CHARACTER(256), intent(in)     :: work_dir     !
  INTEGER, intent(in)            :: ik
  !

  !
  ! output variables
  !
  INTEGER, intent(out)           :: maxnpw

  !
  ! local vars
  !
  CHARACTER(7) :: subname='example'
  !
  CHARACTER(19)  :: sysname, workdir
  CHARACTER(256)           :: dirname, filename, str_units
  !
  INTEGER      :: ierr


!
!----------------------------------------
! main Body
!----------------------------------------
!

  WRITE(stdout, "(/,'< Read Number of Planewaves > ')" )
  !
  open( unit=97, file="../../Common/prefix", form='formatted', status='unknown')
  read( 97, *) sysname
  close( 97 )
  !
  open( unit=97, file="../../Common/work_dir", form='formatted', status='unknown')
  read( 97, *) workdir
  close( 97 )
  !
  ! read stdin
  !


  !
  !==========================
  ! init the qexml library
  !==========================
  !
  WRITE(stdout, "(/, 'Init QEXML library...')" )
  !
  dirname = TRIM(workdir) // '/' // TRIM(sysname) // '.save/'
  CALL qexml_init( iunit, DIR=dirname )

  filename = TRIM(dirname) // "data-file.xml"
  !
  write(6,*) " dirname = ", dirname
  write(6,*) " filename = ", filename
  !
  CALL qexml_openfile( filename, "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'opening dftdata file',ABS(ierr))






  !
  ! now read data specific to a given kpt
  !
  WRITE(stdout, "(/, 'Read ik-specific dims...')" )
  !
  CALL qexml_read_gk( ik, NPWK=maxnpw, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ik dims',ABS(ierr))
  !





  !
  !==========================
  ! finalize the qexml read
  !==========================
  !
  WRITE(stdout, "(/,'Finalize QEXML...')" )
 !
  CALL qexml_closefile ( "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'closing dftdata file',ABS(ierr))
  

  !
  ! local cleanup
  !

  write(6,*) "DONE WITH GETNPW"

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

END SUBROUTINE qenpw2

