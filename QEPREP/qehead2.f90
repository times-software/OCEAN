!
! Copyright (C) 2010 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================
  SUBROUTINE qehead2(prefix,work_dir,maxband,nsppol,nspinor,nkpt)
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
  CHARACTER(6), intent(in)     :: prefix       ! used to locate the .save directory
  CHARACTER(256), intent(in)     :: work_dir     !
  !
  INTEGER            :: ik           ! indexes of the kpt and band wfc
  INTEGER            :: ib           ! to be read

  !
  ! output variables
  !
  INTEGER, intent(out)            :: maxband, nsppol, nspinor, nkpt

!  NAMELIST /INPUT/ prefix, work_dir, ik, ib
!  NAMELIST /INPUT/ prefix, work_dir

  !
  ! local vars
  !
  CHARACTER(19)  :: sysname, workdir
  CHARACTER(7)   :: subname='example'
  !
  REAL(DP)     :: avec(3,3), alat
  INTEGER      :: npw
  !
  CHARACTER(256)           :: dirname, filename, str_units
  INTEGER,     ALLOCATABLE :: igv(:,:), igk(:)
  COMPLEX(DP), ALLOCATABLE :: wfc(:,:)
  !
  INTEGER      :: ierr, ig

  CHARACTER(14)   :: wfname
  INTEGER         :: ic
  INTEGER         :: nk1, nk2, nk3, k1, k2, k3

  CHARACTER(13)            :: frmt
  REAL(DP), ALLOCATABLE    :: eig(:), occ(:)


  INTEGER         :: npwx, ngm, ngms, npol, nbnd, nkpts, nspin


!
!----------------------------------------
! main Body
!----------------------------------------
!

  WRITE(stdout, "(/,'< Read File Header > ')" )
  !
  open( unit=97, file="../../Common/prefix", form='formatted', status='unknown')
  read( 97, *) sysname
  close( 97 )
  !
  open( unit=97, file="../../Common/work_dir", form='formatted', status='unknown')
  read( 97, *) workdir
  close( 97 )
  !
  write(6,*) " Prefix = ", sysname
  write(6,*) " Prefix = ", TRIM(sysname)
  write(6,*) " WorkDir = ", TRIM(workdir)


  !
  !==========================
  ! init the qexml library
  !==========================
  !
  WRITE(stdout, "(/, 'Init QEXML library...')" )
  !
  dirname = TRIM(workdir) // '/' // TRIM(sysname) // '.save/'
  write(6,*) " dirname = ", dirname
  CALL qexml_init( iunit, DIR=dirname )

  filename = TRIM(dirname) // "data-file.xml"
  !
  write(6,*) " filename = ", filename
  !
  CALL qexml_openfile( filename, "read", IERR=ierr )
  IF ( ierr/=0) CALL errore(subname,'opening dftdata file',ABS(ierr))



  !
  !==========================
  ! read band info
  !==========================
  !
  ! num_k_points
  ! nbnd
  ! nspin
  !
  write(6,*) ":Reading the band info"
  !
  CALL qexml_read_bands_info( NBND=nbnd, NUM_K_POINTS=nkpts, NSPIN=nspin, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  write(6,*) ' nbnd  = ', nbnd
  write(6,*) ' nkpts = ', nkpts
  write(6,*) ' nspin = ', nspin


  !
  !==========================
  ! read planewave info
  !==========================
  !
  ! npwx
  !
  write(6,*) ":Reading the planewave info"
  !
  WRITE(stdout, "(/, 'Read the planewave info ...')" )
  !
  CALL qexml_read_planewaves( NPWX=npwx, NGM=ngm, NGMS=ngms, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  write(6,*) ' npwx = ', npwx
  write(6,*) ' ngm  = ', ngm
  write(6,*) ' ngms = ', ngms
  !


  !
  !==========================
  ! read spin info    
  !==========================
  !
  ! npol
  !
  write(6,*) ":Reading the spin info"
  !
  CALL qexml_read_spin( NPOL=npol, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  write(6,*) ' npol = ', npol



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
!  maxnpw = npwx  !KG! this worked for first C attempt
!  maxnpw = ngm  !KG! changed to this
  maxband = nbnd
  nsppol = npol
  nspinor = nspin
  nkpt = nkpts



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

END SUBROUTINE qehead2

