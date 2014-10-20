!
! Copyright (C) 2010 A. Ferretti
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!   from grabwf.f90
!
!#!     subroutine grabwf(filenum, maxband, maxnpw, kg_unshift, kg_shift, &
!#!     &     eigen_un, eigen_sh, occ_un, occ_sh, cg_un, cg_sh, occ_max,      &
!#!     &     unocc_max, nband, un_npw, sh_npw, noshift)


!   from wfconvert2.f90
!
!#!           call grabwf(wfkinfile, maxband, maxnpw, kg_unshift, kg_shift, &
!#!     &     eigen_un, eigen_sh, occ_un, occ_sh, cg_un, cg_sh, brange(2) &
!#!     &     brange(4), nband, un_npw, sh_npw, noshift)



!! IN
!!      prefix, work_dir, maxband, maxnpw, occ_max, unocc_max, noshift 


!! OUT
!!     kg_unshift, eigen_un, occ_un, cg_un, un_npw, nband(1) 
!!     kg_shift,   eigen_sh, occ_sh, cg_sh, sh_npw, nband(2)

!! CONVERSION

!!     occ_max = brange(2)
!!     unocc_max = brange(4)


!========================
  SUBROUTINE grabwf(prefix,work_dir,ik,maxband,maxnpw,kg_unshift,kg_shift, &
     &  eigen_un,eigen_sh,occ_un,occ_sh,cg_un,cg_sh,occ_max,unocc_max,  &
     &  nband,un_npw,sh_npw,noshift)
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
  LOGICAL, intent(in)            :: noshift
  INTEGER, intent(in)            :: ik
  INTEGER, intent(in)            :: maxband, maxnpw, occ_max, unocc_max
  !
!  INTEGER            :: ik           ! indexes of the kpt and band wfc
!  INTEGER            :: ib           ! to be read

  !
  ! output variables
  !
  INTEGER, intent(out)           :: nband(2), un_npw, sh_npw
  INTEGER, intent(out)           :: kg_unshift(3,maxnpw), kg_shift(3,maxnpw)
  REAL(DP), intent(out)          :: eigen_un(maxband), eigen_sh(maxband)
  REAL(DP), intent(out)          :: occ_un(maxband), occ_sh(maxband)
  REAL(DP), intent(out)          :: cg_un(maxband,2*maxnpw), cg_sh(maxband,2*maxnpw)

!  NAMELIST /INPUT/ prefix, work_dir, ik, ib
!  NAMELIST /INPUT/ prefix, work_dir


!#!      integer :: filenum, maxband, maxnpw, nband(2),                    &
!#!     &  iband,ii,un_npw,sh_npw,nspinor, occ_max, unocc_max
!#!      double precision :: eigen_un(maxband), eigen_sh(maxband),         &
!#!     &  occ_un(maxband), occ_sh(maxband), cg_un(maxband,2*maxnpw),      &
!#!     &  cg_sh(maxband,2*maxnpw)
!#!      integer :: kg_unshift(3,maxnpw), kg_shift(3,maxnpw)
!#!      logical :: noshift



  !
  ! local vars
  !
  CHARACTER(7) :: subname='example'
  !
  REAL(DP)     :: avec(3,3), alat
  INTEGER      :: npw
  !
  CHARACTER(19)  :: sysname, workdir
  CHARACTER(256)           :: dirname, filename, str_units
  INTEGER,     ALLOCATABLE :: igv(:,:), igk(:)
  COMPLEX(DP), ALLOCATABLE :: wfc(:,:)
  !
  INTEGER      :: ierr, ig, ngm_test

  CHARACTER(14)   :: wfname
  INTEGER         :: ib,ic
  INTEGER         :: nk1, nk2, nk3, k1, k2, k3

  CHARACTER(13)            :: frmt
  REAL(DP), ALLOCATABLE    :: eig(:), occ(:), temp_eigen_un(:,:), temp_occ_un(:,:)

  COMPLEX(DP)     :: cg_un_tmp(maxnpw,maxband), cg_sh_tmp(maxnpw,maxband)
  INTEGER         :: npwx, ngm, ngms, npol, nbnd, nkpts, nspin

  LOGICAL         :: lsda, lpoint_dir


!
!----------------------------------------
! main Body
!----------------------------------------
!

  WRITE(stdout, "(/,'< Wavefunction Parser > ')" )
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
!  WRITE(stdout, "(/,'Reading input namelist...')" )
  !
  !
!  READ( stdin, INPUT, IOSTAT=ierr)
!  IF ( ierr/=0 ) CALL errore(subname,'reading INPUT namelist',ABS(ierr))
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
  !==========================
  ! read planewave info
  !==========================
  !
  ! npwk
  ! kg
  !
  write(6,*) ":Reading the planewave info"
  !
  CALL qexml_read_gk( IK=ik, NPWK=un_npw, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
!  write(6,*) ' npw  = ', un_npw
!  write(6,*) ' kg_un = ', kg_unshift


  !
  !==========================
  ! read g-vectors
  !==========================
  !
  ! kg
  !
  write(6,*) ":Reading the g-vectors"
  !
  CALL qexml_read_gk( IK=ik, NPWK=ngm_test, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  write(6,*) "kg  is size 3 x ", maxnpw
  write(6,*) "igv is size 3 x ", ngm_test
  !
  CALL qexml_read_gk( IK=ik, IGK=kg_unshift, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  !  write(6,*) ' kg_un = ', kg_unshift



  !
  !==========================
  ! read wavevectors
  !==========================
  !
  ! wf
  !
  write(6,*) ":Reading the wavefunctions"
  !
  write(6,*) " Sending CG in as ", maxnpw, " by ", maxband
  CALL qexml_read_wfc( IBNDS=1, IBNDE=maxband, IK=ik, WF=cg_un_tmp, &
       &  IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  !
  do ib = 1, maxband
    do ig = 1, maxnpw
      cg_un( ib, 2*ig -1 ) = REAL ( cg_un_tmp(ig,ib) )
      cg_un( ib, 2*ig    ) = AIMAG( cg_un_tmp(ig,ib) )
      cg_sh( ib, 2*ig -1 ) = REAL ( cg_un_tmp(ig,ib) )
      cg_sh( ib, 2*ig    ) = AIMAG( cg_un_tmp(ig,ib) )
!      cg_un( ib, 2*(ig-1)   ) = REAL ( cg_un_tmp(ig,ib) )
!      cg_un( ib, 2*(ig-1)+1 ) = AIMAG( cg_un_tmp(ig,ib) )
!      cg_sh( ib, 2*(ig-1)   ) = REAL ( cg_un_tmp(ig,ib) )
!      cg_sh( ib, 2*(ig-1)+1 ) = AIMAG( cg_un_tmp(ig,ib) )
!      cg_sh( ib, 2*(ig-1)   ) = REAL ( cg_sh_tmp(ig,ib) )
!      cg_sh( ib, 2*(ig-1)+1 ) = AIMAG( cg_sh_tmp(ig,ib) )
    end do
  end do
  !
!  write(6,*) ' cg_un = ', cg_un
!  write(6,*) ' cg_sh = ', cg_sh


  !
  !==========================
  ! read occupancies
  !==========================
  !
  ! nband
  !
  write(6,*) ":Reading occupancies"
  !
  CALL qexml_read_occ( NSTATES_UP=nband(1), IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ', ABS(ierr))
  !
!  write(6,*) ' nband  = ', nband(1)




  !
  !==========================
  ! read bands
  !==========================
  !
  ! eigen
  ! occ
  !
  write(6,*) ":Reading the bands"
  !
#ifdef __QE51
  !=======================
  ! This is the worst way to do this, but should work for now
  CALL qexml_read_bands_info( NUM_K_POINTS=nkpts, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML read info ',ABS(ierr))

  allocate( temp_eigen_un( maxband, nkpts ), temp_occ_un( maxband, nkpts ) )
  
  CALL qexml_read_bands_pw( nkpts, maxband, nkpts, .false., .true., filename, ET=temp_eigen_un, WG=temp_occ_un, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
  eigen_un(:) = temp_eigen_un( :, ik )
  occ_un(:)   = temp_occ_un( :, ik )
  deallocate( temp_eigen_un, temp_occ_un )
#else
  CALL qexml_read_bands( IK=ik, EIG=eigen_un, OCC=occ_un, IERR=ierr )
  IF ( ierr/=0 ) CALL errore(subname,'QEXML reading ',ABS(ierr))
#endif
  !
!  write(6,*) ' eigen_un = ', eigen_un
!  write(6,*) ' occ_un   = ', occ_un



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

!  sh_npw = un_npw
!  nband(2) = nband(1)
!  eigen_sh = eigen_un
!  occ_sh = occ_un
!  kg_unshift = 1
!  kg_shift(:,:) = kg_unshift(:,:)
!  cg_un = 1.0
!  cg_sh = cg_un

  write(6,*) "DONE WITH GRABWF"

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

END SUBROUTINE grabwf

