MODULE SelfEnergyMod
  USE IOMod
  USE ErrorMod
  USE constants

  TYPE SEInput
     ! User Input:
     ! Run self-energy > 0
     ! 1 use many-pole model
     ! 2 use broadened pole model.
     INTEGER RunMode
     ! User specified EpsInv
     DOUBLE PRECISION Eps0
     ! User specified loss function.
     CHARACTER(30) LossFile
     ! User specified pole representation.
     CHARACTER(30) PoleFile
     ! EGridFile - File that contains energy grid for output.
     CHARACTER(30) EGridFile
  END TYPE SEInput

  TYPE SEData
     ! Data used by SelfEnergy subroutine
     ! Fermi energy, Wigner-Sietz radius
     DOUBLE PRECISION Rs, EGap
     ! Number of energy points.
     INTEGER NEPts
     ! Energy grid relative to vacuum.
     COMPLEX*16,ALLOCATABLE :: EGrid(:)
     ! Number of poles.
     INTEGER NPoles
     ! Pole data. Position, width, and amplitude.
     DOUBLE PRECISION,ALLOCATABLE :: Omega(:), Width(:), g(:)
     ! Self energy
     COMPLEX*16,ALLOCATABLE :: SigDat(:), Z(:)
  END TYPE SEData
     
CONTAINS

  SUBROUTINE ReadSEInp(InData)
    TYPE(SEInput),INTENT(OUT) :: InData

    CHARACTER(14) SEFile
    CHARACTER(200) ErrorMessage(5)

    SEFile = 'selfenergy.inp'
    ErrorMessage(1) = "File: selfenergy.inp"
    ErrorMessage(2) = "Comment lines begin with #*!cC."
    ErrorMessage(3) = "Data lines are ordered as follows:"
    ErrorMessage(4) = "Integer RunMode: 0 (Many Pole) or 1 (Broadened Pole)."
    ErrorMessage(5) = "Strings: LossFile PoleFile EgridFile"
    
    
    CALL ReadData(SEFile, Int1 = InData%RunMode,ErrorMessage = ErrorMessage)
    CALL ReadData(SEFile, String1 = InData%LossFile, String2 = &
         & InData%PoleFile, String3 = InData%EGridFile)
    CALL ReadData(SEFile, Double1 = InData%Eps0)
  END SUBROUTINE ReadSEInp

  SUBROUTINE ReadSEData(InData,SEData1)
    TYPE(SEData),INTENT(OUT) :: SEData1
    TYPE(SEInput),INTENT(IN) :: InData
    DOUBLE PRECISION Energy
    INTEGER iError

    ! Read NEPts from egrid.dat
    !CALL ReadData(InData%EGridFile, Int1 = SEData1%NEPts)    
    !Allocate space for energy grid.
    !ALLOCATE(SEData1%EGrid(SEData1%NEPts), STAT = iError)
    !CALL CheckAllocation(iError,'Allocation error occured during allocation of energy grid array.')
    ! Read energy points from egrid.dat
    !CALL ReadArrayData(InData%EGridFile, DComplex1 = SEData1%EGrid)
    ! Close file.
    !CALL CloseFl(InData%EGridFile)

    ! Read Rs and EGap from PoleFile.
    CALL ReadData(InData%PoleFile, Double1 = SEData1%Rs, Double2 = SEData1%EGap)
    ! Read NPoles from PoleFile
    CALL ReadData(InData%PoleFile, Int1 = SEData1%NPoles)
    !Allocate space for pole data.
    ALLOCATE(SEData1%Omega(SEData1%NPoles), SEData1%Width(SEData1%NPoles), & 
             SEData1%g(SEData1%NPoles), STAT = iError)    
    CALL CheckAllocation(iError,'Allocation error occured during allocation of pole data arrays.')
    ! Read pole data.
    CALL ReadArrayData(InData%PoleFile, Double1 = SEData1%Omega, Double2 = SEData1%Width, &
                       Double3 = SEData1%g)    
    
    ! Close file.
    CALL CloseFl(InData%PoleFile)

    ! Unit conversion from eV to Hartree.
    !SEData1%EGrid(:) = SEData1%EGrid(:)/hart
    SEData1%EGap = SEData1%EGap/hart
    SEData1%Omega(:) = SEData1%Omega(:)/hart
    SEData1%Width(:) = SEData1%Width(:)/hart

    ! Allocate space
    !ALLOCATE(SEData1%SigDat(SEData1%NEPts))
    !ALLOCATE(SEData1%Z(SEData1%NEPts))
  END SUBROUTINE ReadSEData
  
  SUBROUTINE SEnergy(SEData1, InData)
    TYPE(SEData),INTENT(INOUT) :: SEData1
    TYPE(SEInput),INTENT(IN) :: InData
    COMPLEX*16 SigmaF, Sigma0, Z, En
    INTEGER iE
    LOGICAL UseBP, OnShll
    
    OnShll = .TRUE.
    IF(SEData1%EGap.gt.0.d0) OnShll = .FALSE.
    SELECT CASE(InData%RunMode)
       CASE(0)
          UseBP = .FALSE.
       CASE(1)
          UseBP = .TRUE.
       CASE DEFAULT
          UseBP = .FALSE.
    END SELECT

    ! Find Sigma(EFermi)
    En = 0.001d0
    CALL CSigZ(En, 0.d0, SEData1%Rs, SigmaF, Z, SEData1%Omega, SEData1%Width, SEData1%g, &
               SEData1%NPoles, .TRUE., UseBP, 0.d0)

    ! Now calculate Sigma(E) and Z(E)
    PRINT '(A)', 'Beginning calculation of the self energy.'
    DO iE = 1, SEData1%NEPts
       IF(MOD(iE,100).eq.0) THEN
          PRINT '(A9,i4,A4,i5)', '   Point ', iE, ' of ', SEData1%NEPts
       END IF
       En = SEData1%EGrid(iE)
       CALL CSigZ(En, 0.d0, SEData1%Rs, Sigma0, Z, SEData1%Omega, SEData1%Width, SEData1%g, & 
                  SEData1%NPoles, OnShll, UseBP, SEData1%EGap)
       SEData1%SigDat(iE) = Z*(Sigma0 - SigmaF)
       SEData1%Z = Z
    END DO
  END SUBROUTINE SEnergy

!  SUBROUTINE SEConv(file, 
END MODULE SelfEnergyMod
