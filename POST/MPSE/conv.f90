MODULE Convolution
  USE IOMod
  USE Constants

  IMPLICIT NONE

  TYPE SpecData
     ! Columns to read.
     INTEGER iX, nY
     INTEGER,ALLOCATABLE :: iY(:)
     ! Spectrum file.
     CHARACTER(300) SpectrumFile, OutFile
     INTEGER NData, nDatInt
     REAL(KIND=8) FermiLevel, EGap, ConstBr
     REAL(KIND=8),ALLOCATABLE :: Data(:,:), DataInt(:,:)
  END TYPE SpecData

CONTAINS
    
  ! Read input file.
  SUBROUTINE ReadConvInp(SpecDat)
    TYPE(SpecData),INTENT(OUT) :: SpecDat
    INTEGER iAllocateError

    ! Read spectrum file (file with data to convolve).
    CALL ReadData('seconv.inp',String1=SpecDat%SpectrumFile, CommentCharacters = '#!*')
    ! Read output file
    CALL ReadData('seconv.inp',String1=SpecDat%OutFile, CommentCharacters = '#!*')
    ! Read number of y columns.
    CALL ReadData('seconv.inp',Int1=SpecDat%nY, CommentCharacters = '#!*')
    ALLOCATE(SpecDat%iY(SpecDat%nY), STAT = iAllocateError)
    CALL CheckAllocation(iAllocateError, 'Error: Failed to allocate space for SpecDat%iY' //&
         & ' in subroutine ReadConvInp.')

    ! Read columns to use for convolution. One x column
    ! and up to 9 y columns.
    CALL ReadData('seconv.inp',Int1=SpecDat%iX, CommentCharacters = '#!*')
    CALL ReadArrayData('seconv.inp',Int1=SpecDat%iY, CommentCharacters = '#!*')
    CALL ReadData('seconv.inp', Double1 = SpecDat%FermiLevel, Double2 = SpecDat%EGap, Double3 = SpecDat%ConstBr,&
         & CommentCharacters = '#!*')
  END SUBROUTINE ReadConvInp

  SUBROUTINE ReadSpectrum(SpecDat)
    USE IOFiles
    TYPE(SpecData),INTENT(INOUT) :: SpecDat
    INTEGER iUnit, iAllocateError
    ! Read data from file
    CALL OpenFl(SpecDat%SpectrumFile)
    CALL GetIOFileInfo(SpecDat%SpectrumFile,UnitNumber=iUnit)
    SpecDat%nData = NumberOfLines(SpecDat%SpectrumFile)

    ALLOCATE(SpecDat%Data(SpecDat%nData,SpecDat%iY(SpecDat%nY)), STAT=iAllocateError)

    CALL CheckAllocation(iAllocateError, 'Error: Failed to allocate space for SpecDat%Data' //&
         & ' in subroutine ReadSpectrum.')
    
    CALL Read2D(SpecDat%SpectrumFile, SpecDat%Data, SpecDat%nData, SpecDat%iY(SpecDat%nY))
  END SUBROUTINE ReadSpectrum

  SUBROUTINE InterpolateSpectrum(SpecDat,Delta0)
    TYPE(SpecData),INTENT(INOUT) :: SpecDat
    REAL(KIND=8),INTENT(IN) :: Delta0
    REAL(KIND=8) Delta, DeltaMin
    INTEGER i1, i2, i3, nExtra, iAllocateError

    ! Find smallest interval in data.
    DeltaMin = ABS(SpecDat%Data(2,SpecDat%iX)-SpecDat%Data(1,SpecDat%iX))
    DO i1 = 2, SpecDat%nData - 1
       Delta = ABS(SpecDat%Data(i1+1,SpecDat%iX)-SpecDat%Data(i1,SpecDat%iX))
       IF(Delta.lt.DeltaMin) DeltaMin = Delta
    END DO
    IF(Delta0.lt.DeltaMin) DeltaMin = Delta0
    IF(DeltaMin.lt.0.001) DeltaMin = Delta0 
    DeltaMin = DeltaMin/5
    SpecDat%nDatInt = INT((SpecDat%Data(SpecDat%nData,SpecDat%iX) - SpecDat%Data(1,SpecDat%iX))/DeltaMin) + 1

    ! Allocate space for interpolated array.
    ALLOCATE(SpecDat%DataInt(SpecDat%nDatInt, SpecDat%nY + 1), STAT=iAllocateError)
    CALL CheckAllocation(iAllocateError, 'Error: Failed to allocate space for SpecDat%DataInt' //&
         & ' in subroutine InterpolateSpectrum.' )

    ! Create the grid.
    DO i1 = 1, SpecDat%nDatInt
       SpecDat%DataInt(i1,1) = SpecDat%Data(1,SpecDat%iX) + (i1-1)*DeltaMin
    END DO
    IF(.FALSE.) THEN
      DO i1 = 1, SpecDat%nDatInt
         SpecDat%DataInt(i1,1) = SpecDat%Data(1,SpecDat%iX) + (i1-INT(SpecDat%nDatInt/3))*DeltaMin
      END DO

    ! Pad the beginning of the data. 
      DO i1 = 1, SpecDat%nY
         DO i2 = 1, INT(SpecDat%nDatInt/3)
            SpecDat%DataInt(i2,i1+1) = SpecDat%Data(i1,SpecDat%iY(i1))
         END DO
      END DO

    ! Interpolate.
      Loop1: DO i1 = INT(SpecDat%nDatInt/3) + 1, SpecDat%nDatInt
         Loop2: DO i2 = 1, SpecDat%nY
            IF(SpecDat%DataInt(i1,1).gt.SpecDat%Data(SpecDat%nData,SpecDat%iX)) EXIT Loop1
            CALL terp(SpecDat%Data(:,SpecDat%iX), SpecDat%Data(:,SpecDat%iY(i2)), &
                 & SpecDat%nData, 1, SpecDat%DataInt(i1,1), SpecDat%DataInt(i1,i2+1))
         END DO Loop2
      END DO Loop1
    
    ! Pad the end of the data.
      DO i2 = i1, SpecDat%nDatInt
         DO i3 = 1, SpecDat%nY
            SpecDat%DataInt(i2,i3+1) = SpecDat%Data(SpecDat%nData,SpecDat%iY(i3))
         END DO
      END DO
    ELSE 
      Loop3: DO i1 = 1, SpecDat%nDatInt
         Loop4: DO i2 = 1, SpecDat%nY
            IF(SpecDat%DataInt(i1,1).gt.SpecDat%Data(SpecDat%nData,SpecDat%iX)) EXIT Loop3
            CALL terp(SpecDat%Data(:,SpecDat%iX), SpecDat%Data(:,SpecDat%iY(i2)), &
                 & SpecDat%nData, 1, SpecDat%DataInt(i1,1), SpecDat%DataInt(i1,i2+1))
         END DO Loop4
      END DO Loop3
    END IF
  END SUBROUTINE InterpolateSpectrum

  SUBROUTINE ConvolveSpectrum(SpecDat,Gam,Dele)
    TYPE(SpecData),INTENT(INOUT) :: SpecDat
    REAL(KIND=8) TmpData(SpecDat%nData,SpecDat%nY), TmpIntDat(SpecDat%nDatInt,SpecDat%nY)
    REAL(KIND=8) ConvArray(SpecDat%nDatInt,SpecDat%nY), Norm, Gam(SpecDat%nData), Sum1, Sum2, &
         & DelE(SpecDat%nData), tol, Gam0, Del0
    INTEGER i1, i2, i3, ifirst
    ifirst=0.d0
    tol = 2.d0*(SpecDat%DataInt(2,1)-SpecDat%DataInt(1,1))  
    PRINT '(A)', 'Beginning convolution.'
    DO i1 = 1, SpecDat%nY
       DO i2 = 1, SpecDat%nData
          !IF(SpecDat%Data(i2,1).lt.SpecDat%FermiLevel) CYCLE
          IF(MOD(i2,100).eq.0) THEN
             PRINT '(A9,I5,A4,I5)', '   Point ', (i1-1)*SpecDat%nData + i2, ' of ', SpecDat%nY*SpecDat%nData
          END IF
          ! Fill ConvArray with broadening function.
          IF(Gam(i2).lt.tol) THEN
             ! Use only real part
             IF((SpecDat%Data(i2,1)-Dele(i2)).lt.SpecDat%DataInt(1,1)) THEN
                SpecDat%Data(i2,SpecDat%iY(i1)) = SpecDat%Data(1,SpecDat%iY(i1))
             ELSE
                CALL terp(SpecDat%DataInt(:,1), SpecDat%DataInt(:,i1+1), &
                   & SpecDat%nDatInt, 1, SpecDat%Data(i2,1)-Dele(i2), SpecDat%Data(i2,SpecDat%iY(i1)))
             END IF
             CYCLE
          END IF

          IF(.FALSE.) THEN
             ! Gaussian broadening.
             Norm = 1.d0/(SQRT(pi)*2.d0*Gam(i2))
             DO i3 = 1, SpecDat%nDatInt
                ConvArray(i3,i1) = Norm*EXP(-(SpecDat%Data(i2,SpecDat%iX)-SpecDat%DataInt(i3,1))**2/(2.d0*Gam(i2))**2) * &
                     & SpecDat%DataInt(i3,i1+1)
             END DO
          ELSE
             ! Lorenzian broadening.
             ! Endpoint corrections
!             Sum1 = (0.5d0 - 1.d0/pi*ATAN((SpecDat%Data(i2,SpecDat%iX)-SpecDat%DataInt(1,1)-DelE(i2))/Gam(i2)))*SpecDat%DataInt(1,i1+1)
!             Sum1 = Sum1 + (0.5d0 + 1.d0/pi*ATAN((SpecDat%Data(i2,SpecDat%iX)-SpecDat%DataInt(SpecDat%nDatInt,1)-DelE(i2))/Gam(i2)))* &
!                  & SpecDat%DataInt(SpecDat%nDatInt,i1+1)
             
!             Norm = Gam(i2)/pi
             DO i3 = 1, SpecDat%nDatInt
                !IF(SpecDat%DataInt(i3,1).le.SpecDat%FermiLevel) THEN
                !   ConvArray(i3,i1) = 0.d0
                !   CYCLE
                !END IF
                ! Interpolate self-energy
                CALL terp(SpecDat%Data(1,SpecDat%iX), Gam, &
                     & SpecDat%nData, 1, SpecDat%DataInt(i3,1), Gam0)
                IF(Gam0.lt.tol) Gam0 = tol
                CALL terp(SpecDat%Data(1,SpecDat%iX), Dele, &
                     & SpecDat%nData, 1, SpecDat%DataInt(i3,1), Del0)

                !Norm = Gam0/pi
                Norm = Gam0/(-ATAN((SpecDat%FermiLevel-SpecDat%DataInt(i3,1)-Del0)/Gam0) + &
                     & ATAN((SpecDat%DataInt(i3,1)+Del0)/Gam0))
                IF(.FALSE.) THEN
                   CALL WriteData('SpectFn.dat', Double1 = SpecDat%DataInt(i3,1), Double2 = &
                     & Norm*(1.d0/((SpecDat%Data(i2,SpecDat%iX)-SpecDat%DataInt(i3,1)-Del0)**2 + Gam0**2) - &
                     &                   1.d0/((SpecDat%Data(i2,SpecDat%iX)+SpecDat%DataInt(i3,1)+Del0-SpecDat%FermiLevel)**2 &
                     & + Gam0**2)))
                END IF
                IF(gam0.gt.0.d0) THEN
                   ConvArray(i3,i1) = Norm*(1.d0/((SpecDat%Data(i2,SpecDat%iX)-SpecDat%DataInt(i3,1)-Del0)**2 + Gam0**2) - &
                        &                   1.d0/((SpecDat%Data(i2,SpecDat%iX)+SpecDat%DataInt(i3,1)+Del0-SpecDat%FermiLevel)**2 &
                        & + Gam0**2)) * SpecDat%DataInt(i3,i1+1)
                ELSE IF((SpecDat%Data(i2,SpecDat%iX)-SpecDat%DataInt(i3,1)-Del0).ne.0.d0) THEN
                   ConvArray(i3,i1) = 0.d0
                ELSE
                   ConvArray(i3,i1) = SpecDat%DataInt(i3,i1+1)/(SpecDat%DataInt(2,1)-SpecDat%DataInt(1,1))
                END IF
                  
             END DO
          END IF
          IF(.FALSE.) CALL WriteData('SpectFn.dat', String1 = ' ')
          ifirst=1
          ! Do integration
          CALL trap(SpecDat%DataInt(:,1),ConvArray(:,i1),SpecDat%nDatInt,Sum2)
          SpecDat%Data(i2,SpecDat%iY(i1)) = Sum2
       END DO
    END DO

    ! Shift by real part.
!    SpecDat%Data(:,SpecDat%iX) = SpecDat%Data(:,SpecDat%iX) + DelE(:)             
       
    ! Interpolate onto original grid.
!     DO i1 = 1, SpecDat%nY
!        DO i2 = 1, SpecDat%nData
!           CALL terp(SpecDat%DataInt(:,1), ConvArray(:,i1), &
!                & SpecDat%nDatInt, 1, SpecDat%Data(i2,SpecDat%iX), SpecDat%Data(i2,i1+1))
!        END DO
!     END DO
  END SUBROUTINE ConvolveSpectrum
          
END MODULE Convolution

