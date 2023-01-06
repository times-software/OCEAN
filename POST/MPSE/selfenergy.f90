PROGRAM selfenergy
  USE Convolution
  USE SelfEnergyMod

  TYPE(SpecData) SpecDat
  TYPE(SEInput)  SEInput1
  TYPE(SEData)   SEData1
  REAL(KIND=8),ALLOCATABLE :: Gam(:), Dele(:), EGrid(:)
  REAL(KIND=8) Gam0, Del0, v1(3), v2(3), v3(3)
  INTEGER iE, iSC, nSC
  CALL ReadConvInp(SpecDat)
  CALL ReadSpectrum(SpecDat)
  ! Set SE input.
  SEInput1%RunMode = 0
  SEInput1%Eps0 = 0.d0

  ! Read Rs out of exc.dat
  CALL ReadData('exc.dat', Double1=SEData1%Rs)
  SEData1%EGap = SpecDat%EGap

  ! Read NPoles out of exc.dat and allocate space
  CALL ReadData('exc.dat', Int1 = SEData1%NPoles)
  ALLOCATE(SEData1%Omega(SEData1%NPoles), SEData1%Width(SEData1%NPoles), SEData1%g(SEData1%NPoles), STAT = iError)    
  CALL CheckAllocation(iError,'Allocation error occured during allocation of pole data arrays.')

  ! Read pole data.
  CALL ReadArrayData('exc.dat', Double1 = SEData1%Omega, Double2 = SEData1%Width, Double3 = SEData1%g)    
  
  ! Close file.
  CALL CloseFl('exc.dat')
  
  ! Unit conversion from eV to Hartree.
  SEData1%EGap = SEData1%EGap/hart
  SEData1%Omega(:) = SEData1%Omega(:)/hart
  SEData1%Width(:) = SEData1%Width(:)/hart

  ! Set up energy grid for self energy calculation.
  SEData1%NEPts = SpecDat%nData

  ! Deallocate egrid, reallocate and set to correct grid.
  ALLOCATE(EGrid(SEData1%NEPts),SEData1%EGrid(SEData1%NEPts),SEData1%SigDat(SEData1%NEPts),SEData1%Z(SEData1%NEPts))
  ALLOCATE(Gam(SpecDat%nData),Dele(SpecDat%nData))
  Gam(:) = 0.d0
  Dele(:) = 0.d0
  nSC = 1
  DO iSC = 1, nSC
     SEData1%EGrid(:) = (SpecDat%Data(:,SpecDat%iX) - SpecDat%FermiLevel + Dele(:))/hart
     
     ! Calculate the self energy.
     CALL SEnergy(SEData1, SEInput1)

     ! Set Gam and Dele and convolve.
     DO iE = 1, SpecDat%nData
        Gam(iE)  = -1.d0*DIMAG(SEData1%SigDat(iE))*hart + SpecDat%ConstBr
        Dele(ie) = 1.d0*DBLE(SEData1%SigDat(iE))*hart
!        IF(EGrid(iE).lt.0.d0) THEN
!           Gam(iE) = 0.d0
!           Dele(ie) = 0.d0
!        END IF
!        IF(Gam(iE).lt.0.d0) Gam(iE) = 0.d0
     END DO

     ! Change EGrid to eV
     DO iE = 1, SEData1%NEPts
        EGrid(iE) = DBLE(SEData1%EGrid(iE)*hart)
     END DO

     ! Interpolate to find Sigma(EGap)
     CALL terp(EGrid,Dele,SEData1%NEPts,1,-0.1d0,Del0)
     DelE(:) = DelE(:) - Del0
     ! Find first point above E Fermi.
     DO iE = 1, SEData1%NEPts
        IF(EGrid(iE).gt.0.d0) THEN
           DelE(1:iE) = DelE(iE)
           Gam(1:iE) = Gam(iE)
           EXIT
        END IF
     END DO
!     CALL terp(EGrid,Dele,SEData1%NEPts,1,0.2d0,Del0)
!     DO iE = 1, SEData1%NEPts 
!        IF(EGrid(iE).lt.0.d0) THEN
!           DelE(iE) = Del0
!           Gam(iE) = SpecDat%ConstBr
!        END IF
!     END DO
     !PRINT*, Del0
  END DO
! Test convolution with constant broadening.
!   DO iE = 1, SpecDat%nData
!     IF(DBLE(SEData1%EGrid(iE)).ge.0.d0) THEN
!      IF(.TRUE.) THEN
!         Gam(iE) = 2.5d0
!         Dele(iE) = 0.d0
!      ELSE
!         Gam(iE) = 0.d0
!         DelE(iE) = 0.d0
!      END IF     
!   END DO
  ! Interpolate the spectrum onto an even grid with padding at the ends.
!  CALL InterpolateSpectrum(SpecDat, 0.1d0)

  CALL WriteArrayData('SE.dat', Double1 = SpecDat%Data(:,SpecDat%iX), Double2 = EGrid(:), Double3 = Dele(:), Double4 = Gam(:))

!  CALL ConvolveSpectrum(SpecDat,Gam,Dele)
!  CALL Write2D(SpecDat%OutFile,SpecDat%Data)

END PROGRAM selfenergy
