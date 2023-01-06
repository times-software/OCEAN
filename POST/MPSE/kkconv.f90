! kkconv takes a spectrum and self-energy as input and calculates
! the self-energy convolution and KK-transform simultaneously.
! Input:
! spect - i.e., Im[eps]
! omega - energy grid
! sigma - self-energy
! ne - number of energy points

! Output:
! cspect - complex spectrum with self-energy effects
SUBROUTINE kkconv(spect, omega, sigma, ne, cspect)
  IMPLICIT NONE
  ! Input
  REAL(8), INTENT(IN) :: spect(ne), omega(ne)
  COMPLEX*16, INTENT(IN) :: sigma(ne)
  INTEGER ne
  ! Output
  COMPLEX*16 cspect(ne)

  ! Local
  REAL(8), ALLOCATABLE :: omegaInt(:), spectInt(:)
  COMPLEX*16, ALLOCATABLE :: sigmaInt(:)
  REAL(8) gam0, dele
  INTEGER ie, ieInt, neInt
  COMPLEX*16, PARAMETER :: coni = (0.d0,1.d0)
  REAL(8), PARAMETER :: pi = 3.14159265358979323846264d0
  
  ! Get interpolation grid.
  gam0 = 1.d30
  DO ie = 1, ne
     IF(ABS(DIMAG(sigma(ie))).LT.gam0) gam0 = ABS(DIMAG(sigma(ie)))
  END DO
  gam0 = MAX(gam0,0.005d0)
  dele = gam0/20.d0
  neInt = FLOOR((omega(ne) - omega(1))/dele) + 1
  ALLOCATE(omegaInt(2*neInt), spectInt(2*neInt), sigmaInt(2*neInt))

  DO ie = 1, 2*neInt
     omegaInt(ie) = DBLE(ie - 1)*dele
     IF(ie.LE.neInt) THEN
        CALL terp(omega, spect, ne, 3, omegaInt(ie), spectInt(ie))
        CALL terpc(omega, sigma, ne, 3, omegaInt(ie), sigmaInt(ie))
     ELSE
        spectInt(ie) = spectInt(neInt)*omegaInt(neInt)**2/omegaInt(ie)**2
        sigmaInt(ie) = sigmaInt(neInt)
     END IF
     IF(gam0.GT.ABS(DIMAG(sigmaInt(ie)))) sigmaInt(ie) = DBLE(sigmaInt(ie)) + coni*gam0
  END DO


  DO ie = 1, ne
     IF(MOD(ie,100).EQ.0) PRINT*, 'Point ', ie, ' of ', ne
     cspect(ie) = 0.d0
     DO ieInt = 1, 2*neInt - 1
        cspect(ie) = cspect(ie) + 0.5d0*( spectInt(ieInt) / &
             & ( omega(ie) - omegaInt(ieInt) - sigmaInt(ieInt) ) + &
             & spectInt(ieInt+1) / ( omega(ie) - omegaInt(ieInt+1) - sigmaInt(ieInt+1) ))*dele
        cspect(ie) = cspect(ie) - 0.5d0*( spectInt(ieInt) / &
             & ( omega(ie) + omegaInt(ieInt) + CONJG(sigmaInt(ieInt)) ) + &
             & spectInt(ieInt+1) / ( omega(ie) + omegaInt(ieInt+1) + CONJG(sigmaInt(ieInt+1)) ))*dele
     END DO
  END DO
  cspect(:) = -CONJG(cspect(:))*1.d0/pi
END SUBROUTINE kkconv
  
     
     
