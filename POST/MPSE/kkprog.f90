PROGRAM kkprog
  implicit none
  
  REAL(8),ALLOCATABLE :: respect(:), spect(:), omega(:)
  REAL(8) omega0, delta, gam, tmp, eps1, eps2, loss, rad, theta, &
       &  indref, indabs, ref, reflct, mu, eloss
  COMPLEX*16 refrac
  !     fine structure alpha
  REAL(8), PARAMETER :: alpinv = 137.03598956d0
  REAL(8), PARAMETER :: alphfs = 1 / alpinv
  REAL(8), PARAMETER :: bohr = 0.529177249d0
  REAL(8), PARAMETER :: ryd = 13.605698d0
  REAL(8), PARAMETER :: hart = 2 * ryd
  COMPLEX*16, ALLOCATABLE :: sigma(:), cspect(:)
  
  INTEGER ie, ne, nsig, epsunit, sigunit, epsout, indsout, opout, refout, eps1out, eps2out, lossout
  
  epsunit = 13
  sigunit = 14
  eps1out = 15
  eps2out = 16
  lossout = 17
  indsout = 18
  opout   = 19
  refout  = 20
  gam = 2.d0
  ! Read spectrum
  OPEN(epsunit, FILE='eps2', STATUS='OLD')
  ie = 0
  CALL RdCmt(epsunit, '#!cC') 
  DO 
     READ(epsunit,*,END=5)
     ie = ie + 1
  END DO
5 CONTINUE
  ne = ie - 1
  ALLOCATE(omega(ne),spect(ne),cspect(ne),sigma(ne),respect(ne))
  respect(:) = 0.d0
  REWIND(epsunit)
  CALL RdCmt(epsunit, '#!cC') 
  DO ie = 1, ne
     READ(epsunit,*) omega(ie), spect(ie)     
  END DO
  CLOSE(epsunit)

  ! Read self-energy.
  OPEN(sigunit, FILE='SE.dat', STATUS='OLD')
  CALL RdCmt(sigunit,'#!cC')
  ie = 0
  DO 
     READ(sigunit,*,END=10)
     ie = ie + 1
  END DO
10 CONTINUE
  IF(ne.NE.(ie - 1)) THEN
     PRINT*, 'ERROR: Number of energy ' // &
          & 'points in SE.dat do not match those in eps2.', ne, ie-1
     STOP
  END IF
  REWIND(sigunit)
  CALL RdCmt(sigunit,'#!cC')
  DO ie = 1, ne
     READ(sigunit,*) omega0, tmp, delta, gam
     gam = MAX(gam,0.01)
     IF(omega0.EQ.omega(ie)) THEN
        sigma(ie) = delta + (0.d0,1.d0)*gam
     ELSE
        PRINT*, 'ERROR: Energy grid in SE.dat does not match ' // &
             & 'that in eps2'
        STOP
     END IF
  END DO
  CLOSE(sigunit)

  IF(.FALSE.) THEN
     ! Test with lorentzian
     gam = 0.5d0
     delta = 0.0d0
     omega0 = 10.d0
     DO ie = 1, ne
        omega(ie) = 0.1d0*(ie-1)
        spect(ie) = gam/((omega(ie) - omega0)**2 + gam**2) - gam/((omega(ie) + omega0)**2 + gam**2)
        respect(ie) = -(omega(ie)-omega0)/((omega(ie) - omega0)**2 + gam**2) + &
                      (omega(ie) + omega0)/((omega(ie) + omega0)**2 + gam**2)
        sigma(ie) = 8.d0 + (0.d0,1.d0) * 0.01d0
     END DO
  END IF

  CALL kkconv(spect, omega, sigma, ne, cspect)

  OPEN(eps1out,FILE='eps1_SE',STATUS='REPLACE')
  OPEN(eps2out,FILE='eps2_SE',STATUS='REPLACE')
  OPEN(lossout,FILE='loss_SE',STATUS='REPLACE')
  OPEN(refout,FILE='refl_SE',STATUS='REPLACE')
  OPEN(indsout,FILE='inds_SE',STATUS='REPLACE')
  OPEN(opout,FILE='opcons_SE.dat',STATUS='REPLACE')
  WRITE(opout,'(a)') "#   omega (eV)      epsilon_1       epsilon_2       n"// &
       & "               kappa           mu (cm^(-1))    R"// &
       & "               epsinv"
  DO ie = 1, ne
     eps1 = DBLE(cspect(ie)) + 1.d0
     eps2 = DIMAG(cspect(ie))
     loss = eps2/(eps1**2 + eps2**2)
     rad = SQRT( eps1**2 + eps2**2 )
     theta = ACOS( eps1/rad )/2.d0
     indref = SQRT(rad)*COS(theta)
     indabs = SQRT(rad)*SIN(theta)
     ref = ( (indref - 1)**2 + indabs**2 ) / &
          & ( (indref + 1)**2 + indabs**2) 
     refrac = SQRT(cspect(ie) + 1.d0)
     reflct = ABS((refrac - 1)/(refrac + 1))**2
     mu = 2*omega(ie)/hart*alphfs*DIMAG(refrac)/bohr*1000.d0
     WRITE(eps1out,'(20f20.10)'), omega(ie), eps1
     WRITE(eps2out,'(20f20.10)'), omega(ie), eps2
     WRITE(lossout,'(20f20.10)'), omega(ie), loss
     WRITE(refout,'(20f20.10)'), omega(ie), ref
     WRITE(indsout,'(20e20.10)'), omega(ie), indref, indabs, &
          & indref**2 - indabs**2, 2*indref*indabs     
     WRITE(opout,"(19e16.6)") omega(ie), cspect(ie), refrac-1.d0, mu, reflct, loss
  END DO
  
  CLOSE(eps1out)
  CLOSE(eps2out)
  CLOSE(lossout)
  CLOSE(refout)
  CLOSE(indsout)
  CLOSE(opout)

END PROGRAM kkprog
