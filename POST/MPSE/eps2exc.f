      PROGRAM eps2exc
      DOUBLE PRECISION Dat(50000,2), g(1000), gamma, omi(1000),
     &     Delta(1000),eps0, sumrl, xNElec, csumrl
      INTEGER i,i1,ios1,ios2,imax, NPoles,NPts,NFile
      CHARACTER*80 infl, outfl, comment
      CHARACTER ch
      EXTERNAL IntGrl
      NPts = 1000
      infl='loss.dat'
      outfl='exc.dat'
      PRINT*, '# Enter number of poles:'
      READ*, NPoles
      PRINT*, '# Enter eps^-1 sumrule and N_el:'
      READ*, sumrl, xNElec
      PRINT*, '# Enter dielectric constant: '
      READ*, eps0
      
      gamma = 0.01
      csumrl= xNElec/sumrl
      
      OPEN(unit=12,file=infl,status='old')
      OPEN(unit=13,file=outfl,status='replace')
      CALL rdloss(12,Dat,NPts)

      DO i1 = 1, NPts
         Dat(i1,2) = Dat(i1,2)*csumrl
      END DO
      CALL getomi(Dat(1,1), Dat(1,2), NPts, NPoles, omi, g, Delta, eps0)
      WRITE(13,'(A33,I4,A6)') '# Loss function represented with ',
     $     NPoles, ' poles'
      WRITE(13,'(A23,f8.4)') '# Dielectric constant: ', eps0
      DO i = 1, NPoles
         WRITE(13,'(4f20.10)'), omi(i), omi(i)*gamma, g(i), Delta(i)
      END DO
      
      END
