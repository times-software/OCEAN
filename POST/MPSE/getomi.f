      SUBROUTINE getomi(wgrid, epsInv, NPts, NPoles, omi, gi, Deltai,
     &     eps0)
      IMPLICIT NONE
!     Input variables
!     wGrid(NPts)  - grid of omega
!     epsInv(NPts) - loss function
!     NPoles       - Number of poles to represent epsInv
      
!     Output
!     omi          - omega position of poles
!     gi           - pole strengths
!     Deltai       - energry region about omi where
      INTEGER NPoles, NPts
      DOUBLE PRECISION wGrid(NPts), epsInv(NPts), omi(NPoles),
     &     gi(NPoles), Deltai(NPoles), eps0

!     Local variables
      INTEGER maxPts, MaxPoles, iMax, i1, i2
      PARAMETER (maxPts = 50000)
      DOUBLE PRECISION M0, M0i, M1i, MM1i, M1Old, M1New, omNew, omOld,
     &     pi, Tmp1(maxPts), Tmp2(maxPts), IntGrl, a, b, c, frac, egap,
     &     M0test
      PARAMETER (pi = 3.14159 26535 89793 23846 26433d0)
      EXTERNAL IntGrl
!     Find zeroth moment
      frac = 1.d0-0.01d0/NPoles
      M0 = 0.5d0*epsInv(1)*wGrid(1)
      DO i1 = 2, NPts
         M0 = M0 + 0.5d0*(epsInv(i1)+epsInv(i1-1))*
     &        (wGrid(i1)-wGrid(i1-1))
         Tmp1(i1) = epsInv(i1)/wGrid(i1)
         Tmp2(i1) = epsInv(i1)*wGrid(i1)
      END DO
!     Find gap energy
!      DO i1 = 1, NPts
      
!     Now find energy regions such that the integral of epsinv in each region is
!     equal.
      omOld = 0.d0
      omNew = 0.d0
      M1Old = 0.d0
      M1New = 0.d0
     
      DO i1 = 1, NPoles
         M0i = 0.5d0*epsInv(1)*wGrid(1)
         M1New = 0.5d0*epsInv(1)*wGrid(1)**2
         M0test = DBLE(i1)*frac*M0/DBLE(NPoles)
         DO i2 = 2, NPts
            IF(M0i.ge.M0test) THEN
               omOld = omNew
               omNew = wGrid(i2-1)
               IF(i2.eq.1) THEN
                  M0i = 0.d0
               ELSEIF(i2.eq.2) THEN
                  M0i = 0.5d0*epsInv(1)*wGrid(1)
               ELSE
                  M0i = M0i - 0.5d0*(epsInv(i2-1)+ epsInv(i2-2))*
     $              (wGrid(i2-1)-wGrid(i2-2))
               END IF
               GOTO 5
            ELSE
               M0i = M0i + 0.5d0*(epsInv(i2)+epsInv(i2-1))*
     $              (wGrid(i2)-wGrid(i2-1))
               M1New = M1New + 0.5d0*(epsInv(i2)*wGrid(i2)+ epsInv(i2-1)
     $              *wGrid(i2-1))*(wGrid(i2)-wGrid(i2-1))
            END IF
         END DO
 5       CONTINUE

!        The region has been found up to the wGrid accuracy. Now interpolate
!        to find omi
         IF(i2.eq.2) THEN
            a = epsInv(i2-1)/omNew
         ELSE
            a = (epsInv(i2-1) + epsInv(i2-2))/(omNew - wGrid(i2-2))
            b = epsInv(i2-2)
         END IF

         c = -2.d0/a*(M0i - frac*M0*(DBLE(i1)/DBLE(NPoles)))
         omNew = wGrid(i2-2) + c
!         M0i = DBLE(i1)*frac*M0/DBLE(NPoles)
         IF((omNew.lt.wGrid(i2-2)).or.(omNew.gt.wGrid(i2-1))) THEN
            omNew = (wGrid(i2-2) + wGrid(i2-1))/2.d0
         END IF

!        Now set omi and g s.t. the inverse and first moments are preserved for
!        this region.
         MM1i = 0.d0
         M1i  = 0.d0
         IF(omOld.eq.0.d0) THEN
            MM1i = 0.5d0*epsInv(1)
            M1i  = 0.5d0*epsInv(1)*wGrid(1)**2
            omOld = wGrid(1)
         END IF
         MM1i = MM1i + IntGrl(omOld, omNew, wGrid, Tmp1, NPts)
         M1i  = M1i + IntGrl(omOld, omNew, wGrid, Tmp2, NPts)
         gi(i1)  = 2.d0/pi*MM1i
         omi(i1) = SQRT(M1i/MM1i)
         Deltai(i1) = omNew-omOld
         M1Old  = M1New
      END DO

!     Shift and scale the poles to set inverse moment to 1-1/eps0, without changing
!     first moment. (eps0 is input by user)
      CALL getdom(gi,omi,NPoles,eps0)
      
      RETURN
      END
