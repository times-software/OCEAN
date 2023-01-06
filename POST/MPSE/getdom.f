      SUBROUTINE getdom(gi,omi,NPoles,eps0)
      IMPLICIT NONE

!     Input/output
      INTEGER NPoles
      DOUBLE PRECISION gi(NPoles), omi(NPoles), eps0

!     Local
      INTEGER i1, i2, MxIter
      DOUBLE PRECISION MM1, tol, dom, dom1, dom2, sumg1, sumg2
      PARAMETER(MxIter = 10000, tol = 1.d-6)
      IF(eps0.lt.-1.d0) RETURN
      IF(eps0.lt.0.d0) THEN
         MM1 = 1.d0
      ELSE
         MM1   = 1.d0 - 1.d0/eps0
      END IF
      dom1  = -omi(1)*0.95d0
      dom2  = omi(1)*0.95d0
!     Bracket dom between dom1 and dom2
      DO i1 = 1, MxIter
         sumg1 = 0.d0
         sumg2 = 0.d0
         DO i2 = 1, NPoles
            sumg1 = sumg1 + gi(i2)*(omi(i2)/(omi(i2)-dom1))**2
            sumg2 = sumg2 + gi(i2)*(omi(i2)/(omi(i2)-dom2))**2
         END DO
         IF(sumg1.gt.MM1) dom1 = dom1*2.d0
         IF(sumg2.lt.MM1) dom2 = dom2/2.d0
         IF((sumg1.lt.MM1).and.(sumg2.gt.MM1)) GOTO 5
      END DO
 5    CONTINUE                  ! dom bracketed

      dom = 0.5d0*(dom1+dom2)
      DO i1 = 1, MxIter
         sumg1 = 0.d0
         DO i2 = 1, NPoles
            sumg1 = sumg1 + gi(i2)*(omi(i2)/(omi(i2)-dom))**2
         END DO
         IF(abs((sumg1-MM1)/MM1).lt.tol) GOTO 10
         IF((sumg1-MM1).gt.0.d0) THEN
            dom2 = dom
         ELSE
            dom1 = dom
         END IF
         dom = 0.5d0*(dom1+dom2)
      END DO
 10   CONTINUE

      DO i1 = 1, NPoles
         gi(i1)  = gi(i1)*(omi(i1)/(omi(i1)-dom))**2
         omi(i1) = omi(i1)-dom
      END DO
      RETURN
      END
      
