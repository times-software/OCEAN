      DOUBLE PRECISION FUNCTION IntGrl(a, b, x, y, n)
      INTEGER n, i, iStart, iLast
      DOUBLE PRECISION x(n), y(n), dx, a, b, aEnd, bEnd,
     &     ya, yb

      iStart = -1
      iLast = -1
      IntGrl = 0.d0
      dx = x(2) - x(1)
      
      DO i = 1, n
         IF(x(i).ge.a.and.x(i).lt.b) THEN
            iStart = i
            goto 5
         ELSE
            iStart = -1
         END IF
      END DO
 5    CONTINUE
      
      DO i = n, 1, -1
         IF(x(i).le.b.and.x(i).gt.a) THEN
            iLast = i
            goto 10
         ELSE
            iLast = -1
         END IF
      END DO
 10   CONTINUE

      IF (iStart.lt.0.and.iLast.lt.0) THEN
         CALL terp(x, y, n, 1, a, ya)
         CALL terp(x, y, n, 1, b, yb)
         IntGrl = IntGrl + 0.5d0*(yb+ya)*(b-a)
      ELSE
         DO i = iStart, iLast-1
            dx=x(i+1)-x(i)
            IntGrl = IntGrl + 0.5d0*(y(i)+y(i+1))*dx
         END DO
         
         CALL terp(x, y, n, 1, a, ya)
         CALL terp(x, y, n, 1, b, yb)
         
         IntGrl = IntGrl + 0.5d0*(ya+y(iStart))*(x(iStart) - a)
         IntGrl = IntGrl + 0.5d0*(yb+y(iLast))*(b-x(iLast))
      END IF
      RETURN
      END
