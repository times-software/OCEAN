PROGRAM CellVolume
  ! Find the volume of a unit cell (parallelopiped).
  DOUBLE PRECISION v1(3), v2(3), v3(3), CellVol
  INTEGER i1, i2, i3

  READ*, v1
  READ*, v2
  READ*, v3

  ! (v1 cross v2) dot v3
  CellVol = (v1(2)*v2(3) - v2(2)*v1(3))*v3(1) + &
       &       (v1(3)*v2(1) - v2(3)*v1(1))*v3(2) + &
       &       (v1(1)*v2(2) - v2(1)*v1(2))*v3(3) 

  PRINT*, ABS(CellVol)
END PROGRAM CellVolume 
