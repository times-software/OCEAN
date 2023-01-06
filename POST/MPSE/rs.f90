PROGRAM CalcRs
  DOUBLE PRECISION Rs, CellVol, Pi, NEl, rho
  PARAMETER (Pi = 3.1415926535897932384626433d0)
  READ*, NEl
  READ*, CellVol
  rho = NEl/CellVol
  
  Rs = (3.d0/(4.d0*Pi*rho))**(1.d0/3.d0)
  PRINT*, Rs
END PROGRAM CalcRs
