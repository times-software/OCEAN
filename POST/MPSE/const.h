      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
!     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.03598956d0)
!     fine structure alpha
      parameter (alphfs = 1 / alpinv)
