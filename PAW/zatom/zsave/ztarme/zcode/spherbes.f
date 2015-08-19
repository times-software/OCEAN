      function spherbes( arg, l )
      implicit none
      integer l
      double precision spherbes, arg
      double precision, external :: j0, j1, j2, j3, j4
      double precision, external :: js, jp, jd, jf, jg
      if ( arg .lt. 0.01d0 ) then
         select case( l )
         case( 0 )
           spherbes = js( arg )
         case( 1 )
           spherbes = jp( arg )
         case( 2 )
           spherbes = jd( arg )
         case( 3 )
           spherbes = jf( arg )
         case( 4 )
           spherbes = jg( arg )
         end select
      else
         select case( l )
         case( 0 )
           spherbes = j0( arg )
         case( 1 )
           spherbes = j1( arg )
         case( 2 )
           spherbes = j2( arg )
         case( 3 )
           spherbes = j3( arg )
         case( 4 )
           spherbes = j4( arg )
         end select
      end if
      end
c---------------------------------------
      function j0(x) result(f)
      implicit none
      double precision x, f
      f = dSin(x)/x
      end
c_____________________________________________________________________
      function j1(x) result(f)
      implicit none
      double precision x, f
      f = (dSin(x)/x**2) - dCos(x)/x
      end
c_____________________________________________________________________
      function j2(x) result(f)
      implicit none
      double precision x, f
      f = (3.d0/x**3 - 1.d0/x)*dSin(x) - (3.d0/x**2)*dCos(x)
      end
c_____________________________________________________________________
      function j3(x) result(f)
      implicit none
      double precision x, f
      f = (15.d0/x**4 - 6.d0/x**2)*dSin(x) -
     &    (15.d0/x**3 - 1.d0/x)*dCos(x)
      end
c_____________________________________________________________________
      function j4(x) result(f)
      implicit none
      double precision x, f
      f = (105.d0/x**5 - 45.d0/x**3 + 1.d0/x)*dSin(x) - 
     &    (105.d0/x**4 - 10.d0/x**2)*dCos(x)
      end
c_____________________________________________________________________
      function js(x) result(f)
      implicit none
      double precision x, f
      f = 1.d0 - x**2/6.d0 + x**4/120.d0 - x**6/5040.d0
      end
c_____________________________________________________________________
      function jp(x) result(f)
      implicit none
      double precision x, f
      f = x/3.d0 - x**3/30.d0 + x**5/840.d0 - x**7/45360.d0 +
     &    x**9/3991680.d0 - x**11/51891840.d0
      end
c_____________________________________________________________________
      function jd(x) result(f)
      implicit none
      double precision x, f
      f = x**2.d0/15.d0 - x**4.d0/210.d0 + x**6.d0/7500.d0 -
     &    x**8.d0/488960.d0 + x**10.d0/51891840.d0
      end
c_____________________________________________________________________
      function jf(x) result(f)
      implicit none
      double precision x, f
      f = x**3/105.d0 - x**5/1890.d0 + x**7/83160.d0 - 
     &    x**9/6486480.d0 + x**11/778377600.d0
      end
c____________________________________________________________________
      function jg(x) result(f)
      implicit none
      double precision x, f
      f = x**4/945.d0 - x**6/20790.d0 + x**8/1081080.d0 -
     &    x**10/97297200.d0
      end
