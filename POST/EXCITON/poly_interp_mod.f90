module poly_interp_mod
!   use iso_fortran_env,  only : dp => real64
   implicit none
   private
   public  :: interp_poly_cut

   integer, parameter :: DP = kind(1.0d0 )
   real(dp), parameter :: tol = 1.0d-12        ! numerical tolerance

   !>>>>>>>>>  local helper routines  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
   pure function cross(a,b) result(c)
      real(dp), intent(in) :: a(3), b(3)
      real(dp)             :: c(3)
      c = [ a(2)*b(3)-a(3)*b(2),  &
             a(3)*b(1)-a(1)*b(3), &
             a(1)*b(2)-a(2)*b(1) ]
   end function cross
   !--------------------------------------------------------------------
   pure subroutine barycentric(p, v, w, inside)
      !! p … point (size 3) in param-space
      !! v … vertices of the tetra  v(4,3)
      !! w … barycentric weights
      real(dp), intent(in)  :: p(3), v(4,3)
      real(dp), intent(out) :: w(4)
      logical,  intent(out) :: inside
      real(dp) :: a1(3), a2(3), a3(3), rhs(3), detA

      a1  = v(2,:) - v(1,:)
      a2  = v(3,:) - v(1,:)
      a3  = v(4,:) - v(1,:)
      rhs = p       - v(1,:)

      detA = dot_product(a1, cross(a2,a3))
      if (abs(detA) < tol) then   ! degenerate tetra – should not occur
         inside = .false.;  w = 0._dp
         return
      end if

      ! Cramer’s rule
      w(2) =  dot_product(rhs, cross(a2,a3)) / detA
      w(3) =  dot_product(a1, cross(rhs,a3)) / detA
      w(4) =  dot_product(a1, cross(a2,rhs)) / detA
      w(1) = 1._dp - (w(2)+w(3)+w(4))

      inside = all(w >= -tol) .and. all(w <= 1._dp+tol)
   end subroutine barycentric
   !--------------------------------------------------------------------
   pure function interp_poly_cut(fff,u,v,w,cutoff) result(val)
      real(dp), intent(in) :: fff(8), u, v, w, cutoff
      real(dp)             :: val
      logical              :: active(8)
      integer              :: i,j,k,l
      real(dp), parameter  :: vtx(8,3)=reshape([ &
          0,0,0 , 0,0,1 , 0,1,0 , 0,1,1 , 1,0,0 , 1,0,1 , 1,1,0 , 1,1,1 ], [8,3])
      real(dp) :: p(3), verts(4,3), wts(4)
      logical  :: inside, found

      !--- always compute ordinary trilinear interpolation --------------
      val = tri_linear(fff,u,v,w)
      !--- check if the point lies in any all-active tetra --------------
      active = (fff > cutoff)
      if (.not. any(active .eqv. .false.)) return    ! all active ⇒ done


      p     = [u,v,w]
      found = .false.

      ! try every 4-tuple of the 8 vertices  (70 in total)
      do i=1,5
         do j=i+1,6
            do k=j+1,7
               do l=k+1,8
                  if (.not.(active(i) .and. active(j) .and. active(k) .and. active(l))) cycle
                  verts = vtx([i,j,k,l],:)
                  call barycentric(p,verts,wts,inside)
                  if (inside) then
                     found = .true.
                     exit
                  end if
               end do
               if (found) exit
            end do
            if (found) exit
         end do
         if (found) exit
      end do

      if (.not. found) val = 0._dp                    ! outside the active polyhedron
   end function interp_poly_cut

!---------------- ordinary trilinear interpolation --------------------
   pure function tri_linear(fff,u,v,w) result(val)
      real(dp), intent(in) :: fff(8), u, v, w      ! z-fast, y-mid, x-slow
      real(dp)             :: val
      real(dp) :: omu, omv, omw                    ! (1-u) etc.
      omu = 1.0_dp - u ;  omv = 1.0_dp - v ;  omw = 1.0_dp - w
      val =                                                           &
           fff(1)*omu*omv*omw + fff(2)*omu*omv*w  +                   &
           fff(3)*omu*v  *omw + fff(4)*omu*v  *w  +                   &
           fff(5)*u  *omv*omw + fff(6)*u  *omv*w  +                   &
           fff(7)*u  *v  *omw + fff(8)*u  *v  *w
   end function tri_linear
end module poly_interp_mod

