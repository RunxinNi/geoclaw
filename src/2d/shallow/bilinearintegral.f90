!======================================================================
function bilinearintegral(domain, topoparam, corners)
!======================================================================
      
   use topo_module, only: intersection

   !i/o
   real(kind=8) :: domain(1:4), topoparam(1:4), corners(1:2,1:2)
   real(kind=8) :: bilinearintegral

   !local
   integer :: indicator
   real(kind=8) :: area, sumxi, sumeta, a, b, c, d, dx, dy
   real(kind=8) :: bound(1:4)

!#######################################################################

   !bilinearintegral integrates the bilinear with values z##
   !over the rectangular region xim <= x <= xip, and
                                     !yjm <= y <= yjp

                                          !written by David L George
                                          !Vancouver, WA April 2010
!#######################################################################



   !integrate the portion of the bilinear intersected with the
   !rectangular cell analytically.

   !find limits of integral (this should already be true?)
   call intersection(indicator, area, bound, domain, topoparam)
      

   !find the area of integration
   dx = topoparam(2) - topoparam(1)
   dy = topoparam(4) - topoparam(3)
   sumxi = (bound(1) + bound(2) - 2.d0*topoparam(1))/dx
   sumeta = (bound(3) + bound(4) - 2.d0*topoparam(3))/dy

   !find coefficients of bilinear a*xi + b*eta + c*xi*eta + d
   a = corners(2,1) - corners(1,1)
   b = corners(1,2) - corners(1,1)
   c = corners(2,2) - corners(2,1) - corners(1,2) + corners(1,1)
   d = corners(1,1)

   bilinearintegral = (0.5d0*(a*sumxi + b*sumeta) + 0.25d0*c*sumxi*sumeta + d)*area

   return

end function


!======================================================================
function bilinearintegral_s(domain, topoparam, corners)
!======================================================================

   use topo_module, only: intersection
   use geoclaw_module, only: r2d => rad2deg, d2r => deg2rad
   use geoclaw_module, only: Rearth => earth_radius

   !i/o
   real(kind=8) :: domain(1:4), topoparam(1:4), corners(1:2,1:2)
   real(kind=8) :: bilinearintegral_s

   !local
   integer :: indicator
   real(kind=8) :: bound(1:4), delx, dx, dy, a, b, c, d, area
   real(kind=8) :: xdiffhi, xdifflow, ydiffhi, ydifflow, xdiff2
   real(kind=8) :: adsinint, cbsinint

!#######################################################################

         !bilinearintegral integrates the bilinear with values z##
         !over the rectangular region xim <= x <= xip, and
                                     !yjm <= y <= yjp
         !integration is actually done on the surface of a sphere

                                          !written by David L George
                                          !Vancouver, WA April 2010
!#######################################################################



   !integrate the portion of the bilinear intersected with the
   !rectangular cell analytically.
   !find limits of integral (this should already be true?)
   call intersection(indicator, area, bound, domain, topoparam)

   !find terms for the integration
   xdiffhi = bound(2) - topoparam(1)
   xdifflow = bound(1) - topoparam(1)
   ydiffhi = bound(4) - topoparam(3)
   ydifflow = bound(3) - topoparam(3)
   xdiff2 = 0.5d0*(xdiffhi**2 - xdifflow**2)

   cbsinint = (r2d*cos(d2r*bound(4)) + ydiffhi*sin(d2r*bound(4))) - &
              (r2d*cos(d2r*bound(3)) + ydifflow*sin(d2r*bound(3)))

   adsinint = r2d*(sin(d2r*bound(4)) - sin(d2r*bound(3)))


   !find coefficients of bilinear a*xi + b*eta + c*xi*eta + d
   dx = topoparam(2) - topoparam(1)
   dy = topoparam(4) - topoparam(3)
   delx = bound(2) - bound(1)
   a = (corners(2,1) - corners(1,1))/dx
   b = (corners(1,2) - corners(1,1))/dy
   c = (corners(2,2) - corners(2,1) - corners(1,2) + corners(1,1))/ &
       (dx*dy)
   d = corners(1,1)

   bilinearintegral_s = ((a*xdiff2 + d*delx)*adsinint + &
                         r2d*(c*xdiff2 + b*delx)*cbsinint)*(Rearth*d2r)**2

   return

end function
