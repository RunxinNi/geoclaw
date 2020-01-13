!===========================================================================
function topointegral(domain, topoparam, mx, my, z, intmethod)
!===========================================================================

!###########################################################################
      !topointegral integrates a surface over a rectangular region
      !that is the intersection with a Cartesion grid (of topography data)
      !the surface integrated is defined by a piecewise bilinear through the
      !nodes of the Cartesian grid.

      !The rectangular intersection has coords:
      !xim <= x <= xip, yjm <= y <= yjp

      !The Cartesian grid has coords:
      !xxlow <= x <= xxhi, yylow <= y <= yyhi, with grid cell size dxx by dyy
      !and mxx by myy cells.

                                               !written by David L. George
                                                !Seattle, WA 7/16/08
!###########################################################################

   use geoclaw_module
   use topo_module, only: intersection
   
   !i/o   
   real(kind=8) :: domain(1:4), topoparam(1:6), z(1:mx,1:my)
   real(kind=8) :: topointegral
   
   !local
   real(kind=8) :: area, xim, xip, yjm, yjp
   real(kind=8) :: xlow, ylow, xhi, yhi, dx, dy
   real(kind=8) :: x1, x2, y1, y2, distart, djstart, diend, djend
   real(kind=8) :: theintegral, bilinearintegral, bilinearintegral_s
   real(kind=8) :: inttopoparam(1:4), intcorners(1:2,1:2), rect(1:4)

   !initialize:
   theintegral = 0.d0

   xlow = topoparam(1); xhi = topoparam(2)
   ylow = topoparam(3); yhi = topoparam(4)
   dx = topoparam(5); dy = topoparam(6)

!========TEST FOR SMALL ROUNDING ERROR==========
   call intersection(indicator, area, domain, domain, topoparam(1:4))

   xim = domain(1); xip = domain(2)
   yjm = domain(3); yjp = domain(4)

!=============INTEGRATE PIECEWISE BILINEAR OVER RECTANGULAR REGION====
   if (intmethod .eq. 1) then !use bilinear method

      !don't waste time looping through the entire grid
      !just find indices that include the rectangular region

      distart = (xim - xlow)/dx
      istart = idint(distart) + 1
      
      djstart = (yjm - ylow)/dy
      jstart = idint(djstart) + 1

      diend = (xip - xlow)/dx
      iend = ceiling(diend) + 1

      djend = (yjp-ylow)/dy
      jend = ceiling(djend) + 1

      istart = max(istart, 1)
      jstart = max(jstart, 1)
      iend = min(mx, iend)
      jend = min(my, jend)


      do j = jstart, jend-1
         y1 = ylow + (j - 1.d0)*dy
         y2 = ylow + j*dy
         ! the array zz is indexed from north to south: jjz is the actual index
         ! of interest in the array zz
         jz1 = my - j + 1
         jz2 = jz1 - 1

         do i = istart, iend-1
            x1 = xlow + (i - 1.d0)*dx
            x2 = xlow + i*dx

            intcorners(1,1) = z(i,jz1)
            intcorners(1,2) = z(i,jz2)
            intcorners(2,1) = z(i+1,jz1)
            intcorners(2,2) = z(i+1,jz2)
               
            inttopoparam(1) = x1; inttopoparam(2) = x2
            inttopoparam(3) = y1; inttopoparam(4) = y2

            if (coordinate_system .eq. 1) then !cartesian rectangle
               theintegral = theintegral + bilinearintegral(domain, &
                                                      inttopoparam, &
                                                      intcorners)
            elseif (coordinate_system .eq. 2) then !integrate on surface of sphere
               theintegral = theintegral + bilinearintegral_s(domain, &
                                                        inttopoparam, &
                                                        intcorners)
            else
               write(*,*)  'TOPOINTEGRAL: coordinate_system error'
            endif
         enddo
      enddo

   else
      write(*,*) 'TOPOINTEGRAL: only intmethod = 1,2 is supported'
   endif

   topointegral = theintegral
   return
end function
