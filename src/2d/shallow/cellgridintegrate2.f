!====================================================================
subroutine cellgridintegrate(cell, topoint)
!=====================================================================

! *** Note: xcell and ycell are no longer needed -- should be removed.

   use topo_module, only: topo_data, rectintegral, intersection

   ! arguments
   real(kind=8) , intent(in) :: cell(1:4), topoint
   
   ! local
   integer :: im, mfid, indicator
   real(kind=8), external :: topointegral
   real(kind=8) :: cellintdomain(1:4), topoparam(1:6)
   real(kind=8) :: cellarea, area


!     ##############################################################################
!     cellgridintegrate integrates a unique surface, over a rectangular cell
!     defined from data from multiple regular Cartesian grids
!     (using the finest data available in any region)

!     The rectangle has coords:
!     xim <= x <= xip, yjm <= y <= yjp, with center (x,y) = (xcell, ycell)

!     The intersection (with one particular grid has coords:
!     xintlo <= x <= xinthi, yintlo <= y <= yinthi

!     The _set_ version uses a recursive strategy using the formulas for
!     intersections of sets.

   !initialize the integral of the surface
   topoint = 0.d0

   !determine the type of integration
   im = 1

   !first see if the grid cell is entirely in a fine topofile
   do m = 1, num_topo_files
      !look at topofiles, from fine to coarse
      mfid = topo_order(m)
      !check for intersection of grid cell and this topofile
      cellarea = (cell(2) - cell(1))*(cell(4) - cell(3))
         
      topoparam(1) = topo_data(mfid)%x_lower
      topoparam(2) = topo_data(mfid)%x_upper
      topoparam(3) = topo_data(mfid)%y_lower
      topoparam(4) = topo_data(mfid)%y_upper
      topoparam(5) = topo_data(mfid)%dx
      topoparam(6) = topo_data(mfid)%dy
         
      call intersection(indicator, area, cellintdomain, cell, & 
                        topoparam(1:4))
      if (indicator .eq. 1) then !cell overlaps grid
         if (area .eq. cellarea) then !cell is entirely in grid
            ! (should we check if they agree to some tolerance??)
            !integrate surface and get out of here
            topoint = topoint + topointegral(cellintdomain, &
                     topoparam, topo_data(mfid)%num_cells(1), &
                     topo_data(mfid)%num_cells(2), topo_data(mfid)%z, im)
            return
         else
            call rectintegral(cell, m, topoint)
         endif
      endif
   enddo

   write(6,601) cell(1),cell(2),cell(3),cell(4)
601  format('*** Error, grid cell does not overlap any topo grid',/,
     &     '  xim = ',e24.14,'  xip = ',e24.14,      
     &   /,'  yjm = ',e24.14,'  yjp = ',e24.14)
   stop
end subroutine cellgridintegrate
