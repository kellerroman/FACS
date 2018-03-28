module grid_metrics
use types
implicit none
contains
subroutine calc_center(cells,pnts,nCells)
implicit none
type(tCell)       , intent (inout)            :: cells(:)
real(kind = 8)    , intent (in)               :: pnts(:,:)
integer           , intent (in)               :: nCells


integer                                   :: i

do i = 1, nCells
   cells(i) % center(:) = 0.25d0 * ( pnts(:,cells(i) % pnts(1)) & 
                                   + pnts(:,cells(i) % pnts(2)) & 
                                   + pnts(:,cells(i) % pnts(3)) & 
                                   + pnts(:,cells(i) % pnts(4)) ) 
end do
end subroutine calc_center
end module grid_metrics
