module types
use const
implicit none
type :: tCell
! **************************************************************************************************
!  Point ordering              Neighbors
!
!   4 _____ 3                __4__ 
!    |     |                |     |
!    |     |              1 |     | 2
!    |_____|                |_____|
!   1       2                  3
! 
! Neighbor always points on the Upper / left cell if the adjecent cell is more refined
! **************************************************************************************************
   integer                                :: pnts(4)
   integer                                :: neigh(4)
   integer                                :: refineLevel(2)
   integer                                :: ref
   real(kind = 8)                         :: Q(Q_DIM)
   real(kind = 8)                         :: QC(Q_DIM)
   real(kind = 8)                         :: aux(4)
   real(kind = 8)                         :: center(2)
   real(kind = 8)                         :: grad(2,Q_DIM)

end type tCell
type :: tParentCell
! **************************************************************************************************
!  Point ordering              Neighbors         Child
!
!   4 _____ 3                __4__               _____
!    |     |                |     |             | 1|4 |
!    |     |              1 |     | 2           |-----|
!    |_____|                |_____|             | 2|3 |
!   1       2                  3                 
! 
! Neighbor always points on the Upper / left cell if the adjecent cell is more refined
! **************************************************************************************************
   integer                                :: parent
   integer                                :: child(4)
   integer                                :: pnts(4)
   integer                                :: refineLevel(2)
   integer                                :: ref
   integer                                :: neigh(4)
   integer                                :: cut_type ! 0 = Not Cut
                                                      ! 1 = Cut in i
                                                      ! 2 = Cut in j
                                                      ! 3 = Cut in i and j 
   integer                                :: pos_CanCoarse   ! Position in the CanCoarse List
   !real(kind = 8)                         :: center(2)
   !real(kind = 8)                         :: grad(2)

end type tParentCell
end module types
