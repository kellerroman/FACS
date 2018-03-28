module types
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
   real(kind = 8)                         :: var
   real(kind = 8)                         :: center(2)
   real(kind = 8)                         :: grad(2)

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
   integer                                :: cut_type ! 1 = Cut in i and j 
                                                      ! 2 = Cut in i
                                                      ! 3 = Cut in j
   !real(kind = 8)                         :: var
   !real(kind = 8)                         :: center(2)
   !real(kind = 8)                         :: grad(2)

end type tParentCell
end module types
