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
   real(kind = 8)                         :: q(Q_DIM)
   real(kind = 8)                         :: qc(Q_DIM)
   real(kind = 8)                         :: aux(4)
   real(kind = 8)                         :: center(2)
   real(kind = 8)                         :: grad(2,Q_DIM)
   real(kind = 8)                         :: vol                     ! Inverse of the Volume/Area of the Cell
   integer                                :: nFace                   ! Number of Fluxes the Cell has 
   integer                                :: faces(GIT_DIM*2*2)      ! Dimension of grid * 2 (for 2 neighbor Cells per Side)
                                                                     !                   * 2 (for both sides)
   real(kind = 8)                         :: f_sign(GIT_DIM*2*2)     ! Sign of each Flux of the cell face

   integer                                :: faces_dir(GIT_DIM*2,2)
   logical                                :: faces_dir_ref(GIT_DIM*2)

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
type :: tFace
   integer                                :: stencil(1,2)
   real(kind = 8)                         :: q_const(Q_DIM,2)     ! Reconstructed Value at Cell Interface
   real(kind = 8)                         :: n(GIT_DIM)           ! Normal Vector of the Face
   real(kind = 8)                         :: area                 ! length / Area of the Face
   real(kind = 8)                         :: flux(Q_DIM)          ! Flux in the Face
end type tFace
end module types
