program test
use types
use init
use mod_create_block
use refinement

implicit none
integer         , parameter         :: MAXCELLS = 12
integer         , parameter         :: MAXPNTS  = 50
integer         , parameter         :: MAXLIST = MAXCELLS

type(tCell), allocatable            :: cells(:)
integer                             :: nCells

type(tParentCell), allocatable      :: parentCells(:)
integer                             :: nParentCells

real(kind = 8), allocatable         :: pnts(:,:)
integer                             :: nPnts

integer, allocatable                :: refineType(:)   ! Type of Refinement for each Cell (Global Cell index)
                                                       ! 1 = i-Ref, 2=j-Ref, 3 = both Ref

integer, allocatable                :: refineList(:)   ! List with all the Cells to refine (nRefine Values)
integer                             :: nRefine

integer, allocatable                :: canCoarseList(:)   ! List all the Cells(parentCell) which can be coarsed
integer                             :: nCanCoarse

integer, allocatable                :: doCoarseList(:)   ! List all Cells(Cell[first child]) to coarse
integer                             :: nDoCoarse

integer, allocatable                :: holesParentCells(:)   ! List of Holes in ParentCellsArray
integer                             :: nHolesParentCell

integer, allocatable                :: holesPnts(:)   ! List of Holes in Point Array
integer                             :: nHolesPnt

type(tFace), allocatable            :: faces(:)
integer                             :: nFace

integer, allocatable                :: holesFaces(:)   ! List of Holes in Point Array
integer                             :: nHolesFace

write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'( 3("="),10X,A)') "UNIT TEST NO REFINEMENT"
write(*,'(90("="))') 
write(*,'(90("="))') 

allocate (cells(MAXCELLS))
allocate (parentCells(MAXCELLS))
allocate (pnts(2,MAXPNTS))
allocate (faces(MAXCELLS ))

allocate (refineList       (MAXLIST))
allocate (refineType       (MAXCELLS))
allocate (canCoarseList    (MAXCELLS))
allocate (doCoarseList     (MAXLIST))
allocate (holesParentCells (MAXCELLS))
allocate (holesPnts        (MAXPNTS))
allocate (holesFaces       (MAXCELLS))
nFace             = 0
nRefine           = 0
nDoCoarse         = 0
nCanCoarse        = 0
nHolesParentCell  = 0
nHolesPnt         = 0
nHolesFace        = 0
call create_block(4,4,cells,pnts,nCells,nPnts)
call init_sol(cells,parentCells,pnts,faces,nCells,nParentCells,nFace)
end program test
