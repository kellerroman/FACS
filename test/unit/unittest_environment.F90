module unittest_environment
use types
use init
use mod_create_block
use refinement
use file_io
use array_holes
use choose
implicit none
integer         , parameter         :: MAXCELLS = 5000
integer         , parameter         :: MAXPNTS  = 5000
integer         , parameter         :: MAXLIST = MAXCELLS

type(tCell), allocatable            :: cells(:)
type(holes)                         :: nCells

type(tParentCell), allocatable      :: parentCells(:)
type(holes)                         :: nParentCells

real(kind = 8), allocatable         :: pnts(:,:)
type(holes)                         :: nPnts

integer, allocatable                :: refineType(:)   ! Type of Refinement for each Cell (Global Cell index)
                                                       ! 1 = i-Ref, 2=j-Ref, 3 = both Ref

integer, allocatable                :: refineList(:)   ! List with all the Cells to refine (nRefine Values)
integer                             :: nRefine

integer, allocatable                :: canCoarseList(:)   ! List all the Cells(parentCell) which can be coarsed
integer                             :: nCanCoarse

integer, allocatable                :: doCoarseList(:)   ! List all Cells(Cell[first child]) to coarse
integer                             :: nDoCoarse


type(tFace), allocatable            :: faces(:)
type(holes)                         :: nFace
contains
subroutine ut_init()
implicit none
allocate (cells(MAXCELLS))
allocate (parentCells(MAXCELLS))
allocate (pnts(2,MAXPNTS))
allocate (faces(MAXCELLS*2))

allocate (refineList       (MAXLIST))
allocate (refineType       (MAXCELLS))
allocate (canCoarseList    (MAXCELLS))
allocate (doCoarseList     (MAXLIST))
nRefine           = 0
nDoCoarse         = 0
nCanCoarse        = 0
refineType        = 0
end subroutine ut_init

subroutine ut_cleanup()
implicit none
deallocate (cells)
deallocate (parentCells)
deallocate (pnts)
deallocate (faces)

deallocate (refineList)
deallocate (refineType)
deallocate (canCoarseList)
deallocate (doCoarseList)
end subroutine ut_cleanup
end module
