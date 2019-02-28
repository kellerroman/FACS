program test
use types
use init
use mod_create_block
use refinement
use file_io

implicit none
integer         , parameter         :: MAXCELLS = 50
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

integer :: i

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
nRefine = 1
refineList(nRefine) = 5
refineType(5) = 3

call doRefinement (cells,parentCells,pnts,faces                      &
                  ,nCells,nParentCells,nPnts,nFace                   &
                  ,refineType,refineList,nRefine                     &
                  ,canCoarseList,nCanCoarse                          &
                  ,doCoarseList,nDoCoarse                            &
                  ,holesParentCells,nHolesParentCell                 &
                  ,holesPnts,nHolesPnt                               &
                  ,holesFAces,nHolesFace                             &
                  ,.true.)

call write_sol(cells,pnts,nCells,nPnts,"sol.dat")

do i = 1, nCells
    write(*,*) i, cells(i) % refineLevel, cells(i) % pnts, cells(i) % neigh
end do
end program test
