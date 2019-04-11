program test
use types
use init
use mod_create_block
use refinement
use file_io
use array_holes
use choose

implicit none
integer         , parameter         :: MAXCELLS = 150
integer         , parameter         :: MAXPNTS  = 250
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

integer                             :: i,j

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
nRefine           = 0
nDoCoarse         = 0
nCanCoarse        = 0

call create_block(4,4,cells,pnts,nCells,nPnts)

call init_sol(cells,parentCells,pnts,faces,nCells,nParentCells,nFace)
do j = 1,7
    nRefine = 0
    !do i = 1, 1 ! nCells % nEntry 
    do i = 5, 5
        nRefine       = nRefine + 1
        refineList(nRefine) = i
        refineType(i) = 3
    end do

    !i = 5
    !nRefine       = 1
    !refineList(nRefine) = i
    !refineType(i) = 3

    call smooth_refinement(cells,refineType,refineList,nRefine)

    call doRefinement (cells,parentCells,pnts,faces                      &
                    ,nCells,nParentCells,nPnts,nFace                   &
                    ,refineType,refineList,nRefine                     &
                    ,canCoarseList,nCanCoarse                          &
                    ,doCoarseList,nDoCoarse                            &
                    ,.true.)

    call check_points(pnts,nPnts)
end do


write(*,*) "Coarseable  CELLS ",canCoarseList(1:nCanCoarse)
call write_sol(cells,pnts,nCells,nPnts,"sol.dat")

do i = 1, nCells % nEntry
    !write(*,*) i, cells(i) % refineLevel, cells(i) % pnts, cells(i) % neigh
end do
end program test
