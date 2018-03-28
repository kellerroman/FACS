program adapt
! **************************************************************************************************
! ***                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRle                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   19.01.2018
! Last changes: 19.01.2018
! Version:      V0.0.1
! --------------------------------------------------------------------------------------------------
! Description:
!   Adapts a given Grid file
!   
! --------------------------------------------------------------------------------------------------
! Comments and Notes:
!   
! --------------------------------------------------------------------------------------------------
! References:
!
! --------------------------------------------------------------------------------------------------
! Author and Change History:
!   - 2018-01-19,RK : Started of Project
! **************************************************************************************************
use types
use file_io
use refinement
use grid_metrics
use choose
implicit none
integer         , parameter         :: MAXCELLS = 100000
integer         , parameter         :: MAXPNTS  = 200000
character(len=*),parameter          :: FILENAME_IN = "sol.dat"
integer                             :: nCells
integer                             :: nParentCells
integer                             :: nPnts
type(tCell), allocatable            :: cells(:)
type(tParentCell), allocatable      :: parentCells(:)
real(kind = 8), allocatable         :: pnts(:,:)
integer, allocatable                :: refineType(:)   ! Type of Refinement for each Cell (Global Cell index)
                                                       ! 1 = i-Ref, 2=j-Ref, 3 = both Ref
integer, allocatable                :: refineList(:)   ! List with all the Cells to refine (nRefine Values)
integer                             :: nRefine
integer, allocatable                :: canCoarseList(:)   ! List all the Cells(parentCell) which can be coarsed
integer                             :: nCanCoarse
integer, allocatable                :: doCoarseList(:)   ! List all Cells(Cell[first child]) to coarse
integer                             :: nDoCoarse

integer                             :: iter,n
real(kind = 8)                      :: r

write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'(3("="),10X,A)') "Grid Adapter"

write(*,'(90("="))') 
write(*,'(90("="))') 

allocate (cells(MAXCELLS))
allocate (parentCells(MAXCELLS))
allocate (pnts(2,MAXPNTS))

call read_sol(cells,parentCells,pnts,nCells,nParentCells,nPnts,FILENAME_IN)
call calc_center(cells,pnts,nCells)
allocate (refineList(MAXCELLS))
allocate (refineType(MAXCELLS))
allocate (canCoarseList(MAXCELLS))
allocate (doCoarseList(MAXCELLS))
nCanCoarse = 0

do iter = 1,8
   write(*,'(10("="),3X,I5.5,3X,10("="))') iter
   do n = 1, nCells
      r = sqrt(cells(n) % center(1)**2 + cells(n) % center(2) **2)
      cells(n) % var = tanh(4.0d0* (2.0d0*r-1.0d0))*0.5d0 + 0.5d0
!      if (r > 0.3d0 + iter*0.02d0) then
!         cells(n) % var = 0
!      else
!         cells(n) % var = 1
!      end if
      if (cells(n) % center(1) > 0.1d0 + iter*0.02d0) then
         cells(n) % var = 0
      else
         cells(n) % var = 1
      end if
   end do
   call calc_gradient(cells,nCells)
   call choose_cells(cells,nCells,refineType,refineList,nRefine)
   call smooth_refinement(cells,refineType,refineList,nRefine)
   call choose_coarse(cells,parentCells,refineType,canCoarseList,nCanCoarse,doCoarseList,nDoCoarse)

   call doRefinement (cells,parentCells,pnts,nCells,nParentCells,nPnts  &
                     ,refineType,refineList,nRefine                     &
                     ,canCoarseList,nCanCoarse                          &
                     ,doCoarseList,nDoCoarse                            &
                     ,.false.)
   call write_sol(cells,pnts,nCells,nPnts,FILENAME_IN,iter)
   write(*,*) "Coarseable Cells:",canCoarseList(1:nCanCoarse)
end do
do n = 1, nCells
   r = sqrt(cells(n) % center(1)**2 + cells(n) % center(2) **2)
   cells(n) % var = tanh(4.0d0* (2.0d0*r-1.0d0))*0.5d0 + 0.5d0
end do
call check_neighbors(cells,nCells,pnts)
!call write_sol(cells,pnts,nCells,nPnts,FILENAME_IN)

do n = 1, nCells
   write(666,*) n, cells(n) % neigh
end do
write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'(3("="),10X,A)') "Grid Adapter done!"
write(*,'(90("="))') 
write(*,'(90("="))') 
end program adapt
