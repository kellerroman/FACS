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
use fluxes
implicit none
include "git_version.h"
integer         , parameter         :: MAXCELLS = 200000
integer         , parameter         :: MAXPNTS  = 400000
integer         , parameter         :: MAXLIST = MAXCELLS / 10
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

integer, allocatable                :: holesParentCells(:)   ! List of Holes in ParentCellsArray
integer                             :: nHolesParentCell

integer, allocatable                :: holesPnts(:)   ! List of Holes in Point Array
integer                             :: nHolesPnt

integer                             :: iter,n

integer                             :: niter

real(kind=8) :: fl(Q_DIM) ,fr(Q_DIM),dt
real(kind=8),allocatable :: qt(:,:)

write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'( 3("="),10X,A)') "Grid Adapter"
write(*,'( 3("="),10X,A)') GIT_VERSION
write(*,'( 3("="),10X,A)') GIT_DATE
write(*,'(90("="))') 
write(*,'(90("="))') 

niter = 200
dt = 0.2 / dble(niter)
! dt = dt / dx

allocate (cells(MAXCELLS))
allocate (parentCells(MAXCELLS))
allocate (pnts(2,MAXPNTS))

allocate (refineList       (MAXLIST))
allocate (refineType       (MAXCELLS))
allocate (canCoarseList    (MAXCELLS))
allocate (doCoarseList     (MAXLIST))
allocate (holesParentCells (MAXCELLS))
allocate (holesPnts        (MAXPNTS))
nRefine           = 0
nDoCoarse         = 0
nCanCoarse        = 0
nHolesParentCell  = 0
nHolesPnt         = 0

call read_sol(cells,parentCells,pnts,nCells,nParentCells,nPnts,FILENAME_IN)
call calc_center(cells,pnts,nCells)

write(*,*) "dx", 1.0d0 / dble(ncells), "dt", dt
dt = dt / (1.0d0 / dble(ncells)) 
do n = 1, nCells
   cells(n) % QC(1)       = cells(n) % Q(1)
   cells(n) % QC(2:Q_DIM) = cells(n) % Q(1) * cells(n) % Q(2:Q_DIM)
end do
do iter = 1,niter
   write(*,'(10("="),3X,I5.5,3X,10("="))') iter
   allocate(qt(Q_DIM,nCells))
   do n = 1, nCells
      if (cells(n) % neigh(1) /= NO_CELL) then
         call inv_flux(cells(cells(n) % neigh(1)) % QC, cells(n) % QC,dt, fl)
      else 
         call inv_flux(cells(n) % QC, cells(n) % QC, dt, fl)
         !fl = 0.0d0
      end if
      if (cells(n) % neigh(2) /= NO_CELL) then
         call inv_flux(cells(n) % QC, cells(cells(n) % neigh(2)) % QC,dt, fr)
      else 
         call inv_flux(cells(n) % QC, cells(n) % QC, dt, fr)
         !fr = 0.0d0
      end if
      qt(:,n)  = cells(n) % QC + dt * (fl-fr)
   end do
   do n = 2, nCells-1
      cells(n) % qc(:)      = qt(:,n)
      cells(n) % q(1)       = qt(1,n)
      cells(n) % Q(2:Q_DIM) = qt(2:Q_DIM,n) / qt(1,n) 
   end do
   deallocate(qt)
!   call calc_gradient(cells,nCells)
!   call choose_cells(cells,nCells,refineType,refineList,nRefine)
!   call smooth_refinement(cells,refineType,refineList,nRefine)
!   call choose_coarse(cells,parentCells,refineType,canCoarseList,nCanCoarse,doCoarseList,nDoCoarse)
!
!   call doRefinement (cells,parentCells,pnts,nCells,nParentCells,nPnts  &
!                     ,refineType,refineList,nRefine                     &
!                     ,canCoarseList,nCanCoarse                          &
!                     ,doCoarseList,nDoCoarse                            &
!                     ,holesParentCells,nHolesParentCell                 &
!                     ,holesPnts,nHolesPnt                               &
!                     ,.false.)
   if (mod(iter,10) == 0) &
   call write_sol(cells,pnts,nCells,nPnts,FILENAME_IN,iter)
end do
call check_neighbors(cells,nCells,pnts)
call write_sol(cells,pnts,nCells,nPnts,FILENAME_IN)

write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'( 3("="),10X,A)') "Grid Adapter done!"
write(*,'(90("="))') 
write(*,'(90("="))') 
end program adapt
