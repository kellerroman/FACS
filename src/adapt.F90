program adapt
! **************************************************************************************************
! ***               Fully Adaptive Combustion Solver (FACS)                                      ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   19.01.2018
! Last changes: 30.03.2018
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
!   - 2018-03-28,RK : Grid Refining and Coarseing based on Gradient working for Refining in both 
!                     directions at the same time
!   - 2018-03-30,RK : 1D Euler equation working
! **************************************************************************************************
use types
use file_io
use refinement
use init
use choose
use fluxes
use git_version_module
use array_holes
implicit none
integer         , parameter         :: MAXCELLS = 200000
integer         , parameter         :: MAXPNTS  = 400000
integer         , parameter         :: MAXLIST = MAXCELLS / 10
character(len=*),parameter          :: FILENAME_IN = "sol.dat"

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

integer                             :: iter,n,f

integer                             :: niter

real(kind=8) :: dt,qt(Q_DIM)

write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'( 3("="),10X,A)') "Grid Adapter"
write(*,'( 3("="),10X,A)') GIT_VERSION
write(*,'( 3("="),10X,A)') GIT_DATE
write(*,'(90("="))') 
write(*,'(90("="))') 

niter = 200
dt = 0.2 / dble(niter)

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

call read_sol(FILENAME_IN,cells,pnts,nCells,nPnts)

call init_sol(cells,parentCells,pnts,faces,nCells,nParentCells,nFace)

do iter = 1,niter
   write(*,'(10("="),3X,I5.5,3X,10("="))') iter
   do n = 1, nFace % nEntry
      call inv_flux(cells(faces(n) % stencil(1,1)) % QC, cells(faces(n) % stencil(1,2)) % QC,faces(n) % n, faces(n) % flux)
      faces(n) % flux = faces(n) % flux * faces(n) % area
   end do
   do n = 1, nCells % nEntry
      qt = 0.0d0
      do f = 1, cells(n) % nFace
         qt = qt + cells(n) % f_sign(f) * faces(cells(n) % faces(f) ) % flux
      end do
      qt = dt * cells(n) % vol * qt + cells(n) % qc

      cells(n) % qc(:)      = qt
      cells(n) % q(1)       = qt(1)
      cells(n) % q(2:Q_DIM) = qt(2:Q_DIM) / qt(1) 
   end do
!   call calc_gradient(cells,nCells)
!   call choose_cells(cells,nCells,refineType,refineList,nRefine)
!   call smooth_refinement(cells,refineType,refineList,nRefine)
!   call choose_coarse(cells,parentCells,refineType,canCoarseList,nCanCoarse,doCoarseList,nDoCoarse)
!
!   call doRefinement (cells,parentCells,pnts,faces                      &
!                     ,nCells,nParentCells,nPnts,nFace                   &
!                     ,refineType,refineList,nRefine                     &
!                     ,canCoarseList,nCanCoarse                          &
!                     ,doCoarseList,nDoCoarse                            &
!                     ,.false.)
   if (mod(iter,10) == 0) &
   call write_sol(cells,pnts,nCells,nPnts,FILENAME_IN,iter)
end do
call check_neighbors(cells,nCells%nEntry,pnts)
call write_sol(cells,pnts,nCells,nPnts,FILENAME_IN)

write(*,'(90("="))') 
write(*,'(90("="))') 
write(*,'( 3("="),10X,A)') "Grid Adapter done!"
write(*,'(90("="))') 
write(*,'(90("="))') 
end program adapt
