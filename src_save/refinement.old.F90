module refinement
contains
subroutine doRefinement (cells,parentCells,pnts,nCells,nParentCells,nPnts,refine,debug_in)
use types
implicit none
type(tCell)   , intent (inout)            :: cells(:)
type(tParentCell)   , intent (inout)      :: parentCells(:)
real(kind = 8), intent (inout)            :: pnts(:,:)

integer       , intent (inout)            :: nCells
integer       , intent (inout)            :: nParentCells
integer       , intent (inout)            :: nPnts
integer       , intent (inout)            :: refine(:,:)
logical       , intent (in), optional     :: debug_in


integer                                   :: nRef
integer                                   :: nNewCells
integer                                   :: nNewParentCells
integer                                   :: nNewPnts

integer                                   :: nMaxCells
integer                                   :: nMaxParentCells
integer                                   :: nMaxPnts

integer, allocatable                      :: CellId_old2new(:)
integer, allocatable                      :: CellMove(:)

integer                                   :: i,r,ctr,nc,oc,nc1,nc2,nc3
integer                                   :: opc,npc1,npc2,npc3,npc4              ! New Parent cells
integer                                   :: rfl,nrfl
real(kind = 8)                            :: dvar

logical                                   :: debug


if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
end if
nMaxCells = ubound(cells,1)
nMaxParentCells = ubound(parentCells,1)
nMaxPnts  = ubound(pnts,2)
nRef = 0
nNewCells = nCells
nNewParentCells = nParentCells
nNewPnts  = nPnts
do i = 1, nCells
   if (refine(1,i) > 0 .and. refine(2,i) == 0) then
      nRef      = nRef + 1
      nNewPnts  = nNewPnts  + 2
      nNewCells = nNewCells + 1
      nNewParentCells = nNewParentCells + 2
   else if (refine(1,i) == 0 .and. refine(2,i) > 0) then
      nRef      = nRef + 1
      nNewPnts  = nNewPnts  + 2
      nNewCells = nNewCells + 1
      nNewParentCells = nNewParentCells + 2
   else if (refine(1,i) > 0 .and. refine(2,i) > 0) then
      nRef      = nRef + 1
      nNewPnts  = nNewPnts  + 5
      nNewCells = nNewCells + 3
      nNewParentCells = nNewParentCells + 4
   end if
end do

if (nNewCells > nMaxCells) then
   write(*,*) "Maximum number of Cells reached"
   stop 1
end if

if (nNewParentCells > nMaxParentCells) then
   write(*,*) "Maximum number of ParentCells reached"
   stop 1
end if

if (nNewPnts  > nMaxPnts) then
   write(*,*) "Maximum number of Points reached"
   stop 1
end if

allocate (CellId_old2new(0:nNewCells))
allocate (CellMove(nCells))

CellMove(1) = 0
CellId_old2new(1) = 1
CellId_old2new(0) = 0

do i = 2, nCells
   if (refine(1,i-1) > 0 .and. refine(2,i-1) == 0) then
      CellMove(i) =  CellMove(i-1) + 1
   else if (refine(1,i-1) == 0 .and. refine(2,i-1) > 0) then
      CellMove(i) =  CellMove(i-1) + 1
   else if (refine(1,i-1) > 0 .and. refine(2,i-1) > 0) then
      CellMove(i) =  CellMove(i-1) + 3
   else
      CellMove(i) =  CellMove(i-1)
   end if
   if (debug)write(*,*) i,CellMove(i), refine(:,i)
   CellId_old2new(i) = i + CellMove(i)
end do

!! Movement of the old Cells to their new location
do i = nCells, 2,-1
   ! stop if the last cell to move has been reached
   if (CellMove(i) == 0) then
      exit
   end if
   nc = i+CellMove(i)
   cells(nc) = cells(i) 
end do

! Updated neighbour information to the new Cells structure
do i = 1, nNewCells
   do r = 1,4
      cells(i) % neigh(r) = CellId_old2new(cells(i) % neigh(r))
   end do
end do
write(*,'("Number of Refinements : ",I0)') nRef
write(*,'("Number of New Cells   : ",I0)') nNewCells
write(*,'("Number of New ParCells: ",I0)') nNewParentCells
write(*,'("Number of New Points  : ",I0)') nNewPnts
ctr = 0
do r = 1, nRef
   do
      ctr = ctr + 1
      oc = CellId_old2new(ctr)
      if (debug) write(*,*) "checking Cell:",ctr,oc

! **************************************************************************************************
!           
!                REFINEMENT IN I DIRECTION
!           
! **************************************************************************************************

      if (refine(1,ctr) > 0 .and. refine(2,ctr) == 0) then
         !
         !          ______________
         !         |       |      |
         !         |  OC   |  NC  |          OC Original Cell
         !         |       |      |          NC NEW Cell
         !         |_______|______|
         !        
         !        
         if (debug)write(*,*) "Refining Cell in i-Dir:", ctr, oc
         !!!!  ADD NEW POINTS
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(2)) ) * 0.5d0
         if (debug)write(*,*) "New Point:", nPnts, Pnts(:,nPnts)
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(4)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
         if (debug)write(*,*) "New Point:", nPnts,Pnts(:,nPnts)

         nc = oc + 1
         !! UPDATE NEW CELL
         cells(nc) % pnts(1) = nPnts - 1    
         cells(nc) % pnts(2) =  cells(oc) % pnts(2)
         cells(nc) % pnts(3) =  cells(oc) % pnts(3)
         cells(nc) % pnts(4) = nPnts
         cells(nc) % center(:) = 0.25d0 * ( pnts(:,cells(nc) % pnts(1)) & 
                                          + pnts(:,cells(nc) % pnts(2)) & 
                                          + pnts(:,cells(nc) % pnts(3)) & 
                                          + pnts(:,cells(nc) % pnts(4)) ) 

         cells(nc) % refineLevel(1) = cells(oc) % refineLevel(1) + 1
         cells(nc) % refineLevel(2) = cells(oc) % refineLevel(2)
         cells(nc) % neigh(1) = oc
         cells(nc) % neigh(2) = cells(oc) % neigh(2)
         cells(nc) % neigh(3) = cells(oc) % neigh(3)
         cells(nc) % neigh(4) = cells(oc) % neigh(4)
         ! Update Neighboor of right adjecent cell
         nc1 = cells(nc) % neigh(2) 
         cells(nc1) % neigh(1) = nc
         
         !!!! UPDATE OLD CELL
         cells(oc) % pnts(2) =  cells(nc) % pnts(1)
         cells(oc) % pnts(3) =  cells(nc) % pnts(4)
         cells(oc) % center(:) = 0.25d0 * ( pnts(:,cells(oc) % pnts(1)) & 
                                          + pnts(:,cells(oc) % pnts(2)) & 
                                          + pnts(:,cells(oc) % pnts(3)) & 
                                          + pnts(:,cells(oc) % pnts(4)) ) 
         cells(oc) % refineLevel(1) = cells(oc) % refineLevel(1) + 1
         cells(oc) % neigh(2) = nc

         dvar = cells(oc) % grad(1) * (cells(nc) % center(1) - cells(oc) % center(1)) * 0.5d0
         cells(nc) % var = cells(oc) % var ! + dvar
         !cells(oc) % var = cells(oc) % var - dvar
         cells(nc) % grad = cells(oc) % grad

         ! Parent Cell

         npc1 = nParentCells + 1
         npc2 = nParentCells + 2
         nParentCells = npc2
         opc = cells(oc) % ref                  ! old parent cell
         parentCells(opc) % ref = -1               ! Cell has childs, thus no tin cells array anymore
         parentCells(opc) % cut_type = 2
         parentCells(opc) % child(1) = npc1
         parentCells(opc) % child(2) = npc2

         parentCells(npc1) % parent   = opc
         parentCells(npc1) % ref      = oc                                 ! referencing
         cells(oc)         % ref      = npc1
         parentCells(npc1) % child    = -1
         parentCells(npc1) % neigh(1) = parentCells(opc) % neigh(1)
         parentCells(npc1) % neigh(2) = npc2
         parentCells(npc1) % neigh(3) = parentCells(opc) % neigh(3)
         parentCells(npc1) % neigh(4) = parentCells(opc) % neigh(4)
         parentCells(npc1) % refineLevel = cells(oc) % refineLevel

         parentCells(npc2) % parent = opc
         parentCells(npc2) % ref    = nc
         cells(nc)         % ref      = npc2
         parentCells(npc2) % child  = -1
         parentCells(npc2) % neigh(1) = npc1
         parentCells(npc2) % neigh(2) = parentCells(opc) % neigh(2)
         parentCells(npc2) % neigh(3) = parentCells(opc) % neigh(3)
         parentCells(npc2) % neigh(4) = parentCells(opc) % neigh(4)
         parentCells(npc2) % refineLevel = cells(nc) % refineLevel
         exit

! **************************************************************************************************
!           
!                REFINEMENT IN Y DIRECTION
!           
! **************************************************************************************************
      else if (refine(2,ctr) > 0 .and. refine(1,ctr) == 0) then
         !
         !          ________
         !         |       |
         !         |  OC   |          OC Original Cell
         !         |  NPC1 |
         !         |_______|
         !         |       |
         !         |  NC   |
         !         |  NPC2 |          NC NEW Cell
         !         |_______|
         !        
         !        
         if (debug) write(*,*) "Refining Cell in j-Dir:", ctr, oc
         !!!!  ADD NEW POINTS
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(4)) ) * 0.5d0
         if (debug) write(*,*) "New Point:", nPnts, Pnts(:,nPnts)
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(2)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
         if (debug) write(*,*) "New Point:", nPnts,Pnts(:,nPnts)

         nc = oc + 1
         !! UPDATE NEW CELL
         cells(nc) % pnts(1) = cells(oc) % pnts(1)
         cells(nc) % pnts(2) = cells(oc) % pnts(2)
         cells(nc) % pnts(3) = nPnts
         cells(nc) % pnts(4) = nPnts - 1    
         cells(nc) % center(:) = 0.25d0 * ( pnts(:,cells(nc) % pnts(1)) & 
                                          + pnts(:,cells(nc) % pnts(2)) & 
                                          + pnts(:,cells(nc) % pnts(3)) & 
                                          + pnts(:,cells(nc) % pnts(4)) ) 
         cells(nc) % refineLevel(1) = cells(oc) % refineLevel(1)
         cells(nc) % refineLevel(2) = cells(oc) % refineLevel(2) + 1
         cells(nc) % neigh(1) = cells(oc) % neigh(1)
         cells(nc) % neigh(2) = cells(oc) % neigh(2)
         cells(nc) % neigh(3) = cells(oc) % neigh(3)
         cells(nc) % neigh(4) = oc
         ! Update Neighboor of lower adjecent cell
         nc1 = cells(nc) % neigh(3) 
         cells(nc1) % neigh(4) = nc
         !!!! UPDATE OLD CELL
         cells(oc) % pnts(1) =  cells(nc) % pnts(4)
         cells(oc) % pnts(2) =  cells(nc) % pnts(3)
         cells(oc) % center(:) = 0.25d0 * ( pnts(:,cells(oc) % pnts(1)) & 
                                          + pnts(:,cells(oc) % pnts(2)) & 
                                          + pnts(:,cells(oc) % pnts(3)) & 
                                          + pnts(:,cells(oc) % pnts(4)) ) 
         cells(oc) % refineLevel(2) = cells(oc) % refineLevel(2) + 1
         cells(oc) % neigh(3) = nc

         cells(nc) % var = cells(oc) % var ! + dvar
         cells(nc) % grad = cells(oc) % grad

         ! Parent Cell

         npc1 = nParentCells + 1
         npc2 = nParentCells + 2
         nParentCells = npc2
         opc = cells(oc) % ref                  ! old parent cell
         parentCells(opc) % ref = -1               ! Cell has childs, thus no tin cells array anymore
         parentCells(opc) % cut_type = 3
         parentCells(opc) % child(1) = npc1
         parentCells(opc) % child(2) = npc2

         parentCells(npc1) % parent   = opc
         parentCells(npc1) % ref      = oc                                 ! referencing
         cells(oc)         % ref      = npc1
         parentCells(npc1) % child    = -1
         parentCells(npc1) % neigh(1) = parentCells(opc) % neigh(1)
         parentCells(npc1) % neigh(2) = parentCells(opc) % neigh(2)
         parentCells(npc1) % neigh(3) = npc2
         parentCells(npc1) % neigh(4) = parentCells(opc) % neigh(4)
         parentCells(npc1) % refineLevel = cells(oc) % refineLevel

         parentCells(npc2) % parent   = opc
         parentCells(npc2) % ref      = nc
         cells(nc)         % ref      = npc2
         parentCells(npc2) % child    = -1
         parentCells(npc2) % neigh(1) = parentCells(opc) % neigh(1)
         parentCells(npc2) % neigh(2) = parentCells(opc) % neigh(2)
         parentCells(npc2) % neigh(3) = parentCells(opc) % neigh(3)
         parentCells(npc2) % neigh(4) = npc1
         parentCells(npc2) % refineLevel = cells(nc) % refineLevel
         exit

! **************************************************************************************************
!           
!                REFINEMENT IN Y DIRECTION
!           
! **************************************************************************************************
      else if (refine(2,ctr) > 0 .and. refine(1,ctr) > 0) then
         !
         !          ______-1______
         !         |       |      |
         !         |  OC   |  NC3 |          OC Original Cell
         !         |       |      |          NC1 NEW Cell
         !      -4 |_______|______|-3        NC2 NEW Cell
         !         |      0|      |          NC3 NEW Cell
         !         |  NC1  |  NC2 |
         !         |       |      |
         !         |_______|______|
         !                -2 
         !        
         if (debug)write(*,*) "Refining Cell in both Dirs:", ctr, oc
         !!!!  ADD NEW POINTS
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(4)) ) * 0.5d0
         if (debug)write(*,*) "New Point:", nPnts, Pnts(:,nPnts)
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(2)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
         if (debug)write(*,*) "New Point:", nPnts,Pnts(:,nPnts)
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(2)) ) * 0.5d0
         if (debug)write(*,*) "New Point:", nPnts, Pnts(:,nPnts)
         nPnts = nPnts + 1
         Pnts(:,nPnts) = ( Pnts(:,cells(oc) % pnts(4)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
         if (debug)write(*,*) "New Point:", nPnts,Pnts(:,nPnts)
         nPnts = nPnts + 1
         Pnts(:,nPnts) = cells(oc) % center(:)
         if (debug)write(*,*) "New Point:", nPnts,Pnts(:,nPnts)
         nc1 = oc + 1
         nc2 = oc + 2
         nc3 = oc + 3
         cells(oc) % refineLevel(:) = cells(oc) % refineLevel(:) + 1
         !! UPDATE NEW CELL
         cells(nc1) % pnts(1) = cells(oc) % pnts(1)
         cells(nc1) % pnts(2) = nPnts - 2
         cells(nc1) % pnts(3) = nPnts
         cells(nc1) % pnts(4) = nPnts - 4
         cells(nc1) % refineLevel(:) = cells(oc) % refineLevel(:)
         nc   = cells(oc ) % neigh(1)
         if (nc > 0) then
           rfl  = cells(nc1) % refineLevel(2)
            nrfl = cells(nc ) % refineLevel(2)
            if ( nrfl < rfl) then
               cells(nc1) % neigh(1) = nc
            else if (nrfl == rfl) then
               nc   = cells(nc) % neigh(3) !! get south neighbor
               nrfl = cells(nc) % refineLevel(2)
               if (nrfl == rfl) then
                  cells(nc1) % neigh(1) = nc
                  cells(nc ) % neigh(2) = nc1
               else
                  write(*,*) "Neighbor West: Second Neighbor: refinement level not supported",__LINE__,__FILE__
                  stop 1
               end if
            else if (nrfl == rfl+1) then
               write(*,*) nc1, cells(nc1) % refineLevel, cells(oc) % neigh
               write(*,*) nc , cells(nc ) % refineLevel
               stop 1
            else
               write(*,*) "Neighbor WEST with higher refinment level not implemented yet",__LINE__,__FILE__
               write(*,*) nc1, cells(nc1) % refineLevel
               write(*,*) nc , cells(nc ) % refineLevel
               stop 1
            end if
         else
            cells(nc1) % neigh(1) = nc
         end if

         cells(nc1) % neigh(2) = nc2
         cells(nc1) % neigh(3) = cells(oc) % neigh(3)
         cells(nc1) % neigh(4) = oc
         !! UPDATE NEW CELL
         cells(nc2) % pnts(1) = nPnts - 2
         cells(nc2) % pnts(2) = cells(oc) % pnts(2)
         cells(nc2) % pnts(3) = nPnts - 3
         cells(nc2) % pnts(4) = nPnts
         cells(nc2) % refineLevel(:) = cells(oc) % refineLevel(:)
         cells(nc2) % neigh(1) = nc1
         cells(nc2) % neigh(2) = cells(oc) % neigh(2)
         nc   = cells(oc ) % neigh(3)
         if (nc > 0) then
            rfl  = cells(nc2) % refineLevel(1)
            nrfl = cells(nc ) % refineLevel(1)
            if ( nrfl < rfl) then
               cells(nc2) % neigh(3) = nc
            else if (nrfl == rfl) then
               nc   = cells(nc) % neigh(2) !! get east neighbor
               nrfl = cells(nc) % refineLevel(1)
               if (nrfl == rfl) then
                  cells(nc2) % neigh(3) = nc
                  cells(nc ) % neigh(4) = nc2
               else
                  write(*,*) "Neighbor SOUTH: Second Neighbor: refinement level not supported",__LINE__,__FILE__
                  write(*,*) nc2, cells(nc2) % refineLevel, cells(nc2) % neigh
                  write(*,*) nc, cells(nc) % refineLevel, cells(nc) % neigh
                  write(*,*) cells(oc ) % neigh(3)
                  stop 1
               end if
            else
               write(*,*) "Neighbor SOUTH with higher refinment level not implemented yet",__LINE__,__FILE__
               stop 1
            end if
         else
            cells(nc2) % neigh(3) = nc
         end if
         cells(nc2) % neigh(4) = nc3
         !! UPDATE NEW CELL
         cells(nc3) % pnts(1) = nPnts
         cells(nc3) % pnts(2) = nPnts - 3
         cells(nc3) % pnts(3) = cells(oc) % pnts(3)
         cells(nc3) % pnts(4) = nPnts - 1
         cells(nc3) % refineLevel(:) = cells(oc) % refineLevel(:)
         cells(nc3) % neigh(1) = oc
         cells(nc3) % neigh(2) = cells(oc) % neigh(2)
         cells(nc3) % neigh(3) = nc2
         cells(nc3) % neigh(4) = cells(oc) % neigh(4)
         !!! UPDATE REFINED CELL
         cells(oc) % pnts(1) = nPnts - 4
         cells(oc) % pnts(2) = nPnts
         cells(oc) % pnts(3) = nPnts - 1
         cells(oc) % neigh(2) = nc3
         cells(oc) % neigh(3) = nc1

         ! Update Neighboor right
         nc = cells(nc3) % neigh(2) 
         if (nc /= 0) &
         cells(nc) % neigh(1) = nc3

         ! Update Neighboor bottom
         nc = cells(nc1) % neigh(3) 
         if (nc /= 0) &
         cells(nc) % neigh(4) = nc1

         do nc = oc, nc3
            cells(nc) % center(:) = 0.25d0 * ( pnts(:,cells(nc) % pnts(1)) & 
                                             + pnts(:,cells(nc) % pnts(2)) & 
                                             + pnts(:,cells(nc) % pnts(3)) & 
                                             + pnts(:,cells(nc) % pnts(4)) ) 
            cells(nc) % var       = cells(oc) % var
            cells(nc) % grad      = cells(oc) % grad
         end do

         ! Parent Cell

         npc1 = nParentCells + 1
         npc2 = nParentCells + 2
         npc3 = nParentCells + 3
         npc4 = nParentCells + 4
         nParentCells = npc4
         opc = cells(oc) % ref                  ! old parent cell
         parentCells(opc) % ref = -1               ! Cell has childs, thus no tin cells array anymore
         parentCells(opc) % cut_type = 1
         parentCells(opc) % child(1) = npc1
         parentCells(opc) % child(2) = npc2
         parentCells(opc) % child(3) = npc3
         parentCells(opc) % child(4) = npc4

         ! 1 Child upper left
         parentCells(npc1) % parent   = opc
         parentCells(npc1) % ref      = oc                                 ! referencing
         cells(oc)         % ref      = npc1
         parentCells(npc1) % child    = -1
         parentCells(npc1) % refineLevel = cells(oc) % refineLevel
         parentCells(npc1) % neigh(1) = parentCells(opc) % neigh(1)
         parentCells(npc1) % neigh(2) = npc4
         parentCells(npc1) % neigh(3) = npc2
         parentCells(npc1) % neigh(4) = parentCells(opc) % neigh(4)

         parentCells(npc2) % parent   = opc
         parentCells(npc2) % ref      = nc1
         cells(nc1)        % ref      = npc2
         parentCells(npc2) % child    = -1
         parentCells(npc2) % neigh(1) = parentCells(opc) % neigh(1)
         parentCells(npc2) % neigh(2) = npc3
         parentCells(npc2) % neigh(3) = parentCells(opc) % neigh(3)
         parentCells(npc2) % neigh(4) = npc1
         parentCells(npc2) % refineLevel = cells(nc) % refineLevel

         parentCells(npc3) % parent   = opc
         parentCells(npc3) % ref      = nc2
         cells(nc2)        % ref      = npc3
         parentCells(npc3) % child    = -1
         parentCells(npc3) % neigh(1) = npc2
         parentCells(npc3) % neigh(2) = parentCells(opc) % neigh(2)
         parentCells(npc3) % neigh(3) = parentCells(opc) % neigh(3)
         parentCells(npc3) % neigh(4) = npc4
         parentCells(npc3) % refineLevel = cells(nc) % refineLevel

         parentCells(npc4) % parent   = opc
         parentCells(npc4) % ref      = nc3
         cells(nc3)        % ref      = npc4
         parentCells(npc4) % child    = -1
         parentCells(npc4) % neigh(1) = npc1
         parentCells(npc4) % neigh(2) = parentCells(opc) % neigh(2)
         parentCells(npc4) % neigh(3) = npc3
         parentCells(npc4) % neigh(4) = parentCells(opc) % neigh(4)
         parentCells(npc4) % refineLevel = cells(nc) % refineLevel
         exit
         write(*,*) "refinment in both direction at same time not supported yet"
         stop 1
      end if
   end do
end do

nCells = nNewCells
end subroutine doRefinement

subroutine check_neighbors (cells,parentCells,nCells,nParentCells,debug_in)
use types
implicit none
type(tCell)   , intent (inout)            :: cells(:)
type(tParentCell)   , intent (inout)      :: parentCells(:)

integer       , intent (inout)            :: nCells
integer       , intent (inout)            :: nParentCells
logical       , intent (in), optional     :: debug_in

logical                                   :: debug

integer                                   :: i,n,ni,nn
if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
end if
do i = 1, nCells
   do n = 1, 4
      if (cells(i) % neigh(n) == 0) cycle
      ni = cells(i) % neigh(n)
      nn = n + mod(n,2)*2-1
      if ( i /= cells(ni) % neigh(nn)) then
         if (n <= 2 .and. cells(i) % refineLevel(2) > cells(ni) % refineLevel(2)) then
         else
            write(*,*) "Cell:",i, "Neighbors:",n,"points to", ni,"which points to", nn,cells(ni) % neigh(nn)
            write(*,*) "Cell:",i ,"Neighbors:",cells(i) % neigh, cells(i) % refineLevel
            write(*,*) "Cell:",ni,"Neighbors:",cells(ni) % neigh, cells(ni) % refineLevel
            write(*,*) parentCells(i), nParentCells
            stop 1
         end if
      end if
   end do
end do
end subroutine check_neighbors

end module refinement
