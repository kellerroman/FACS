module refinement
use const
use types
use face, only: add_face
implicit none

contains
subroutine doRefinement (cells,parentCells,pnts,faces                      &
                        ,nCells,nParentCells,nPnts,nFace                   &
                        ,refineType,refineList,nRefine                     &
                        ,canCoarseList,nCanCoarse                          &
                        ,doCoarseList,nDoCoarse                            &
                        ,holesParentCells,nHolesParentCell                 &
                        ,holesPnts       ,nHolesPnt                        &
                        ,holesFaces      ,nHolesFace                       &
                        ,debug_in)
implicit none
type(tCell)   , intent (inout)            :: cells(:)
type(tParentCell)   , intent (inout)      :: parentCells(:)
real(kind = 8), intent (inout)            :: pnts(:,:)
type(tFace      )   , intent (inout)      :: faces(:)

integer       , intent (inout)            :: nCells
integer       , intent (inout)            :: nParentCells
integer       , intent (inout)            :: nPnts
integer       , intent (inout)            :: nFace

integer       , intent (in)               :: refineType(:)
integer       , intent (inout)            :: refineList(:)
integer       , intent (in)               :: nRefine

integer       , intent (inout)            :: canCoarseList(:)
integer       , intent (inout)            :: nCanCoarse

integer       , intent (inout)            :: doCoarseList(:)
integer       , intent (in)               :: nDoCoarse

integer       , intent (inout)            :: holesParentCells(:)
integer       , intent (inout)            :: nHolesParentCell             

integer       , intent (inout)            :: holesPnts(:)
integer       , intent (inout)            :: nHolesPnt             

integer       , intent (inout)            :: holesFaces(:)
integer       , intent (inout)            :: nHolesFace             

logical       , intent (in), optional     :: debug_in


integer                                   :: nNewCells
integer                                   :: nNewParentCells
integer                                   :: nNewPnts

integer                                   :: nMaxCells
integer                                   :: nMaxParentCells
integer                                   :: nMaxPnts

integer, allocatable                      :: CellId_old2new(:)

! Cell Movement
integer                                   :: ci,ri,dm, nIntervall
integer, allocatable                      :: intervall_starts(:)
integer, allocatable                      :: intervall_ends(:)
integer, allocatable                      :: intervall_dirs(:)
integer, allocatable                      :: intervall_dms(:)


integer                                   :: i,r,ctr,nc,oc,nc1,nc2,nc3
integer                                   :: opc,npc1,npc2,npc3,npc4              ! New Parent cells
integer                                   :: rfl,nrfl
integer                                   :: p(5)

integer                                   :: mf,myFace
logical                                   :: found
logical                                   :: debug


if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
   write(*,*) nHolesFace, holesFaces(1)
end if
nIntervall = nRefine + nDoCoarse
if (nDoCoarse > 0 .and. nRefine > 0) then
   if (max(doCoarseList(nDoCoarse) , refineList(nRefine)) < nCells) nIntervall = nIntervall + 1
else if (nRefine > 0) then
   if (refineList(nRefine) < nCells) nIntervall = nIntervall + 1
else if (nDoCoarse > 0) then
   if (doCoarseList(nDoCoarse) < nCells) nIntervall = nIntervall + 1
else
   return
end if
nMaxCells = ubound(cells,1)
nMaxParentCells = ubound(parentCells,1)
nMaxPnts  = ubound(pnts,2)

nNewPnts        = nPnts        + 5 * (nRefine - nDoCoarse)
nNewCells       = nCells       + 3 * (nRefine - nDoCoarse)
nNewParentCells = nParentCells + 4 * (nRefine - nDoCoarse)


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

allocate (CellId_old2new(NO_CELL:nCells))
! Sort of Refinemnet and Coarseing List
nc1 = nDoCoarse
do while (nc1 > 0)
   nc2 = 0
   do i = 2,nc1 
      if (doCoarseList(i-1) > doCoarseList(i)) then
         nc3 = doCoarseList(i)
         doCoarseList(i) = doCoarseList(i-1)
         doCoarseList(i-1) = nc3
         nc2 = i 
      end if
   end do
   nc1 = nc2
end do

nc1 = nRefine
do while (nc1 > 0)
   nc2 = 0
   do i = 2,nc1 
      if (refineList(i-1) > refineList(i)) then
         nc3 = refineList(i)
         refineList(i) = refineList(i-1)
         refineList(i-1) = nc3
         nc2 = i 
      end if
   end do
   nc1 = nc2
end do

allocate(intervall_starts (nIntervall))
allocate(intervall_ends   (nIntervall))
allocate(intervall_dirs   (nIntervall))
allocate(intervall_dms    (nIntervall))

! to avaoid error List get a additonal element which prevents their furter
! selection
doCoarseList(nDoCoarse+1) = nCells+1
refineList(nRefine+1)     = nCells+1

ci = 1
ri = 1
intervall_starts(1) = 1
intervall_dirs(:) = 1
intervall_ends(nIntervall) = nCells
intervall_dms(1) = 0
do i = 2, nIntervall
   intervall_ends(i-1) = min(doCoarseList(ci),refineList(ri))
   intervall_starts(i) = intervall_ends(i-1) + 1
   if (doCoarseList(ci) > refineList(ri)) then
      ri = ri + 1
      intervall_dms(i) = intervall_dms(i-1) + 3
   else if (doCoarseList(ci) < refineList(ri)) then
      ci = ci + 1
      intervall_dms(i) = intervall_dms(i-1) - 3
      intervall_starts(i) = intervall_starts(i) + 3
      nc = intervall_ends(i - 1)
      !write(*,*) "Cell gets combined",nc,cells(nc) % neigh

      ! Since Cell is overwritten after copy, the information must be copied
      ! before hand.
      cells(nc) % neigh(2) = cells(nc+3) % neigh(2)
      cells(nc) % neigh(3) = cells(nc+1) % neigh(3)
      cells(nc) % pnts(1)  = cells(nc+1) % pnts(1) 
      cells(nc) % pnts(2)  = cells(nc+2) % pnts(2) 
      cells(nc) % pnts(3)  = cells(nc+3) % pnts(3) 
      cells(nc) % refineLevel = cells(nc) % refineLevel - 1
      cells(nc) % Q = 0.25d0 * ( cells(nc  ) % Q &
                               + cells(nc+1) % Q &
                               + cells(nc+2) % Q &
                               + cells(nc+3) % Q )

      cells(nc) % QC = 0.25d0 * ( cells(nc  ) % QC &
                                + cells(nc+1) % QC &
                                + cells(nc+2) % QC &
                                + cells(nc+3) % QC )
      nHolesPnt = nHolesPnt + 5
      holesPnts(nHolesPnt-4)= cells(nc+1) % pnts(2) 
      holesPnts(nHolesPnt-3)= cells(nc+1) % pnts(3) 
      holesPnts(nHolesPnt-2)= cells(nc+1) % pnts(4) 
      holesPnts(nHolesPnt-1)= cells(nc+3) % pnts(2) 
      holesPnts(nHolesPnt  )= cells(nc+3) % pnts(4) 

      ! Since Cells are ignored in the intervalls, the newId must be set here to
      ! the surviving child new position
      CellId_old2new(nc+1) = nc+intervall_dms(i-1)
      CellId_old2new(nc+2) = nc+intervall_dms(i-1)
      CellId_old2new(nc+3) = nc+intervall_dms(i-1)
      !write(*,*) "Cell new neigh ",cells(nc) % neigh
      
   end if
end do
do i = 2, nIntervall
   if (intervall_dms(i) > 0) then
      dm                  = intervall_starts(i)
      intervall_starts(i) = intervall_ends(i)
      intervall_ends(i)   = dm                
      intervall_dirs(i)   = -1
   else
      intervall_dirs(i)   = 1
   end if
end do

nc1 = nIntervall
do while (nc1 > 0)
   nc2 = 0
   do i = 3,nc1 
      if (intervall_starts(i-1) < intervall_starts(i) .and. &
          intervall_dirs(i-1) == -1 .and. &
          intervall_dirs(i)   == -1) then
         nc3 = intervall_starts(i)
         intervall_starts(i) = intervall_starts(i-1)
         intervall_starts(i-1) = nc3

         nc3 = intervall_ends(i)
         intervall_ends(i) = intervall_ends(i-1)
         intervall_ends(i-1) = nc3

         nc3 = intervall_dms(i)
         intervall_dms(i) = intervall_dms(i-1)
         intervall_dms(i-1) = nc3

         nc2 = i 
      end if
   end do
   nc1 = nc2
end do

CellId_old2new(1) = 1
CellId_old2new(NO_CELL) = NO_CELL

do ci = 1, nIntervall
   dm = intervall_dms(ci)
   do i = intervall_starts(ci),intervall_ends(ci),intervall_dirs(ci)
      nc = i + dm 
      if (debug) write(*,*) "Moving Cell from", i, nc , cells(i) % neigh
      CellId_old2new(i) = nc
      cells(nc) = cells(i)
      parentCells(cells(nc) % ref) % ref = nc
   end do
end do

! Updated neighbour information to the new Cells structure
do ci = 1, nIntervall
   dm = intervall_dms(ci)
   do i = intervall_starts(ci),intervall_ends(ci),intervall_dirs(ci)
      nc = CellId_old2new(i)
      !write(*,*) "Updating Cell Neighbors", nc, nCells , cells(nc) % neigh
      do r = 1,4
         cells(nc) % neigh(r) = CellId_old2new(cells(nc) % neigh(r))
      end do
   end do
end do

! Updating Face stencil information
do i = 1, nFace
   faces(i) % stencil(1,1) = CellId_old2new(faces(i) % stencil(1,1))
   faces(i) % stencil(1,2) = CellId_old2new(faces(i) % stencil(1,2))
end do
! **************************************************************************************************
!           
!                COARSENING IN BOTH DIRECTIONS
!           
! **************************************************************************************************
do r = 1, nDoCoarse
   ctr = doCoarseList(r)
   oc = CellId_old2new(ctr)
   cells(oc) % center(:) = 0.25d0 * ( pnts(:,cells(oc) % pnts(1)) & 
                                    + pnts(:,cells(oc) % pnts(2)) & 
                                    + pnts(:,cells(oc) % pnts(3)) & 
                                    + pnts(:,cells(oc) % pnts(4)) ) 
   npc1 = cells(oc) % ref  ! ParentCell Equivalent
   opc  = parentCells(npc1) % parent ! Parent of the Cell, which will be new reference
   npc2 = parentCells(opc) % parent
   do i = 1, 4
      nHolesParentCell = nHolesParentCell + 1
      holesParentCells(nHolesParentCell) = parentCells(opc) % child(i)
   end do
   parentCells(opc) % child = NO_CELL
   parentCells(opc) % ref   = oc
   cells(oc) % ref          = opc
    
   ! delete old parent cell from refinement
   ! If parent exist, the parent maybe can be added to the list avoiding
   ! deleting the entry
   if (npc2 /= NO_CELL) then
      found = .true.
      do i = 1,4
         if (parentCells(parentCells(npc2) % child(i)) % child(1) /= NO_CELL) then
            found = .false.
            exit
         end if
      end do
   else
      found = .false.
   end if

   ! if parent has other children with children or doesnt exist, entry must be
   ! deleted
   ci = parentCells(opc) % pos_CanCoarse
   parentCells(opc) % pos_CanCoarse = NO_CELL
   if (found) then
      canCoarseList(ci) = npc2
      parentCells(npc2) % pos_CanCoarse = ci
   else
      nCanCoarse = nCanCoarse - 1
      do i = ci, nCanCoarse
         canCoarseList(i) = canCoarseList(i+1)
         parentCells(canCoarseList(i)) % pos_CanCoarse = i
      end do
   end if
   
end do
do r = 1, nRefine
   ctr = refineList(r)
   oc = CellId_old2new(ctr)
   if (debug) write(*,*) "Refining Cell: old:",ctr,"new:",oc

! **************************************************************************************************
!           
!                REFINEMENT IN BOTH DIRECTIONS
!           
! **************************************************************************************************
   if (refineType(ctr) == 3) then
      nc1 = oc + 1
      nc2 = oc + 2
      nc3 = oc + 3
      if (debug) then
         write(*,'(30("="),9X,"Working in Cells: ",I0,"->",I0," RefLevel:",2(1X,I0))') oc,nc3,cells(oc) % refineLevel
         write(*,'("Old Neigh:",4(1X,I0,"(",I0,1X,I0,")"))') (cells(oc) % neigh(nc),  &
            cells(cells(oc) % neigh(nc)) % refineLevel(1), &
            cells(cells(oc) % neigh(nc)) % refineLevel(2),nc=1,4)
      end if
      !
      !          ______P4______
      !         |       |      |
      !         |  OC   |  NC3 |          OC Original Cell
      !         |       |      |          NC1 NEW Cell
      !      P1 |_______|______|P2        NC2 NEW Cell
      !         |      P5      |          NC3 NEW Cell
      !         |  NC1  |  NC2 |
      !         |       |      |
      !         |_______|______|
      !                P3 
      !        
      if (debug) write(*,*) "Refining Cell in both Dirs:", ctr, oc
      !!!!  ADD NEW POINTS
      do i = 1, 5
         if (nHolesPnt > 0) then
            p(i) = holesPnts(nHolesPnt) 
            nHolesPnt = nHolesPnt - 1
         else
            nPnts = nPnts + 1
            p(i) = nPnts
         end if
         select case(i)
            case(1)
               Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(4)) ) * 0.5d0
            case(2)
               Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(2)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
            case(3)
               Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(2)) ) * 0.5d0
            case(4)
               Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(4)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
            case(5)
               Pnts(:,p(i)) = cells(oc) % center(:)
         end select
         if (debug) write(*,*) "New Point:", p(i), Pnts(:,p(i))
      end do

      !do nc = 1,4
      !   if (cells(cells(oc) % neigh(nc) ) % refineLevel(1) < cells(oc) % refineLevel(1)) then
      !      write(*,*) "CEll",oc,"is connected to",cells(oc) % neigh(nc), nc
      !      !stop 1
      !   end if
      !end do
      cells(oc) % refineLevel(:) = cells(oc) % refineLevel(:) + 1
      !! UPDATE NEW CELL
      cells(nc1) % pnts(1) = cells(oc) % pnts(1)
      cells(nc1) % pnts(2) = p(3)
      cells(nc1) % pnts(3) = p(5)
      cells(nc1) % pnts(4) = p(1)
      cells(nc1) % refineLevel(:) = cells(oc) % refineLevel(:)

      ! neighbor LEFT
      nc = cells(oc) % neigh(1)
      if (nc /= NO_CELL) then
         rfl  = cells(nc1) % refineLevel(2)
         nrfl = cells(nc ) % refineLevel(2)
         if ( nrfl < rfl) then    ! CASE 1
!
!      ----------------|-----------|-----
!                      |           |
!                      |   oc      |
!                      |           |
!                      |           |
!         nc           [===========]-----
!                      [           ]
!                      [   nc1     ]
!                      [           ]
!                      [           ]
!      ----------------[===========]-----
!
            cells(nc1) % neigh(1) = nc
            if (debug) write(*,*) nc1,"setting left neighbor to", nc
            call add_face(nc1,nc,LEFT,cells,pnts,faces,holesFaces,nHolesFace,nFace,debug)
         else if (nrfl == rfl) then   
            nc   = cells(nc) % neigh(3) !! get south neighbor
            nrfl = cells(nc) % refineLevel(2)
            if (nrfl == rfl) then    ! CASE 2
!
!              ---|---------|--------|-----
!                 |         |        |
!                 |    *    |   oc   |
!                 |   \|/   |        |
!                 |---------[========]-----
!                 |         [        ]
!                 |    nc   [   nc1  ]
!                 |         [        ]
!              ---|---------[========]-----
!
               !modify Right face 1 face of  LEFT neighbor cell
               myFace = cells(nc) % faces_dir(RIGHT,1)
               faces(myFace) % stencil(1,2) = nc1

               ! number of faces in new cell 
               mf = cells(nc1) % nFace + 1
               cells(nc1) % nFace = mf
               cells(nc1) % faces(mf) = myFace
               cells(nc1) % faces_dir_ref(LEFT)   = .false.
               cells(nc1) % faces_dir    (LEFT,1) = myFace
               cells(nc1) % f_sign(mf) = -1.0d0
            
               cells(nc1) % neigh(1) = nc
               if (debug) write(*,*) nc1,"setting left neighbor to", nc
               cells(nc ) % neigh(2) = nc1
               if (debug) write(*,*) nc,"setting right neighbor to", nc1
            else if (nrfl == rfl+1) then!!! CASE 4
!
!              ---|---------|---------|-----
!                 |         |         |
!                 | *       |   oc    |
!                 |\|/      |         |
!                 |----|----[=========]-----
!                 | -> | nc [         ]
!                 |---------[   nc1   ]
!                 |    |    [         ]
!              ---|----|----[=========]-----
!
               nc   = cells(nc) % neigh(2) !! get RIGHT neighbor
               cells(nc1) % neigh(1) = nc
               if (debug) write(*,*) nc1,"setting left neighbor to", nc
               cells(nc ) % neigh(2) = nc1
               if (debug) write(*,*) nc,"setting right neighbor to", nc1
               nc   = cells(nc) % neigh(3) !! get south neighbor
               cells(nc ) % neigh(2) = nc1
               if (debug) write(*,*) nc,"setting right neighbor to", nc1

               !modify Right face 1 face of  LEFT neighbor cell
               myFace = cells(nc) % faces_dir(RIGHT,1)
               faces(myFace) % stencil(1,2) = nc1

               mf = cells(nc1) % nFace + 1
               cells(nc1) % faces_dir_ref(LEFT)   = .true.
               cells(nc1) % faces_dir    (LEFT,1) = myFace
               cells(nc1) % f_sign(mf) = -1.0d0

               !modify Right face 1 face of LEFT LOWER neighbor cell
               nc   = cells(nc) % neigh(SOUTH) !! get RIGHT neighbor
               myFace = cells(nc) % faces_dir(RIGHT,1)
               faces(myFace) % stencil(1,2) = nc1

               ! number of faces in new cell 
               mf = cells(nc1) % nFace + 1
               cells(nc1) % nFace = mf
               cells(nc1) % faces(mf) = myFace
               cells(nc1) % faces_dir    (LEFT,2) = myFace
               cells(nc1) % f_sign(mf) = -1.0d0
            else
               write(*,*) "Neighbor LEFT: Second Neighbor: refinement level not supported",__LINE__,__FILE__
               write(*,*) "old:",ctr,"Cell",oc,"Neighbor",cells(oc) % neigh,"REF:",rfl,nrfl
               nc = cells(oc) % neigh(1)
               write(*,*) "Neighbor from OC:", nc, cells(nc) % refineLevel, cells(nc) % neigh
               nc   = cells(nc) % neigh(3) !! get south neighbor
               write(*,*) "South Neighbor", nc, cells(nc) % refineLevel, cells(nc) % neigh

               write(*,'("neighbor:",I0,"(",I0,1X,I0,") Neigh:",4(1X,I0))') nc, cells(nc) % refineLevel,cells(nc) % neigh
               stop 1
            end if
         else if (nrfl == rfl+1) then
!
!              ---|----|----|--------|-----
!                 |    |\/ *|        |
!                 |----|----|   oc   |
!                 |    |\|/ |        |
!                 |----|----[========]-----
!                 |         [        ]
!                 |    ?    [   nc1  ]
!                 |         [        ]
!              ---|---------[========]-----
!
            nc   = cells(nc) % neigh(3) !! get south neighbor
            nc   = cells(nc) % neigh(3) !! get south neighbor
            nrfl = cells(nc) % refineLevel(2)
            if (nrfl == rfl) then    !  CASE 3
!
!              ---|----|----|--------|-----
!                 |    |\/ *|        |
!                 |----|----|   oc   |
!                 |    |\|/ |        |
!                 |----|----[========]-----
!                 |         [        ]
!                 |    nc   [   nc1  ]
!                 |         [        ]
!              ---|---------[========]-----
!
               cells(nc1) % neigh(1) = nc
               if (debug) write(*,*) nc1,"setting left neighbor to", nc
               cells(nc ) % neigh(2) = nc1
               if (debug) write(*,*) nc,"setting right neighbor to", nc1
               !modify Right face 1 face of  LEFT neighbor cell
               myFace = cells(nc) % faces_dir(RIGHT,1)
               faces(myFace) % stencil(1,2) = nc1

               mf = cells(nc1) % nFace + 1
               cells(nc1) % nFace = mf
               cells(nc1) % faces(mf) = myFace
               cells(nc1) % faces_dir_ref(LEFT)   = .true.
               cells(nc1) % faces_dir    (LEFT,1) = myFace
               cells(nc1) % f_sign(mf) = -1.0d0

            else if (nrfl == rfl+1) then ! CASE 5
!
!              ---|----|----|--------|-----
!                 |    |\/ *|        |
!                 |----|----|   oc   |
!                 |    |\|/ |        |
!                 |----|----[========]-----
!                 |    | nc [        ]
!                 |----|----[   nc1  ]
!                 |    |    [        ]
!              ---|----|----[========]-----
!
               cells(nc1) % neigh(1) = nc
               if (debug) write(*,*) nc1,"setting left neighbor to", nc
               cells(nc ) % neigh(2) = nc1
               if (debug) write(*,*) nc,"setting right neighbor to", nc1
               nc = cells(nc) % neigh(3) !! get south neighbor
               cells(nc ) % neigh(2) = nc1
               if (debug) write(*,*) nc,"setting right neighbor to", nc1
               !modify Right face 1 face of  LEFT neighbor cell
               myFace = cells(nc) % faces_dir(RIGHT,1)
               faces(myFace) % stencil(1,2) = nc1

               mf = cells(nc1) % nFace + 1
               cells(nc1) % nFace = mf
               cells(nc1) % faces(mf) = myFace
               cells(nc1) % faces_dir_ref(LEFT)   = .true.
               cells(nc1) % faces_dir    (LEFT,1) = myFace
               cells(nc1) % f_sign(mf) = -1.0d0

               !modify Right face 1 face of LEFT LOWER neighbor cell
               nc   = cells(nc) % neigh(SOUTH) !! get SOUTH neighbor
               myFace = cells(nc) % faces_dir(RIGHT,1)
               faces(myFace) % stencil(1,2) = nc1

               ! number of faces in new cell 
               mf = cells(nc1) % nFace + 1
               cells(nc1) % nFace = mf
               cells(nc1) % faces(mf) = myFace
               cells(nc1) % faces_dir    (LEFT,2) = myFace
               cells(nc1) % f_sign(mf) = -1.0d0
            else
               write(*,*) "Neighbor West: Second Neighbor RFL/=: refinement level not supported",__LINE__,__FILE__
               write(*,*) nc,nrfl,rfl
               stop 1
            end if
         else
            write(*,*) "Neighbor WEST with higher refinment level not implemented yet",__LINE__,__FILE__
            write(*,*) nc1, cells(nc1) % refineLevel
            write(*,*) nc , cells(nc ) % refineLevel
            stop 1
         end if
      else
         cells(nc1) % neigh(1) = nc
         if (debug) write(*,*) nc1,"setting left neighbor to", nc
         myFace = cells(oc) % faces_dir(RIGHT,1)

         mf = cells(nc1) % nFace + 1
         cells(nc1) % nFace = mf
         cells(nc1) % faces(mf) = myFace
         cells(nc1) % faces_dir_ref(LEFT)   = .true.
         cells(nc1) % faces_dir    (LEFT,1) = myFace
         cells(nc1) % f_sign(mf) = -1.0d0
      end if

      ! NEIGHBOR RIGHT
      cells(nc1) % neigh(RIGHT) = nc2
      if (debug) write(*,*) nc1,"setting right neighbor to", nc2
      call add_face(nc1,nc2,RIGHT,cells,pnts,faces,holesFaces,nHolesFace,nFace,debug)

      !Neighbor SOUTH
      nc = cells(oc) % neigh(3)
      cells(nc1) % neigh(3) = nc
      if (debug) write(*,*) nc1,"setting south neighbor to", nc
      ! FACE
      myFace = cells(oc) % faces_dir(SOUTH,1)
      ! set normal vector of face
      faces(myFace) % n = [0.0d0,1.0d0]
      if (nc /= NO_CELL) then
         if (cells(nc) % neigh(4) == oc) then
            cells(nc) % neigh(4) = nc1
            if (debug) write(*,*) nc,"setting north neighbor to", nc1
            rfl  = cells(nc1) % refineLevel(1)
            nrfl = cells(nc ) % refineLevel(1)
            if (nrfl == rfl + 1) then
               nc = cells(nc) % neigh(2)
               cells(nc) % neigh(4) = nc1
               if (debug) write(*,*) nc,"setting north neighbor to", nc1
            else if (nrfl == rfl) then
            else if (nrfl == rfl -1) then
            else if (nrfl == rfl -2) then
            else
               write(*,*) "SOUTH:, only RFL + 1 supported"
               write(*,*) nc1, cells(nc1) % refineLevel, cells(nc1) % neigh
               write(*,*) nc , cells(nc ) % refineLevel, cells(nc ) % neigh
               stop 1
            end if
         end if
      end if
      cells(nc1) % neigh(4) = oc
      if (debug) write(*,*) nc1,"setting north neighbor to", oc

      !! UPDATE NEW CELL

      cells(nc2) % pnts(1) = p(3)
      cells(nc2) % pnts(2) = cells(oc) % pnts(2)
      cells(nc2) % pnts(3) = p(2)
      cells(nc2) % pnts(4) = p(5)
      cells(nc2) % refineLevel(:) = cells(oc) % refineLevel(:)
      !Neighbor LEFT
      cells(nc2) % neigh(1) = nc1
      if (debug) write(*,*) nc2,"setting left neighbor to", nc1
      !Neighbor RIGHT
      nc = cells(oc) % neigh(2)
      if (nc /= NO_CELL) then
         rfl  = cells(nc2) % refineLevel(2)
         nrfl = cells(nc ) % refineLevel(2)
         if ( nrfl < rfl) then
            cells(nc2) % neigh(2) = nc
            if (debug) write(*,*) nc2,"setting right neighbor to", nc
         else if (nrfl == rfl) then
            nc   = cells(nc) % neigh(3) !! get south neighbor
            nrfl = cells(nc) % refineLevel(2)
            cells(nc2) % neigh(2) = nc
            if (debug) write(*,*) nc2,"setting right neighbor to", nc
            cells(nc ) % neigh(1) = nc2
            if (debug) write(*,*) nc,"setting left neighbor to", nc2
            if (nrfl == rfl) then
            else if (nrfl == rfl+1) then
               nc   = cells(nc) % neigh(3) !! get south neighbor
               cells(nc ) % neigh(1) = nc2
               if (debug) write(*,*) nc,"setting left neighbor to", nc2
            else
               write(*,*) "Neighbor RIGHT: Second Neighbor: refinement level not supported",__LINE__,__FILE__
               write(*,*) nc2, cells(nc2) % refineLevel, cells(nc2) % neigh
               write(*,*) nc, cells(nc) % refineLevel, cells(nc) % neigh
               write(*,*) cells(oc ) % neigh(2)
               stop 1
            end if
         else if (nrfl == rfl + 1) then
            nc   = cells(nc) % neigh(3) !! get south neighbor
            nc   = cells(nc) % neigh(3) !! get south neighbor
            nrfl = cells(nc) % refineLevel(2)
            cells(nc2) % neigh(2) = nc
            if (debug) write(*,*) nc2,"setting right neighbor to", nc
            cells(nc ) % neigh(1) = nc2
            if (debug) write(*,*) nc,"setting left neighbor to", nc2
            if (nrfl == rfl) then
            else if (nrfl == rfl+1) then
               nc = cells(nc) % neigh(3) !! get south neighbor
               cells(nc ) % neigh(1) = nc2
               if (debug) write(*,*) nc,"setting left neighbor to", nc2
            else
               write(*,*) "Neighbor RIGHT: Second Neighbor RFL/=: refinement level not supported",__LINE__,__FILE__
               write(*,*) nc,nrfl,rfl
               stop 1
            end if
            
         else
            write(*,*) "Neighbor RIGHT with higher refinment level not implemented yet",__LINE__,__FILE__
            stop 1
         end if
      else
         cells(nc2) % neigh(2) = nc
         if (debug) write(*,*) nc2,"setting right neighbor to", nc
      end if

      ! Neighbor SOUTH
      nc = cells(oc) % neigh(3)
      if (nc /= NO_CELL) then
         rfl  = cells(nc2) % refineLevel(1)
         nrfl = cells(nc ) % refineLevel(1)
         if ( nrfl < rfl) then
            cells(nc2) % neigh(3) = nc
            if (debug) write(*,*) nc2,"setting south neighbor to", nc
         else if (nrfl == rfl) then
            nc   = cells(nc) % neigh(2) !! get right neighbor
            nrfl = cells(nc) % refineLevel(1)
            cells(nc2) % neigh(3) = nc
            if (debug) write(*,*) nc2,"setting south neighbor to", nc
            cells(nc ) % neigh(4) = nc2
            if (debug) write(*,*) nc,"setting north neighbor to", nc2
            if (nrfl == rfl) then
            else if (nrfl == rfl+1) then
               nc   = cells(nc) % neigh(2) !! get right neighbor
               cells(nc ) % neigh(4) = nc2
               if (debug) write(*,*) nc,"setting north neighbor to", nc2
            else
               write(*,*) "Neighbor SOUTH: Second Neighbor: refinement level not supported",__LINE__,__FILE__
               write(*,*) nc2, cells(nc2) % refineLevel, cells(nc2) % neigh
               write(*,*) nc, cells(nc) % refineLevel, cells(nc) % neigh
               write(*,*) cells(oc ) % neigh(3)
               stop 1
            end if
         else if (nrfl == rfl + 1) then
            nc   = cells(nc) % neigh(2) !! get right neighbor
            nc   = cells(nc) % neigh(2) !! get right neighbor
            nrfl = cells(nc) % refineLevel(1)
            cells(nc2) % neigh(3) = nc
            if (debug) write(*,*) nc2,"setting south neighbor to", nc
            cells(nc ) % neigh(4) = nc2
            if (debug) write(*,*) nc,"setting north neighbor to", nc2
            if (nrfl == rfl) then
            else if (nrfl == rfl+1) then
               nc = cells(nc) % neigh(2) !! get right neighbor
               cells(nc ) % neigh(4) = nc2
               if (debug) write(*,*) nc,"setting north neighbor to", nc2
            else
               write(*,*) "Neighbor SOUTH: Second Neighbor RFL/=: refinement level not supported",__LINE__,__FILE__
               write(*,*) nc,nrfl,rfl
               stop 1
            end if
            
         else
            write(*,*) "Neighbor SOUTH with higher refinment level not implemented yet",__LINE__,__FILE__
            stop 1
         end if
      else
         cells(nc2) % neigh(3) = nc
         if (debug) write(*,*) nc2,"setting south neighbor to", nc
      end if
      ! Neighbor NORTH
      cells(nc2) % neigh(4) = nc3
      if (debug) write(*,*) nc2,"setting north neighbor to", nc3

      !! UPDATE NEW CELL NORTH RIGHT
      cells(nc3) % pnts(1) = p(5)
      cells(nc3) % pnts(2) = p(2)
      cells(nc3) % pnts(3) = cells(oc) % pnts(3)
      cells(nc3) % pnts(4) = p(4)
      cells(nc3) % refineLevel(:) = cells(oc) % refineLevel(:)
      cells(nc3) % neigh(1) = oc
      if (debug) write(*,*) nc3,"setting left neighbor to", oc

      ! Neighbor RIGHT
      nc = cells(oc) % neigh(2)
      cells(nc3) % neigh(2) = nc
      if (debug) write(*,*) nc3,"setting right neighbor to", nc
      if (nc /= NO_CELL) then
         if (cells(nc) % neigh(1) == oc) then
            cells(nc) % neigh(1) = nc3
            if (debug) write(*,*) nc,"setting left neighbor to", nc3
            rfl  = cells(nc3) % refineLevel(2)
            nrfl = cells(nc ) % refineLevel(2)
            if (nrfl > rfl) then
               nc = cells(nc) % neigh(3)
               cells(nc) % neigh(1) = nc3
               if (debug) write(*,*) nc,"setting left neighbor to", nc3
            end if
         end if
      end if
      ! Neighbor SOUTH
      cells(nc3) % neigh(3) = nc2
      if (debug) write(*,*) nc3,"setting south neighbor to", nc2
      ! Neighbor NORTH
      nc   = cells(oc ) % neigh(4)
      if (nc /= NO_CELL) then
         rfl  = cells(nc3) % refineLevel(1)
         nrfl = cells(nc ) % refineLevel(1)
         if (debug) write(*,*) "cell",nc3,"ref:",rfl,"Neigh N",nc,"Ref:",nrfl
         if ( nrfl < rfl) then
            cells(nc3) % neigh(4) = nc
            if (debug) write(*,*) nc3,"setting north neighbor to", nc
         else if (nrfl == rfl) then
            nc   = cells(nc) % neigh(2) !! get right neighbor
            nrfl = cells(nc) % refineLevel(1)
            if (nrfl == rfl) then
               cells(nc3) % neigh(4) = nc
               if (debug) write(*,*) nc3,"setting north neighbor to", nc
               cells(nc ) % neigh(3) = nc3
               if (debug) write(*,*) nc,"setting south neighbor to", nc3
            else if (nrfl == rfl + 1) then
               nc   = cells(nc) % neigh(3) !! get south neighbor
               cells(nc3) % neigh(4) = nc
               if (debug) write(*,*) nc3,"setting north neighbor to", nc
               cells(nc ) % neigh(3) = nc3
               if (debug) write(*,*) nc,"setting south neighbor to", nc3
               nc   = cells(nc) % neigh(2) !! get right neighbor
               cells(nc ) % neigh(3) = nc3
               if (debug) write(*,*) nc,"setting south neighbor to", nc3
            else
               write(*,*) "Neighbor NORTH: Second Neighbor: refinement level not supported",__LINE__,__FILE__
               write(*,*) nc3, cells(nc3) % refineLevel, cells(nc3) % neigh
               write(*,*) nc, cells(nc) % refineLevel, cells(nc) % neigh
               write(*,*) cells(oc ) % neigh(4)
               stop 1
            end if
         else if (nrfl == rfl + 1) then
            nc   = cells(nc) % neigh(2) !! get right neighbor
            if (debug) write(*,*) "looking at",nc
            nc   = cells(nc) % neigh(2) !! get right neighbor
            nrfl = cells(nc) % refineLevel(1)
            if (debug) write(*,*) "looking at",nc,nrfl
            if (nrfl == rfl) then
               cells(nc3) % neigh(4) = nc
               if (debug) write(*,*) nc3,"setting north neighbor to", nc
               cells(nc ) % neigh(3) = nc3
               if (debug) write(*,*) nc,"setting south neighbor to", nc3
            else if (nrfl == rfl + 1) then
               cells(nc3) % neigh(4) = nc
               if (debug) write(*,*) nc3,"setting north neighbor to", nc
               cells(nc ) % neigh(3) = nc3
               if (debug) write(*,*) nc,"setting south neighbor to", nc3
               nc   = cells(nc) % neigh(2) !! get right neighbor
               cells(nc ) % neigh(3) = nc3
               if (debug) write(*,*) nc,"setting south neighbor to", nc3
            else 
               write(*,*) "somthing going wrong"
               stop 1
            end if
         else
            write(*,*) "Neighbor NORTH with higher refinment level not implemented yet",__LINE__,__FILE__
            stop 1
         end if
      else
         cells(nc3) % neigh(4) = nc
         if (debug) write(*,*) nc3,"setting north neighbor to", nc
      end if
      !!! UPDATE REFINED CELL
      cells(oc) % pnts(1) = p(1)
      cells(oc) % pnts(2) = p(5)
      cells(oc) % pnts(3) = p(4)
      cells(oc) % neigh(2) = nc3
      cells(oc) % neigh(3) = nc1

      do nc = oc, nc3
         cells(nc) % center(:) = 0.25d0 * ( pnts(:,cells(nc) % pnts(1)) & 
                                          + pnts(:,cells(nc) % pnts(2)) & 
                                          + pnts(:,cells(nc) % pnts(3)) & 
                                          + pnts(:,cells(nc) % pnts(4)) ) 
         cells(nc) % Q         = cells(oc) % Q  
         cells(nc) % QC        = cells(oc) % QC
         cells(nc) % grad      = cells(oc) % grad
      end do

      ! Parent Cell

      if (nHolesParentCell > 0) then
         npc1 = holesParentCells(nHolesParentCell) 
         nHolesParentCell = nHolesParentCell - 1
      else
         nParentCells = nParentCells + 1
         npc1 = nParentCells
      end if
         
      if (nHolesParentCell > 0) then
         npc2 = holesParentCells(nHolesParentCell) 
         nHolesParentCell = nHolesParentCell - 1
      else
         nParentCells = nParentCells + 1
         npc2 = nParentCells
      end if
      if (nHolesParentCell > 0) then
         npc3 = holesParentCells(nHolesParentCell) 
         nHolesParentCell = nHolesParentCell - 1
      else
         nParentCells = nParentCells + 1
         npc3 = nParentCells
      end if
      if (nHolesParentCell > 0) then
         npc4 = holesParentCells(nHolesParentCell) 
         nHolesParentCell = nHolesParentCell - 1
      else
         nParentCells = nParentCells + 1
         npc4 = nParentCells
      end if
      opc = cells(oc) % ref                      ! old parent cell
      parentCells(opc) % ref = NO_CELL                ! Cell has childs, thus no tin cells array anymore
      parentCells(opc) % cut_type = 3
      parentCells(opc) % child(1) = npc1
      parentCells(opc) % child(2) = npc2
      parentCells(opc) % child(3) = npc3
      parentCells(opc) % child(4) = npc4
      if (parentCells(opc) % parent == NO_CELL) then  ! Cell on first level, no need to delet an old cell 
                                                      ! from the list (the parent)
         found = .false.
      else !parent can be in the list
         nc = parentCells(parentCells(opc) % parent) % pos_CanCoarse
         if (nc /= NO_CELL) then
            canCoarseList(nc) = opc
            parentCells(opc) % pos_CanCoarse = nc
            parentCells(parentCells(opc) % parent) % pos_CanCoarse = NO_CELL
            found = .true.
         else
            found = .false.
         end if
      end if
      if (.not. found) then
         nCanCoarse = nCanCoarse + 1 
         canCoarseList(nCanCoarse) = opc
         parentCells(opc) % pos_CanCoarse = nCanCoarse
      end if

      ! 1 Child upper left
      parentCells(npc1) % parent   = opc
      parentCells(npc1) % ref      = oc                                 ! referencing
      cells(oc)         % ref      = npc1
      parentCells(npc1) % child    = NO_CELL
      parentCells(npc1) % refineLevel = cells(oc) % refineLevel
      parentCells(npc1) % neigh(1) = parentCells(opc) % neigh(1)
      parentCells(npc1) % neigh(2) = npc4
      parentCells(npc1) % neigh(3) = npc2
      parentCells(npc1) % neigh(4) = parentCells(opc) % neigh(4)

      parentCells(npc2) % parent   = opc
      parentCells(npc2) % ref      = nc1
      cells(nc1)        % ref      = npc2
      parentCells(npc2) % child    = NO_CELL
      parentCells(npc2) % neigh(1) = parentCells(opc) % neigh(1)
      parentCells(npc2) % neigh(2) = npc3
      parentCells(npc2) % neigh(3) = parentCells(opc) % neigh(3)
      parentCells(npc2) % neigh(4) = npc1
      parentCells(npc2) % refineLevel = cells(oc) % refineLevel

      parentCells(npc3) % parent   = opc
      parentCells(npc3) % ref      = nc2
      cells(nc2)        % ref      = npc3
      parentCells(npc3) % child    = NO_CELL
      parentCells(npc3) % neigh(1) = npc2
      parentCells(npc3) % neigh(2) = parentCells(opc) % neigh(2)
      parentCells(npc3) % neigh(3) = parentCells(opc) % neigh(3)
      parentCells(npc3) % neigh(4) = npc4
      parentCells(npc3) % refineLevel = cells(oc) % refineLevel

      parentCells(npc4) % parent   = opc
      parentCells(npc4) % ref      = nc3
      cells(nc3)        % ref      = npc4
      parentCells(npc4) % child    = NO_CELL
      parentCells(npc4) % neigh(1) = npc1
      parentCells(npc4) % neigh(2) = parentCells(opc) % neigh(2)
      parentCells(npc4) % neigh(3) = npc3
      parentCells(npc4) % neigh(4) = parentCells(opc) % neigh(4)
      parentCells(npc4) % refineLevel = cells(oc) % refineLevel

      !call check_neighbors(cells,nNewCells,pnts)
   else
      write(*,*) "Refining Type not supported",__LINE__,__FILE__
      stop 1
   end if
end do

write(*,'("Number of Refine/Coarse              : ",I0,"/",I0)') nRefine,nDoCoarse
write(*,'("Number of Cells/Parents/Points       : ",I0,2("/",I0))') nNewCells &
                                                                  , nParentCells - nHolesParentCell &
                                                                  , nPnts - nHolesPnt
write(*,'("Arraysize req Cells/Parents/Points   : ",I0,2("/",I0))') nNewCells,nParentCells, nPnts
write(*,'("Number of Holes Parent/Pnt           : ",I0,"/",I0)') nHolesParentCell,nHolesPnt
nCells = nNewCells
end subroutine doRefinement

subroutine check_neighbors (cells,nCells,pnts,debug_in)
implicit none
real(kind = 8), parameter                 :: EPSI = 1.0D-10
type(tCell)   , intent (inout)            :: cells(:)
real(kind = 8), intent (inout)            :: pnts(:,:)

integer       , intent (inout)            :: nCells
logical       , intent (in), optional     :: debug_in

logical                                   :: debug

integer                                   :: i,n,ni,nn,ip1,ip2,inp1,inp2
real(kind = 8)                            :: p1(2), np1(2), p2(2) , np2(2)
if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
end if
do i = 1, nCells
   do n = 1, 4                ! over all dircetions (neighbors)
      if (cells(i) % neigh(n) == NO_CELL) cycle
      ni = cells(i) % neigh(n)    ! neighbor index
      nn = n + mod(n,2)*2-1       ! neighbor Direction
      if (n == 1) then
         ip1 = cells(i) % pnts(1)
         ip2 = cells(i) % pnts(4)
      else if (n == 2) then
         ip1 = cells(i) % pnts(2)
         ip2 = cells(i) % pnts(3)
      else if (n == 3) then
         ip1 = cells(i) % pnts(2)
         ip2 = cells(i) % pnts(1)
      else
         ip1 = cells(i) % pnts(3)
         ip2 = cells(i) % pnts(4)
      end if
      if (nn == 1) then
         inp1 = cells(ni) % pnts(1)
         inp2 = cells(ni) % pnts(4)
      else if (nn == 2) then
         inp1 = cells(ni) % pnts(2)
         inp2 = cells(ni) % pnts(3)
      else if (nn == 3) then
         inp1 = cells(ni) % pnts(2)
         inp2 = cells(ni) % pnts(1)
      else
         inp1 = cells(ni) % pnts(3)
         inp2 = cells(ni) % pnts(4)
      end if
      p2(:)  = pnts(:,ip2)
      np2(:) = pnts(:,inp2)
      if (cells(ni) % refineLevel(1) == cells(i) % refineLevel(1) .and. &
          cells(ni) % refineLevel(2) == cells(i) % refineLevel(2) ) then
         p1(:)  = pnts(:,ip1 )
         np1(:) = pnts(:,inp1)
      else if (cells(ni) % refineLevel(1) == 1 + cells(i) % refineLevel(1) .and. &
               cells(ni) % refineLevel(2) == 1 + cells(i) % refineLevel(2) ) then
               ! neighbor cell is more refined
         p1(:)  = (pnts(:,ip1 ) + pnts(:,ip2)  ) * 0.5d0
         np1(:) =  pnts(:,inp1)
      else if (cells(ni) % refineLevel(1) + 1 == cells(i) % refineLevel(1) .and. &
               cells(ni) % refineLevel(2) + 1 == cells(i) % refineLevel(2) ) then
               ! cell is more refined
         p1(:)  =  pnts(:,ip1 )
         if (abs(p2(1) - np2(1)) > EPSI .or. &
             abs(p2(2) - np2(2)) > EPSI ) then
            np2(:) = (pnts(:,inp1) + pnts(:,inp2) ) * 0.5d0
            np1(:) = pnts(:,inp1)
         else
            np1(:) = (pnts(:,inp1) + pnts(:,inp2) ) * 0.5d0
         end if
      else
         !write(*,*) "Different ref Level not  yet supp", i,ni
         !write(*,*) "Cell:", i, "NEIGHBOR", ni, "Direction", n
         !write(*,*) cells(i) % refineLevel, cells(ni) % refineLevel
         !write(*,*) pnts(:,i)
         !stop 1
         cycle
      end if
      if (abs(p1(1) - np1(1)) > EPSI .or. &
          abs(p1(2) - np1(2)) > EPSI .or. &
          abs(p2(1) - np2(1)) > EPSI .or. &
          abs(p2(2) - np2(2)) > EPSI ) then
         write(*,*) "Neighbor Connection is wrong"
         write(*,*) "Connection from",i,"to",ni,"Direction", n
         write(*,*) "Refinement Level",cells(i) % refineLevel, cells(ni) % refineLevel
         write(*,*) "      MyCell:",i, "Neighbor:",cells(i) % neigh
         write(*,*) p1, p2
         write(*,*) "NeighborCell:",ni,"Neighbor:",cells(ni) % neigh
         write(*,*) np1,np2
         write(*,*) nCells
         stop 1
      end if

   end do
end do
write(*,*) "Neighbor cells checked", nCells
end subroutine check_neighbors

end module refinement
