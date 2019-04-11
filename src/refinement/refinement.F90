module refinement
use const
use types
use face, only: add_face
use array_holes
implicit none

contains
subroutine doRefinement (cells,parentCells,pnts,faces                      &
                        ,nCells,nParentCells,nPnts,nFace                   &
                        ,refineType,refineList,nRefine                     &
                        ,canCoarseList,nCanCoarse                          &
                        ,doCoarseList,nDoCoarse                            &
                        ,debug_in)
implicit none
type(tCell)   , intent (inout)            :: cells(:)
type(tParentCell)   , intent (inout)      :: parentCells(:)
real(kind = 8), intent (inout)            :: pnts(:,:)
type(tFace      )   , intent (inout)      :: faces(:)

type(holes)   , intent (inout)            :: nCells
type(holes)   , intent (inout)            :: nParentCells
type(holes)   , intent (inout)            :: nPnts
type(holes)   , intent (inout)            :: nFace

integer       , intent (in)               :: refineType(:)
integer       , intent (inout)            :: refineList(:)
integer       , intent (in)               :: nRefine

integer       , intent (inout)            :: canCoarseList(:)
integer       , intent (inout)            :: nCanCoarse

integer       , intent (inout)            :: doCoarseList(:)
integer       , intent (in)               :: nDoCoarse

logical       , intent (in), optional     :: debug_in


integer                                   :: nNewCells
integer                                   :: nNewParentCells
integer                                   :: nNewPnts

integer                                   :: nMaxCells
integer                                   :: nMaxParentCells
integer                                   :: nMaxPnts

integer                                   :: i,r,nc,oc,nc1,nc2,nc3,ci
integer                                   :: opc,npc1,npc2,npc3,npc4              ! New Parent cells
integer                                   :: rfl,nrfl
integer                                   :: p(5)

integer                                   :: mf,myFace
integer                                   :: new_Cells(4)
logical                                   :: found
logical                                   :: debug


if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
   !write(*,*) nHolesFace, holesFaces(1)
end if

nMaxCells       = ubound(cells,1)
nMaxParentCells = ubound(parentCells,1)
nMaxPnts        = ubound(pnts,2)

nNewPnts        = nPnts%nEntry        + 5 * (nRefine - nDoCoarse)
nNewCells       = nCells%nEntry       + 3 * (nRefine - nDoCoarse)
nNewParentCells = nParentCells%nEntry + 4 * (nRefine - nDoCoarse)


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

! **************************************************************************************************
!           
!                COARSENING IN BOTH DIRECTIONS
!           
! **************************************************************************************************
do r = 1, nDoCoarse
   oc = doCoarseList(r)
   cells(oc) % center(:) = 0.25d0 * ( pnts(:,cells(oc) % pnts(1)) & 
                                    + pnts(:,cells(oc) % pnts(2)) & 
                                    + pnts(:,cells(oc) % pnts(3)) & 
                                    + pnts(:,cells(oc) % pnts(4)) ) 
   npc1 = cells(oc) % ref  ! ParentCell Equivalent
   opc  = parentCells(npc1) % parent ! Parent of the Cell, which will be new reference
   npc2 = parentCells(opc) % parent
   do i = 1, 4
      !nHolesParentCell = nHolesParentCell + 1
      !holesParentCells(nHolesParentCell) = parentCells(opc) % child(i)
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

call sort_refine_list(nRefine,refineList,cells)



do r = 1, nRefine
   oc = refineList(r)

   if (debug) write(*,*) "Refining Cell: ",oc

! **************************************************************************************************
!           
!                REFINEMENT IN BOTH DIRECTIONS
!           
! **************************************************************************************************
   if (refineType(oc) == 3) then
      nc1 = nCells % newEntry()
      nc2 = nCells % newEntry()
      nc3 = nCells % newEntry()
      if (debug) then
         write(*,*) "Refining Cell in both Dirs",cells(oc) % refineLevel

         write(*,'(30("="),9X,"New Cells: ",3(I0,1X))') nc1,nc2,nc3
         write(*,'("Old Neigh:",4(1X,I0,"(",I0,1X,I0,")"))') (cells(oc) % neigh(nc),  &
            cells(cells(oc) % neigh(nc)) % refineLevel(1), &
            cells(cells(oc) % neigh(nc)) % refineLevel(2),nc=1,4)
      end if
      !
      !          ______P4______
      !         |       |      |
      !         |  OC   |  NC3 |          OC Original Cell
      !         |       |      |          NC1 NEW Cell
      !      P1 |_______|______| P2       NC2 NEW Cell
      !         |      P5      |          NC3 NEW Cell
      !         |  NC1  |  NC2 |
      !         |       |      |
      !         |_______|______|
      !                P3 
      !        
      !!!!  ADD NEW POINTS
      npc1 = cells(oc) % refineLevel(1)
      npc2 = cells(oc) % refineLevel(2)
      do i = 1, 5
         select case(i)
            case(1)
               nc = cells(oc) % neigh(1)
               if (nc == NO_CELL) then
                   npc3 = 0
               else
                   npc3 = cells(nc) % refineLevel(2)
               end if
               write(*,*) "POINTS:",oc,nc,npc2,npc3
               if (npc3 <= npc2) then
                  p(i) = nPnts % newEntry()
                  Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(4)) ) * 0.5d0
               else
                   !nc = cells(nc) % neigh(3)
                   p(i) = cells(nc) % pnts(2)
               end if
            case(2)
               nc = cells(oc) % neigh(2)
               if (nc == NO_CELL) then
                   npc3 = 0
               else
                   npc3 = cells(nc) % refineLevel(2)
               end if
               write(*,*) "POINTS:",oc,nc,npc2,npc3
               if (npc3 <= npc2) then
                  p(i) = nPnts % newEntry()
                  Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(2)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
               else
                   !nc = cells(nc) % neigh(3)
                   p(i) = cells(nc) % pnts(1)
               end if
            case(3)
               nc = cells(oc) % neigh(3)
               if (nc == NO_CELL) then
                   npc3 = 0
               else
                   npc3 = cells(nc) % refineLevel(1)
               end if
               write(*,*) "POINTS:",oc,nc,npc1,npc3
               if (npc3 <= npc1) then
                  p(i) = nPnts % newEntry()
                  Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(1)) + Pnts(:,cells(oc) % pnts(2)) ) * 0.5d0
               else
                   !nc = cells(nc) % neigh(2)
                   p(i) = cells(nc) % pnts(3)
               end if
            case(4)
               nc = cells(oc) % neigh(4)
               if (nc == NO_CELL) then
                   npc3 = 0
               else
                   npc3 = cells(nc) % refineLevel(1)
               end if
               write(*,*) "POINTS:",oc,nc,npc1,npc3
               if (npc3 <= npc1) then
                  p(i) = nPnts % newEntry()
                  Pnts(:,p(i)) = ( Pnts(:,cells(oc) % pnts(4)) + Pnts(:,cells(oc) % pnts(3)) ) * 0.5d0
               else
                   !nc = cells(nc) % neigh(2)
                   p(i) = cells(nc) % pnts(2)
               end if
            case(5)
               p(i) = nPnts % newEntry()
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

      !! UPDATING POINTS OF CELLS
      cells(nc1) % pnts(1) = cells(oc) % pnts(1)
      cells(nc1) % pnts(2) = p(3)
      cells(nc1) % pnts(3) = p(5)
      cells(nc1) % pnts(4) = p(1)
      cells(nc1) % refineLevel(:) = cells(oc) % refineLevel(:) + 1

      cells(nc2) % pnts(1) = p(3)
      cells(nc2) % pnts(2) = cells(oc) % pnts(2)
      cells(nc2) % pnts(3) = p(2)
      cells(nc2) % pnts(4) = p(5)
      cells(nc2) % refineLevel(:) = cells(oc) % refineLevel(:) + 1

      cells(nc3) % pnts(1) = p(5)
      cells(nc3) % pnts(2) = p(2)
      cells(nc3) % pnts(3) = cells(oc) % pnts(3)
      cells(nc3) % pnts(4) = p(4)
      cells(nc3) % refineLevel(:) = cells(oc) % refineLevel(:) + 1

      cells(oc ) % pnts(1) = p(1)
      cells(oc ) % pnts(2) = p(5)
      cells(oc ) % pnts(3) = p(4)
      cells(oc ) % refineLevel(:) = cells(oc) % refineLevel(:) + 1
      ! neighbor LEFT

      nc = cells(oc) % neigh(1)
      if (debug) write(*,*) "=== NEIGHBOR LEFT"
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
            !call add_face(nc1,nc,LEFT,cells,pnts,faces,holesFaces,nHolesFace,nFace,debug)
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
              ! !modify Right face 1 face of  LEFT neighbor cell
              ! myFace = cells(nc) % faces_dir(RIGHT,1)
              ! faces(myFace) % stencil(1,2) = nc1

              ! ! number of faces in new cell 
              ! mf = cells(nc1) % nFace + 1
              ! cells(nc1) % nFace = mf
              ! cells(nc1) % faces(mf) = myFace
              ! cells(nc1) % faces_dir_ref(LEFT)   = .false.
              ! cells(nc1) % faces_dir    (LEFT,1) = myFace
              ! cells(nc1) % f_sign(mf) = -1.0d0
            
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
               write(*,*) "Neighbor LEFT: Second Neighbor: refinement level not supported" &
                   ,__LINE__,__FILE__
               write(*,*) "Cell",oc,"Neighbor",cells(oc) % neigh,"REF:",rfl,nrfl
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
               write(*,*) "Neighbor West: Second Neighbor RFL/=: refinement level not supported" &
                   ,__LINE__,__FILE__
               write(*,*) nc,nrfl,rfl
               stop 1
            end if
         else
            write(*,*) "Neighbor WEST with higher refinment level not implemented yet" &
                ,__LINE__,__FILE__
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
      if (debug) write(*,*) "=== NEIGHBOR RIGHT"
      cells(nc1) % neigh(RIGHT) = nc2
      if (debug) write(*,*) nc1,"setting right neighbor to", nc2
      !call add_face(nc1,nc2,RIGHT,cells,pnts,faces,holesFaces,nHolesFace,nFace,debug)

      !Neighbor SOUTH
      if (debug) write(*,*) "=== NEIGHBOR SOUTH"
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
               write(*,*) "Neighbor RIGHT: Second Neighbor: refinement level not supported" &
                   ,__LINE__,__FILE__
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
               write(*,*) "Neighbor RIGHT: Second Neighbor RFL/=: refinement level not supported" &
                   ,__LINE__,__FILE__
               write(*,*) nc,nrfl,rfl
               stop 1
            end if
            
         else
            write(*,*) "Neighbor RIGHT with higher refinment level not implemented yet" &
                ,__LINE__,__FILE__
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
               write(*,*) "Neighbor SOUTH: Second Neighbor: refinement level not supported" &
                   ,__LINE__,__FILE__
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
               write(*,*) "Neighbor SOUTH: Second Neighbor RFL/=: refinement level not supported" &
                   ,__LINE__,__FILE__
               write(*,*) nc,nrfl,rfl
               stop 1
            end if
            
         else
            write(*,*) "Neighbor SOUTH with higher refinment level not implemented yet" &
                ,__LINE__,__FILE__
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
      if (debug) write(*,*) nc3," setting south neighbor to", nc2
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
               write(*,*) "Neighbor NORTH: Second Neighbor: refinement level not supported" &
                   ,__LINE__,__FILE__
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
            write(*,*) "Neighbor NORTH with higher refinment level not implemented yet" &
                ,__LINE__,__FILE__
            stop 1
         end if
      else
         cells(nc3) % neigh(4) = nc
         if (debug) write(*,*) nc3,"setting north neighbor to", nc
      end if
      !!! UPDATE REFINED CELL
      if (debug) write(*,*) "Updating Refined Cells"
      cells(oc) % neigh(2) = nc3
      cells(oc) % neigh(3) = nc1

      new_Cells = [oc,nc1,nc2,nc3]
      do i = 1,4
         nc = new_cells(i)
         if (debug) write(*,*) "Updating Refined Cells", nc
         cells(nc) % center(:) = 0.25d0 * ( pnts(:,cells(nc) % pnts(1)) & 
                                          + pnts(:,cells(nc) % pnts(2)) & 
                                          + pnts(:,cells(nc) % pnts(3)) & 
                                          + pnts(:,cells(nc) % pnts(4)) ) 
         cells(nc) % Q         = cells(oc) % Q  
         cells(nc) % QC        = cells(oc) % QC
         cells(nc) % grad      = cells(oc) % grad
      end do

      ! Parent Cell

      npc1 = nParentCells % newEntry()
      npc2 = nParentCells % newEntry()
      npc3 = nParentCells % newEntry()
      npc4 = nParentCells % newEntry()

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
                                                                  , nParentCells%nEntry &
                                                                  , nPnts % nEntry
!write(*,'("Number of Holes Parent/Pnt           : ",I0,"/",I0)') nHolesParentCell,nHolesPnt
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

subroutine check_points (pnts,nPnts,debug_in)
implicit none
real(kind = 8), intent (in)               :: pnts(:,:)

type(holes)   , intent (in)               :: nPnts

logical       , intent (in), optional     :: debug_in

logical                                   :: debug

integer                                   :: i,j, np

if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
end if
np = nPnts % nEntry
do i = 1, np
    if (debug) write(*,*) "= checking",i,pnts(:,i)
    do j = i+1, np
        if (debug) write(*,*) "==  with",j,pnts(:,j)
        if (pnts(1,i) == pnts(1,j) .and. pnts(2,i) == pnts(2,j)) then
            write(*,*) "POINTS ALREADY EXISTS"
            write(*,*) i,j
            write(*,*) pnts(:,i)
            !stop 1
        end if
    end do
end do

end subroutine check_points

subroutine sort_refine_list(nRefine,refineList,cells)
integer, intent(in) :: nRefine
integer, intent(inout) :: refineList(:)
type(tCell), intent(in) :: cells(:)
integer :: temp
integer :: i, j
integer :: nc
logical :: swapped

do j = nRefine-1, 1, -1
  swapped = .false.
  do i = 1, j
    if (cells(refineList(i)) % refineLevel(1) > cells(refineList(i+1)) % refineLevel(1)) then
      temp = refineList(i)
      refineList(i) = refineList(i+1)
      refineList(i+1) = temp
      swapped = .true.
    end if
  end do
  if (.not. swapped) exit
end do

end subroutine sort_refine_list

end module refinement
