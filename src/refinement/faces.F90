module face
use const
use types
implicit none

contains
   subroutine add_face(cell1,cell2,dir1                                 &
                      ,cells, pnts, faces                               &
                      ,holesFaces,nHolesFace                            &
                      ,nFace,debug_in )
implicit none

integer        , intent(in)               :: cell1
integer        , intent(in)               :: cell2
integer        , intent(in)               :: dir1


type(tCell)   , intent (inout)            :: cells(:)
real(kind = 8), intent (inout)            :: pnts(:,:)
type(tFace      )   , intent (inout)      :: faces(:)

integer       , intent (inout)            :: holesFaces(:)
integer       , intent (inout)            :: nHolesFace             

integer       , intent (inout)            :: nFace
logical       , intent (in), optional     :: debug_in

logical                                   :: debug

integer :: myFace
integer :: p1, p2
integer :: mf,nf

logical :: pos_dir
real(kind = 8) :: dx,dy

if (present(debug_in)) then
   debug = debug_in
else
   debug = .false.
end if
if (dir1 == LEFT .or. dir1 == SOUTH) then
   pos_dir = .false.
else
   pos_dir = .true.
end if
! FACE
!create new face
if (nHolesFace > 0 ) then
   myFace = holesFaces(nHolesFace)
   nHolesFace = nHolesFace - 1
else
   nFace = nFace + 1
   myFace = nFace
end if
! set normal vector of face

! set length of face
if (dir1 == LEFT) then
   p1 = 4
   p2 = 1
else if (dir1 == RIGHT) then
   p1 = 3
   p2 = 2
else
   stop 1
end if

dx = pnts(1,cells(cell1) % pnts(p1)) - pnts(1,cells(cell1) % pnts(p2))
dy = pnts(2,cells(cell1) % pnts(p1)) - pnts(2,cells(cell1) % pnts(p2))
faces(myFace) % area = sqrt(dx*dx+dy*dy)
faces(myFace) % n(1) =  dy
faces(myFace) % n(2) = -dx

! add right cell to stencil
if (pos_dir) then
   faces(myFace) % stencil(1,1) = cell1
   faces(myFace) % stencil(1,2) = cell2
else
   faces(myFace) % stencil(1,1) = cell2
   faces(myFace) % stencil(1,2) = cell1
end if
! number of faces in new cell 
mf = cells(cell1) % nFace + 1
cells(cell1) % nFace = mf
! add face to cell face array
cells(cell1) % faces(mf) = myFace
! add face to cell directinal face array and mark side as one phase only
cells(cell1) % faces_dir_ref(dir1)   = .false.
cells(cell1) % faces_dir    (dir1,1) = myFace
! set face sign
! number of faces for LEFT neighbor cell
nf = cells(cell2) % nFace + 1
cells(cell2) % nFace = nf
! add face to LEFT neighbor's face array
cells(cell2) % faces(nf) = myFace
! add face to LEFT neighbor's directional face array
! RIGHT side of cell has two faces now
! since the new cell is the LOWER one, the new face is number 2 in the RIGHT direction
if (cells(cell2) % refineLevel(1) < cells(cell1) % refineLevel(1) ) then
   cells(cell2) % faces_dir_ref(RIGHT)   = .true.
   cells(cell2) % faces_dir    (RIGHT,2) = myFace
else if (cells(cell2) % refineLevel(1) == cells(cell1) % refineLevel(1) ) then 
   cells(cell2) % faces_dir_ref(RIGHT)   = .false.
   cells(cell2) % faces_dir    (RIGHT,1) = myFace
end if

! set face sign
if (pos_dir) then
   cells(cell1) % f_sign(mf) =  1.0d0
   cells(cell2) % f_sign(nf) = -1.0d0
else
   cells(cell1) % f_sign(mf) =  -1.0d0
   cells(cell2) % f_sign(nf) = 1.0d0
end if
end subroutine add_face
end module
