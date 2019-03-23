module init
use const
use types
use array_holes


implicit none
contains

subroutine init_sol(cells,parentCells,pnts,faces,nCells,nParentCells,nFace)
implicit none
type(tCell)      , intent(inout), allocatable      :: cells(:)
type(tParentCell), intent(inout), allocatable      :: parentCells(:)
real(kind = 8)   , intent(inout), allocatable      :: pnts(:,:)
type(tFace)      , intent(inout), allocatable      :: faces(:)
type(holes)      , intent(in)                      :: nCells
type(holes)      , intent(out)                     :: nParentCells, nFace

integer                                            :: i,f,n,mf,nf
integer                                            :: ngf, nc, npc


nc = nCells % nEntry

call nParentCells % init(nc)

call nFace % Init(0)

npc = nParentCells % nEntry

do i = 1, nPC
   parentCells(i) % pnts          = cells(i) % pnts
   parentCells(i) % neigh         = cells(i) % neigh
   parentCells(i) % refineLevel   = cells(i) % refineLevel
   parentCells(i) % parent        = NO_CELL 
   parentCells(i) % child         = NO_CELL
   cells(i) % ref                 = i
   parentCells(i) % ref           = i
   parentCells(i) % pos_CanCoarse = NO_CELL
end do

do i = 1, nc
   cells(i) % center(:) = 0.25d0 * ( pnts(:,cells(i) % pnts(1)) & 
                                   + pnts(:,cells(i) % pnts(2)) & 
                                   + pnts(:,cells(i) % pnts(3)) & 
                                   + pnts(:,cells(i) % pnts(4)) ) 
   cells(i) % QC(1)       = cells(i) % Q(1)
   cells(i) % QC(2:Q_DIM) = cells(i) % Q(1) * cells(i) % Q(2:Q_DIM)
   cells(i) % vol = (pnts(1,cells(i) % pnts(2))-pnts(1,cells(i) % pnts(1))) & 
                  * (pnts(2,cells(i) % pnts(3))-pnts(2,cells(i) % pnts(2)))
   cells(i) % vol = 1.0d0 / cells(i) % vol
   cells(i) % nFace = 0
end do
do i = 1, nc
   do f = LEFT,RIGHT
      n = cells(i) % neigh(f)
      if (n > i .and. n /= NO_CELL) then
         ngf = nFace % newEntry()
         faces(ngf) % n = [1.0d0,0.0d0]

         cells(i) % nFace = cells(i) % nFace + 1
         mf = cells(i) % nFace   ! MyFace
         cells(i) % faces(mf)  = ngf
         cells(i) % faces_dir_ref(f)   = .false.
         cells(i) % faces_dir    (f,1) = ngf

         cells(n) % nFace = cells(n) % nFace + 1
         nf = cells(n) % nFace   ! NeighborFace
         cells(n) % faces(nf)  = ngf
         cells(n) % faces_dir_ref(3-f)   = .false.
         cells(n) % faces_dir    (3-f,1) = ngf

         if (f == LEFT) then
            cells(i) % f_sign(mf) =  1.0d0
            cells(n) % f_sign(nf) = -1.0d0
            faces(ngf) % area = pnts(2,cells(i) % pnts(4)) - pnts(2,cells(i) % pnts(1))
            faces(ngf) % stencil(1,1) = n
            faces(ngf) % stencil(1,2) = i
         else
            cells(i) % f_sign(mf) = -1.0d0
            cells(n) % f_sign(nf) =  1.0d0
            faces(ngf) % area = pnts(2,cells(i) % pnts(3)) - pnts(2,cells(i) % pnts(2))
            faces(ngf) % stencil(1,1) = i
            faces(ngf) % stencil(1,2) = n
         end if
      else if (n == NO_CELL) then
         ngf = nFace % newEntry()
         faces(ngf) % n = [1.0d0,0.0d0]

         cells(i) % nFace = cells(i) % nFace + 1
         mf = cells(i) % nFace   ! MyFace
         cells(i) % faces(mf)  = ngf
         cells(i) % faces_dir_ref(f)   = .false.
         cells(i) % faces_dir    (f,1) = ngf

         faces(ngf) % stencil(1,1) = i
         faces(ngf) % stencil(1,2) = i
         if (f == LEFT) then
            cells(i) % f_sign(mf) =  1.0d0
            faces(ngf) % area = pnts(2,cells(i) % pnts(4)) - pnts(2,cells(i) % pnts(1))
         else
            cells(i) % f_sign(mf) = -1.0d0
            faces(ngf) % area = pnts(2,cells(i) % pnts(3)) - pnts(2,cells(i) % pnts(2))
         end if
      end if
   end do
   do f = SOUTH,NORTH
      n = cells(i) % neigh(f)
      if (n > i .and. n /= NO_CELL) then
         ngf = nFace % newEntry()
         faces(ngf) % n = [0.0d0,1.0d0]

         cells(i) % nFace = cells(i) % nFace + 1
         mf = cells(i) % nFace   ! MyFace
         cells(i) % faces(mf)  = ngf
         cells(i) % faces_dir_ref(f)   = .false.
         cells(i) % faces_dir    (f,1) = ngf

         cells(n) % nFace = cells(n) % nFace + 1
         nf = cells(n) % nFace   ! NeighborFace
         cells(n) % faces(nf)  = ngf
         cells(n) % faces_dir_ref(7-f)   = .false.
         cells(n) % faces_dir    (7-f,1) = ngf
         if (f == SOUTH) then
            cells(i) % f_sign(mf) =  1.0d0
            cells(n) % f_sign(nf) = -1.0d0
            faces(ngf) % area = pnts(1,cells(i) % pnts(2)) - pnts(1,cells(i) % pnts(1))
            faces(ngf) % stencil(1,1) = n
            faces(ngf) % stencil(1,2) = i
         else
            cells(i) % f_sign(mf) = -1.0d0
            cells(n) % f_sign(nf) =  1.0d0
            faces(ngf) % area = pnts(1,cells(i) % pnts(3)) - pnts(1,cells(i) % pnts(4))
            faces(ngf) % stencil(1,1) = i
            faces(ngf) % stencil(1,2) = n
         end if
      else if (n == NO_CELL) then
         ngf = nFace % newEntry()
         faces(ngf) % n = [0.0d0,1.0d0]

         cells(i) % nFace = cells(i) % nFace + 1
         mf = cells(i) % nFace   ! MyFace
         cells(i) % faces(mf)  = ngf
         cells(i) % faces_dir_ref(f)   = .false.
         cells(i) % faces_dir    (f,1) = ngf

         faces(ngf) % stencil(1,1) = i
         faces(ngf) % stencil(1,2) = i
         if (f == SOUTH) then
            cells(i) % f_sign(mf) =  1.0d0
            faces(ngf) % area = pnts(1,cells(i) % pnts(2)) - pnts(1,cells(i) % pnts(1))
         else
            cells(i) % f_sign(mf) = -1.0d0
            faces(ngf) % area = pnts(1,cells(i) % pnts(3)) - pnts(1,cells(i) % pnts(4))
         end if
      end if
   end do
end do
end subroutine init_sol
end module init
