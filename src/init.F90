module init
use const
use types


implicit none
contains

subroutine init_sol(cells,parentCells,pnts,faces,nCells,nParentCells,nFace)
implicit none
type(tCell)      , intent(inout), allocatable      :: cells(:)
type(tParentCell), intent(inout), allocatable      :: parentCells(:)
real(kind = 8)   , intent(inout), allocatable      :: pnts(:,:)
type(tFace)      , intent(inout), allocatable      :: faces(:)
integer          , intent(out)                     :: nCells, nParentCells, nFace

integer                                            :: i,f,n,mf,nf


nParentCells = nCells

do i = 1, nParentCells
   parentCells(i) % pnts          = cells(i) % pnts
   parentCells(i) % neigh         = cells(i) % neigh
   parentCells(i) % refineLevel   = cells(i) % refineLevel
   parentCells(i) % parent        = NO_CELL 
   parentCells(i) % child         = NO_CELL
   cells(i) % ref                 = i
   parentCells(i) % ref           = i
   parentCells(i) % pos_CanCoarse = NO_CELL
end do

do i = 1, nCells
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
do i = 1, nCells
   do f = 1,2
      n = cells(i) % neigh(f)
      if (n > i .and. n/= NO_CELL) then
         nFace = nFace + 1
         faces(nFace) % n = [1.0d0,0.0d0]

         cells(i) % nFace = cells(i) % nFace + 1
         mf = cells(i) % nFace   ! MyFace
         cells(i) % faces(mf)  = nFace

         cells(n) % nFace = cells(n) % nFace + 1
         nf = cells(n) % nFace   ! NeighborFace
         cells(n) % faces(nf)  = nFace
         if (f == LEFT) then
            cells(i) % f_sign(mf) =  1.0d0
            cells(n) % f_sign(nf) = -1.0d0
            faces(nFace) % area = pnts(2,cells(i) % pnts(4)) - pnts(2,cells(i) % pnts(1))
            faces(nFace) % stencil(1,1) = n
            faces(nFace) % stencil(1,2) = i
         else
            cells(i) % f_sign(mf) = -1.0d0
            cells(n) % f_sign(nf) =  1.0d0
            faces(nFace) % area = pnts(2,cells(i) % pnts(3)) - pnts(2,cells(i) % pnts(2))
            faces(nFace) % stencil(1,1) = i
            faces(nFace) % stencil(1,2) = n
         end if
      end if
   end do
end do
end subroutine init_sol
end module init
