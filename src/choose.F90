module choose
   use const
   use types
implicit none
logical, parameter :: DO_MULTI_DIM_REFINEMENT = .false.
real(kind = 8), parameter :: GRAD_MIN_COARSE = 1.00D-10
real(kind = 8) :: grad_ref = 110.0d-2
real(kind = 8) :: grad_coa =  10.0d-2
real(kind = 8) :: grad_avg(Q_DIM)
real(kind = 8) :: grad_max(Q_DIM)
real(kind = 8) :: grad_min(Q_DIM)

contains
subroutine choose_cells(cells,nCells,refineType,refineList,nRefine)
implicit none
type(tCell)   , intent (inout)            :: cells(:)
integer       , intent (inout)            :: nCells
integer       , intent (out)              :: refineType(:)
integer       , intent (out)              :: refineList(:)
integer       , intent (out)              :: nRefine
integer                                   :: i
real(kind = 8)                            :: min_grad, max_grad

nRefine = 0
refineType = 0 
max_grad =                     grad_ref * grad_avg(1)
min_grad = max(GRAD_MIN_COARSE,grad_ref * grad_avg(1))
do i = 1, nCells
   if (DO_MULTI_DIM_REFINEMENT) then
      if ( abs(Cells(i) % grad(1,1)) > max_grad) then
         if (cells(i) % refineLevel(1) < MAX_REF_LEVEL) then
            nRefine = nRefine + 1
            refineList(nRefine) = i
            refineType(i) = 1
         end if
      end if
      if ( abs(Cells(i) % grad(2,1)) > max_grad) then
         if (cells(i) % refineLevel(2) < MAX_REF_LEVEL) then
            if (refineType(i) == 0) then
               nRefine = nRefine + 1
               refineList(nRefine) = i
               refineType(i) = 2
            else
               refineType(i) = 3
            end if
         end if
      end if
   else
      if ( abs(cells(i) % grad(1,1)) > max_grad .or. &
           abs(cells(i) % grad(2,1)) > max_grad) then
         if (cells(i) % refineLevel(1) < MAX_REF_LEVEL) then
            nRefine = nRefine + 1
            refineList(nRefine) = i
            refineType(i) = 3
         end if
      else if ( abs(cells(i) % grad(1,1)) < min_grad .and. &
                abs(cells(i) % grad(2,1)) < min_grad) then
         refineType(i) = -3
      end if
   end if

end do
end subroutine choose_cells

!!!!! This routine ensures  that there is no bigger difference in refinment level than 1
subroutine smooth_refinement(cells,refineType,refineList,nRefine)
implicit none
type(tCell)   , intent (inout)            :: cells(:)
integer       , intent (inout)            :: refineType(:)
integer       , intent (inout)            :: refineList(:)
integer       , intent (inout)            :: nRefine
integer                                   :: i,d,n,ni,ci

ci = 1
do while (ci <= nRefine)
   i = refineList(ci)
   do n = 1, 4
      if (cells(i) % neigh(n) == 0) cycle
      ni = cells(i) % neigh(n)
      do d = 1, 2
         if (cells(i) % refineLevel(d) > cells(ni) % refineLevel(d)) then
            if (refineType(ni) <= 0) then
               refineType(ni) = d
               nRefine = nRefine + 1
               refineList(nRefine) = ni
            else if (refineType(ni) == 3 - d ) then
               refineType(ni) = 3
            end if
         end if
      end do
   end do
   ci = ci + 1
end do
end subroutine smooth_refinement

subroutine choose_coarse(cells,parentCells,refineType,canCoarseList,nCanCoarse,doCoarseList,nDoCoarse)

use types
implicit none
type(tCell)   , intent (in)               :: cells(:)
type(tParentCell)   , intent (inout)      :: parentCells(:)

integer       , intent (in)               :: refineType(:)
integer       , intent (in)               :: canCoarseList(:)
integer       , intent (in)               :: nCanCoarse
integer       , intent (out)              ::  doCoarseList(:)
integer       , intent (out)              :: nDoCoarse

integer     :: i, pc,cc,c,n,nc
logical     :: doCoarse
nDoCoarse = 0
do i = 1, nCanCoarse
   ! go thru every possible cell which could be coarsed
   pc = canCoarseList(i)
   doCoarse = .true.
   do c = 1, 4
      ! look for all of its childern if their can be coarsed
      cc = parentCells(pc) % child(c)
      cc = parentCells(cc) % ref
      if (refineType(cc) /= -3) then
         doCoarse = .false.
         exit
      end if
      ! check if coarsing is not possible due to:
      do n = 1, 4
         nc = cells(cc) % neigh(n)
         if (nc == NO_CELL) cycle
         !  1 Neighbor cell is finer
         if (cells(nc) % refineLevel(1) > cells(cc) % refineLevel(1) .or. &
             cells(nc) % refineLevel(2) > cells(cc) % refineLevel(2)) then
            doCoarse = .false.
            exit
         end if
         !  2 Neighbor cell on same refinement level is getting refined 
         if ((cells(nc) % refineLevel(1) == cells(cc) % refineLevel(1) .and. &
              ( refineType(nc) == 1 .or. refinetype(nc) == 3)) .or. &
             (cells(nc) % refineLevel(2) == cells(cc) % refineLevel(2) .and. &
              ( refineType(nc) == 2 .or. refineType(nc) == 3))   ) then
            doCoarse = .false.
            exit
         end if
      end do
      if (.not. doCoarse) exit
   end do
   if (doCoarse) then
      cc = parentCells(pc) % child(1)
      c = parentCells(cc) % ref
      nDoCoarse = nDoCoarse + 1
      doCoarseList(nDoCoarse) = c
      !write(*,*) "Combined:",c,"ParentCell:",cc,"ParentsParent:",pc, nDoCoarse, parentCells(pc) % pos_CanCoarse,i
      !write(*,*) cells(c) % neigh
   end if
end do
end subroutine

subroutine calc_gradient(cells,nCells)
implicit none
type(tCell)   , intent (inout)            :: cells(:)
integer       , intent (in)            :: nCells


integer                                   :: i
integer                                   :: nl,nh    ! Heighbor LOW and HIGH
real(kind = 8) :: grad(Q_DIM)

grad_avg = 0.0d0
grad_min = 1.0D20
grad_max = 1.0D-20
do i = 1, nCells
   !I
   nl = cells(i) % neigh(1)
   nh = cells(i) % neigh(2)
   if (nl /= 0 .and. nh /= 0) then
      grad = (cells(nh) % Q   - cells(nl) % Q  ) * (cells(nh) % center(1) - cells(nl) % center(1))
   else if (nl == 0 .and. nh /= 0) then
      grad = (cells(nh) % Q   - cells(i) % Q  ) * (cells(nh) % center(1) - cells(i) % center(1))
   else if (nl /= 0 .and. nh == 0) then
      grad = (cells(i) % Q   - cells(nl) % Q  ) * (cells(i) % center(1) - cells(nl) % center(1))
   else
      write(*,*) "I-Derivative not yet supported"
      write(*,*) "cell:",i,"Neigh-Left:",nl,"Neigh-Right:",nh
      write(*,*) "neighbors:", cells(i) % neigh
      stop 1
   end if
   cells(i) % grad(1,:) = grad
   grad = abs(grad)
   grad_min = min(grad_min,grad)
   grad_max = max(grad_max,grad)
   grad_avg = grad_avg + grad
   
   !J
   nl = cells(i) % neigh(3)
   nh = cells(i) % neigh(4)
   if (nl /= 0 .and. nh /= 0) then
      grad = (cells(nh) % Q - cells(nl) % Q) * (cells(nh) % center(2) - cells(nl) % center(2))
   else if (nl == 0 .and. nh /= 0) then
      grad = (cells(nh) % Q - cells(i) % Q) * (cells(nh) % center(2) - cells(i) % center(2))
   else if (nl /= 0 .and. nh == 0) then
      grad = (cells(i) % Q - cells(nl) % Q) * (cells(i) % center(2) - cells(nl) % center(2))
   else
      !write(*,*) "J-Derivative not yet supported",i,nl,nh
      grad = 0
!      stop 1
   end if
   cells(i) % grad(2,:) = grad
   grad = abs(grad)
   grad_min = min(grad_min,grad)
   grad_max = max(grad_max,grad)
   grad_avg = grad_avg + grad
end do
grad_avg = grad_avg / dble(2*nCells)
write(*,*) grad_avg,grad_min,grad_max

end subroutine calc_gradient

end module choose
