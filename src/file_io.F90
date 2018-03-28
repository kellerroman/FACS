module file_io
use const
use types


implicit none
contains

subroutine read_sol(cells,parentCells,pnts,nCells,nParentCells,nPnts,filename)
implicit none
type(tCell)      , intent(inout), allocatable      :: cells(:)
type(tParentCell), intent(inout), allocatable      :: parentCells(:)
real(kind = 8)   , intent(inout), allocatable      :: pnts(:,:)
integer          , intent(out)                     :: nCells, nParentCells, nPnts
character(len=*) , intent(in)                      :: filename

integer                                            :: nMaxCells
integer                                            :: nMaxPnts
integer                                            :: i
integer                                            :: fi


open(newunit = fi,file=trim(filename))

read(fi,'(I9,1X,I9)') nCells, nPnts
write(*,'("Grid contains ",I0," Cells and ",I0," Points")') nCells, nPnts

if (.not. allocated(cells)) then
   allocate(cells(nCells))
else
   nMaxCells = ubound(cells,1)
   if (nCells > nMaxCells) then
      write(*,*) "Maximum number of Cells reached"
      stop 1
   end if
end if
if (.not. allocated(parentCells)) then
   allocate(parentCells(nCells))
else
   nMaxCells = ubound(parentCells,1)
   if (nCells > nMaxCells) then
      write(*,*) "Maximum number of ParentCells reached"
      stop 1
   end if
end if

if (.not. allocated(pnts)) then
   allocate(pnts(2,nPnts))
else
   nMaxPnts  = ubound(pnts,2)
   if (nPnts  > nMaxPnts) then
      write(*,*) "Maximum number of Points reached"
      stop 1
   end if
end if

do i = 1, nPnts
   read (fi,'(10X,3(1X,ES12.5))') pnts(:,i)
end do

do i = 1, nCells
   read (fi,'(10X,4(1X,I9))') cells(i) % pnts
end do

do i = 1, nCells
   read (fi,'(10X,4(1X,I9))') cells(i) % neigh
end do

do i = 1, nCells
   read (fi,'(10X,3(1X,I9))') cells(i) % refineLevel
end do

do i = 1, nCells
   read (fi,'(10X,1(1X,ES12.5))') cells(i) % var
end do

do i = 1, nCells
   read (fi,'(10X,2(1X,ES12.5))') cells(i) % grad
end do
close(fi)
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
end subroutine read_sol

subroutine write_sol(cells,pnts,nCells,nPnts,filename,iter)
implicit none
type(tCell)      , intent(in), allocatable         :: cells(:)
real(kind = 8)   , intent(in), allocatable         :: pnts(:,:)
integer          , intent(in)                      :: nCells, nPnts
character(len=*) , intent(in)                      :: filename
integer          , intent(in), optional            :: iter

integer                                            :: i
integer                                            :: fi
character(len=20)                                  :: filen
if (present(iter)) then
   write(filen,'("sol_",I5.5,".dat")') iter
   open(newunit = fi,file=trim(filen))
else
   open(newunit = fi,file=trim(filename))
end if

write(fi,'(I9,",",I9)') nCells, nPnts

do i = 1, nPnts
      write (fi,'(I9,":",3(1X,ES12.5))') i , pnts(:,i)
end do
do i = 1, nCells
   write (fi,'(I9,":",4(1X,I9))') i, cells(i) % pnts
end do
do i = 1, nCells
   write (fi,'(I9,":",4(1X,I9))') i,cells(i) % neigh
end do
do i = 1, nCells
   write (fi,'(I9,":",2(1X,I9))') i,cells(i) % refineLevel
end do
do i = 1, nCells
   write (fi,'(I9,":",1(1X,ES12.5))') i,cells(i) % var
end do
do i = 1, nCells
   write (fi,'(I9,":",2(1X,ES12.5))') i,cells(i) % grad
end do
close(fi)
end subroutine write_sol
end module file_io
