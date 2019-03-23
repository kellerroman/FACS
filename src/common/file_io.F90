module file_io
use const
use types
use array_holes


implicit none
contains

subroutine read_sol(filename,cells,pnts,nCells,nPnts)
implicit none
type(tCell)      , intent(inout), allocatable      :: cells(:)
real(kind = 8)   , intent(inout), allocatable      :: pnts(:,:)
type(holes)      , intent(out)                     :: nCells, nPnts
character(len=*) , intent(in)                      :: filename

integer                                            :: nMaxCells
integer                                            :: nMaxPnts
integer                                            :: i, read_cells, read_pnts
integer                                            :: fi



open(newunit = fi,file=trim(filename))

read(fi,'(I9,1X,I9,1X,I3)') read_cells, read_pnts, i
call nCells % init(read_cells)
call nPnts  % init(read_pnts)
write(*,'("Grid contains ",I0," Cells and ",I0," Points")') read_cells, read_pnts

if (i /= Q_DIM) then
   write(*,*) "Solution File contains different StateVector Size"
   stop 1
end if
if (.not. allocated(cells)) then
   allocate(cells(read_cells))
else
   nMaxCells = ubound(cells,1)
   if (read_cells > nMaxCells) then
      write(*,*) "Maximum number of Cells reached"
      stop 1
   end if
end if

if (.not. allocated(pnts)) then
   allocate(pnts(2,read_pnts))
else
   nMaxPnts  = ubound(pnts,2)
   if (read_pnts  > nMaxPnts) then
      write(*,*) "Maximum number of Points reached"
      stop 1
   end if
end if

do i = 1, read_pnts
   read (fi,'(10X,3(1X,ES12.5))') pnts(:,i)
end do

do i = 1, read_cells
   read (fi,'(10X,4(1X,I9))') cells(i) % pnts
end do

do i = 1, read_cells
   read (fi,'(10X,4(1X,I9))') cells(i) % neigh
end do

do i = 1, read_cells
   read (fi,'(10X,3(1X,I9))') cells(i) % refineLevel
end do

do i = 1, read_cells
   read (fi,'(10X,1(1X,ES12.5))') cells(i) % Q
end do

do i = 1, read_cells
   read (fi,'(10X,2(1X,ES12.5))') cells(i) % grad
end do
close(fi)
end subroutine read_sol

subroutine write_sol(cells,pnts,nCells,nPnts,filename,iter)
implicit none
type(tCell)      , intent(in), allocatable         :: cells(:)
real(kind = 8)   , intent(in), allocatable         :: pnts(:,:)
type(holes)      , intent(in)                      :: nCells, nPnts
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

write(fi,'(I9,1X,I9,1X,I3)') nCells%nEntry, nPnts%nEntry, Q_DIM

do i = 1, nPnts % nEntry
      write (fi,'(I9,":",3(1X,ES12.5))') i , pnts(:,i)
end do
do i = 1, nCells % nEntry
   write (fi,'(I9,":",4(1X,I9))') i, cells(i) % pnts
end do
do i = 1, nCells % nEntry
   write (fi,'(I9,":",4(1X,I9))') i,cells(i) % neigh
end do
do i = 1, nCells % nEntry
   write (fi,'(I9,":",2(1X,I9))') i,cells(i) % refineLevel
end do
do i = 1, nCells % nEntry
   write (fi,'(I9,":",1(1X,ES12.5))') i,cells(i) % Q
end do
do i = 1, nCells % nEntry
   write (fi,'(I9,":",2(1X,ES12.5))') i,cells(i) % grad
end do
close(fi)
end subroutine write_sol
end module file_io
