program sol2para
! **************************************************************************************************
! ***                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRle                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   19.01.2018
! Last changes: 19.01.2018
! Version:      V0.0.1
! --------------------------------------------------------------------------------------------------
! Description:
!   Generates a Paraview output File from a Solution file
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
! **************************************************************************************************
use file_io
implicit none
character(len=100)                  :: filename_in = "sol.dat"
character(len=100)                  :: filename_pv = "paraview.vtk"

type(tCell), allocatable            :: cells(:)
type(tParentCell), allocatable      :: parentCells(:)
real(kind = 8), allocatable         :: pnts(:,:)
integer                             :: nCells
integer                             :: nParentCells
integer                             :: nPnts

integer                             :: fo
integer                             :: i

character(len = 10) :: arg

i = 1
DO
   CALL get_command_argument(i, arg)
   IF (LEN_TRIM(arg) == 0) EXIT
   if( i == 1) then
      CALL get_command_argument(i, filename_in)
   else if (i == 2) then
      CALL get_command_argument(i, filename_pv)
   end if   
   i = i+1
END DO

call read_sol(cells,parentCells,pnts,ncells,nParentCells,npnts,filename_in)

open(newunit = fo,file=trim(filename_pv))
!write(*,'("Grid contains ",I0," Cells and ",I0," Points")') nCells, nPnts

write(fo,"(a)") '# vtk DataFile Version 2.0'
write(fo,"(a)") 'GRID-ADAPTION unstructured_grid'
write(fo,"(a)") 'ASCII'
write(fo,"(a)") ""
write(fo,"(a)") "DATASET UNSTRUCTURED_GRID"
write(fo,"(a,1X,i0,1X,a)") "POINTS",nPnts,"float"

do i = 1, nPnts
      write(fo,'(3(F20.13,1X))') pnts(:,i), 0.0d0
end do
write(fo,*)

write(fo,"(a,1X,i0,1X,i0)") "CELLS", nCells, nCells*5
do i = 1, nCells
   write(fo,'(5(i0,1X))') 4,cells(i) % pnts-1
end do
write(fo,*)

write(fo,"(a,1X,i0)") "CELL_TYPES", nCells
do i = 1, nCells
   write(fo,'(i0)') 9
end do
write(fo,*)

!write(fo,"(A,1X,I0)") "POINT_DATA",git % nPoint
!write(fo,"(A)") 'SCALARS Movement_Restriction double'
!write(fo,"(A)") 'LOOKUP_TABLE Default'
!do i = 1, git % nPoint
!   write(fo,*) temp
!end do

!write(fo,"(A)") 'VECTORS MOVEMENT_VECTOR double'
!do i = 1, git % nPoint
!   !write(fo,'(3(D20.13,1X))') git % point_move_rest_vector(:,i)
!   write(fo,*) git % point_move_rest_vector(:,i)
!end do

write(fo,"(A,1X,I0)") "CELL_DATA",nCells
write(fo,"(A)") 'SCALARS Refinement_Level_I int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % refineLevel(1)
end do

write(fo,"(A)") 'SCALARS Refinement_Level_J int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % refineLevel(2)
end do

write(fo,"(A)") 'SCALARS Cell_Number int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) i
end do

write(fo,"(A)") 'SCALARS Neighbor_1_LEFT int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % neigh(1)
end do

write(fo,"(A)") 'SCALARS Neighbor_2_RIGHT int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % neigh(2)
end do

write(fo,"(A)") 'SCALARS Neighbor_3_SOUTH int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % neigh(3)
end do

write(fo,"(A)") 'SCALARS Neighbor_4_NORTH int'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % neigh(4)
end do

write(fo,"(A)") 'SCALARS Var double'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % var
end do

write(fo,"(A)") 'SCALARS GradX double'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % grad(1)
end do

write(fo,"(A)") 'SCALARS GradY double'
write(fo,"(A)") 'LOOKUP_TABLE Default'
do i = 1, nCells
   write(fo,*) cells(i) % grad(2)
end do

close(fo)
end program sol2para
