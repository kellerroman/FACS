program sol2csv
! **************************************************************************************************
! ***                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRle                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   19.01.2018
! Last changes: 19.01.2018
! Version:      V0.0.1
! --------------------------------------------------------------------------------------------------
! Description:
!   Generates a CSV output File from a Solution file (must be 1D)
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
character(len=*),parameter          :: FILENAME_IN = "sol.dat"
character(len=*),parameter          :: FILENAME_PV = "sol.csv"

type(tCell), allocatable            :: cells(:)
type(tParentCell), allocatable            :: parentCells(:)
real(kind = 8), allocatable         :: pnts(:,:)
integer                             :: nCells
integer                             :: nParentCells
integer                             :: nPnts

integer                             :: fo
integer                             :: i
real(kind = 8)                      :: center


call read_sol(cells,parentCells,pnts,nCells,nParentCells,nPnts,FILENAME_IN)

open(newunit = fo,file=trim(filename_pv))
write(*,'("Grid contains ",I0," Cells and ",I0," Points")') nCells, nPnts

do i = 1, nCells
   center = 0.25d0 * ( pnts(1,cells(i) % pnts(1)) + pnts(1,cells(i) % pnts(2)) &
                     + pnts(1,cells(i) % pnts(3)) + pnts(1,cells(i) % pnts(4)))
   write(fo,*) center,cells(i) % refineLevel(1), cells(i) % Q, cells(i) % grad
end do

close(fo)
end program sol2csv
