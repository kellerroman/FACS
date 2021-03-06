program gridgen
! **************************************************************************************************
! ***                          RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRle                   ***
! **************************************************************************************************
! Author:       Roman Keller(RK)
! Start date:   19.01.2018
! Last changes: 19.01.2018
! Version:      V0.0.1
! --------------------------------------------------------------------------------------------------
! Description:
!   Gridgen for Sod Shock Tube
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

! **************************************************************************************************
!  Point ordering              Neighbors
!
!   4 _____ 3                __4__ 
!    |     |                |     |
!    |     |              1 |     | 2
!    |_____|                |_____|
!   1       2                  3
! 
! **************************************************************************************************
use file_io
use mod_create_block
use array_holes
implicit none

character(len=*), parameter         :: FILENAME = "sol.dat"


integer                             :: ni = 4
integer                             :: nj = 4

type(tCell), allocatable            :: cells(:)
real(kind = 8), allocatable         :: pnts(:,:)
type(holes)                         :: nCells
type(holes)                         :: nPnts
integer                             :: nc,np

write(*,'(A)') "Gridgen for SOD SHOCK TUBE"

nc = (ni-1) * (nj-1)
np  =  ni    *  nj
allocate (cells(nc))
allocate (pnts(2,np))
call create_block(ni,nj,cells,pnts,nCells,nPnts)
write(*,'("Grid contains ",I0," Cells and ",I0," Points")') nCells % nEntry, nPnts % nEntry
call write_sol(cells,pnts,nCells,nPnts,FILENAME)
write(*,'(A)') "Gridgen done!"
end program gridgen
