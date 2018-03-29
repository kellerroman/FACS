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
implicit none

character(len=*), parameter         :: FILENAME = "sol.dat"


integer                             :: ni = 5
integer                             :: nj = 5

real(kind = 8)                      :: xmin = -1.0d0
real(kind = 8)                      :: xmax = 1.0D0
real(kind = 8)                      :: ymin = -1.0d0
real(kind = 8)                      :: ymax =  1.0d0

type(tCell), allocatable            :: cells(:)
real(kind = 8), allocatable         :: pnts(:,:)
integer                             :: nCells
integer                             :: nPnts


integer                             :: i,j
integer                             :: n1
integer                             :: n,i1
real(kind = 8)                      :: dx,dy

write(*,'(A)') "Gridgen for SOD SHOCK TUBE"

nCells = (ni-1) * (nj-1)
nPnts  =  ni    *  nj
 
write(*,'("Grid contains ",I0," Cells and ",I0," Points")') nCells, nPnts

allocate (cells(nCells))
allocate (pnts(2,nPnts))

dx = (xmax - xmin) / dble(ni-1)
dy = (ymax - ymin) / dble(nj-1)

n = 0
do j = 0, nj-1
   do i = 0, ni-1
      n = n + 1
      pnts(1,n) = dx * dble(i) + xmin
      pnts(2,n) = dy * dble(j) + ymin
   end do
end do
n = 0
do j = 1, nj-1
   do i = 1, ni-1
      n = n + 1
      i1 = i + (j-1) * ni
      cells(n) % pnts(1) = i1
      cells(n) % pnts(2) = i1+1
      cells(n) % pnts(3) = i1+ni+1
      cells(n) % pnts(4) = i1+ni
      cells(n) % refineLevel = 0
      cells(n) % Q       = 1.0d0
   end do
end do

n = 0
do j = 1, nj-1
   do i = 1, ni-1
      n = n + 1
      if (i > 1) then
         n1 = n-1
      else
         n1 = 0
      end if
      cells(n) % neigh(1) = n1

      if (i < ni-1) then
         n1 = n + 1
      else 
         n1 = 0
      end if
      cells(n) % neigh(2) = n1

      if (j > 1) then
         n1 = n -(ni-1)
      else
         n1 = 0
      end if
      cells(n) % neigh(3) = n1

      if (j < nj -1) then
         n1 = n + (ni-1)
      else 
         n1 = 0
      end if
      cells(n) % neigh(4) = n1
   end do
end do

call write_sol(cells,pnts,nCells,nPnts,FILENAME)
write(*,'(A)') "Gridgen done!"
end program gridgen
