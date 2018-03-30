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
character(len=100)                  :: filename_in = "sol.dat"
character(len=100)                  :: filename_pv = "sol.csv"

type(tCell), allocatable            :: cells(:)
real(kind = 8), allocatable         :: pnts(:,:)
integer                             :: nCells
integer                             :: nPnts

integer                             :: fo
integer                             :: i
real(kind = 8)                      :: center

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

call read_sol(filename_in,cells,pnts,nCells,nPnts)

open(newunit = fo,file=trim(filename_pv))
write(fo,'(A15)',ADVANCE="NO") "X"
write(fo,'(",",A15)',ADVANCE="NO") "rho"
write(fo,'(",",A15)',ADVANCE="NO") "vel"
write(fo,'(",",A15)',ADVANCE="NO") "E"
write(fo,'(",",A15)',ADVANCE="NO") "p"
write(fo,*)

do i = 1, nCells
   center = 0.25d0 * ( pnts(1,cells(i) % pnts(1)) + pnts(1,cells(i) % pnts(2)) &
                     + pnts(1,cells(i) % pnts(3)) + pnts(1,cells(i) % pnts(4)))
   write(fo,'(E15.5)'    , ADVANCE="NO") center
   write(fo,'(3(",",E15.5))', ADVANCE="NO") cells(i) % Q
   write(fo,'(1(",",E15.5))', ADVANCE="NO") (cells(i) % Q(3) - 0.5*cells(i) % Q(2)**2)*cells(i) % Q(1) * 0.4d0
   write(fo,*)
end do

close(fo)
end program sol2csv
