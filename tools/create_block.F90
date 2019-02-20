module mod_create_block
    implicit none 
    contains
    subroutine create_block(ni,nj,cells,pnts,nCells,nPnts)
    use types
    implicit none
    integer,intent(in)                  :: ni
    integer,intent(in)                  :: nj
    type(tCell), allocatable            :: cells(:)
    real(kind = 8), allocatable         :: pnts(:,:)
    integer,intent(out)                 :: nCells
    integer,intent(out)                 :: nPnts


    real(kind = 8)                      :: xmin = 0.0d0
    real(kind = 8)                      :: xmax = 1.0D0
    real(kind = 8)                      :: ymin =  0.0d+0
    real(kind = 8)                      :: ymax =  1.0d+0

    integer                             :: i,j
    integer                             :: n1
    integer                             :: n,i1
    real(kind = 8)                      :: dx,dy

    nCells = (ni-1) * (nj-1)
    nPnts  =  ni    *  nj
 
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
            cells(n) % Q = 1.0d0
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
    end subroutine
end module
