module array_holes
implicit none
integer, parameter :: N_HOLES_MAX = 1000

type :: holes
    private 
    !integer :: MaxArray
    integer, public :: lastEntry                          ! Position where to add when no more holes 
    integer, allocatable :: holes(:)
    integer, public :: nEntry                             ! Number of Entries (not necessary continues
    integer, public :: nHoles                             ! Number of Holes in the Array
    contains
        procedure :: init       => holes_Init
        procedure :: newEntry => new_Entry
        procedure :: removeEntry => remove_Entry
        procedure :: removeLast => remove_Last
end type

contains
function new_Entry(this) result(nId)
    implicit none
    Class(holes), intent(inout) :: this
    integer :: nId

    this % nEntry = this % nEntry + 1
    if (this % nHoles == 0) then
       this % lastEntry = this % lastEntry + 1
       nId = this % lastEntry
       return
   else
       nId = this % holes(this % nHoles)
       this % nHoles = this % nHoles- 1
   end if
end function new_Entry

function remove_Last(this) result (le)
    implicit none
    Class(holes), intent(inout) :: this
    integer :: le
    le = this % lastEntry
    call remove_Entry(this,this % lastEntry)
end function remove_Last

subroutine remove_entry(this,id)
    implicit none
    Class(holes), intent(inout) :: this
    integer, intent(in) :: id

    integer :: mypos, t
    this % nEntry = this % nEntry - 1
    if (id == this % lastEntry) then
        ! write(*,*) "Deleting last Entry: ",id
        mypos = 0
        t = id - 1
        this % lastEntry = this % lastEntry - 1
        ! check if end of array is all holes
        do 
            if (this % holes(mypos+1) == t .and. mypos < this % nHoles) then
                ! write(*,*) "Deleting Hole : ",t,"@",mypos+1
                t = t - 1
                mypos = mypos + 1
            else 
                this % lastEntry = this % lastEntry - mypos 
                this % holes(1:this % nholes-mypos) = this % holes(mypos+1:this % nholes)
                this % nHoles = this % nHoles - mypos
                exit
            end if
        end do
    else if (id > this % lastEntry) then
        write(*,*) "Array Value does not exist"
        stop 1
    else
        this % nHoles = this % nHoles + 1
        this % holes(this % nHoles) = id
        mypos = this % nHoles
         ! sort the holes by larg at last
        do 
            if (mypos < 2) exit
            t = this % holes(mypos - 1)
            if (id > t ) then
                this % holes(mypos) =  t
                this % holes(mypos - 1) = id
                mypos = mypos - 1
            else 
                exit
            end if
        end do
    end if

end subroutine remove_entry

subroutine holes_Init(this, nEntry)
    implicit none
    Class(holes), intent(inout) :: this
    integer, intent(in) :: nEntry
    this % lastEntry = nEntry
    this % nEntry    = nEntry
    allocate( this % holes( N_HOLES_MAX) ) 
    this % nHoles = 0
    !write(*,*) "Initial Size of Array:", nEntry
end subroutine holes_Init

end module array_holes
