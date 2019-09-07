module array_holes
implicit none
integer, parameter :: N_HOLES_MAX = 1000

type :: holes
    private 
    !integer :: MaxArray
    integer :: lastEntry                          ! Position where to add when no more holes 
    integer, allocatable :: holes(:)
    integer, public :: nEntry                             ! Number of Entries (not necessary continues
    integer, public :: nHoles                             ! Number of Holes in the Array
    contains
        procedure :: init       => holes_Init
        procedure :: newEntry => new_Entry
        procedure :: removeEntry => remove_Entry
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

subroutine remove_entry(this,id)
    implicit none
    Class(holes), intent(inout) :: this
    integer, intent(in) :: id

    this % nEntry = this % nEntry - 1
    if (id == this % lastEntry) then
        this % lastEntry = this % lastEntry - 1
    else if (id > this % lastEntry) then
        write(*,*) "Array does not exist"
        stop 1
    else
        this % nHoles = this % nHoles + 1
        this % holes(this % nHoles) = id
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
