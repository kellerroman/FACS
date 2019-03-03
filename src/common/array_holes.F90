module array_holes
implicit none
 integer, parameter :: HOLES_SIZE = 1000

type :: holes
    private 
    integer :: MaxArray
    integer :: nArray
    integer :: nHoles
    integer, allocatable :: holes(:)
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

    if (this % nHoles == 0) then
       this % nArray = this % nArray + 1
       nId = this % nArray
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

    if (id == this % nArray) then
        this % nArray = this % nArray - 1
    else if (id > this % nArray) then
        write(*,*) "Array does not exist"
        stop 1
    else
        this % nHoles = this % nHoles + 1
        this % holes(this % nHoles) = id
    end if

end subroutine remove_entry



subroutine holes_Init(this, nArray)
    implicit none
    Class(holes), intent(inout) :: this
    integer, intent(in) :: nArray
    this % nArray = nArray
    allocate( this % holes( HOLES_SIZE) ) 
    this % nHoles = 0
end subroutine holes_Init

end module array_holes
