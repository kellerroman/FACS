program test
use array_holes

implicit none

type(holes) :: ar

call ar%init(100)
write(*,*) ar%nEntry

call ar%removeEntry(100)
write(*,*) ar%nEntry
write(*,*) ar%newEntry()
call ar%removeEntry(60)
write(*,*) ar%nEntry
write(*,*) ar%newEntry()

! leads to an error, but 100% code coverage
call ar%removeEntry(101)
end program test
