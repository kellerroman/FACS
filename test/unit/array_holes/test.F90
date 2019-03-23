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
end program test
