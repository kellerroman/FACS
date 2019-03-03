program test
use array_holes

implicit none

type(holes) :: ar

call ar%init(100)

call ar%removeEntry(100)
write(*,*) ar%newEntry()
call ar%removeEntry(60)
write(*,*) ar%newEntry()
end program test
