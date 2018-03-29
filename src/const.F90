module const
implicit none

integer, parameter                        :: MAX_REF_LEVEL = 6
character(len=1), parameter :: FACES(6) = ["L","R","S","N","F","B"]

integer         , parameter :: LEFT     = 1
integer         , parameter :: RIGHT    = 2
integer         , parameter :: SOUTH    = 3
integer         , parameter :: NORTH    = 4
integer         , parameter :: FRONT    = 5
integer         , parameter :: BACK     = 6

integer         , parameter :: NO_CELL  = 0

end module const
