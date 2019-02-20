module const
implicit none

integer, parameter                        :: GIT_DIM = 2
integer, parameter                        :: MAX_REF_LEVEL = 3
integer, parameter                        :: Q_DIM = 4
character(len=1), parameter               :: FACE_NAME_SHORT(6) = ["L","R","S","N","F","B"]

enum, bind(C)
   enumerator :: LEFT=1, RIGHT, SOUTH, NORTH, FRONT, BACK 
end enum

integer         , parameter               :: NO_CELL  = 0

integer        ,parameter                 :: PNTS_ON_FACE(2,4) = reshape([4,1,3,2,2,1,3,4],[2,4])

end module const
