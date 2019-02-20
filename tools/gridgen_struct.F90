module gridgen_struct
implicit none
integer         , parameter :: MAX_BLOCK             = 50 !Maximale Anzahl Bloecke (Dimension des temp block arrays)
contains
subroutine add_block(ni,nj,nk) 
implicit none
integer, intent(in) :: ni,nj,nk
end subroutine add_block
end module gridgen_struct
