
## RESTART CAPABILTY
- include parentCells into file
- HDF5

## REFINEMENT FACES

- fully refinement:

- 4 internal faces

- lower-left:
*     new face with cell NORTH (sibling)
*     new face with cell RIGHT (sibling)
*     new face                if left neighbor has lower level
*     change stencil of face  if left neighbor has same level

- lower-right:
*     new face                if lower neighbor has lower level
*     change stencil of face  if lower neighbor has same level
*     new face                if right neighbor has lower level
*     change stencil of face  if right neighbor has same level

- upper-right:

## Code Coverage gcov & lgoc with pfUnit
## include pfUnit into repo with subrepo
