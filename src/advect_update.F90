!!------------------------------------------------------------------
!! A Fortran example code for solving advection-diffusion function.
!!
!! This code is written by Huazhi Ge for AMS 209 final project.
!!
!! * Three parts of solving the PDEs.
!! 1) 1D advection solution
!! 2) 1D diffusion solution
!! 3) 1D advection & diffusion
!!
!! This program is used to update the cell from upwind method for
!! advection.
!!
!!------------------------------------------------------------------

module advect_update

#include "definition.h"

  use advect_init,    only: setAdvect
  use setup_module,   only: gridN, maxNT, methodType, nghost
  !! subroutines for discretization methods.
  use setup_grids,    only: deltaX, dT

contains

  subroutine advectUpdate(cells, newCell, nIter)

    implicit none
    !! Cells with boundary condition & initial condition.
    real, allocatable, dimension(:,:), intent(INOUT) :: cells
    !! Updated new cell.
    real, allocatable, dimension(:), intent(INOUT) :: newCell
    real, allocatable, dimension(:) :: oldCell
    integer :: nIter

    allocate( oldCell(gridN+2*nghost) )

    end if
    !! Free the memory.
    deallocate(oldCell)

  end subroutine advectUpdate

end module advect_update
