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
!! This program is used to set the upwind method. This method will
!!  calculate the value of a next time step (Ui,n+1) by using (Ui,n)
!! and (Ui-1,n).
!!
!!------------------------------------------------------------------

module upwind

#include "definition.h"

  use setup_module, only: gridN, velocity, nghost, CFL_Ca, runName
  use setup_grids,  only: deltaX, dT
  use boundary,     only: setBound

contains

  subroutine upwind_advect(oldCell, newCell)

    implicit none
    !! Past cell is the Un array, new cell is the Un+1 array.
    real, allocatable, dimension(:), intent(IN) :: oldCell
    real, allocatable, dimension(:), intent(INOUT) :: newCell
    !! Loop index
    integer :: nIter
    !!stabilizer = 1/sqrt( 1 + CFL_Ca**2 * sin(1./real(gridN)) )
    !! Upwind method: Ui,n+1 = Ui,n - velocity * (Ui,n - Ui,n-1) * dT/dX
    do nIter = nghost + 1, gridN + nghost, 1
      newCell(nIter) = oldCell(nIter) - CFL_Ca * &
                      (oldCell(nIter) - oldCell(nIter - 1))
    end do
    !! Boundary conditions.
    call setBound(oldCell, newCell)

  end subroutine upwind_advect

end module upwind
