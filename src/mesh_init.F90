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
!! This program is used to set the grid.
!!
!!------------------------------------------------------------------

module mesh_init

#include "definition.h"

  use setup_module, only: nx1, nghost, x1Beg, x1End
  real, save :: dx1
  real, allocatable, dimension(:), save :: dx1_array

contains

  subroutine setGrids()

    implicit none
    integer :: i
    !! allocate the space for Ghost cells
    allocate(dx1_array(nx1 + 2 * nghost))
    !! Calculate the dx, dy, dz from the boundary.
    dx1 = (x1End - x1Beg) / nx1
    !! set up the coordinate for x,y,z
    do i = 1, nx1 + 2 * nghost
      dx1_array(i) = x1Beg + (i - (0.5 + nghost)) * dx1
    end do

  end subroutine setGrids

end module mesh_init
