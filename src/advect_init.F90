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
!! This program is used to set the initial condition for the diffusion
!! part.
!!
!! Problem Generator.
!!
!!------------------------------------------------------------------

module advect_init

#include "definition.h"
  !! Import the runtime parameters for the PDE.
  use setup_module, only: nx1, x1Beg, x1End, left_speed, right_speed, &
                          limit_dt, nghost, problem
  use mesh_init,  only: setGrids, dx1_array, dx1

contains
  subroutine advection_init(u,w)

    !! Conservative variables
    real, allocatable, dimension(:), intent(INOUT):: u
    !! Primitive variables
    real, allocatable, dimension(:), intent(INOUT) :: w
    !! Define constants.
    real :: pi
    !! Define iterators
    integer :: i
    !! Define Pi
    pi = acos(-1.)

    !! Initialize the variatble: cells
    u = 0.
    w = 0.
    flux = 0.
    !! Get the x grids & t grids.
    call setGrids()

    !! Problem Generator
    if ( problem .eq. 'sin_wave' ) then
      !! Set the initial condition for the propagation of sinusoidal wave.
      print *, dx1_array
      w = SIN(2 * pi * dx1_array)
    else if ( problem .eq. "shock" ) then
      !! set the initial condition for the propagation of shock.
      w = right_speed
      do i = 1, int(nx1/2.) + nghost, 1
        w(i) = left_speed
      end do
    end if


  end subroutine advection_init

end module advect_init
