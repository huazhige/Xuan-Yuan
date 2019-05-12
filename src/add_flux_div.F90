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

module add_flux_divergence

#include "definition.h"
  !! Import the runtime parameters for the PDE.
  use setup_module,   only: nx1, nghost
  use mesh_init,      only: dx1
contains

  subroutine add_flux_div(u,flux,dt)

    !! Primitive variables
    real, allocatable, dimension(:), intent(IN) :: flux
    !! Fluxes at cell interfaces
    real, allocatable, dimension(:), intent(INOUT) :: u
    !! New time step
    real, intent(IN) :: dt

    !! Calculate flux divergence and update conservative variables.
    do i = nghost + 1, nx1 + nghost, 1
      u(i) = u(i) - dt * ( flux( i - nghost + 1) - flux( i - nghost ) ) / dx1
    end do

  end subroutine add_flux_div

end module add_flux_divergence
