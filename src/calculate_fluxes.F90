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

module calculate_fluxes

#include "definition.h"
  !! Import the runtime parameters for the PDE.
  use setup_module,   only: reconstruction, method
  use riemann_solver, only: first_godunov, second_non_TVD, second_TVD

contains

  subroutine calculate_flux(w,flux,dt)

    !! Primitive variables
    real, allocatable, dimension(:), intent(IN) :: w
    !! Fluxes at cell interfaces
    real, allocatable, dimension(:), intent(INOUT) :: flux
    !! New time step
    real, intent(IN) :: dt

    !! Flux Reconstruction
    if ( reconstruction .eq. 1 ) then
      call first_godunov(w, flux, dt)
    else if ( reconstruction .eq. 2 ) then
      if ( method .eq. "non_TVD" ) then
        call second_non_TVD(w, flux, dt)
      else if ( method .eq. "TVD" ) then
        call second_TVD(w, flux, dt)
      end if
    end if

  end subroutine calculate_flux

end module calculate_fluxes
