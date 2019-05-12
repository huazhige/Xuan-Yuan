!!------------------------------------------------------------------
!! A Fortran example code for solving advection-diffusion function.
!!
!! This code is written by Huazhi Ge for AMS 209 final project.
!!
!! * Three parts of solving nonlinear Euler equations.
!! 1) 1D advection solution
!!
!! * This is a driver routine which calls subroutines to solve PDEs:
!!
!!   advecct_diff:
!!         | -
!!         | - Read inputs (setup_module.F90)
!!         | - Initialize grids (grid_init.F90)
!!         | - Initialize advection
!!         | -
!!         | - Initialize time step (cfl.F90)
!!         | - Flux reconstruction (reconstruction.F90)
!!         | - Add flux divergence (advect_update.F90)
!!         | - Setup boundary condition (bc.F90)
!!         | -
!!         | - output_write (from output_module)
!!
!!
!!------------------------------------------------------------------

program advection

  !! include a C-type header file:
  !! this is why the file extensions are .F90, instead of .f90
#include "definition.h"
  !! Begining of the real implementation of the driver
  !! define usages of module variables and subroutines

  use mesh_init,           only: setGrids, dx1, dx1_array
  use advect_init,         only: advection_init
  use new_time_step,       only: new_dt
  use calculate_fluxes,    only: calculate_flux
  use add_flux_divergence, only: add_flux_div
  use boundary,            only: set_boundary
  use write_data,          only: write_out

  use setup_module,        only: setup_init, nx1, limit_dt, limit_time, &
                           printScale, nghost

  implicit none

  !! Conservative variables
  real, allocatable, dimension(:) :: u
  !! Primitive variables
  real, allocatable, dimension(:) :: w
  !! Fluxes at cell boundary
  real, allocatable, dimension(:) :: flux
  !! Define iterator
  integer :: nIter
  !! Define time step and simulation time.
  real :: dt, time

  !! Setup the runtime parameters at here.
  call setup_init()
  !! Allocate space for the old and new cell.
  allocate(   u(nx1 + 2 * nghost))
  allocate(   w(nx1 + 2 * nghost))
  allocate(flux(nx1 + 1))
  !! Initialize the cells & the other runtime parameters.
  nIter = 1
  u = 0.
  w = 0.
  flux = 0.
  time = 0.


  !! Initialize the initial condition and boundary conditions for the ghost cell.
  call advection_init(u,w)
  call set_boundary(u)
  w = u
  !! Task List
  do while ( ( nIter .lt. limit_dt ) .and. ( time .lt. limit_time ) )
    !! Step I: Compute the new time step.
    call new_dt(u,dx1,dt)
    !! Calculate the simulation time.
    time = time + dt
    !! Step II: Calculate flux at cell interfces.
    call calculate_flux(w,flux,dt)
    !! Step III: Add flux divergence.
    call add_flux_div(u,flux,dt)
    !! Step IV: Setup ghost cells for boundary condition.
    call set_boundary(u)
    !! Step V: Conservative variables to Primitive variables.
    w = u
    !! Step VI: Output variable of each time step.
    !! Step VII: Update iterator.
    nIter = nIter + 1
  end do

  !! Output the data.
  call write_out(nIter, dx1_array, u)

end program advection
