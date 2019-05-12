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
!! This program is used to use the CFL condition to choose a suitable
!! foot step.
!!
!!------------------------------------------------------------------

module new_time_step

#include "definition.h"

  use setup_module, only: CFL, nx1, nghost, x1Beg, x1End

contains
  !! Use the CDF value to calculate the time foot step.
  subroutine new_dt(u,dx1,dt)
    implicit none

    real, intent(IN) :: dx1
    real, allocatable, dimension(:), intent(IN) :: u
    real, intent(OUT) :: dt
    !! Set the temporary time step.
    real :: temp_dt
    !! Set the iterator.
    integer :: i
    !! Set the new time step.
    dt = CFL * dx1 / u(1)

    do i = 1, nx1, 1
     temp_dt = CFL * dx1 / u(i)
     if ( temp_dt .lt. dt ) then
       dt = temp_dt
     end if
    end do

  end subroutine new_dt

end module new_time_step
