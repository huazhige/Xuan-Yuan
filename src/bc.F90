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
!! This program is used to set the boundary condition for each method.
!!
!!------------------------------------------------------------------

module boundary

  use setup_module, only: nx1, nghost, problem, x1_inner, x1_outer

contains

  subroutine set_boundary(u)
    implicit none
    integer :: ighost
    !! Primitive variable and conservative variables.
    real, allocatable, dimension(:), intent(INOUT) :: u
    !! Set the boundary condition

    if ( (problem .eq. 'advection') .or. (problem .eq. 'advect_diff') ) then
      !! periodic boudnary condition for inner boundary of x1
      if (x1_inner .eq. 'periodic') then
        do ighost = 1, nghost, 1
          u(ighost) = u(nx1 + ighost)
        end do
      !! periodic boudnary condition for outer boundary of x1
      else if (x1_outer .eq. 'periodic') then
        do ighost = 1, nghost, 1
          u(nx1 + ighost + nghost) = u(nghost + ighost)
        end do
      !! outflow boudnary condition for inner boundary of x1
      else if (x1_inner .eq. 'outflow') then
        do ighost = 1, nghost, 1
          u(ighost) = u(nghost + 1)
        end do
      !! outflow boudnary condition for outer boundary of x1
      else if (x1_outer .eq. 'outflow') then
        do ighost = 1, nghost, 1
          u(nx1 + ighost + nghost) = u(nx1 + nghost)
        end do
      end if
    end if

  end subroutine set_boundary

end module boundary
