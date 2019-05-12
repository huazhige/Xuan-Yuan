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

module riemann_solver

#include "definition.h"
  !! Import the runtime parameters for the PDE.
  use setup_module,   only: nx1, nghost, reconstruction, slope_type
  use mesh_init,      only: setGrids, dx1

contains

  !! Subroutine for getting the sign of the variable.
  subroutine get_sign(a, sign)
    real, intent(IN) :: a
    real, intent(OUT) :: sign

    !! Return +1 if a >= 0, return -1 if a < 0.
    if ( a .ge. 0 ) then
      sign = 1.
    else if ( a .lt. 0 ) then
      sign = -1.
    end if

  end subroutine get_sign

  !! First order Gudnov reconstruction (constant reconstruction).
  subroutine first_godunov(w,flux,dt)

    !! Primitive variables
    real, allocatable, dimension(:), intent(IN) :: w
    !! Fluxes at cell interface
    real, allocatable, dimension(:), intent(INOUT) :: flux
    !! Define the current time step.
    real, intent(IN) :: dt
    !! Define shock speed and characteristic speed at each cell interface
    !! from Burgers' equation
    real :: shock_s, character_s
    !! Define iterator
    integer :: i
    !! Initialize the shock speed.
    shock_s = 0.
    character_s = 0.

    !! Calculate fluxes with Godunov's method.
    do i = nghost, nghost + nx1, 1
      !! Caculate shock speed at interface.
      shock_s = 0.5 * ( w(i) + w(i+1) )

      !! U_{i}^{n} > U_{i+1}^{n}
      if ( w(i) .ge. w(i+1) ) then
        !! S > 0
        if ( shock_s .gt. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * w(i)**2.
        !! S < 0
        else if ( shock_s .le. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * w(i+1)**2.
        end if
      !! U_{i}^{n} < U_{i+1}^{n}

      else if ( w(i) .lt. w(i+1) ) then
        !! 0 < U_{i}^{n}
        if ( w(i) .ge. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * w(i)**2.
        !! U_{i+1}^{n} > 0
        else if ( w(i+1) .le. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * w(i+1)**2.
        !! U_{i}^{n} < 0 < U_{i+1}^{n}
        else if ( (w(i) .lt. character_s) .and. &
                (w(i+1) .gt. character_s) ) then
          flux( i - nghost + 1 ) = 0.
        end if
      end if
    end do

  end subroutine first_godunov

  !! Second order reconstruction with non-TVD slope
  !! (Piecewise Linear Method with Non-TVD slope, PLM).
  subroutine second_non_TVD(w,flux,dt)

    !! Primitive variables
    real, allocatable, dimension(:), intent(IN) :: w
    !! Fluxes at cell interface
    real, allocatable, dimension(:), intent(INOUT) :: flux
    !! Define the current time step.
    real, intent(IN) :: dt
    !! Define the primitive variables at the left and right boundary.
    real, allocatable, dimension(:) :: wl
    real, allocatable, dimension(:) :: wr
    real, allocatable, dimension(:) :: slope
    !! Define shock speed.
    real :: shock_s, character_s
    !! Allocate memory space for variables.
    allocate(    wl(nx1 + 1) )
    allocate(    wr(nx1 + 1) )
    allocate( slope(nx1 + 2) )
    !! Initialize shock speed, wr, wl, and slope.
    shock_s = 0.
    wr = 0.
    wl = 0.
    slope = 0.
    character_s = 0.

    !! Calculate the slope for non-TVD methods: upwind, centered & downwind.
    !! Index starts from (nghost) to (nx1 + 2 nghost - 1).
    do i = 1, nx1 + 2, 1
      !! Compute PLM slope for each cell, linear interpolation.
      !! Upwind method
      if (slope_type .eq. "upwind") then
        slope(i) = ( w(i + nghost - 1) - w(i + nghost - 2) )/dx1
      !! centered method
      else if (slope_type .eq. "centered") then
        slope(i) = ( w(i + nghost) - w(i - 2 + nghost) )/dx1
      !! Downwind method
      else if (slope_type .eq. "downwind") then
        slope(i) = ( w(i + nghost) - w(i + nghost - 1) )/dx1
      end if
    end do

    !! Compute the primitive variable at the left cell interface.
    do i = 1, nx1 + 1, 1
      !! Computer the shock speed.
      shock_s  = 0.5 * ( w(i + nghost - 1) + w(i + nghost) )
      !! Compute primitive variables at left cell interfaces.
      wl(i) = w(i + nghost) - 0.5 * dx1 * slope(i + 1) * &
              (1 + shock_s * dt/dx1)
    end do

    !! Compute the primitive variable at the right cell interface.
    do i = 1, nx1 + 1, 1
      !! Computer the shock speed.
      shock_s  = 0.5 * ( w(i + nghost - 1) + w(i + nghost) )
      !! Compute primitive variables at left cell interfaces.
      wr(i) = w(i + nghost - 1) + 0.5 * dx1 * slope(i) * &
              (1 - shock_s * dt/dx1)
    end do

    !! Compute fluxes at cell interfaces.
    do i = nghost, nx1 + nghost, 1
      shock_s = 0.5 * (w(i) + w(i+1))
      !! U_{i}^{n} > U_{i+1}^{n}
      if ( w(i) .ge. w(i+1) ) then
        !! S > 0
        if ( shock_s .gt. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wr(i - nghost + 1)**2.
        !! S < 0
        else if ( shock_s .le. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wl(i - nghost + 1)**2.
        end if

      !! U_{i}^{n} < U_{i+1}^{n}
      else if ( w(i) .lt. w(i+1) ) then
        !! 0 < U_{i}^{n}
        if ( w(i) .ge. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wr(i - nghost + 1)**2.
        !! U_{i+1}^{n} > 0
        else if ( w(i+1) .le. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wl(i - nghost + 1)**2.
        !! U_{i}^{n} < 0 < U_{i+1}^{n}
        else if ( (w(i) .lt. character_s) .and. &
                (w(i+1) .gt. character_s) ) then
          flux( i - nghost + 1 ) = 0.
        end if
      end if
    end do

  end subroutine second_non_TVD

  !! Second order reconstruction with TVD slope
  !! (Piecewise Linear Method with TVD slope, PLM).
  subroutine second_TVD(w,flux,dt)

    !! Primitive variables
    real, allocatable, dimension(:), intent(IN) :: w
    !! Fluxes at cell interface
    real, allocatable, dimension(:), intent(INOUT) :: flux
    !! Define the current time step.
    real, intent(IN) :: dt

    !! Define the primitive variables at the left and right boundary.
    real, allocatable, dimension(:) :: wl
    real, allocatable, dimension(:) :: wr
    real, allocatable, dimension(:) :: slope, upwind_slope, centered_slope, &
                                       downwind_slope
    !! Define shock speed.
    real :: shock_s, character_s
    !! Define signs of the slope for minmod function.
    real :: sign_a, sign_b
    !! Allocate memory space for primitive variables.
    allocate( wl(nx1 + 1) )
    allocate( wr(nx1 + 1) )
    !! Allocate memory space for slopes.
    allocate(   upwind_slope(nx1 + 2) )
    allocate( centered_slope(nx1 + 2) )
    allocate( downwind_slope(nx1 + 2) )
    allocate(          slope(nx1 + 2) )

    !! Initialize variables.
    shock_s = 0.
    wr = 0.
    wl = 0.
    upwind_slope = 0.
    centered_slope = 0.
    downwind_slope = 0.
    character_s = 0.
    sign_a = 0.
    sign_b = 0.

    !! Calculate the slope for TVD methods: upwind, centered & downwind.
    !! Index starts from (nghost) to (nx1 + 2 nghost - 1).
    do i = 1, nx1 + 2, 1
      !! Compute PLM slope for each cell, linear interpolation.
      !! Upwind method
      upwind_slope(i) = ( w(i + nghost - 1) - w(i + nghost - 2) )/dx1
      !! centered method
      centered_slope(i) = ( w(i + nghost) - w(i + nghost - 2) )/dx1
      !! Downwind method
      downwind_slope(i) = ( w(i + nghost) - w(i + nghost - 1) )/dx1
      !! Get the sign of the slope
      call get_sign(upwind_slope(i), sign_a)
      call get_sign(downwind_slope(i), sign_b)

      !! Minmod method.
      if (slope_type .eq. "minmod") then
        slope(i) = 0.5 * (sign_a + sign_b) * &
                   MIN(ABS( upwind_slope(i)), downwind_slope(i) )
      !! Van Leer method.
      else if (slope_type .eq. "van_leer") then
        if ( centered_slope(i) .ne. 0. ) then
          slope(i) = 0.5 * (sign_a + sign_b) * &
                     upwind_slope(i) * downwind_slope(i) / &
                     ABS(upwind_slope(i) + downwind_slope(i))
        else
          slope(i) = 0.
        end if
      !! Mc limiter method.
      else if (slope_type .eq. "mc_limiter") then
        slope(i) = 0.5 * (sign_a + sign_b) * &
                   MIN( MIN(ABS( upwind_slope(i)), downwind_slope(i) ), &
                       0.25 * ABS( upwind_slope(i)) + downwind_slope(i) )
      end if
    end do

    !! Compute the primitive variable at the left cell interface.
    do i = 1, nx1 + 1, 1
      !! Computer the shock speed.
      shock_s  = 0.5 * ( w(i + nghost - 1) + w(i + nghost) )
      !! Compute primitive variables at left cell interfaces.
      wl(i) = w(i + nghost) - 0.5 * dx1 * slope(i + 1) * &
              (1 + shock_s * dt/dx1)
    end do

    !! Compute the primitive variable at the right cell interface.
    do i = 1, nx1 + 1, 1
      !! Computer the shock speed.
      shock_s  = 0.5 * ( w(i + nghost - 1) + w(i + nghost) )
      !! Compute primitive variables at left cell interfaces.
      wr(i) = w(i + nghost - 1) + 0.5 * dx1 * slope(i) * &
              (1 - shock_s * dt/dx1)
    end do

    !! Compute fluxes at cell interfaces.
    do i = nghost, nx1 + nghost, 1
      shock_s = 0.5 * (w(i) + w(i+1))
      !! U_{i}^{n} > U_{i+1}^{n}
      if ( w(i) .ge. w(i+1) ) then
        !! S > 0
        if ( shock_s .gt. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wr(i - nghost + 1)**2.
        !! S < 0
        else if ( shock_s .le. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wl(i - nghost + 1)**2.
        end if

      !! U_{i}^{n} < U_{i+1}^{n}
      else if ( w(i) .lt. w(i+1) ) then
        !! 0 < U_{i}^{n}
        if ( w(i) .ge. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wr(i - nghost + 1)**2.
        !! U_{i+1}^{n} > 0
        else if ( w(i+1) .le. character_s ) then
          flux( i - nghost + 1 ) = 0.5 * wl(i - nghost + 1)**2.
        !! U_{i}^{n} < 0 < U_{i+1}^{n}
        else if ( (w(i) .lt. character_s) .and. &
                (w(i+1) .gt. character_s) ) then
          flux( i - nghost + 1 ) = 0.
        end if
      end if
    end do

  end subroutine second_TVD

end module riemann_solver
