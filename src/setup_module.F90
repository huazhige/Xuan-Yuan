!!------------------------------------------------------------------
!! A Fortran example code for solving advection-diffusion function.
!!
!! This code is written by Huazhi Ge for AMS 209 final project.
!!
!! * Three parts of solving the PDEs.
!! 1) 1D advection solution
!!
!! This program is used to use the CFL condition to choose a suitable
!! foot step.
!!
!!------------------------------------------------------------------

module setup_module

#include "definition.h"

  use read_initFile_module

  implicit none

  real, save :: left_speed, right_speed
  real, save :: x1Beg, x1End, x3End, CFL
  real, save :: limit_time, threshold, printScale
  integer, save :: nx1, limit_dt, nghost, reconstruction
  character(len=80), save :: x1_inner, x1_outer, problem, slope_type, method
  integer, save :: multiplicity

contains

  subroutine setup_init()

    implicit none

    !! Import keywords and parameters about dynamical core.
    call read_initFileChar('advection.init','problem',problem)
    PRINT 100, 'The problem name is: ', problem
    call read_initFileInt('advection.init','reconstruction',reconstruction)
    PRINT 200, 'The reconstruction method is: ', reconstruction
    call read_initFileChar('advection.init','slope_type',slope_type)
    PRINT 100, 'Reconstruction method: ', slope_type
    call read_initFileChar('advection.init','method',method)
    PRINT 100, 'TVD method: ', method
    !! Import parameters of time and convergence.
    call read_initFileInt ('advection.init', 'limit_dt', limit_dt)
    PRINT 200, 'Maximum number of time steps: ', limit_dt
    call read_initFileReal ('advection.init', 'limit_time', limit_time)
    PRINT 300, 'Time limitation: ', limit_time
    call read_initFileReal ('advection.init', 'printScale', printScale)
    PRINT 300, 'Print Scale: ', printScale
    call read_initFileReal('advection.init', 'CFL', CFL)
    PRINT 300, 'CFL condition parameter: ', CFL
    !!call read_initFileReal('advection.init', 'threshold', threshold)
    !!PRINT 400, 'Threshold value for solution convergence: ', threshold
    !! Import parameters of problem.
    call read_initFileReal ('advection.init', 'left_speed', left_speed)
    PRINT 300, 'Velocity of the advection of the left part: ', left_speed
    call read_initFileReal ('advection.init', 'right_speed', right_speed)
    PRINT 300, 'Velocity of the advection of the right part: ', right_speed
    !! Import parameters of mesh.
    call read_initFileInt ('advection.init', 'nx1', nx1)
    PRINT 200, 'Number of cells in x1 direction: ', nx1
    !!call read_initFileInt ('advection.init', 'nx2', nx2)
    !!PRINT 200, 'Number of cells in x2 direction: ', nx2
    !!call read_initFileInt ('advection.init', 'nx3', nx3)
    !!PRINT 200, 'Number of cells in x3 direction: ', nx3
    call read_initFileInt ('advection.init', 'nghost', nghost)
    PRINT 200, 'Number of ghost cells: ', nghost
    call read_initFileReal('advection.init', 'x1_beg', x1Beg)
    PRINT 300, 'Setting up the boundary: ', x1Beg
    call read_initFileReal('advection.init', 'x1_end', x1End)
    PRINT 300, 'Setting up the boundary: ', x1End
    !!call read_initFileReal('advection.init', 'x2_beg', x2Beg)
    !!PRINT 300, 'Setting up the boundary: ', x2Beg
    !!call read_initFileReal('advection.init', 'x2_end', x2End)
    !!PRINT 300, 'Setting up the boundary: ', x2End
    !!call read_initFileReal('advection.init', 'x3_beg', x3Beg)
    !!PRINT 300, 'Setting up the boundary: ', x3Beg
    !!call read_initFileReal('advection.init', 'x3_end', x3End)
    !!PRINT 300, 'Setting up the boundary: ', x3End
    call read_initFileChar('advection.init','x1_inner_bc', x1_inner)
    PRINT 100, 'The inner boundary condition of x1 direction is: ', x1_inner
    call read_initFileChar('advection.init','x1_outer_bc', x1_outer)
    PRINT 100, 'The inner boundary condition of x1 direction is: ', x1_outer
    !!call read_initFileChar('advection.init','x2_inner_bc', x2_inner)
    !!PRINT 100, 'The inner boundary condition of x2 direction is: ', x2_inner
    !!call read_initFileChar('advection.init','x2_outer_bc', x2_outer)
    !!PRINT 100, 'The inner boundary condition of x2 direction is: ', x2_outer
    !!call read_initFileChar('advection.init','x3_inner_bc', x3_inner)
    !!PRINT 100, 'The inner boundary condition of x3 direction is: ', x3_inner
    !!call read_initFileChar('advection.init','x3_outer_bc', x3_outer)
    !!PRINT 100, 'The inner boundary condition of x3 direction is: ', x3_outer

100 FORMAT(A70,A15)
200 FORMAT(A70,I7)
300 FORMAT(A70,F7.3)
400 FORMAT(A70,F10.9)
500 FORMAT(A70,I1)

  end subroutine setup_init


end module setup_module
