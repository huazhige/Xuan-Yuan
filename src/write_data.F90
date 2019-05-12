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
!! This program is used write data into the .dat files.
!!
!!------------------------------------------------------------------


!! Don't change anything in this file.
module write_data

  use setup_module, only : problem, method, slope_type, nx1, nghost
  implicit none

contains

  subroutine write_out(nIter, dx1_array, u)
    implicit none
    integer, intent(IN) :: nIter
    real, dimension(:), intent(IN) :: dx1_array, u
    integer :: i
    character(len=35) :: ofile

    !! File name for ascii output
    ofile = "./output/"//trim(problem)//"_"//&
            trim(slope_type)//'.dat'

    !! File open
    open(unit = 20, file = ofile, status='unknown')

    !! Write into a file:
    !! iteration number, search result x, function value f(x), and residual
    do i = nghost + 1, nx1 + nghost , 1
       write(20,920) i , dx1_array(i), u(i)
    end do

    !! Output format specifiers
920 format(1x, i5, f16.8, 1x, f16.8)

    !! File close
    close(20)

  end subroutine write_out

end module write_data
