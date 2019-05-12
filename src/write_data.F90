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

  use setup_module, only : runName
  implicit none

contains

  subroutine writeDiff(nIter, printCell)
    implicit none
    integer, intent(IN) :: length
    real, dimension(:), intent(IN) :: x,f,residual
    integer :: i
    character(len=35) :: ofile

    !! File name for ascii output
    ofile = 'rootFinder_'//trim(runName)//'.dat'

    !! File open
    open(unit=20,file=ofile,status='unknown')

    !! Write into a file:
    !!   iteration number, search result x, function value f(x), and residual
    do i=1,length
       write(20,920)i,x(i),f(i),residual(i)
    end do

    !! Output format specifiers
920 format(1x, i5, f16.8, 1x, f16.8, 1x, f16.12)

    !! File close
    close(20)

  end subroutine output_write

end module output_module
