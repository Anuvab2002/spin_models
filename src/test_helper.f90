!> @file test_helper.f90
!> @brief provides assisting routines for unit testing
!> @author ap
module test_helper_m
contains
!> @brief povides test_banner_header
!> @param[in]     testname        name of the test
  subroutine test_banner_header(testname)
    implicit none
    character(len=256), intent(in)    :: testname
    !
    write(*,*) ""
    write(*,*) "---------------- Test: ", trim(adjustl(testname)), " ----------------"
  end subroutine test_banner_header
!> @brief povides test_banner_footer
!> @param[in]     testname        name of the test
  subroutine test_banner_footer(testname)
    implicit none
    character(len=256), intent(in)    :: testname
    !
    write(*,*) "---------------- ----- ----- ----- ----- ----- ----- ----------------"
  end subroutine test_banner_footer
end module test_helper_m
