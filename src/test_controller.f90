!> @file test_controller.f90
!> @brief controlls tests
!> @author ap
module test_controller_m
  use test_driver_m
  use test_helper_m
contains
  subroutine testing_list()
    implicit none
    logical             :: test_stat
    character(len=256)  :: testname
    !
    write(*,*) "==============================================================================="
    write(*,*) "Hi Po! Reporting from test controller. :-)>"
    write(*,*) "==============================================================================="
    !--------------------------------------------
    testname = "if_null_c"
    call test_banner_header(testname)
    call test_if_null_c(test_stat)
    if (test_stat) then
      write(*,*) "IF_NULL_C passed."
    else
      write(*,*) "IF_NULL_C failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "calculate_trace_c"
    call test_banner_header(testname)
    call test_calculate_trace_c(test_stat)
    if (test_stat) then
      write(*,*) "calculate_trace_c passed."
    else
      write(*,*) "calculate_trace_c failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "kron_product"
    call test_banner_header(testname)
    call test_kron_product(test_stat)
    if (test_stat) then
      write(*,*) "kron_product passed."
    else
      write(*,*) "kron_product failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "diagonalize_matrix"
    call test_banner_header(testname)
    call test_diagonalize_matrix(test_stat)
    if (test_stat) then
      write(*,*) "diagonalize_matrix passed."
    else
      write(*,*) "diagonalize_matrix failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "calculate_commutator_c"
    call test_banner_header(testname)
    call test_calculate_commutator_c(test_stat)
    if (test_stat) then
      write(*,*) "calculate_commutator_c passed."
    else
      write(*,*) "calculate_commutator_c failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "check_hermiticity"
    call test_banner_header(testname)
    call test_check_hermiticity(test_stat)
    if (test_stat) then
      write(*,*) "check_hermiticity passed."
    else
      write(*,*) "check_hermiticity failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "unitarity_check"
    call test_banner_header(testname)
    call test_unitarity_check(test_stat)
    if (test_stat) then
      write(*,*) "unitarity_check passed."
    else
      write(*,*) "unitarity_check failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    testname = "inner_product_dis"
    call test_banner_header(testname)
    call test_inner_product_dis(test_stat)
    if (test_stat) then
      write(*,*) "inner_product_dis passed."
    else
      write(*,*) "inner_product_dis failed."
    end if
    call test_banner_footer(testname)
    !---------------------------------------------
    write(*,*) "==============================================================================="
  end subroutine testing_list
end module test_controller_m
