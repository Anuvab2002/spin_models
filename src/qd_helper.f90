!> @file qd_helper.f90
!> @brief provides different helper routines for performing quantum dynamics
!> @author ap
module qd_helper_m
contains
!> @brief function for calculating inner poducts of two state vectors
!> @param[in]     v_1         vector 1
!> @param[in]     v_2         vector 2
!> @param[in]     weight      weight for the sum
!> @return        val         the inner product value
!> @todo unit testing
  function inner_product(v_1,v_2, weight)result(val)
    implicit none
    ! io variables
    complex(8), dimension(:), intent(in)            :: v_1
    complex(8), dimension(:), intent(in)            :: v_2
    double precision, dimension(:), intent(in)      :: weight
    complex(8)                                      :: val
    ! internal variables
    integer                                         :: idx
    integer                                         :: n
    complex(8), allocatable, dimension(:)           :: bra_v
    complex(8), allocatable, dimension(:)           :: ket_v
    !
    n = size(ket_v)
    if (size(weight).ne.n .or. size(bra_v).ne.n) then
      write(*,*) "Execution error! Size inconsitency in inner product calculation."
      val = 0.d0
      return
    end if
    !
    allocate(bra_v(n), ket_v(n))
    do idx =1,n
      ket_v(idx) = v_1(idx)*weight(idx)
    end do
    bra_v = conjg(v_2)
    !
    val = dot_product(bra_v,ket_v)
    deallocate(bra_v, ket_v)
  end function inner_product
!> @brief function for calculating the expectation value of some operator in certain quantum state
!> @param[in]     state_v       state vector
!> @param[in]     op_mat        operator matrix
!> @param[in]     weight        weight for summation
!> @return        val           expectation value
!> @todo unit testing
  function expectation_value(state_v, op_mat, weight)result(val)
    implicit none
    ! io variable
    complex(8), dimension(:), intent(in)          :: state_v
    complex(8), dimension(:,:), intent(in)        :: op_mat
    double precision, dimension(:), intent(in)    :: weight
    complex(8)                                    :: val
    ! internal variable
    integer                                       :: idx
    integer                                       :: n
    complex(8), allocatable, dimension(:)         :: bra_v
    complex(8), allocatable, dimension(:)         :: ket_v
    !
    n = size(state_v)
    if (size(op_mat,1).ne.n .or. size(op_mat,2).ne.n .or. size(weight).ne.n) then
      write(*,*) "Execution error! Size inconsitency in expectation value calculation."
      val = 0.d0
      return
    end if
    !
    allocate(bra_v(n), ket_v(n))
    bra_v = conjg(state_v)
    ket_v = matmul(op_mat, state_v)
    do idx = 1,n
      ket_v(idx) = ket_v(idx)*weight(idx)
    end do
    !
    val = dot_product(bra_v,ket_v)
    deallocate(bra_v, ket_v)
  end function expectation_value
end module qd_helper_m
