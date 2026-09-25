!> @file qd_helper.f90
!> @brief provides different helper routines for quantum dynamics calculations
!> @author ap
module qd_helper_m
contains
!> @brief function for calculating inner poducts of two state vectors in discrete basis
!> @param[in]     v_1         vector 1
!> @param[in]     v_2         vector 2
!> @return        val         the inner product value
!> @todo unit testing
  function inner_product_dis(v_1, v_2)result(val)
    implicit none
    ! io variables
    complex(8), dimension(:), intent(in)            :: v_1
    complex(8), dimension(:), intent(in)            :: v_2
    complex(8)                                      :: val
    ! internal variables
    integer                                         :: idx
    integer                                         :: n
    complex(8), allocatable, dimension(:)           :: bra_v
    complex(8), allocatable, dimension(:)           :: ket_v
    !
    n = size(v_1)
    if (size(v_2).ne.n) then
      write(*,*) "Execution error! Size inconsitency in inner product calculation."
      val = 0.d0
      return
    end if
    !
    allocate(bra_v(n), ket_v(n))
    ket_v = v_2
    bra_v = v_1
    !
!> @note dot_product takes care of taking conjugate of the bra vector.
    val = dot_product(bra_v,ket_v)
    deallocate(bra_v, ket_v)
  end function inner_product_dis
!> @brief function for calculating the expectation value of some operator in certain quantum state in discrete basis
!> @param[in]     state_v       state vector
!> @param[in]     op_mat        operator matrix
!> @return        val           expectation value
!> @todo unit testing
  function expectation_value_dis(state_v, op_mat)result(val)
    implicit none
    ! io variable
    complex(8), dimension(:), intent(in)          :: state_v
    complex(8), dimension(:,:), intent(in)        :: op_mat
    complex(8)                                    :: val
    ! internal variable
    integer                                       :: idx
    integer                                       :: n
    complex(8), allocatable, dimension(:)         :: bra_v
    complex(8), allocatable, dimension(:)         :: ket_v
    !
    n = size(state_v)
    if (size(op_mat,1).ne.n .or. size(op_mat,2).ne.n) then
      write(*,*) "Execution error! Size inconsitency in expectation value calculation."
      val = 0.d0
      return
    end if
    !
    allocate(bra_v(n), ket_v(n))
    bra_v = state_v
    ket_v = matmul(op_mat, state_v)
    !
    val = dot_product(bra_v,ket_v)
    deallocate(bra_v, ket_v)
  end function expectation_value_dis
!> @brief function for calculating operators in interaction picture
!> @param[in]          op_s        operator in Schrodinger picture
!> @param[in]          ham_0       time independent part of the hamiltonian
!> @param[in]          t           time
!> return              op_i        operator in interaction picture
!> @unit testing to be done
  function operator_interaction(op_s, ham_0, t, bch_order)result(op_i)
    use math_helper_m
    implicit none
    ! io variables
    complex(8), intent(in), dimension(:,:)      :: op_s
    complex(8), intent(in), dimension(:,:)      :: ham_0
    double precision, intent(in)                :: t
    integer, intent(in)                         :: bch_order
    complex(8), allocatable, dimension(:,:)     :: op_i
    ! internal variable
    !
    call bch_c(ham_0*(iota*t/hbar), op_s, bch_order, op_i)
  end function
end module qd_helper_m
