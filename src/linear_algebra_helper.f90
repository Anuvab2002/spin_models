module linear_algebra_helper_m
  use global_m
contains
!> @brief subroutine for checking whether a complex matrix is null or not
!> @param[in]     a_matrix        input matrix
!> @param[out]    null_status     status of whether the matrix is null or not
!> @todo unit testing
  subroutine if_null_c(a_matrix, null_status)
    implicit none
    ! io variables
    complex(8), intent(in)    :: a_matrix(:,:)
    logical, intent(out)      :: null_status
    ! internal variables
    complex(8)                :: sum_elements
    double precision          :: sum_elements_abs
    !
    sum_elements = sum(abs(a_matrix))
    sum_elements_abs = abs(sum_elements)
    !
    if (sum_elements_abs .lt. zero) then
      null_status = .true.
    else
      null_status = .false.
    end if
  end subroutine if_null_c
!> @brief subroutine for calculating the trace of a complex matrix
!> @param[in]     a_matrix      input matrix
!> @param[out]    a_trace       trace of the matrix
!> @todo unit testing
  subroutine calculate_trace_c(a_matrix, a_trace)
    implicit none
    ! io variables
    complex(8), intent(in)    :: a_matrix(:,:)
    complex(8), intent(out)   :: a_trace
    ! internal variables
    integer                   :: i
    integer                   :: j
    integer                   :: a_r
    integer                   :: a_c
    !
    a_r = size(a_matrix,1)
    a_c = size(a_matrix,2)
    !
    if (a_r .ne. a_c) then
      write(*,*) "Execution error! Trace can't be calculated, not a square matrix!"
      return !EXITS the subroutine
    end if
    !
    a_trace = cmplx(0.0d0, 0.0d0)
    do i = 1,a_r
       a_trace = a_trace + a_matrix(i,i)
    end do
  end subroutine calculate_trace_c
!> @brief subroutine for calculating kronecker product of two complex matrices
!> @param[in]   a_matrix      one input matrix
!> @param[in]   b_matrix      another input matrix
!> @param[out]  c_matrix      the kronecker product matrix
!> @todo this subroutine is a bit tidious as it contains for nested loops, and therefore needs more effective implimentation.
!> @todo unit testing
  subroutine kron_product(a_matrix, b_matrix, c_matrix)
    implicit none
    ! io variables
    complex(8), intent(in)                :: a_matrix(:,:)
    complex(8), intent(in)                :: b_matrix(:,:)
    complex(8), intent(out), allocatable  :: c_matrix(:,:)
    integer                               :: i
    integer                               :: j
    integer                               :: k
    integer                               :: l
    integer                               :: m
    integer                               :: n
    !
    allocate(c_matrix(size(a_matrix, 1) * size(b_matrix, 1), &
       size(a_matrix, 2) * size(b_matrix, 2))) !C(ma*mb,na*nb)
    !
    do i = 1, size(a_matrix, 1)
      do j = 1, size(a_matrix, 2)
        do k = 1, size(b_matrix, 1)
          do l = 1, size(b_matrix, 2)
            m = (i - 1) * size(b_matrix, 1) + k! row indx of an element of the product matrix
            n = (j - 1) * size(b_matrix, 2) + l! col indx of an element of the product matrix
            c_matrix(m, n) = a_matrix(i, j) * b_matrix(k, l) !calculates the particlar element
          end do
        end do
      end do
    end do
  end subroutine kron_product
!> @brief subroutine for diagonalizing a square complex matrix
!> @param[in]    n_dim      dimension of the matrix (n_dim X n_dim)
!> @param[in]    mat        input matrix
!> @param[out]   eig_vect   the eigen vector matrix
!> @param[out]   eig_val    the diagonal eigen value matrix
!> @note LAPACK has been used for this subroutine
!> @todo unit testing
  subroutine diagonalize_matrix(n_dim, mat, eig_vect, eig_vals)
    implicit none
    ! io variables
    integer, intent(in)                :: n_dim
    complex(8), intent(in)             :: mat(n_dim,n_dim)
    complex(8), intent(out)            :: eig_vect(n_dim,n_dim)
    double precision, intent(out)      :: eig_vals(n_dim)
    ! internal variables
    integer                            :: lda
    integer                            :: lwork
    integer                            :: info
    complex(8), allocatable            :: work(:)
    double precision, allocatable      :: rwork(:)
    !
    ! Derived dimensions
    lda = n_dim
    lwork = 2 * n_dim - 1
    !
    ! Allocate workspace
    allocate(work(lwork), rwork(3 * n_dim - 2))
    !
    ! Copy input matrix to eig_vect for in-place diagonalization
    eig_vect = mat
    !
    ! Call LAPACK's ZHEEV routine
    call zheev('V', 'U', n_dim, eig_vect, lda, eig_vals, work, lwork, rwork, info)
    !
    ! Check for errors
    if (info /= 0) then
      print *, "Error: LAPACK ZHEEV routine failed with INFO =", info
      deallocate(work, rwork)
      return
    end if
    ! Deallocate workspace
    deallocate(work, rwork)
  end subroutine diagonalize_matrix
!> @brief subroutine for calculating cmmutator between two complex matrices
!> @param[in]     a_matrix      one input matrix
!> @param[in]     b_matrix      another input matrix
!> @param[out]    c_matrix      the commutator matrix
!> @todo unit testing
  subroutine calculate_commutator_c(a_matrix, b_matrix, commutator)
    implicit none
    ! io variables
    complex(8), intent(in) :: a_matrix(:,:)
    complex(8), intent(in) :: b_matrix(:,:)
    complex(8), intent(out), allocatable :: commutator(:,:)
    !
    allocate(commutator(size(a_matrix,1),size(a_matrix,2)))
    commutator = matmul(a_matrix, b_matrix) - matmul(b_matrix, a_matrix)
  end subroutine calculate_commutator_c
!> @brief subroutine for calculating the hermitian conjugate of a complex matrix
!> @param[in]    a_matrix       input matrix
!> @param[out]  a_herm_conjg    hermitian conjugate of the matrix
!> @todo unit testing
  subroutine hermitian_conjugate(a_matrix, a_herm_conjg)
    implicit none
    ! io variables
    complex(8), intent(in)                :: a_matrix(:,:)
    complex(8), intent(out), allocatable  :: a_herm_conjg(:,:)
    ! internal variables
    complex(8), allocatable               :: a_transpose(:,:)
    integer                               :: a_r
    integer                               :: a_c
    !
    a_r = size(a_matrix,1)
    a_c = size(a_matrix,2)
    !
    allocate(a_transpose(a_c, a_r))
    a_transpose = transpose(a_matrix)
    allocate(a_herm_conjg(a_c,a_r))
    a_herm_conjg = conjg(a_transpose)
    deallocate(a_transpose)
  end subroutine hermitian_conjugate
!> @brief subroutine for checking whether a complex matrix is hermitian or not
!> @param[in]     a_matrix      input matrix
!> @param[out]    herm_stat     status for hermiticity
!> @todo unit testing
  subroutine check_hermiticity(a_matrix, herm_stat)
    implicit none
    ! io variables
    complex(8), intent(in)        :: a_matrix(:,:)
    logical, intent(out)          :: herm_stat
    ! internal variables
    complex(8), allocatable       :: a_herm_conjg(:,:)
    complex(8), allocatable       :: diff(:,:)
    integer                       :: a_r
    integer                       :: a_c
    integer                       :: i
    integer                       :: j
!
    a_r = size(a_matrix,1)
    a_c = size(a_matrix,2)
    if (a_r .ne. a_c) then
      write(*,*) "Execution error! This is not a square matrix so can't be hermitian!"
      herm_stat = .false.
      return
    end if
!
    allocate(a_herm_conjg(a_r,a_r)) !square matrix
    call hermitian_conjugate(a_matrix, a_herm_conjg)
!
    allocate(diff(a_r,a_r)) !square matrix
    diff = a_matrix - a_herm_conjg
!
    call if_null(diff, herm_stat)
  end subroutine check_hermiticity
!> @brief subroutine for checking whether a matrix is right unitary or not
!> @param[in]     a_matrix      input matrix
!> @param[out]    ru_stat       status for right unitarity
!> @todo unit testing
  subroutine right_unitarity_check(a_matrix, ru_stat)
    use matrix_generator_m
    implicit none
    ! io variable
    complex(8), intent(in)     :: a_matrix(:,:)
    logical, intent(out)       :: ru_stat
    ! internal variable
    complex(8), allocatable    :: a_herm_conjg(:,:)
    complex(8), allocatable    :: mult(:,:)
    complex(8), allocatable    :: diff(:,:)
    complex(8), allocatable    :: identity_matrix(:,:)
    integer                    :: a_r
    integer                    :: a_c
    !
    a_r = size(a_matrix, 1)
    a_c = size(a_matrix, 2)
    if (a_r .ne. a_c) then
      print*, "Execution error! It's not a square matrix, therefore can't be unitary!"
      ru_stat = .false.
      return
    end if
    !
    allocate(a_herm_conjg(a_r,a_c)) !square matrix
    call hermitian_conjugate(a_matrix, a_herm_conjg)
    !
    allocate(mult(a_r,a_c))
    mult = matmul(a_matrix, a_herm_conjg)
    !
    allocate(diff(a_r,a_c))
    identity_matrix = identity_matrix_complex(a_r)
    diff = identity_matrix - mult
    call if_null(diff, ru_stat)
    deallocate(a_herm_conjg, mult, diff)
  end subroutine right_unitarity_check
!> @brief subroutine for checking whether a matrix is left unitary or not
!> @param[in]     a_matrix      input matrix
!> @param[out]    lu_stat       status for left unitarity
!> @todo unit testing
  subroutine left_unitarity_check(a_matrix, lu_stat)
    use matrix_generator_m
    implicit none
    ! io variable
    complex(8), intent(in)     :: a_matrix(:,:)
    logical, intent(out)       :: lu_stat
    ! internal variable
    complex(8), allocatable    :: a_herm_conjg(:,:)
    complex(8), allocatable    :: mult(:,:)
    complex(8), allocatable    :: diff(:,:)
    complex(8), allocatable    :: identity_matrix(:,:)
    integer                    :: a_r
    integer                    :: a_c
    !
    a_r = size(a_matrix, 1)
    a_c = size(a_matrix, 2)
    if (a_r .ne. a_c) then
      print*, "Execution error! It's not a square matrix, therefore can't be unitary!"
      lu_stat = .false.
      return
    end if
    !
    allocate(a_herm_conjg(a_r,a_r)) !square matrix
    call hermitian_conjugate(a_matrix, a_herm_conjg)
    !
    allocate(mult(a_r,a_r))
    mult = matmul(a_herm_conjg, a_matrix)
    !
    allocate(diff(a_r,a_r))
    identity_matrix = identity_matrix_complex(a_r)
    diff = abs(identity_matrix - mult)
    call if_null(diff, lu_stat)
    deallocate(a_herm_conjg, mult, diff)
  end subroutine left_unitarity_check
!> @brief subroutine for checking whether a matrix is unitary or not
!> @param[in]     a_matrix        input matrix
!> param[out]     u_stat          status tag for unitarity
!> @todo unit testing
  subroutine unitarity_check(a_matrix, u_stat)
    implicit none
    ! io variables
    complex(8), intent(in)    :: a_matrix(:,:)
    logical, intent(out)      :: u_stat
    ! internal variable
    logical                   :: lu_stat
    logical                   :: ru_stat
    !
    call right_unitarity_check(a_matrix, ru_stat)
    call left_unitarity_check(a_matrix, lu_stat)
    !
    if (ru_stat .and. lu_stat) then
      u_stat = .true.
    else
      u_stat = .false.
    end if
  end subroutine unitarity_check
!> @brief subroutine for checking whether a complex matrix is involutory or not
!> @param[in]      a_matrix       input matrix
!> param[out]      invol_stat      status tag for the involutarity
!> @todo unit testing
  subroutine if_involutory(a_matrix, invol_stat)
    use matrix_generator_m
    implicit none
    ! io variables
    complex(8), intent(in)    :: a_matrix(:,:)
    logical, intent(out)      :: invol_stat
    ! internal variable
    integer                   :: a_dim
    complex(8), allocatable   :: a_matrix_sq(:,:)
    complex(8), allocatable   :: diff(:,:)
    complex(8), allocatable   :: identity(:,:)
    !
    if (size(a_matrix,1) .ne. size(a_matrix,2)) then
      print*, "Execution error! Not a square matrix, so can't be idempotent!"
      invol_stat = .false.
      return
    end if
    a_dim = size(a_matrix,1)! must be a square matrix
    !
    allocate(a_matrix_sq(a_dim,a_dim))
    a_matrix_sq = matmul(a_matrix,a_matrix)
    !
    identity = identity_matrix_complex(a_dim)! idenitity matrix of the same dimensions
    allocate(diff(a_dim,a_dim))
    diff = abs(a_matrix_sq - identity)
    !
    call if_null(diff, invol_stat)
    deallocate(a_matrix_sq, identity, diff)
  end subroutine if_involutory
!> @brief subroutine for checking whether a matrix is idempotent or not
!> @param[in]       a_matrix        input matrix
!> @param[out]      idem_stat       status tag for idempotency
!> @todo unit testing
  subroutine if_idempotent(a_matrix, idem_stat)
    ! io variables
    implicit none
    complex(8), intent(in)   :: a_matrix(:,:)
    logical, intent(out)     :: idem_stat
    ! internal variables
    integer                  :: a_dim
    complex(8), allocatable  :: a_matrix_sq(:,:)
    complex(8), allocatable  :: diff(:,:)
    !
    if (size(a_matrix,1) .ne. size(a_matrix,2)) then
      print*, "Execution error! Not a square matrix, so can't be idempotent!"
      idem_stat = .false.
      return
    end if
    a_dim = size(a_matrix,1)! must be a square matrix
    !
    allocate(a_matrix_sq(a_dim,a_dim))
    a_matrix_sq = matmul(a_matrix,a_matrix)
    !
    allocate(diff(a_dim,a_dim))
    diff = abs(a_matrix_sq - a_matrix)
    !
    call if_null(diff, idem_stat)
    deallocate(a_matrix_sq, diff)
  end subroutine if_idempotent
end module linear_algebra_helper_m
