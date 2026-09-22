!> @file matrix_generator.f90
!> @brief generates different kinds of matrices
!> @author ap
module matrix_generator_m
contains
!> @brief funciton for generating random complex matrix
!> @param[in]    r                number of rows
!> @param[in]    c                number of columns
!> @return       random_matrix    complex random matrix of given dimensions
  function random_complex_matrix(r,c) result(random_matrix)
    implicit none
    ! io variables
    integer, intent(in)             :: r
    integer, intent(in)             :: c
    complex(8), allocatable         :: random_matrix(:,:)
    ! internal variables
    double precision, allocatable   :: real_part(:,:)
    double precision, allocatable   :: imaginary_part(:,:)
    !
    allocate(random_matrix(r,c))
    allocate(real_part(r,c), imaginary_part(r,c))
    call random_number(real_part)
    call random_number(imaginary_part)
    real_part = real_part * 10.5d0
    imaginary_part = imaginary_part * 10.8d0
!
    random_matrix = cmplx(real_part, imaginary_part)
  end function random_complex_matrix
!> @brief funciton for generating complex identity matrix
!> @param[in]     dim                 number of rows and columns
!> @return        identity_matrix     identity matrix of given dimensions
  function identity_matrix_complex(dim) result(identity_matrix)
    implicit none
    ! io variables
    integer, intent(in) :: dim !dimension of the squre matrix (dim X dim)
    complex(8)          :: identity_matrix(dim,dim)
    ! internal variables
    integer             :: i
    integer             :: j
!
    identity_matrix = cmplx(0.0d0,0.0d0)
    do i = 1,dim
      identity_matrix(i,i) = cmplx(1.0d0, 0.0d0)
    end do
  end function identity_matrix_complex
!> @brief funciton for generating a complex square null matrix
!> @param[in]   dim                 number of rows and columns
!> @return      null_matrix         the null matrix of the given dimensions
  function null_matrix_complex(dim) result(null_matrix)
    implicit none
    ! io variable
    integer, intent(in)   :: dim !dimension of the squre matrix (dim X dim)
    complex(8)            :: null_matrix(dim,dim)
    ! internal variables
    integer :: i
    integer :: j
    !
    null_matrix = cmplx(0.0d0,0.0d0)
  end function null_matrix_complex
!> @brief function for constructing Pauli matrices
!> @param[in]     n_pauli         number tag of the Pauli matrices
!> @return        pauli_matrix    Pauli matrix of the given tag
  function pauli_matrices(n_pauli)result(pauli_matrix)
    implicit none
    ! io variables
    integer, intent(in)  :: n_pauli
    complex(8)           :: pauli_matrix(2,2)
    !the order in which fortran reads and assigns the alemets: (1,1), (2,1), (1,2), (2,2)
    select case(n_pauli)
    case(1)
      pauli_matrix = reshape([cmplx(0.0d0, 0.0d0), cmplx(1.0d0, 0.0d0), &
                     cmplx(1.0d0, 0.0d0), cmplx(0.0d0, 0.0d0)], [2,2]) !S_x
    case(2)
      pauli_matrix = reshape([cmplx(0.0d0, 0.0d0), cmplx(0.0d0, 1.0d0), &
                     cmplx(0.0d0, -1.0d0), cmplx(0.0d0, 0.0d0)], [2,2]) !S_y
    case(3)
      pauli_matrix = reshape([cmplx(1.0d0, 0.0d0), cmplx(0.0d0, 0.0d0), &
                     cmplx(0.0d0, 0.0d0), cmplx(-1.0d0, 0.0d0)], [2,2]) !S_z
    end select
  end function pauli_matrices
end module matrix_generator_m
