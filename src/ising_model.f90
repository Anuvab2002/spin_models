!> @file ising_model.f90
!> @brief provides different routines for quantum isning model calculations
!> @author ap
module ising_model_m
contains
!> @brief function for calculating the spin-spin interaction terms of the Hamiltonian
!> @param[in]      s            spin matrix of given direction
!> @param[in]      n            size of the chain
!> @param[in]      bc           boundary condition
!> @param[in]      j            coupling strength
!> @return         spin_int     spin-spin interaction term
!> @todo unit test to be done
  function spin_spin_interaction(s, n, boundary, j)result(spin_int)
    use linear_algebra_helper_m
    use matrix_generator_m
    implicit none
    ! io variables
    complex(8), dimension(:,:), intent(in)    :: s
    integer, intent(in)                       :: n
    character(len=*), intent(in)              :: boundary
    double precision, intent(in)              :: j
    complex(8), allocatable, dimension(:,:)   :: spin_int
    ! internal variable
    integer                                   :: dim
    integer                                   :: dimm
    complex(8), allocatable, dimension(:,:)   :: op_s
    integer                                   :: idx
    complex(8), allocatable, dimension(:,:)   :: i1
    complex(8), allocatable, dimension(:,:)   :: i2
    integer                                   :: p
    integer                                   :: q
    complex(8), allocatable, dimension(:,:)   :: temp1
    complex(8), allocatable, dimension(:,:)   :: temp2
    !
    dim = size(s,1)
    dimm = dim**n
    !
    allocate(spin_int(dimm, dimm))
    spin_int = cmplx(0.d0, 0.d0)
    !
    allocate(op_s(dim**2,dim**2))
    call kron_product(s, s, op_s)
    !
    allocate(temp2(dimm,dimm))
    do idx = 1,n-1
      p = idx-1
      q = n-idx-1
      i1 = identity_matrix_complex(dim**p)
      i2 = identity_matrix_complex(dim**q)
      allocate(temp1(dim**(idx+1),dim**(idx+1)))
      call kron_product(i1, op_s, temp1)
      call kron_product(temp1, i2, temp2)
      deallocate(i1,i2)
      deallocate(temp1)
      spin_int = spin_int + j*temp2
    end do
    deallocate(temp2)
  end function spin_spin_interaction
end module ising_model_m
