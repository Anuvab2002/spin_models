!> @file math_helper.f90
!> @brief provides necessary routines for mathematical operations
!> @author ap
module math_helper_m
contains
!> @brief function for constructing delta function
!> @param[in]     idx      one integer
!> @param[in]     jdx      another integer
!> @return        delta    value of the delta function
  function delta_ij(idx, jdx)result(delta)
    use global_m
    implicit none
    ! io variables
    integer, intent(in)             :: idx
    integer, intent(in)             :: jdx
    double precision                :: delta
    ! internal variables
    !
    if (idx .eq. jdx) then
      delta = 1.d0
    else
      delta = 0.d0
    end if
  end function delta_ij
!> @brief function for constructing three variable levi-civita
!> @param[in]     idx       one integer
!> @param[in]     jdx       another integer
!> @param[in]     kdx       another integer
!> @return        epsilon   value of the levi-civita
  function levi_civita(idx, jdx, kdx)result(epsilon)
    implicit none
    ! io variables
    integer, intent(in)     :: idx
    integer, intent(in)     :: jdx
    integer, intent(in)     :: kdx
    double precision        :: epsilon
    !
    if (idx.eq.jdx .or. jdx.eq.kdx .or. kdx.eq.idx) then
      epsilon = 0.d0
    else if (idx.lt.jdx .and. jdx.lt.kdx) then
      epsilon = 1.d0
    else if (idx.lt.jdx .and. jdx.gt.kdx .and. idx.gt.kdx) then
      epsilon = 1.d0
    else if (idx.gt.jdx .and. jdx.lt.kdx .and. idx.gt.kdx) then
      epsilon = 1.d0
    else
      epsilon = -1.d0
    end if
  end function levi_civita
!> @brief function for calculating factorial of an integer
!> @param[in]         n         input number
!> @return            n_fact    factorial of the number
  recursive function factorial(n)result(n_fact)
    implicit none
    ! io variables
    integer, intent(in)             :: n
    double precision                :: n_fact
    ! internal variables
    integer                         :: idx
    !
    if (n .eq. 0) then
      n_fact = 1.d0
    end if
    if (n .ge. 1) then
      n_fact = n*factorial(n-1)
    end if
  end function factorial
!> @brief subroutine for partial Baker-Campbell-Hausdorff calculations of given order
!> @note e^{A}Be^{-A} = B + [A,B] + 1/2! [A,[A,B]] + 1/3! [A,[A,[A,B]]] + ...
!> @param[in]     m1      one matrix
!> @param[in]     m2      another matrix
!> @param[in]     order   order of the calculation
!> @return        bch_m   calculated matrix
!> @todo unit testing to be done.
  subroutine bch_c(m1, m2, order, bch_m)
    use linear_algebra_helper_m
    implicit none
    ! io variables
    complex(8), dimension(:,:), intent(in)                   :: m1
    complex(8), dimension(:,:), intent(in)                   :: m2
    integer, intent(in)                                      :: order
    complex(8), allocatable, dimension(:,:), intent(out)     :: bch_m
    ! internal variable
    integer                                                  :: idx
    integer                                                  :: r
    integer                                                  :: c
    complex(8), allocatable, dimension(:,:)                  :: temp
    !
    r = size(m1,1)
    c = size(m1,2)
    if (size(m2,1).ne.r .or. size(m2,2).ne.c) then
      write(*,*) "Execution error! BCH calculation can't be done due to size mismatch!"
      bch_m = cmplx(0.d0,0.d0)
      return
    end if
    !
    allocate(temp(r,c))
    allocate(bch_m(r,c))
    temp = m2
    bch_m = m2
    do idx = 1,order
      call calculate_commutator_c(m1, temp, temp)
      bch_m = bch_m + (1.d0/factorial(order))*temp
    end do
    deallocate(temp)
  end subroutine bch_c
end module math_helper_m
