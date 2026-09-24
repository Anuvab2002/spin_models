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
    if (idx.lt.jdx .and. jdx.lt.kdx) then
      epsilon = 1.0d0
    else if (idx.eq.jdx .or. jdx.eq.kdx .or. kdx.eq.idx) then
      epsilon = 0.d0
    else
      epsilon = -1.d0
    end if
  end function levi_civita
end module math_helper_m
