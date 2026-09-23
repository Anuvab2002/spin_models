!> @file global.f90
!> @brief provides global parameters
!> @author ap
module global_m
  double precision, parameter           :: zero=1.0d-8! zero
  double precision, parameter           :: tol=1.0d-5! tolerance
  complex(8), parameter                 :: iota=(0.d0,1.d0)! iota
  double precision, parameter           :: pi=4.d0*atan(1.d0)! pi
end module global_m
