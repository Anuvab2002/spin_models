program main_p
  use test_controller_m
  implicit none
  integer             :: num_args
  character(len=10)   :: arg
  num_args = command_argument_count()
  if (num_args>0) then
    call get_command_argument(1,arg)
    if (trim(adjustl(arg))=="test") then
      call test_driver()
    else
      write(*,*) "Invalid argument! Execution stopped."
    end if
  end if
end program main_p
