! gfortran b.\ large_precision_test.f90 -o b.\ large_precision_test
program Loop_of_log
    implicit none
  
    integer(kind=8) :: n
    integer :: i
    real(kind=4) :: x_4
    real(kind=8) :: x_8
    real(kind=16) :: x_16

    real(kind=8) :: t1_4,t2_4,dt_4
    real(kind=8) :: t1_8,t2_8,dt_8
    real(kind=8) :: t1_16,t2_16,dt_16

  
    ! Number of iterations
    n = 5*10**9
    
    call cpu_time(t1_4)
    do i = 1, n
        x_4 = 1.0 + real(i)
    end do
    call cpu_time(t2_4)
    dt_4 = t2_4-t1_4
    
    call cpu_time(t1_8)
    do i = 1, n
        x_8 = 1.0 + real(i)
    end do
    call cpu_time(t2_8)
    dt_8 = t2_8-t1_8

    call cpu_time(t1_16)
    do i = 1, n
        x_16 = 1.0 + real(i)
    end do
    call cpu_time(t2_16)
    dt_16 = t2_16-t1_16

    print *, "dt = ", dt_4, " (32-bit)"
    print *, "dt = ", dt_8, " (64-bit)"
    print *, "dt = ", dt_16, " (128-bit)"
  
end program Loop_of_log

