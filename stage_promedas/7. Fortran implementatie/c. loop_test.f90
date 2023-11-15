! gfortran b.\ large_precision_test.f90 -o b.\ large_precision_test
program Loop_of_log
    implicit none
  
    integer :: i, n
    real(kind=4), dimension(100,101) :: y_4
    real(kind=4), dimension(100,101) :: x_4
    real(kind=8), dimension(100,101) :: x_8
    real(kind=8), dimension(100,101) :: y_8
    real(kind=16), dimension(100,101) :: x_16
    real(kind=16), dimension(100,101) :: y_16

    real(kind=8) :: t1_4,t2_4,dt_4
    real(kind=8) :: t1_8,t2_8,dt_8
    real(kind=8) :: t1_16,t2_16,dt_16

  
    ! Number of iterations
    n = 2**11
    print *, "log(x_matrix) --> n=2^11"
    
    call cpu_time(t1_4)
    do i = 1, n
        y_4 = log(x_4)
    end do
    call cpu_time(t2_4)
    dt_4 = t2_4-t1_4
    print *, "dt = ", dt_4, " (32-bit)"
    
    call cpu_time(t1_8)
    do i = 1, n
        y_8 = log(x_8)        
    end do
    call cpu_time(t2_8)
    dt_8 = t2_8-t1_8
    print *, "dt = ", dt_8, " (64-bit)"

    call cpu_time(t1_16)
    do i = 1, n
        y_16 = log(x_16)
    end do
    call cpu_time(t2_16)
    dt_16 = t2_16-t1_16
    print *, "dt = ", dt_16, " (128-bit)"
  
end program Loop_of_log

