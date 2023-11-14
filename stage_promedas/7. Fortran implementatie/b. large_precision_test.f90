! gfortran b.\ large_precision_test.f90 -o b.\ large_precision_test

program large_precision
    implicit none
  
    real :: pi_single 
    real :: one_single
    real :: four_single
    real(kind=8) :: pi_double 
    real(kind=8) :: one_double 
    real(kind=8) :: four_double 
    real(kind=16) :: pi_quad 
    real(kind=16) :: one_quad 
    real(kind=16) :: four_quad 
    
    one_single = 1.0
    four_single = 4.0
    pi_single = four_single * atan(one_single)
    one_double = 1.0
    four_double = 4.0
    pi_double = four_double * atan(one_double)
    one_quad = 1.0
    four_quad = 4.0
    pi_quad = four_quad * atan(one_quad)


    ! Display pi with different precisions
    print *, "Single precision pi (32-bit):    ", pi_single
    print *, "Double precision pi (64-bit):    ", pi_double
    print *, "Quadruple precision pi (128-bit):", pi_quad

end program large_precision


  