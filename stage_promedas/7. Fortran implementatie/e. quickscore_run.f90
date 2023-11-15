! gfortran e.\ quickscore_run.f90 -o e.\ quickscore_run

program QuickScoreFortran
  implicit none

  integer, parameter :: n=870, m=9
  integer :: i, j, k, file_previn, rc_previn, file_pfmin, rc_pfmin, file_pfminneg, rc_pfminneg, v(m)
  integer, allocatable :: myset(:)
  real :: t1, t2, dt
  real(kind=16) :: previn(n), pfmin(m,n), pfminneg(n), prev(n+1,n), prevminneg(n+1,n), pfplus(n+1)
  real(kind=16) :: pfmin_matrix(m,n+1,n), term(n+1), P_joint(n), posterior(n)
  

  open (action="read", file="../variables/cpp/patient404_case_3_previn.csv", iostat=rc_previn, unit=file_previn)
  ! if (rc_previn /= 0) stop 'Error: open failed'
  do i=1, n
    read(file_previn,*,iostat=rc_previn) previn(i)
    if (rc_previn /= 0) exit 
  end do 
  close (file_previn)

  open (action="read", file="../variables/cpp/patient404_case_3_pfmin.csv", iostat=rc_pfmin, unit=file_pfmin)
  ! if (rc_pfmin /= 0) stop 'Error: open failed'
  do j = 1, m
    read(file_pfmin, *, iostat=rc_pfmin) pfmin(j,:)
    if (rc_pfmin /= 0) exit 
  end do 
  close (file_pfmin)

  open (action="read", file="../variables/cpp/patient404_case_3_pfminneg.csv", iostat=rc_pfminneg, unit=file_pfminneg)
  ! if (rc_pfminneg /= 0) stop 'Error: open failed'
  do i=1, n
    read(file_pfminneg,*,iostat=rc_pfminneg) pfminneg(i)
    if (rc_pfminneg /= 0) exit 
  end do 
  close (file_pfminneg)

  do i = 1, n+1
    do j = 1, n
      if (i==j+1) then 
        prev(i,j) = 1.0
      else
        prev(i, j) = previn(j)
      end if
    end do
  end do
  
  do i=1, n+1
    prevminneg(i,:) = prev(i,:) * pfminneg(:)
    pfmin_matrix(:,i,:) = pfmin
  end do



  pfplus = 0
  call cpu_time(t1)
  do k=1,2**m
    do j = 1, m
      v(j) = merge(1, 0, iand(k, 2**(j-1)) /= 0)
    end do
    myset = pack([(i,i=1,m)]*v, [(i,i=1,m)]*v /= 0)
    term = (-1)**size(myset)*product(1e-50 + product(pfmin_matrix(myset,:,:),dim=1) * prevminneg + (1-prev),dim=2)
    pfplus = pfplus + term
  end do
  call cpu_time(t2)
  dt = t2 - t1
  
  P_joint = pfplus(2:size(pfplus)) * previn
  posterior = P_joint / pfplus(1)
  
  print *
  do i=1, 5
    print *, "previn:", previn(i)
  end do
  print *
  do i=1, 7
    print *, "pfmin: ", pfmin(i,1:3)
  end do
  print *
  do i=1, 5
    print *, "pfminneg:", pfminneg(i)
  end do
  print *
  print *, "P_joint:        size:",size(P_joint)
  do i=1, 5
    print *, "P_joint: ", P_joint(i)
  end do
  print *
  print *, "posterior:        size:",size(posterior)
  do i=1, 5
    print *, "posterior:", posterior(i)
  end do
  print *
  print *, "Running time:", dt
  print *

end program QuickScoreFortran





!   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
!   REAL(dp), DIMENSION(:,:), ALLOCATABLE :: previn, pfmin, pfminneg, prev, pfplus, P_joint, posterior
!   INTEGER :: m, n, i, j, set_size, num_sets
!   REAL(dp) :: dt
!   CHARACTER(20) :: method


!   ! Set parameters
!   m = SIZE(pfmin, 1)
!   n = SIZE(previn, 2)
!   ALLOCATE(prev(n+1, n), pfplus(n+1, 1), P_joint(n), posterior(n))

!   PRINT *, "m = ", m, " (positive tests)"

!   ! Copy prevalences
!   prev(:, :) = TRANSPOSE(MATMUL(previn, 1.0_dp))

!   ! Set entries for condition on d_i=1
!   DO i = 1, n
!   prev(i+1, i) = 1.0_dp
!   END DO

!   ! Absorb negative findings in prevminneg
!   IF (SIZE(pfminneg, 2) /= 0) THEN
!   prev(:, :) = prev(:, :) * TRANSPOSE(pfminneg)
!   ELSE
!   prev(:, :) = prev(:, :) * 1.0_dp
!   END IF

!   ! Initialize pfplus
!   pfplus(:, 1) = 0.0_dp

!   ! Iterate over 2^m possibilities
!   num_sets = 2**m
!   dt = 0.0_dp
!   method = "exp-sum-log"

!   DO i = 0, num_sets - 1
!   CALL Progress(i, num_sets)

!   ! Create vector of 0's and 1's
!   CALL CreateBinaryVector(i, m, set_size)

!   ! Calculate loop term and update pfplus
!   CALL LoopTerm(set_size, pfmin, prev, pfplus, method)
!   END DO

!   PRINT *, "Running time: ", dt

!   ! Calculate P_joint and posterior
!   P_joint = pfplus(2:, 1) * previn
!   posterior = P_joint / pfplus(1)

!   ! Deallocate memory
!   DEALLOCATE(previn, pfmin, pfminneg, prev, pfplus, P_joint, posterior)

!   CONTAINS

!   SUBROUTINE ReadCSV(filename, data)
!   CHARACTER(*), INTENT(IN) :: filename
!   REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: data
!   INTEGER :: i, j, num_rows, num_cols
!   CHARACTER(200) :: line
!   LOGICAL :: header

!   ! Count the number of rows and columns
!   OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=i)
!   IF (i /= 0) THEN
!     PRINT *, "Error opening file: ", filename
!     STOP
!   END IF

!   header = .TRUE.
!   num_rows = 0
!   DO
!     READ(10, '(A)', IOSTAT=i) line
!     IF (i /= 0) EXIT
!     IF (.NOT. header) THEN
!       num_rows = num_rows + 1
!     END IF
!   END DO

!   num_cols = COUNT(line, DELIM=',') + 1

!   CLOSE(10)

!   ! Allocate memory for the data array
!   ALLOCATE(data(num_rows, num_cols))

!   ! Read data from the CSV file
!   OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=i)
!   IF (i /= 0) THEN
!     PRINT *, "Error opening file: ", filename
!     STOP
!   END IF

!   DO i = 1, num_rows
!     READ(10, *, IOSTAT=i) data(i, :)
!   END DO

!   CLOSE(10)
!   END SUBROUTINE ReadCSV

!   SUBROUTINE CreateBinaryVector(decimal, size, binary_vector)
!   INTEGER, INTENT(IN) :: decimal, size
!   INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: binary_vector
!   INTEGER :: i

!   ALLOCATE(binary_vector(size))

!   DO i = 1, size
!     binary_vector(i) = MOD(decimal, 2)
!     decimal = decimal / 2
!   END DO
!   END SUBROUTINE CreateBinaryVector

!   SUBROUTINE Progress(iteration, total)
!   INTEGER, INTENT(IN) :: iteration, total

!   ! Implement a progress bar here based on iteration and total
!   ! This may vary depending on your Fortran environment
!   PRINT *, "Progress: ", 100 * REAL(iteration) / REAL(total), "%"
!   END SUBROUTINE Progress

!   SUBROUTINE LoopTerm(myset, pfmin, prev, pfplus, method)
!   INTEGER, DIMENSION(:), INTENT(IN) :: myset
!   REAL(dp), DIMENSION(:,:), INTENT(IN) :: pfmin, prev
!   REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: pfplus
!   CHARACTER(20), INTENT(IN) :: method
!   INTEGER :: i

!   IF (TRIM(method) == "exp-sum-log") THEN
!     pfplus(:, 1) = pfplus(:, 1) + ((-1)**SIZE(myset)) * EXP(SUM(LOG(1.0_dp + (PRODUCT(pfmin(myset, :), DIM=1) * prev * (1.0_dp - prev))), DIM=2))
!   ELSEIF (TRIM(method) == "prod") THEN
!     pfplus(:, 1) = pfplus(:, 1) + ((-1)**SIZE(myset)) * PRODUCT(1.0_dp + (PRODUCT(pfmin(myset, :), DIM=1) * prev * (1.0_dp - prev)), DIM=2)
!   END IF
! END SUBROUTINE LoopTerm
  

