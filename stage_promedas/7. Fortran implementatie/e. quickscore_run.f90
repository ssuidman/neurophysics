! gfortran e.\ quickscore_run.f90 -o e.\ quickscore_run

program QuickScoreFortran
  implicit none
  
  ! integer, parameter :: n=773, m=7
  ! integer, parameter :: n=870, m=9
  character(len=1):: case_nr
  character(len=100) :: output_path
  integer :: m, n, i, j, k, file_previn, rc_previn, file_pfmin, rc_pfmin, file_pfminneg, rc_pfminneg, file_output, rc_output
  integer, allocatable :: myset(:)
  real :: t1, t2, dt
  real(kind=16), allocatable :: previn(:), pfmin(:,:), pfminneg(:), prev(:,:), prevminneg(:,:), pfplus(:)
  real(kind=16), allocatable :: pfmin_matrix(:,:,:), term(:), P_joint(:), posterior(:)
  integer, allocatable :: v(:)

  call get_command_argument(1, value=case_nr)

  print *, "Case:",case_nr
  if (case_nr == '1') then
    n = 773
    m = 7 
    open (action="read", file="../variables/cpp/patient404_case_1_previn.csv", iostat=rc_previn, unit=file_previn)
    open (action="read", file="../variables/cpp/patient404_case_1_pfmin.csv", iostat=rc_pfmin, unit=file_pfmin)
    open (action="read", file="../variables/cpp/patient404_case_1_pfminneg.csv", iostat=rc_pfminneg, unit=file_pfminneg)
  else if (case_nr =='2') then
    n = 809
    m = 8
    open (action="read", file="../variables/cpp/patient404_case_2_previn.csv", iostat=rc_previn, unit=file_previn)
    open (action="read", file="../variables/cpp/patient404_case_2_pfmin.csv", iostat=rc_pfmin, unit=file_pfmin)
    open (action="read", file="../variables/cpp/patient404_case_2_pfminneg.csv", iostat=rc_pfminneg, unit=file_pfminneg)
  else if (case_nr =='3') then
    n = 870
    m = 9
    open (action="read", file="../variables/cpp/patient404_case_3_previn.csv", iostat=rc_previn, unit=file_previn)
    open (action="read", file="../variables/cpp/patient404_case_3_pfmin.csv", iostat=rc_pfmin, unit=file_pfmin)
    open (action="read", file="../variables/cpp/patient404_case_3_pfminneg.csv", iostat=rc_pfminneg, unit=file_pfminneg)
  else 
    print *, "No case_nr defined!!!"
  end if
  allocate(previn(n), pfmin(m, n), pfminneg(n), prev(n+1, n), prevminneg(n+1, n), pfplus(n+1))
  allocate(pfmin_matrix(m, n+1, n), term(n+1), P_joint(n), posterior(n), v(m))

  do i=1, n
    read(file_previn,*,iostat=rc_previn) previn(i)
    if (rc_previn /= 0) exit 
  end do 
  do j = 1, m
    read(file_pfmin, *, iostat=rc_pfmin) pfmin(j,:)
    if (rc_pfmin /= 0) exit 
  end do 
  do i=1, n
    read(file_pfminneg,*,iostat=rc_pfminneg) pfminneg(i)
    if (rc_pfminneg /= 0) exit 
  end do
  close (file_previn)
  close (file_pfmin)
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

  write(output_path, "('../variables/fortran_patient404/output_case_', A, '.csv')") trim(case_nr)
  open(newunit=file_output, file=output_path, action="write", status="replace", iostat=rc_output)
  do i = 1, n
    if (i==1) then
      write(file_output, advance='NO') posterior(i), P_joint(i), dt
    else
      write(file_output, *) posterior(i), P_joint(i)
    end if
  end do
  print *, "output_path:", output_path
  close(file_output)

end program QuickScoreFortran




