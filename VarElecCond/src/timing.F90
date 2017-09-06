!==========================================================
! Timing module for monitoring code performance
!
! NOTE: This timing routine may need altering for different 
! serial architectures
!==========================================================
MODULE timing
  use types
  use control

  implicit none

!------------------------------------
! Declarations for "process" measures
!------------------------------------
  real(kind=dp) :: time_start 
  real(kind=dp) :: time_finish
  real(kind=dp) :: time_dx
  real(kind=dp) :: time_d2x
  real(kind=dp) :: time_dy
  real(kind=dp) :: time_d2y
  real(kind=dp) :: time_dz
  real(kind=dp) :: time_d2z
  real(kind=dp) :: time_dzup
  real(kind=dp) :: time_bound
  real(kind=dp) :: time_dealias
  real(kind=dp) :: time_forward
  real(kind=dp) :: time_reverse
  real(kind=dp) :: time_tmp
  real(kind=dp) :: time_tmp2
  real(kind=dp) :: start_iteration
  real(kind=dp) :: finish_iteration

!---------------------------------
! Declarations for "term" measures
!---------------------------------
  real(kind=dp) :: t1_time, t2_time, t3_time, t4_time, t5_time 
  real(kind=dp) :: t6_time, t7_time, t8_time, t9_time, t10_time 
  real(kind=dp) :: t11_time, t12_time, t13_time, t14_time, t15_time 
  real(kind=dp) :: t16_time, t17_time

CONTAINS

!------------------------------------------
! Subroutine to initialise timing variables
!------------------------------------------
subroutine init_timing()
  implicit none

!---------------------
! Calculate start time
!---------------------
  call findtime(time_start)

!--------------------
! "Process" variables
!--------------------
  time_finish = 0.0d0
  time_dx  = 0.0d0
  time_d2x = 0.0d0
  time_dy  = 0.0d0
  time_d2y = 0.0d0
  time_dz  = 0.0d0
  time_d2z = 0.0d0
  time_dzup = 0.0d0
  time_bound = 0.0d0
  time_dealias = 0.0d0
  time_forward = 0.0d0
  time_reverse = 0.0d0

!-----------------
! "Term" variables
!-----------------
  t1_time  = 0.0d0
  t2_time  = 0.0d0
  t3_time  = 0.0d0
  t4_time  = 0.0d0
  t5_time  = 0.0d0
  t6_time  = 0.0d0
  t7_time  = 0.0d0
  t8_time  = 0.0d0
  t9_time  = 0.0d0
  t10_time  = 0.0d0
  t11_time  = 0.0d0
  t12_time  = 0.0d0
  t13_time  = 0.0d0
  t14_time  = 0.0d0
  t15_time  = 0.0d0
  t16_time  = 0.0d0
  t17_time  = 0.0d0

  time_tmp = 0.0d0
  time_tmp2 = 0.0d0
 
end subroutine init_timing

!-----------------------------------------------------
! 
subroutine print_timing()
  implicit none

!----------------------
! Calculate finish time
!----------------------
  call findtime(time_finish)

!------------------------------------
! Print out "process" terms if needed
!------------------------------------
  if (timer) then
    print*,''
    print*,' Dbydx      = ', time_dx
    print*,' Dbydy      = ', time_dy
    print*,' Dbydz      = ', time_dz
    print*,' D2bydx     = ', time_d2x
    print*,' D2bydy     = ', time_d2y
    print*,' D2bydz     = ', time_d2z
    print*,' Dbydz_up   = ', time_dzup
    print*,' Boundaries = ', time_bound
    print*,' Dealias    = ', time_dealias
    print*,' Foward FFT = ', time_forward
    print*,' Invert FFT = ', time_reverse
    print*,''
  end if
  
!---------------------------
! Print out "term" variables
!---------------------------
  print*,' T1  = ', t1_time
  print*,' T2  = ', t2_time
  print*,' T3  = ', t3_time
  print*,' T4  = ', t4_time
  print*,' T5  = ', t5_time
  print*,' T6  = ', t6_time
  print*,' T7  = ', t7_time
  print*,' T8  = ', t8_time
  print*,' T9  = ', t9_time
  print*,' T10 = ', t10_time
  print*,' T11 = ', t11_time
  print*,' T12 = ', t12_time
  print*,' T13 = ', t13_time
  print*,' T14 = ', t14_time
  print*,' T15 = ', t15_time
  print*,' T16 = ', t16_time
  print*,' T17 = ', t17_time

!-----------------------
! Work out code duration
!-----------------------
  print*,'TOTAL time ', time_finish - time_start

end subroutine print_timing

!------------------------------------------------------
! Timing function - note this is architecture dependent
!------------------------------------------------------
subroutine findtime(currenttime) 
  implicit none

!-----------------------
! Define output variable
!-----------------------
  real(kind=dp), intent(out)  :: currenttime

!-------------------
! Local declarations
!-------------------
  real           :: timetemp
  real           :: tarray(2)
  real           :: etime

#ifdef MPI
  include 'mpif.h'
  currenttime = mpi_wtime()
#else
  timetemp = etime(tarray)
  currenttime = timetemp
#endif
  
  return

end subroutine findtime


END MODULE timing
