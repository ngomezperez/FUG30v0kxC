!=========================================
! Makes default data type double precision
!=========================================
MODULE types
  implicit none

!---------------------------------------
! Defines 'dp' to imply double precision
!---------------------------------------
  integer, parameter :: dp = kind(1.0d0)

END MODULE types
