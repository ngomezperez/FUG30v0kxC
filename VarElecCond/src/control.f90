!===========================
! Defines control parameters
!===========================
MODULE control
  implicit none

!-----------------------------------------------------------------
! Time-stepping flags
!-----------------------------------------------------------------
! first: Tells the code to do an Euler step
! second: Tells the code to do a second order Adams-Bashforth step
!-----------------------------------------------------------------
  logical            :: first
  logical            :: second

!----------------------------------------------------------
! Boundary condition flags
!----------------------------------------------------------
! radbc: Radiative upper boundary condition for temperature
! constflux: Constant flux lower boundary condition
! perfect: Perfectly conducting boundaries (condition on B)
!----------------------------------------------------------
  logical            :: radbc
  logical            :: constflux
  logical            :: perfect

!----------------------------------------------------------------
! Magnetic field flags
!----------------------------------------------------------------
! kinematic: Puts the code into kinematic mode
! dynamo: Replaces Bx, By, Bz with a zero-flux boundary condition
!----------------------------------------------------------------
  logical            :: kinematic
  logical            :: dynamo

!----------------------------------
! Other flags
!----------------------------------
! timer: Turns on additional timing
! dalias: Turns on dealiasing
!----------------------------------
  logical            :: timer
  logical            :: dalias

END MODULE control
