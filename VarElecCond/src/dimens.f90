!============================================
! This module defines the dimension variables
!============================================
MODULE dimens
  implicit none

!---------------
! Box dimensions
!---------------
  integer             :: nx
  integer             :: ny
  integer             :: nz

!----------------------------
! Data distribution variables
!----------------------------
  integer             :: izstart
  integer             :: izfinish
  integer             :: nlayer

END MODULE dimens

