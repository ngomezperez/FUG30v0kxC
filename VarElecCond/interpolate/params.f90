!===========================
! Defines program parameters
!===========================
MODULE params
  use types

  implicit none

  real(kind=dp)      :: gamma
  real(kind=dp)      :: sigma
  real(kind=dp)      :: theta
  real(kind=dp)      :: ck
  real(kind=dp)      :: cm
  real(kind=dp)      :: f
  real(kind=dp)      :: zeta
  real(kind=dp)      :: xmax
  real(kind=dp)      :: ymax
  real(kind=dp)      :: sf
  real(kind=dp)      :: chand
  real(kind=dp)      :: rayl
  real(kind=dp)      :: tayl
  real(kind=dp)      :: psi

  integer            :: nframe
  integer            :: ntotal
  integer            :: ndu06
  integer            :: ndu16
  integer            :: ndurestart

END MODULE params
