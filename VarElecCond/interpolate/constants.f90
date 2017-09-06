!============================================================
! This module holds some useful constants used in the program
!============================================================
MODULE constants
  use types
  use params
  use dimens
 
  implicit none

!------------------
! Set up parameters
!------------------
  real(kind=dp), parameter  :: one  = 1.0d0
  real(kind=dp), parameter  :: third  = 1.0d0/3.0d0
  real(kind=dp), parameter  :: elf  = 1.0d0/11.0d0

!--------------------------
! Define constant variables
!--------------------------
  real(kind=dp)      :: gck
  real(kind=dp)      :: gm1
  real(kind=dp)      :: sck
  real(kind=dp)      :: zck
  real(kind=dp)      :: scko3
  real(kind=dp)      :: sckg
  real(kind=dp)      :: fzckg
  real(kind=dp)      :: thmplus1
  real(kind=dp)      :: nusstmp
  real(kind=dp)      :: fo2
  real(kind=dp)      :: dx
  real(kind=dp)      :: dy
  real(kind=dp)      :: dz
  real(kind=dp)      :: dd
  real(kind=dp)      :: xn
  real(kind=dp)      :: xzn
  real(kind=dp)      :: g1
  real(kind=dp)      :: taylt
  real(kind=dp)      :: cor
  real(kind=dp)      :: cpsi
  real(kind=dp)      :: spsi
 
!--------------------------------
! Time-step determining variables
!--------------------------------
  real(kind=dp)      :: dtvisc0x
  real(kind=dp)      :: dtvisc0z
  real(kind=dp)      :: dtdiff0x
  real(kind=dp)      :: dtdiff0z

CONTAINS

!===================================
! This subroutine sets the constants
!===================================
SUBROUTINE set_constants
  implicit none

!-----------------------
! Mesh-related constants
!-----------------------
  dx       = xmax/dble(nx)
  dy       = ymax/dble(ny)
  dz       = 1.0d0/dble(nz-1)
  dd       = min( dx, dy )
  xn       = 1.0d0/dble(nx*ny)
  xzn      = dx*dy*dz

!-------------------
! Rotation constants
!-------------------
  taylt    =  tayl
  cpsi     =  cos(psi)
  spsi     =  sin(psi)
  cor      =  ck*sigma*sqrt(taylt)

!----------------------------
! Diffusion-related constants
!----------------------------
  gck      =  gamma*ck
  sck      =  sigma*ck
  zck      =  zeta*ck
  scko3    =  sck/3.0d0
  sckg     =  sigma*ck*(gamma-1.0d0)
  fzckg    =  f*zeta*ck*(gamma-1.0d0)

!----------------
! Other constants
!----------------
  thmplus1 =  theta*(cm+1.0d0)
  fo2      =  f/2.0d0
  g1       =  sqrt(gamma)
  gm1      =  gamma - one
  nusstmp  =  theta*( cm+1.0d0 )*( 1.0d0 - (1.0d0/gamma) )

!---------------------------------------------------------------------
! Time-step coefficients are set up here
! 
! Safety factors are different for horizontal and vertical diffusion.
! Horizontal coefficeint is 1/pi^2, vertical is 3/16 (in theory).
! Vertical on can be relaxed in highly compressible case (large theta0
! Hence the ad hoc factors of ( 1.0 + theta/7.0).
!---------------------------------------------------------------------
  dtvisc0x = 0.1d0*0.75d0*dd*dd/sigma/ck
  dtvisc0z = 0.2d0*dz*dz/sigma/ck*(1.0d0+theta/7.0d0)
  dtdiff0x = 0.1d0*dd*dd/gamma/ck
  dtdiff0z = 0.2d0*dz*dz/gamma/ck*(1.0d0+theta/7.0d0)

END SUBROUTINE set_constants

END MODULE constants
