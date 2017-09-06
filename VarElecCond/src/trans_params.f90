!====================================
! Defines program transport parameters
!===================================
MODULE trans_params
  use types
  use comms
  use dimens
  use params
  use constants
  use arrays

CONTAINS


SUBROUTINE ini_trans_params() 

  implicit none


!-------------------
! Local declarations
!-------------------
  real(kind=dp), parameter  :: amp = 0.1d0
  integer                   :: i, j, k
  real(kind=dp)             :: z, t0, r0,eps
  real(kind=dp)             :: constA,constB 
  integer                   :: id

 ! real(kind=dp)		    :: kappaZ(1:nlayer)
 ! real(kind=dp)		    :: dkappadZ(1:nlayer)
   
 
!	eps=1.0d-10
!	constA=(1.0d0-(exp(1.0d0)*eps))/(1-exp(1.0d0))
!	constB=eps-constA
	constA=10.0d0
! This part was set up to follow a polytropic EOS
! Test
!    z=dble(izstart-2)*dz
    do k=1,nlayer
	z=dble(k+izstart-2)*dz
	t0=1.0d0+z*theta
	r0=t0**cm
!Test
!    print*, 'Inside: ', dble(k+izstart-2)*dz, z

!=====================================
!definition of kappa(z)/k(z=0)
!lets try a linear change: 
!====================================
!	kappaZ(k)=z
!	dkappadZ(k)=1.0d0		
!=====================================
!definition of kappa(z)/k(z=0)
!lets try an exponential change (avoid kappa =0 and use kappa=eps) 
	
!	kappaZ(k)=constA+constB*exp(z)
!	dkappadZ(k)=constB*exp(z)
	
	constB=constA*(z-5.0d-1)
	kappaZ(k)=5.0d-1*(tanh(constB)+1)
	dkappadZ(k)=5.0d-1*constA*cosh(constB)**(-2.0d0)

!	print*, 'z,T,K',z,t0,kappaZ(k)


    end do

END SUBROUTINE ini_trans_params

END MODULE trans_params
