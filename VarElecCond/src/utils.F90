!==============================================
! This module contains various useful utilities
!==============================================
MODULE utils
  use types
  use dimens
  use params
  use trans_params
  use constants
  use arrays
  use comms
  use fft
  use dbydz

  implicit none

CONTAINS

!=======================================================
! This subroutine calculates various diagnostic averages
!=======================================================
subroutine Globalaverages( ubar, vbar, kenergy, menergy)
  implicit none

!-------------------------
! Declare output variables
!-------------------------
  real(kind=dp), intent(out)   :: ubar
  real(kind=dp), intent(out)   :: vbar
  real(kind=dp), intent(out)   :: kenergy
  real(kind=dp), intent(out)   :: menergy

!-------------------
! Local declarations
!-------------------
  integer               :: iz
  integer               :: i, j, k
  real(kind=dp)         :: urav, vrav, rav, keav, bsqav
  real(kind=dp)         :: rbar, urbar, vrbar
  real(kind=dp)         :: rr, uu, vv, ww 
#ifdef MPI
  integer               :: ier
#endif

!--------------------   
! Zero all quantities
!--------------------
  ubar = 0.0d0
  vbar = 0.0d0
  rbar = 0.0d0
  urbar = 0.0d0
  vrbar = 0.0d0
  menergy = 0.0d0
  kenergy = 0.0d0
  urav = 0.0d0
  vrav = 0.0d0
  rav = 0.0d0
  keav = 0.0d0
  bsqav = 0.0d0
   
!----------------------------------------------
! Work out layer averages within each processor
! Integrals done using the trapezium rule
!----------------------------------------------
  do k=1,nlayer
    urav = 0.0d0
    vrav = 0.0d0
    rav = 0.0d0
    keav= 0.0d0
    bsqav = 0.0d0
    do j=1,ny
      do i=1,nx  
        uu = u(i,j,k)
        vv = v(i,j,k) 
        rr = r(i,j,k)
        ww = w(i,j,k)
        rav = rav + rr
        urav = urav + rr*uu
        vrav = vrav + rr*vv
        keav = keav + 0.5d0*rr*(uu*uu + vv*vv + ww*ww)
        bsqav = bsqav + bx(i,j,k)*bx(i,j,k) + by(i,j,k)*by(i,j,k) + &
                        bz(i,j,k)*bz(i,j,k)
      end do
    end do
    iz = izstart + k - 1
    if ( (iz .eq. 1).or.( iz .eq. nz )) then
      urbar = urbar + 0.5d0*urav 
      vrbar = vrbar + 0.5d0*vrav 
      rbar  = rbar + 0.5d0*rav 
      menergy = menergy + 0.25d0*bsqav
      kenergy = kenergy + 0.5d0*keav
    else 
      urbar = urbar + urav 
      vrbar = vrbar + vrav 
      rbar  = rbar  + rav 
      menergy = menergy + 0.5d0*bsqav
      kenergy = kenergy + keav
    end if
  end do

!----------------------------------
! Now do communication if necessary
!----------------------------------
#ifdef MPI
  sendbuf(1) = urbar
  sendbuf(2) = vrbar
  sendbuf(3) = rbar
  sendbuf(4) = kenergy
  sendbuf(5) = menergy
  call MPI_ALLREDUCE( sendbuf, recvbuf, 5, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, comm, ier )
  urbar = recvbuf(1)
  vrbar = recvbuf(2)
  rbar  = recvbuf(3)
  kenergy = recvbuf(4)
  menergy = recvbuf(5)
#endif

!-----------------------------------------
! Work out CoM velocity and total energies
!-----------------------------------------  
  ubar = urbar/rbar
  vbar = vrbar/rbar
  menergy = menergy*xzn
  kenergy = kenergy*xzn

!-------------------------------------
! This loop adjusts the centre of mass
!-------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx  
        u(i,j,k) = u(i,j,k) - ubar
        v(i,j,k) = v(i,j,k) - vbar
      end do
    end do
  end do

END SUBROUTINE GlobalAverages

!=========================================
! This subroutine calculates the time-step
!=========================================
SUBROUTINE calcdt( dt, outflag)
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=8), intent(out)     :: dt
  integer,      intent(in)      :: outflag

!-------------------
! Local declarations
!-------------------
  integer      :: i, j, k
  real(kind=dp) :: u0, v0, w0, r0, t0, bz0, bx0, by0, bsq0
  real(kind=dp) :: xx, xz
  real(kind=dp) :: velmxx
  real(kind=dp) :: velmxz
  real(kind=dp) :: denmn
  real(kind=dp) :: dtcflx
  real(kind=dp) :: dtcflz
  real(kind=dp) :: dtviscx
  real(kind=dp) :: dtviscz
  real(kind=dp) :: dtdiffx
  real(kind=dp) :: dtdiffz
#ifdef MPI
  integer       :: ier
#endif

!-------------------------
! Pick some initial values
!-------------------------
  denmn =  r(1,1,1)
  velmxx = -1.0d-50
  velmxz = -1.0d-50

!----------------------------------------------------
! Find layer maxima for velocities and density minima
!----------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
         u0 = u(i,j,k)
         v0 = v(i,j,k)
         w0 = w(i,j,k)
         r0 = r(i,j,k)
         t0 = t(i,j,k)
         bz0= bz(i,j,k)
         bx0= bx(i,j,k)
         by0= by(i,j,k)
         bsq0= bx0*bx0 + by0*by0 + bz0*bz0
         xx = sqrt( u0*u0 + v0*v0 ) + g1*sqrt(t0) + sqrt((bsq0*f)/r0)
         xz = sqrt( w0*w0 ) + g1*sqrt(t0) + sqrt((bsq0*f)/r0)
         denmn  = min(denmn,r0)
         velmxx = max(velmxx, xx )
         velmxz = max(velmxz, xz )
      end do
    end do
  end do

!-----------------------------------------------------
! Do communication if necessary to find global extrema
!-----------------------------------------------------
#ifdef MPI
  sendbuf(1) = -denmn
  sendbuf(2) = velmxx
  sendbuf(3) = velmxz
  call MPI_ALLREDUCE( sendbuf, recvbuf, 3, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, comm, ier )

  denmn  = - recvbuf(1)
  velmxx =   recvbuf(2)
  velmxz =   recvbuf(3)
#endif

!----------------------------
! Work out critical timesteps
!----------------------------
  dtcflx = 0.16d0*dd/velmxx
  dtcflz = 0.3d0*dz/velmxz
  dtviscx = dtvisc0x*denmn
  dtviscz = dtvisc0z*denmn
  dtdiffx = dtdiff0x*denmn
  dtdiffz = dtdiff0z*denmn

!----------------------------
! Write some output if needed
!----------------------------
  if ( outflag .eq. 0 ) then
    if ( myrank .eq. 0 ) then
      print*,' velmxx= ', velmxx, ' velmxz= ', velmxz
      print*,' dtcflx= ', dtcflx, ' dtcflz= ', dtcflz
      print*,' dtdiffx= ', dtdiffx, ' dtdiffz= ', dtdiffz
    end if
  end if

!-------------------------------------------------------
! Find minimal value of dt and multiply by safety factor
!-------------------------------------------------------  
  dt = sf*min( dtcflx, dtcflz, dtviscx, dtviscz, dtdiffx, dtdiffz)

END SUBROUTINE calcdt

!=============================
! This subroutine checks div B
!=============================
SUBROUTINE checkdivb( divbmax, divbmin )
  implicit none
  
!------------------------------
! Define input/output variables  
!------------------------------
  real(kind=dp), intent(out) :: divbmax
  real(kind=dp), intent(out) :: divbmin

!-------------------
! Local declarations
!-------------------
  integer       :: i, j, k
  real(kind=dp) :: t1
  real(kind=dp) :: tlmin, tlmax 
#ifdef MPI
  integer       :: ier
#endif

!----------------------------------------
! returns components of div b (un-summed)
!----------------------------------------
  call d1bydx( bx, wk1 )
  call d1bydy( by, wk2 )
  call d1bydz( bz, wk3, 5 )

!--------------------------------------
! Set dummy values of output parameters
!--------------------------------------
  divbmax = -1.0d30
  divbmin =  1.0d00

!-------------------
! Search for max/min
!-------------------
  do k=1,nlayer
    tlmax =-1.0d30
    tlmin = 1.0d30
    do j=1,ny
      do i=1,nx  
        t1 = wk1(i,j,k) + wk2(i,j,k) + wk3(i,j,k)
        tlmax = max(tlmax,t1)
        tlmin = min(tlmin,t1)
      end do
    end do
    divbmax = max(divbmax,tlmax)
    divbmin = min(divbmin,tlmin)
  end do

!---------------------------------------
! Now communicate to find global extrema
!---------------------------------------
#ifdef MPI

  sendbuf(1) = -divbmin
  sendbuf(2) =  divbmax

!--------------------------------------------------------------
! this finds the maximum values of sendbuf across all processes
! ALL implies that the result appears in all receive buffers
!--------------------------------------------------------------
  call MPI_ALLREDUCE( sendbuf, recvbuf, 2, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, comm, ier )

  divbmin = -recvbuf(1)
  divbmax =  recvbuf(2)

#endif

END SUBROUTINE checkdivb

!========================================================
! This subroutine finds global max/min for given function
!========================================================
SUBROUTINE sminmax(z, zmin, zmax )
  implicit none

!------------------------------
! Define input/output variables
!------------------------------
  real(kind=dp), intent(in)    :: z(nx, ny, nlayer )
  real(kind=dp), intent(out)   :: zmin
  real(kind=dp), intent(out)   :: zmax

!-------------------
! Local declarations
!-------------------
  integer               :: i, j, k
#ifdef MPI
  integer               :: ier
#endif

!---------------------------------
! Sets dummy values for zmin, zmax
!---------------------------------
  zmin = z(1,1,1)
  zmax = z(1,1,1)

!--------------------
! Finds layer extrema
!--------------------
  do k=1,nlayer
     do j=1,ny
       do i=1,nx  
         zmin = min( zmin, z(i,j,k))
         zmax = max( zmax, z(i,j,k))
       end do
     end do
  end do

!-------------------------------------------------
! Communicates if necessary to find global extrema
!-------------------------------------------------
#ifdef MPI
  sendbuf(1) =  zmin
  sendbuf(2) = -zmax
  call MPI_ALLREDUCE( sendbuf, recvbuf, 2, MPI_DOUBLE_PRECISION, &
                      MPI_MIN, comm, ier )
  zmin =  recvbuf(1)
  zmax = -recvbuf(2)
#endif

END SUBROUTINE sminmax

END MODULE utils
