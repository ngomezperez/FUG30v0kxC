!-------------------------------------
! This module handles the input to the 
! code (i.e. initial conditions and 
! parameter values)
!-------------------------------------
MODULE initial
  use types
  use comms
  use params
  use dimens
  use arrays
  use fft
  use control
  use constants
  use netcdf 

  implicit none

CONTAINS

!==========================================================
! Reads initial conditions from restart.nc
!==========================================================
SUBROUTINE readparam(iter, ttime, dt1, ier)
  implicit none

!-------------------------
! Declare output variables
!-------------------------               
  integer,       intent(out) ::  iter               
  integer,       intent(out) ::  ier
  real(kind=dp), intent(out) ::  ttime               
  real(kind=dp), intent(out) ::  dt1

!----------------------------------------
! Sets default values of return variables
!----------------------------------------
  ier = 0
  iter = 0
  dt1 = 0.0d0
  ttime = 0.0d0

!-------------------------
! Sets dummy value of ncid
!-------------------------
  ncid = -1

!-------------------
! Set data file name
!-------------------
  restart_file = 'restart.nc'

!-------------------
! Work on myrank = 0
!-------------------
  if (myrank .eq. 0) then

!----------------------------------
! Read in parameters from data file        
!----------------------------------
    call OpenNetCDF( iter, ttime, dt1, ier ) 
    print*,'RESTART ttime ', ttime, ' iter ', iter
    print*,'Stored time-step ', dt1

!----------------
! Check for error
!----------------
    if ( ier .lt. 0 )  return

!-----------------------
! End work on myrank = 0
!-----------------------
  end if

!-------------------------------
! Now do communication if needed
!-------------------------------
#ifdef MPI
  if (myrank .eq. 0) then
     isendbuf(1) = nx
     isendbuf(2) = ny
     isendbuf(3) = nz
     isendbuf(4) = iter
     isendbuf(5) = ncid
     isendbuf(6) = ier

     sendbuf(1)  = gamma
     sendbuf(2)  = sigma
     sendbuf(3)  = theta
     sendbuf(4)  = ck
     sendbuf(5)  = cm
     sendbuf(6)  = f
     sendbuf(7)  = zeta
     sendbuf(8) = xmax
     sendbuf(9) = ymax
     sendbuf(10) = ttime
     sendbuf(11) = dt1
     sendbuf(12) = tayl
     sendbuf(13) = psi
  end if

  ncount =  13
  call MPI_BCAST( sendbuf, ncount, MPI_DOUBLE_PRECISION,  &
                  0, comm, ier )
  ncount =  6
  call MPI_BCAST( isendbuf, ncount, MPI_INTEGER,          &
                  0, comm, ier )

  if (myrank .ne. 0) then
      nx     = isendbuf(1)
      ny     = isendbuf(2)
      nz     = isendbuf(3)
      iter   = isendbuf(4)
      ncid   = isendbuf(5)      
      ier    = isendbuf(6)

      gamma  = sendbuf(1)
      sigma  = sendbuf(2)
      theta  = sendbuf(3)
      ck     = sendbuf(4)
      cm     = sendbuf(5)
      f      = sendbuf(6)
      zeta   = sendbuf(7)
      xmax   = sendbuf(8)
      ymax   = sendbuf(9)
      ttime  = sendbuf(10)
      dt1    = sendbuf(11)
      tayl   = sendbuf(12)
      psi    = sendbuf(13)
   end if
#endif

   return 

END SUBROUTINE readparam

END MODULE initial
