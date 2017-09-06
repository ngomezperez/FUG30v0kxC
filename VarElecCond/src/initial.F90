!-------------------------------------
! This module handles the input to the 
! code (i.e. initial conditions and 
! parameter values)
!-------------------------------------
MODULE initial
  use types
  use comms
  use params
  use trans_params
  use dimens
  use arrays
  use fft
  use control
  use constants
  use netcdf 

  implicit none

CONTAINS

!================================
! Sets up static polytropic layer
!================================
SUBROUTINE static()
  implicit none

!-------------------
! Local declarations
!-------------------
  real(kind=dp), parameter  :: amp = 0.1d0
  integer                   :: i, j, k
  real(kind=dp)             :: z, t0, r0
  integer                   :: id

!---------------
! Sets variables
!---------------
  do k=1,nlayer
    z = dble(k+izstart-2)*dz
    t0 = 1.0d0 + z*theta
    r0 = t0**cm

!----------------------------------------------------------
! Each layer is randomized separately using the layer index
! as the seed. This way the same initial conditions can be
! set up independent of the number of processors
!
! At this stage, t contains the random noise only
!----------------------------------------------------------
    id = k + izstart - 1
    call ranset(id, t(1,1,k), nx*ny, amp )

    do j=1,ny
      do i=1,nx
        u(i,j,k)  = 0.0d0
        v(i,j,k)  = 0.0d0
        w(i,j,k)  = 0.0d0
        r(i,j,k) = r0
        t(i,j,k) = t(i,j,k) + t0
      end do
    end do
    if (perfect) then
      do j=1,ny
        do i=1,nx
          bx(i,j,k) = 1.0d0
          by(i,j,k) = 0.0d0
          bz(i,j,k) = 0.0d0
        end do
      end do 
    else
      do j=1,ny
        do i=1,nx
          bx(i,j,k) = 0.0d0
          by(i,j,k) = 0.0d0
          bz(i,j,k) = 1.0d0
        end do
      end do 
     
    end if

  end do



END SUBROUTINE static

!==========================================================
! Reads initial conditions from restart.nc or INPUT
! nframe = 0 --> clean start
! nframe = 1 --> restart from file
! nframe = 2 --> read data from file, parameters from INPUT
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

!-------------------
! Local declarations
!-------------------
  character(5)               ::  input_file
  integer                    ::  flag               

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

!--------------------
! Set data file names
!--------------------
  input_file = 'INPUT'
  restart_file = 'restart.nc'

!-------------------
! Work on myrank = 0
!-------------------
  if (myrank .eq. 0) then

!-----------------------------------
! Read in parameters from INPUT file
!-----------------------------------

!-----------------------------
! First we need the input flag
!-----------------------------
    OPEN(3,FILE=input_file,err=500)
    READ(3,*)
    READ(3,*)
    READ(3,*) flag
    READ(3,*)
    READ(3,*)

!------------------------------------
! Now read the other input parameters
!------------------------------------
    READ(3,*) GAMMA
    READ(3,*) SIGMA
    READ(3,*) THETA
    READ(3,*) CK
    READ(3,*) CM
    READ(3,*) F
    READ(3,*) TAYL
    READ(3,*) PSI
    READ(3,*) ZETA
    READ(3,*) XMAX
    READ(3,*) YMAX
    READ(3,*) SF
    READ(3,*) NFRAME
    READ(3,*) NTOTAL
    READ(3,*) NDU06
    READ(3,*) NDU16
    READ(3,*) ndurestart
    READ(3,*) nx
    READ(3,*) ny
    READ(3,*) nz 
    CLOSE(3)

    if (flag .eq. 1) then

!-----------------
! Work out R and Q
!-----------------
      chand = F/(ZETA*SIGMA*CK*CK)
      rayl = (((CM+1.0d0)*THETA*THETA)/(SIGMA*CK*CK*GAMMA)) &
      *(CM+1.0d0-( CM*GAMMA ))                              &
      *(1.0d0+(0.5d0*THETA))**((2.0d0*CM)-1.0d0)
    else if (flag .eq. 2) then
    
!-----------------------------------------
! For flag =2 R and Q are input parameters
! We need to relabel these before working 
! out ck and f
!-----------------------------------------
      rayl = ck
      chand = F
      ck = dsqrt((((CM+1.0d0)*THETA*THETA)/(SIGMA*rayl*GAMMA)) &
      *(CM+1.0d0-( CM*GAMMA )) &
      *(1.0d0+(0.5d0*THETA))**((2.0d0*CM)-1.0d0))
      F = chand*(ZETA*SIGMA*ck*ck)

    else 
      print*,'Unknown input flag ',flag
      ier = -4
      return
    end if  

!-----------------------
! Sanity checks on nframe
!-----------------------
  if (nframe .lt. 0) then
    print*,'Unknown value for nframe ',nframe
    ier = -2
    return
  end if

  if (nframe .gt. 2) then
    print*,'Unknown value for nframe ',nframe
    ier = -4
    return
  end if

!-----------------------------------------------------
! Read in parameters from data file if nframe = 1 or 2        
!-----------------------------------------------------
    if (nframe .gt. 0) then
      call OpenNetCDF( iter, ttime, dt1, ier ) 
      print*,'RESTART ttime ', ttime, ' iter ', iter
      print*,'Stored time-step ', dt1

!----------------
! Check for error
!----------------
      if ( ier .lt. 0 )  return
    end if 

!-------------------------------------------------------
! Read in parameters from INPUT file again if nframe = 2
!-------------------------------------------------------
    if (nframe .eq. 2) then
 
!-----------------------------
! First we need the input flag
!-----------------------------
      OPEN(3,FILE=input_file,err=500)
      READ(3,*)
      READ(3,*)
      READ(3,*) flag
      READ(3,*)
      READ(3,*)

!------------------------------------
! Now read the other input parameters
!------------------------------------
      READ(3,*) GAMMA
      READ(3,*) SIGMA
      READ(3,*) THETA
      READ(3,*) CK
      READ(3,*) CM
      READ(3,*) F
      READ(3,*) TAYL
      READ(3,*) PSI
      READ(3,*) ZETA
      READ(3,*) XMAX
      READ(3,*) YMAX
      READ(3,*) SF
      READ(3,*) NFRAME
      READ(3,*) NTOTAL
      READ(3,*) NDU06
      READ(3,*) NDU16
      READ(3,*) ndurestart
      READ(3,*) nx
      READ(3,*) ny
      READ(3,*) nz 
      CLOSE(3)

      if (flag .eq. 1) then

!-----------------
! Work out R and Q
!-----------------
        chand = F/(ZETA*SIGMA*CK*CK)
        rayl = (((CM+1.0d0)*THETA*THETA)/(SIGMA*CK*CK*GAMMA)) &
        *(CM+1.0d0-( CM*GAMMA ))                              &
        *(1.0d0+(0.5d0*THETA))**((2.0d0*CM)-1.0d0)
      else if (flag .eq. 2) then
    
!-----------------------------------------
! For flag =2 R and Q are input parameters
! We need to relabel these before working 
! out ck and f
!-----------------------------------------
        rayl = ck
        chand = F
        ck = dsqrt((((CM+1.0d0)*THETA*THETA)/(SIGMA*rayl*GAMMA)) &
        *(CM+1.0d0-( CM*GAMMA )) &
        *(1.0d0+(0.5d0*THETA))**((2.0d0*CM)-1.0d0))
        F = chand*(ZETA*SIGMA*ck*ck)

      else 
        print*,'Unknown input flag ',flag
        ier = -4
        return
      end if  
    end if 

!------------------------------------------------
! Rescale parameters for dynamo problem if needed
!------------------------------------------------
    if (dynamo) then
      f = 1.0d0
      chand = 0.0d0
      print*,'Rescaling F for the dynamo problem : F = ',F
    end if

!---------------------------------------------------
! Rescale parameters for kinematic problem if needed
!---------------------------------------------------
    if (kinematic) then 
      f = 0.0d0
      chand = 0.0d0
      print*,'Setting F=0 for the kinematic problem'
    end if

!-----------------------
! End work on myrank = 0
!-----------------------
  end if
!-------------------------------
! Now do communication if needed
!-------------------------------
#ifdef MPI
  if (myrank .eq. 0) then
     isendbuf(1) = nframe
     isendbuf(2) = ntotal
     isendbuf(3) = ndu06
     isendbuf(4) = ndu16
     isendbuf(5) = nx
     isendbuf(6) = ny
     isendbuf(7) = nz
     isendbuf(8) = iter
     isendbuf(9) = ndurestart
     isendbuf(10) = ncid
     isendbuf(11) = ier

     sendbuf(1)  = gamma
     sendbuf(2)  = sigma
     sendbuf(3)  = theta
     sendbuf(4)  = ck
     sendbuf(5)  = cm
     sendbuf(6)  = f
     sendbuf(7)  = zeta
     sendbuf(8)  = xmax
     sendbuf(9)  = ymax
     sendbuf(10) = sf
     sendbuf(11) = ttime
     sendbuf(12) = dt1
     sendbuf(13) = tayl
     sendbuf(14) = psi
  end if

  ncount =  14
  call MPI_BCAST( sendbuf, ncount, MPI_DOUBLE_PRECISION,  &
                  0, comm, ier )
  ncount =  11
  call MPI_BCAST( isendbuf, ncount, MPI_INTEGER,          &
                  0, comm, ier )

  if (myrank .ne. 0) then
      nframe = isendbuf(1)
      ntotal = isendbuf(2)
      ndu06  = isendbuf(3)
      ndu16  = isendbuf(4)
      nx     = isendbuf(5)
      ny     = isendbuf(6)
      nz     = isendbuf(7)
      iter   = isendbuf(8)
      ndurestart = isendbuf(9)
      ncid   = isendbuf(10)      
      ier    = isendbuf(11)

      gamma = sendbuf(1)
      sigma = sendbuf(2)
      theta = sendbuf(3)
      ck    = sendbuf(4)
      cm    = sendbuf(5)
      f     = sendbuf(6)
      zeta  = sendbuf(7)
      xmax  = sendbuf(8)
      ymax  = sendbuf(9)
      sf    = sendbuf(10)
      ttime  = sendbuf(11)
      dt1    = sendbuf(12)
      tayl  = sendbuf(13)
      psi   = sendbuf(14)
   end if
#endif

   return 

500 continue
   print*,'Unable to open input file ', input_file
   ier = -1 
   return

END SUBROUTINE readparam

END MODULE initial
