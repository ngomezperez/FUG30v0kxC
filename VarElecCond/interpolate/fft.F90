!=================================================
! This module contains all subroutines 
! that require the use of FFTs
!-------------------------------------------------
! NOTES: nx must be even (best if a power of two)
!
! For the FFTW routines: T^_-1[T](x)=n*x so must 
! divide by n! 
!
! Check documentation to work out how real-complex 
! routines store transformed data. In brief:
!
! FFTW disposes of complex parts of DC and Nyquist 
! so stores all in an n dim array (no padding). 
! First half real increasing, second
! half complex decreasing. 
!=================================================
MODULE fft
  use types
  use dimens
  use comms
  use timing
  use params
  use control

  implicit none

!------------------------------------------
! Include file for FFTW included explicitly
!------------------------------------------
  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)

!--------------------------------
! Constants for the MPI wrappers:
!--------------------------------
  integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
  integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
  parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
  parameter(FFTW_SCRAMBLED_INPUT=8192)
  parameter(FFTW_SCRAMBLED_OUTPUT=16384)

!--------------------
! Declare wavenumbers
!--------------------
  real(kind=dp), allocatable :: akx1(:)
  real(kind=dp), allocatable :: akx2(:)
  real(kind=dp), allocatable :: aky1(:)
  real(kind=dp), allocatable :: aky2(:)

!------------------------
! Declare local tmp array
!------------------------
  real(kind=dp), allocatable :: tmp(:,:,:)
  real(kind=dp), allocatable :: tmpnew(:,:,:)
  real(kind=dp), allocatable :: tmpnew4(:,:,:)
  real(kind=dp), allocatable :: interp(:,:,:)
  
!-----------------------
! Declare Plan variables
!-----------------------
  integer(kind=8)            :: planfx, planrx, planfy, planry
  integer(kind=8)            :: planfxout, planfyout
  integer(kind=8)            :: planfxnew, planrxnew, planfynew
  integer(kind=8)            :: planrynew, planfxoutnew, planfyoutnew

!--------------------------
! Other global declarations
!--------------------------
  integer                       :: icutoffx
  integer                       :: icutoffy
 

CONTAINS

!================================
! Initialisation routine for FFTs
!================================
SUBROUTINE init_fft(iret)
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer, intent(out)        :: iret

!-------------------
! Local declarations
!-------------------
  integer                     :: ier(8)
  integer                     :: i, nxmax, nymax
  real(kind=dp)               :: cx, cy

!-----------------
! Set iret to zero
!-----------------
  iret=0
  
! --------------------------------------------------------
! An odd value of 'nx' messes up the derivative taking for
! FFTW - check that here
!---------------------------------------------------------
  if ( 2*(nx/2) .ne. nx ) then
    iret=-1
    return
  end if
  if ( 2*(ny/2) .ne. ny ) then
    iret=-2
    return
  end if

!---------------------------------
! Create plans for FFTW transforms
! "inplace" by default
!---------------------------------
  call rfftw_f77_create_plan( planfxout, nx, fftw_real_to_complex, &
                             fftw_estimate )
  call rfftw_f77_create_plan( planfyout, ny, fftw_real_to_complex, &
                             fftw_estimate )
  call rfftw_f77_create_plan( planfx, nx, fftw_real_to_complex, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planfy, ny, fftw_real_to_complex, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planrx, nx, fftw_complex_to_real, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planry, ny, fftw_complex_to_real, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planfxoutnew, 2*nx, fftw_real_to_complex, &
                             fftw_estimate )
  call rfftw_f77_create_plan( planfyoutnew, 2*ny, fftw_real_to_complex, &
                             fftw_estimate )
  call rfftw_f77_create_plan( planfxnew, 2*nx, fftw_real_to_complex, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planfynew, 2*ny, fftw_real_to_complex, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planrxnew, 2*nx, fftw_complex_to_real, &
                             fftw_in_place + fftw_estimate )
  call rfftw_f77_create_plan( planrynew, 2*ny, fftw_complex_to_real, &
                             fftw_in_place + fftw_estimate )

!-----------------------------------------
! Allocate wavenumber and temporary arrays
!-----------------------------------------
  allocate(tmp(nx,ny,nlayer), stat=ier(1))
  allocate(tmpnew(2*nx,2*ny,nlayer), stat=ier(2))
  allocate(tmpnew4(2*nx,ny,nlayer), stat=ier(3))
  allocate(interp(2*nx,ny,nlayer), stat=ier(4))
  allocate(akx1(nx), stat=ier(5))
  allocate(akx2(nx), stat=ier(6))
  allocate(aky1(ny), stat=ier(7))
  allocate(aky2(ny), stat=ier(8))

!-------------
! Check result
!-------------
  if ( (ier(1) .gt. 0) .or.    &
       (ier(2) .gt. 0) .or.    &
       (ier(3) .gt. 0) .or.    &
       (ier(4) .gt. 0) .or.    &
       (ier(5) .gt. 0) .or.    &
       (ier(6) .gt. 0) .or.    &
       (ier(7) .gt. 0) .or.    &
       (ier(8) .gt. 0) )   then
     iret = -5
     return
  end if

!------------------
! Initialise arrays
!------------------
  tmp = 0.0d0
  tmpnew = 0.0d0
  tmpnew4 = 0.0d0
  interp = 0.0d0
  akx1 = 0.0d0
  akx2 = 0.0d0
  aky1 = 0.0d0
  aky2 = 0.0d0 

!-----------------------------
! Calculate dealiasing cutoffs
!-----------------------------
  icutoffx = (nx/2)/3
  icutoffy = (ny/2)/3

!----------------------------------------------------------------
! Now calculate the phase factors for calculating the derivatives
! As the convention f(k) = SUM [ exp(-ik.r) f(r) ] is being used
! the derivative is obtained by multipling by "+ik" where 
! k is actually (2*pi/l)k   --- to see what this factor is it is
! best to go to a continuous distribution, take the derivative
! and then transform back to to the discrete
!
! cx and cy are 2Pi/xmax and 2Pi/ymax resp.
!----------------------------------------------------------------
  cx = 8.0d0*atan(1.0d0)/xmax
  cy = 8.0d0*atan(1.0d0)/ymax

!----------------------------------------------
! nxmax and nymax determine maximal wavenumbers
!----------------------------------------------
  nxmax = (nx/2) - 1
  nymax = (ny/2) - 1

!------------------
! Set x-wavenumbers
!------------------
  do i=1, nxmax
     akx1(i+1) = dble(i)*cx
     akx1(nx+1-i) = dble(i)*cx
     akx2(i+1) = -dble(i)*dble(i)*cx*cx
     akx2(nx+1-i) = -dble(i)*dble(i)*cx*cx
  end do

!-------------------------------------------------------------
! The Nyquist frequency contributes to second derivatives only
!-------------------------------------------------------------
  akx2(nxmax+2) = -dble(nx/2)*dble(nx/2)*cx*cx

!------------------
! Set y-wavenumbers
!------------------
  do i=1, nymax
     aky1(i+1) = dble(i)*cy
     aky1(ny+1-i) = dble(i)*cy
     aky2(i+1) = -dble(i)*dble(i)*cy*cy
     aky2(ny+1-i) = -dble(i)*dble(i)*cy*cy
  end do

!-------------------------------------------------------------
! The Nyquist frequency contributes to second derivatives only
!-------------------------------------------------------------
  aky2(nymax+2) = -dble(ny/2)*dble(ny/2)*cy*cy

END SUBROUTINE init_fft

!=======================================
! Deallocate the memory used for the FFT
!=======================================
SUBROUTINE close_fft
  implicit none

!-------------------------------------------
! Deallocate wavenumber and temporary arrays
!-------------------------------------------
  deallocate(akx1)
  deallocate(akx2)
  deallocate(aky1)
  deallocate(aky2)
  deallocate(tmp)

!-------------------
! Destroy FFTW Plans
!-------------------
  call rfftw_f77_destroy_plan(planfx)
  call rfftw_f77_destroy_plan(planfy)
  call rfftw_f77_destroy_plan(planfxout)
  call rfftw_f77_destroy_plan(planfyout)
  call rfftw_f77_destroy_plan(planrx)
  call rfftw_f77_destroy_plan(planry)
  call rfftw_f77_destroy_plan(planfxnew)
  call rfftw_f77_destroy_plan(planfynew)
  call rfftw_f77_destroy_plan(planfxoutnew)
  call rfftw_f77_destroy_plan(planfyoutnew)
  call rfftw_f77_destroy_plan(planrxnew)
  call rfftw_f77_destroy_plan(planrynew)

END SUBROUTINE close_fft

!===================
! First x-derivative
!===================
subroutine d1bydx( fin, fout )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: fout(nx, ny, nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: i, j, k
  real(kind=dp)                 :: fac

!--------------------
! Set timer if needed
!--------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!----------------
! Initialise fout
!----------------
  fout = 0.0d0

!--------------------------------
! Forward, out of place transform
!--------------------------------
  do k=1, nlayer
    call rfftw_f77(planfxout, ny, fin(1,1,k), 1, nx, fout(1,1,k), 1, nx)
  end do

!-----------------------------------
! Multiply by "+ik" in Fourier space
!-----------------------------------
  do k=1,nlayer
    do j=1,ny
      fout(1,j,k)=0.0d0
      fout((nx/2)+1,j,k)=0.0d0
      do i=1,nx/2-1
        fac = fout(i+1,j,k)
        fout(i+1,j,k)   = -fout(nx+1-i,j,k)*akx1(i+1)
        fout(nx+1-i,j,k) =  fac*akx1(nx+1-i)
      end do
    end do
  end do

!------------------------------------
! Inverse, inplace, Fourier transform
!------------------------------------
   do k=1, nlayer
     call rfftw_f77(planrx, ny, fout(1,1,k), 1, nx, tmp, 1, nx)
     do j=1, ny
       do i=1, nx
         fout(i,j,k)=fout(i,j,k)/dble(nx)
       end do
     end do
   end do

!------------------------------
! Increment time_dx if required
!------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_dx = time_dx + (time_tmp2 - time_tmp)
  end if

end subroutine d1bydx

!===================
! First y-derivative
!===================
subroutine d1bydy( fin, fout )
  implicit none

!--------------------------------
! Declare input and output arrays
!--------------------------------
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: fout(nx, ny, nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: i, j, k
  real(kind=dp)                 :: fac

!------------------------
! Start timer if required
!------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!----------------
! Initialise fout
!----------------
  fout = 0.0d0

!------------------------------------------------
! Forward, out of place, transform in y direction
!------------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planfyout, nx, fin(1,1,k), nx, 1, fout(1,1,k), nx, 1)
  end do

!-----------------------------------
! Multiply by "+ik" in Fourier space
!-----------------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        fout(i,1,k)=0.0d0
        fout(i,(ny/2)+1,k)=0.0d0
        fac = fout(i,j+1,k)
        fout(i,j+1,k)   = -fout(i,ny+1-j,k)*aky1(j+1)    !Real part
        fout(i,ny+1-j,k) =  fac*aky1(ny+1-j)   !Imaginary part
      end do
    end do
  end do

!-------------------------------------
! Inverse, in-place, Fourier transform
!-------------------------------------
  do k=1, nlayer
    call rfftw_f77(planry, nx, fout(1,1,k), nx, 1, tmp, nx, 1)
    do j=1, ny
      do i=1, nx
        fout(i,j,k)=fout(i,j,k)/dble(ny)
      end do
    end do
  end do

!------------------------------
! Increment time_dy if required
!------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_dy = time_dy + (time_tmp2 - time_tmp)
  end if

end subroutine d1bydy

!=================
! 2nd X-derivative
!=================
subroutine d2bydx( fin, f1out, f2out )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: f1out(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: f2out(nx, ny, nlayer)

!-------------------
! Local Declarations
!-------------------
  integer                       :: i, j, k

!----------------------
! Start timer if needed
!----------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!---------------------------
! Initialise f1out and f2out
!---------------------------
  f1out = 0.0d0
  f2out = 0.0d0

!---------------------------
! Forward, out of place, FFT
!---------------------------
  do k=1, nlayer
    call rfftw_f77(planfxout, ny, fin(1,1,k), 1, nx, f2out(1,1,k), 1, nx)
  end do

!-----------------------------------------
! For first derivative multiply by "ik" 
! for second derivative multiply by "-k*k" 
!-----------------------------------------
  do k=1,nlayer
    do j=1,ny
      f1out(1,j,k)=0.0d0
      f1out((nx/2)+1,j,k)=0.0d0
      f2out(1,j,k)=0.0d0
      do i=1,nx/2-1
        f1out(i+1,j,k)    = -f2out(nx+1-i,j,k)*akx1(i+1)
        f1out(nx+1-i,j,k) = f2out(i+1,j,k)*akx1(nx+1-i) 
        f2out(i+1,j,k)    = f2out(i+1,j,k)*akx2(i+1) 
        f2out(nx+1-i,j,k) = f2out(nx+1-i,j,k)*akx2(nx+1-i) 
      end do
      f2out((nx/2)+1,j,k)  = f2out((nx/2)+1,j,k)*akx2((nx/2)+1) 
    end do
  end do

!---------------------------------------------
! Inverse, in-place FFT to find 1st derivative
!---------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planrx, ny, f1out(1,1,k), 1, nx, tmp, 1, nx)
    do j=1, ny
      do i=1, nx
        f1out(i,j,k)=f1out(i,j,k)/dble(nx)
      end do
    end do
  end do

!---------------------------------------------
! Inverse, in-place FFT to find 2nd derivative
!---------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planrx, ny, f2out(1,1,k), 1, nx, tmp, 1, nx)
    do j=1, ny
      do i=1, nx
        f2out(i,j,k)=f2out(i,j,k)/dble(nx)
      end do
    end do
  end do

!-------------------------------
! Increment time_d2x if required
!-------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_d2x = time_d2x + (time_tmp2 - time_tmp)
  endif

end subroutine d2bydx


!=================
! 2nd Y derivative
!=================
subroutine d2bydy( fin, f1out, f2out )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: f1out(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: f2out(nx, ny, nlayer)

!-------------------
! Local Declarations
!-------------------
  integer                       :: i, j, k

!----------------------
! Start timer if needed
!----------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!---------------------------
! Initialise f1out and f2out
!---------------------------
  f1out = 0.0d0
  f2out = 0.0d0

!--------------------------
! Forward, out-of-place FFT
!--------------------------
  do k=1, nlayer
    call rfftw_f77(planfyout, nx, fin(1,1,k), nx, 1, tmp(1,1,k), nx, 1)
  end do

!-----------------------------------------
! For first derivative multiply by "ik" 
! for second derivative multiply by "-k*k" 
!-----------------------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        f1out(i,j+1,k)    = -tmp(i,ny+1-j,k)*aky1(j+1)
        f1out(i,ny+1-j,k) = tmp(i,j+1,k)*aky1(ny+1-j)
        f2out(i,j+1,k)    = tmp(i,j+1,k)*aky2(j+1) 
        f2out(i,ny+1-j,k) = tmp(i,ny+1-j,k)*aky2(ny+1-j) 
        f1out(i,1,k) = 0.0d0
        f2out(i,1,k) = 0.0d0
        f1out(i,(ny/2)+1,k) = 0.0d0
        f2out(i,(ny/2)+1,k)  = tmp(i,(ny/2)+1,k)*aky2((ny/2)+1) 
      end do
    end do
  end do

!---------------------------------------------
! Inverse, in-place FFT to find 1st derivative
!---------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planry, nx, f1out(1,1,k), nx, 1, tmp, nx, 1)
    do j=1, ny
      do i=1, nx
        f1out(i,j,k)=f1out(i,j,k)/dble(ny)
      end do
    end do
  end do

!---------------------------------------------
! Inverse, in-place FFT to find 2nd derivative
!---------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planry, nx, f2out(1,1,k), nx, 1, tmp, nx, 1)
    do j=1, ny
      do i=1, nx
        f2out(i,j,k)=f2out(i,j,k)/dble(ny)
      end do
    end do
  end do

!-----------------------------
! Increment time_2dy if needed
!-----------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_d2y = time_d2y + (time_tmp2 - time_tmp)
  end if

end subroutine d2bydy

!===================
! Dealiasing routine
!===================
subroutine dealias(z)
  implicit none

!------------------------------
! Declare input/output variable
!------------------------------
  real(kind=dp), intent(inout)  :: z(nx, ny, nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: i, j, k

!------------------------
! Start timer if required
!------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!--------------------------------------
! Dealias in the x-direction
! 
! Forward, in-place, FFT in x direction
!--------------------------------------
  do k=1, nlayer
    call rfftw_f77(planfx, ny, z(1,1,k), 1, nx, tmp(1,1,k), 1, nx)
  end do

!------------------------------
! Cut out high wavenumber modes
!------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,icutoffx
        z((nx/2)+1+i,j,k) = 0.0d0
        z((nx/2)+1-i,j,k) = 0.0d0
      end do
      z((nx/2)+1,j,k) = 0.0d0
    end do
  end do

!--------------------------------------
! Inverse, in-place, FFT in x direction
!--------------------------------------
  do k=1, nlayer
    call rfftw_f77(planrx, ny, z(1,1,k), 1, nx, tmp(1,1,k), 1, nx)
    do j=1, ny
      do i=1, nx
        z(i,j,k)=z(i,j,k)/dble(nx)
      end do
    end do
  end do

!--------------------------------------
! Dealias in the y-direction
! 
! Forward, in-place, FFT in y direction
!--------------------------------------
  do k=1, nlayer
    call rfftw_f77(planfy, nx, z(1,1,k), nx, 1, tmp(1,1,k), nx, 1)
  end do

!------------------------------
! Cut out high wavenumber modes
!------------------------------
  do k=1,nlayer
    do i=1,nx
      z(i,(ny/2)+1,k) = 0.0d0
      do j=1,icutoffy
        z(i,(ny/2)+1+j,k) = 0.0d0
        z(i,(ny/2)+1-j,k) = 0.0d0
      end do
    end do
  end do

!--------------------------------------
! Inverse, in-place, FFT in y direction
!--------------------------------------
  do k=1, nlayer
    call rfftw_f77(planry, nx, z(1,1,k), nx, 1, tmp(1,1,k), nx, 1)
    do j=1, ny
      do i=1, nx
        z(i,j,k)=z(i,j,k)/dble(ny)
      end do
    end do
  end do

!---------------------------------
! Increment time_dealias if needed
!---------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_dealias = time_dealias + (time_tmp2 - time_tmp)
  end if

end subroutine dealias

!=======================
! Forward (2D) transform
!=======================
subroutine forwardfft( fin, fout )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: fout(nx, ny, nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: k

!--------------------
! Set timer if needed
!--------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!----------------
! Initialise fout
!----------------
  fout = 0.0d0

!---------------------------------------------------
! Forward, out of place transform in the x-direction
!---------------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planfxout, ny, fin(1,1,k), 1, nx, fout(1,1,k), 1, nx)
  end do

!-----------------------------------------------
! Forward, in-place transform in the y-direction
!-----------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planfy, nx, fout(1,1,k), nx, 1, tmp(1,1,k), nx, 1)
  end do

!---------------------------------
! Increment time_forward if needed
!---------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_forward = time_forward + (time_tmp2 - time_tmp)
  end if

end subroutine forwardfft

!-----------------------
! Inverse (2D) transform
!-----------------------
subroutine inversefft( fin, fout )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: fout(nx, ny, nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: i, j, k

!--------------------
! Set timer if needed
!--------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!----------------
! Initialise fout
!----------------
  fout = fin

!-----------------------------------------------
! Inverse, in-place transform in the x-direction
!-----------------------------------------------
   do k=1, nlayer
     call rfftw_f77(planrx, ny, fout(1,1,k), 1, nx, tmp, 1, nx)
     do j=1, ny
       do i=1, nx
         fout(i,j,k)=fout(i,j,k)/dble(nx)
       end do
     end do
   end do

!-----------------------------------------------
! Inverse, in-place transform in the y-direction
!-----------------------------------------------
  do k=1, nlayer
    call rfftw_f77(planry, nx, fout(1,1,k), nx, 1, tmp, nx, 1)
    do j=1, ny
      do i=1, nx
        fout(i,j,k)=fout(i,j,k)/dble(ny)
      end do
    end do
  end do


!---------------------------------
! Increment time_forward if needed
!---------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_reverse = time_reverse + (time_tmp2 - time_tmp)
  end if

end subroutine inversefft

!==================================================================

subroutine expand( fin, fout )

  use types
  use dimens
#ifdef TIMING
  use timing
#endif

  implicit none
  real(kind=dp), intent(in)     :: fin(nx, ny, nlayer)
  real(kind=dp), intent(out)    :: fout(2*nx, 2*ny, nlayer)
  integer                       :: i,j,k


  fout(:,:,:)=0.0d0
  
  do k=1,nlayer
    do j=1,ny
      do i=1,2*nx
        interp(i,j,k)=0.0d0
        tmpnew4(i,j,k)=0.0d0
      end do
    end do
  end do

  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        tmp(i,j,k)=0.0d0
      end do
    end do
  end do

  do k=1,nlayer
    do j=1,2*ny
      do i=1,2*nx
        fout(i,j,k)=0.0d0
        tmpnew(i,j,k)=0.0d0
      end do
    end do
  end do

  do k=1, nlayer
    call rfftw_f77(planfx, ny, fin(1,1,k), 1, nx, tmp(1,1,k), 1, nx)
  end do

  do k=1,nlayer
    do j=1,ny
      interp(1,j,k)=fin(1,j,k)
      do i=1,nx/2-1
        interp(i+1,j,k)=fin(i+1,j,k)     
        interp(2*nx+1-i,j,k)=fin(nx+1-i,j,k)
      end do
      interp(nx/2+1,j,k)=fin(nx/2+1,j,k)
    end do
  end do

  do k=1, nlayer
    call rfftw_f77(planrxnew, ny, interp(1,1,k), 1, 2*nx, tmpnew4(1,1,k), &
                    1, 2*nx)
    do j=1, ny
      do i=1, 2*nx
         interp(i,j,k)=interp(i,j,k)/dble(nx)
      end do
    end do
  end do

  do k=1, nlayer
    call rfftw_f77(planfy, 2*nx, interp(1,1,k), 2*nx, 1, tmpnew4(1,1,k), &
                    2*nx, 1)
  end do

  do k=1,nlayer
    do i=1,2*nx   
      fout(i,1,k)=interp(i,1,k) 
      do j=1,ny/2-1
        fout(i,j+1,k)=interp(i,j+1,k)
        fout(i,2*ny+1-j,k)=interp(i,ny+1-j,k)      
      end do
      fout(i,ny/2+1,k)=interp(i,ny/2+1,k)
    end do
  end do

  do k=1, nlayer
    call rfftw_f77(planrynew, 2*nx, fout(1,1,k), 2*nx, 1, tmpnew(1,1,k), &
                     2*nx, 1)
    do j=1, 2*ny
      do i=1, 2*nx
        fout(i,j,k)=fout(i,j,k)/dble(ny)
      end do
    end do
  end do

end subroutine expand


END MODULE fft
