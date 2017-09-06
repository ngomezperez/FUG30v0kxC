!-----------------------------------------------------------------------
! test_dbydxy: Test the x and y derivative routines
! 
! Note: All output bar final third of dealiased bits should be 
! approximately equal to zero for spectral code
!
! Compact code wont look as good on here as we're measuring it against
! the spectral basis functions. Since the error terms go as the 7th or
! 8th derivative, this gets very large for these special functions for
! large wavenumber eg. (k*PI)^7, so error term cannot remain small.
! (PJB 8/12/03)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
PROGRAM test_dbydxy

  use types
  use comms
  use params
  use trans_params
  use dimens
  use control
  use fft
  use timing
  use constants
  use arrays
  use initial

  implicit none
 
  integer        :: ier

  real(kind=dp)  :: dt1
  real(kind=dp)  :: t1, t2

  integer        :: i, j, k
  real(kind=dp)  :: wspace

  real(kind=dp) :: ttime
  real(kind=dp) :: start_time, finish_time
  real(kind=dp) :: omega
  real(kind=dp) :: diff1, diff1a, diff2, diffd
  integer       :: iw
  integer       :: iter
  
!===================================================================

!------------------------------
! Start communication if needed
!------------------------------
  call init_comms(ier)
  if (myrank .eq. 0) then 
  
  end if

!------------------------------
! Read parameters and constants
!------------------------------ 
  If (myrank .eq. 0) then
  print*,'calling readparam '
  end if
  call readparam(iter, ttime, dt1, ier)
  call set_constants()  

!--------------------------------------------------
! Print short greeting message and array dimensions
!--------------------------------------------------
  if (myrank .eq. 0) then
    print*,' Testing derivatives '
    print*,''
    print*,' nx ', nx, ' ny ', ny, ' nz ', nz
  end if

  izstart = ((myrank)*nz)/nproc + 1
  izfinish = ((myrank+1)*nz)/nproc 
  nlayer = izfinish - izstart + 1

  if ( myrank .eq. 0 ) then
    print*,'DIMENSION of arrays ', nx*ny*nlayer
  end if

!------------------------
! Initialize FFT routines
!------------------------
  call init_fft(ier)

!---------------------
! Allocate work arrays
!---------------------
  call alloc_arrays(ier)
  if ( ier .lt. 0 ) then
    print*,'TROUBLE allocating work arrays '
    call abort_comms(ier)
  end if

!--------------
! Start timings
!--------------
  call findtime(start_time)

!------------
! Start tests
!------------
  if (myrank .eq. 0) then
   print*,'-----------------------------------'
   print*,'Test the x cos derivatives  ', xmax
   print*,'-----------------------------------'
  end if

!========================================================================
! Test derivatives using a single frequency --- NOTE that
! this has to be commensurate with box or else there are
! discontinuities at the boundaries and these show up in
! the derivatives
!=========================================================================== 
!  f(x) = cos( omega*(2*pi/xmax)*x )
!============================================================================
   do iw=0, nx/2
!============================================================================

   omega = dble(iw)*(8.0d0*atan(1.0d0)/xmax)

   do k=1,nlayer
     do j=1,ny
       wspace = 0.0d0
       do i=1,nx
         wk1(i,j,k) =  k*cos( omega*wspace )
         wspace = wspace + dx 
       end do
     end do
   end do


   call d1bydx(wk1, wk2 )
   call d2bydx(wk1, wk3, wk4 )

   call dealias(wk1)

   diffd  = 0.0d0
   diff1  = 0.0d0
   diff2  = 0.0d0
   diff1a = 0.0d0
   dx = xmax/dble(nx)
   do k=1,nlayer
     do j=1,ny
       wspace = 0.0d0
       do i=1,nx
         diff1  = max( diff1,  abs( wk2(i,j,k) + k*omega*sin(omega*wspace) ))
         diff1a = max( diff1a, abs( wk3(i,j,k) + k*omega*sin(omega*wspace) ) )
         diff2  = max( diff2 , abs( wk4(i,j,k) + k*omega*omega*cos(omega*wspace) ))
         diffd  = max( diffd , abs( wk1(i,j,k) - k*cos(omega*wspace) ))
         wspace = wspace + dx 
       end do
     end do
   end do
   if (myrank .eq. 0) then
     print*,iw, diff1, diff1a, diff2, diffd
   end if

!============================================================================
   end do
!============================================================================


  if (myrank .eq. 0) then
   print*,'-----------------------------------'
   print*,'Test the y cos derivatives  ', ymax
   print*,'-----------------------------------'
  end if

!=========================================================================== 
!  f(x) = cos( omega*(2*pi/ymax)*y )
!============================================================================
   do iw=0, ny/2
!============================================================================

   omega = dble(iw)*(8.0d0*atan(1.0d0)/ymax)

   do k=1,nlayer
     wspace = 0.0d0
     do j=1,ny
       t1 = k*cos( omega*dble(j-1)*dy )
       do i=1,nx
         wk1(i,j,k) =  t1 
       end do
       wspace = wspace + dy 
     end do
   end do

   call d1bydy(wk1, wk2 )
   call d2bydy(wk1, wk3, wk4 )

   call dealias(wk1)

   diffd  = 0.0d0
   diff1  = 0.0d0
   diff2  = 0.0d0
   diff1a = 0.0d0
   do k=1,nlayer
     wspace = 0.0d0
     do j=1,ny
       t1 = k*cos( omega*wspace)
       t2 = k*sin( omega*wspace)
       do i=1,nx
         diff1  = max( diff1,  abs( wk2(i,j,k) + omega*t2)  )
         diff1a = max( diff1a, abs( wk3(i,j,k) + omega*t2) )
         diff2  = max( diff2 , abs( wk4(i,j,k) + omega*omega*t1) )
         diffd  = max( diffd , abs( wk1(i,j,k) - t1 ))
       end do
       wspace = wspace + dy 
     end do
   end do
   if (myrank .eq. 0) then
     print*,iw, diff1, diff1a, diff2, diffd
   end if

!============================================================================
   end do
!============================================================================

  if (myrank .eq. 0) then
   print*,'-----------------------------------'
   print*,'Test the x sine derivatives  ', xmax
   print*,'-----------------------------------'
  end if


!
! Test derivatives using a single frequency --- NOTE that
! this has to be commensurate with box or else there are
! discontinuities at the boundaries and these show up in
! the derivatives
!=========================================================================== 
!  f(x) = sin( omega*(2*pi/xmax)*x )
!============================================================================
   do iw=0, nx/2
!============================================================================

   omega = dble(iw)*(8.0d0*atan(1.0d0)/xmax)

   do k=1,nlayer
     do j=1,ny
       wspace = 0.0d0
       do i=1,nx
         wk1(i,j,k) =  k*sin( omega*wspace )
         wspace = wspace + dx 
       end do
     end do
   end do


   call d1bydx(wk1, wk2 )
   call d2bydx(wk1, wk3, wk4 )

   call dealias(wk1)

   diffd  = 0.0d0
   diff1  = 0.0d0
   diff2  = 0.0d0
   diff1a = 0.0d0
   dx = xmax/dble(nx)
   do k=1,nlayer
     do j=1,ny
       wspace = 0.0d0
       do i=1,nx
         diff1  = max( diff1,  abs( wk2(i,j,k) - k*omega*cos(omega*wspace) ))
         diff1a = max( diff1a, abs( wk3(i,j,k) - k*omega*cos(omega*wspace) ) )
         diff2  = max( diff2 , abs( wk4(i,j,k) + k*omega*omega*sin(omega*wspace) ))
         diffd  = max( diffd , abs( wk1(i,j,k) - k*sin(omega*wspace) ))
         wspace = wspace + dx 
       end do
     end do
   end do
   if (myrank .eq. 0) then
     print*,iw, diff1, diff1a, diff2, diffd
   end if
!============================================================================
   end do
!============================================================================


  if (myrank .eq. 0) then
   print*,'-----------------------------------'
   print*,'Test the y sine derivatives  ', ymax
   print*,'-----------------------------------'
  end if

!=========================================================================== 
!  f(x) = sin( omega*(2*pi/ymax)*y )
!============================================================================
   do iw=0, ny/2
!============================================================================

   omega = dble(iw)*(8.0d0*atan(1.0d0)/ymax)

   do k=1,nlayer
     wspace = 0.0d0
     do j=1,ny
       t1 = k*sin( omega*dble(j-1)*dy )
       do i=1,nx
         wk1(i,j,k) =  t1 
       end do
       wspace = wspace + dy 
     end do
   end do

   call d1bydy(wk1, wk2 )
   call d2bydy(wk1, wk3, wk4 )

   call dealias(wk1)

   diffd  = 0.0d0
   diff1  = 0.0d0
   diff2  = 0.0d0
   diff1a = 0.0d0
   do k=1,nlayer
     wspace = 0.0d0
     do j=1,ny
       t1 = k*sin( omega*wspace)
       t2 = k*cos( omega*wspace)
       do i=1,nx
         diff1  = max( diff1,  abs( wk2(i,j,k) - omega*t2)  )
         diff1a = max( diff1a, abs( wk3(i,j,k) - omega*t2) )
         diff2  = max( diff2 , abs( wk4(i,j,k) + omega*omega*t1) )
         diffd  = max( diffd , abs( wk1(i,j,k) - t1 ))
       end do
       wspace = wspace + dy 
     end do
   end do
   if (myrank .eq. 0) then 
     print*,iw, diff1, diff1a, diff2, diffd
   end if

!============================================================================
   end do
!============================================================================


!----------------------------------------------------------------------
! Tidy up and finish
!----------------------------------------------------------------------
   call findtime(finish_time)
   if (myrank .eq. 0) then
     print*,'Time taken = ',finish_time-start_time
   end if

   call close_fft()
   call close_arrays()

#ifdef MPI
   call MPI_Barrier(comm,ier)
#endif
   call exit_comms()


 END PROGRAM test_dbydxy 
