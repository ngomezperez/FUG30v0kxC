!==============================================
! This module holds the large arrays 
!==============================================
MODULE arrays
  use types
  use comms
  use dimens

  implicit none

!------------------
! Declare variables
!------------------
  real(kind=dp), allocatable :: u(:,:,:)
  real(kind=dp), allocatable :: v(:,:,:)
  real(kind=dp), allocatable :: w(:,:,:)
  real(kind=dp), allocatable :: t(:,:,:)
  real(kind=dp), allocatable :: r(:,:,:)
  real(kind=dp), allocatable :: bx(:,:,:)
  real(kind=dp), allocatable :: by(:,:,:)
  real(kind=dp), allocatable :: bz(:,:,:)
  real(kind=dp), allocatable :: pol(:,:,:)
  real(kind=dp), allocatable :: tor(:,:,:)
  real(kind=dp), allocatable :: jsq(:,:,:)
  real(kind=dp), allocatable :: ex(:,:,:)
  real(kind=dp), allocatable :: ey(:,:,:)
  real(kind=dp), allocatable :: ez(:,:,:)
  real(kind=dp), allocatable :: bxbar(:)
  real(kind=dp), allocatable :: bybar(:)
!=============
! NGP: Aug 29, 2017
!=============
  real(kind=dp), allocatable :: kappaZ(:)
  real(kind=dp), allocatable :: dkappadZ(:)

!--------------------
! Declare derivatives
!--------------------
  real(kind=dp), allocatable :: udot(:,:,:,:)
  real(kind=dp), allocatable :: vdot(:,:,:,:)
  real(kind=dp), allocatable :: wdot(:,:,:,:)
  real(kind=dp), allocatable :: tdot(:,:,:,:)
  real(kind=dp), allocatable :: rdot(:,:,:,:)
  real(kind=dp), allocatable :: poldot(:,:,:,:)
  real(kind=dp), allocatable :: tordot(:,:,:,:)
  real(kind=dp), allocatable :: bxbardot(:,:)
  real(kind=dp), allocatable :: bybardot(:,:)  

!--------------------
! Declare work arrays
!--------------------
  real(kind=dp), allocatable :: wk1(:,:,:)
  real(kind=dp), allocatable :: wk2(:,:,:)
  real(kind=dp), allocatable :: wk3(:,:,:)
  real(kind=dp), allocatable :: wk4(:,:,:)

!------------------------
! Declare real work array
!------------------------
  real, allocatable          :: rwk(:)

!-----------------------
! Declare 1D work arrays
!-----------------------
  real(kind=dp), allocatable :: wk1_1D(:)
  real(kind=dp), allocatable :: wk2_1D(:)
!Nat: Aug24 declaring the new variables for thermal conductivity
  real(kind=dp), allocatable :: ka(:)
  real(kind=dp), allocatable :: dkadz(:)

!-----------------------------------
! Best place to put mean value of Bz
!-----------------------------------
  real(kind=dp)              :: bzbar

CONTAINS

!=========================================
! Allocate memory and initialise variables
!=========================================
SUBROUTINE alloc_arrays(iret)
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer, intent(out)  :: iret

!-------------------
! Local declarations
!-------------------
!=============
! NGP: Aug 29, 2017. 
! Changed 32 for 34 on ier for cheking array status 
!=============
  integer               :: ier(34)
  integer               :: nlymax
  integer               :: i, ist, ifn
 
!-----------------------------------
! Set iret to the number of bytes to
! be allocated
!-----------------------------------
  iret = 8*24*nx*ny*nlayer

!-------------------
! Allocate variables
!-------------------
  allocate(t(nx,ny,nlayer), stat=ier(1))
  allocate(r(nx,ny,nlayer), stat=ier(2))
  allocate(u(nx,ny,nlayer), stat=ier(3))
  allocate(v(nx,ny,nlayer), stat=ier(4))
  allocate(w(nx,ny,nlayer), stat=ier(5))
  allocate(bx(nx,ny,nlayer), stat=ier(6))
  allocate(by(nx,ny,nlayer), stat=ier(7))
  allocate(bz(nx,ny,nlayer), stat=ier(8))
  allocate(pol(nx,ny,nlayer), stat=ier(9))
  allocate(tor(nx,ny,nlayer), stat=ier(10))
  allocate(jsq(nx,ny,nlayer), stat=ier(11))
  allocate(ex(nx,ny,nlayer), stat=ier(12))
  allocate(ey(nx,ny,nlayer), stat=ier(13))
  allocate(ez(nx,ny,nlayer), stat=ier(14))
  allocate(bxbar(nlayer), stat=ier(15))
  allocate(bybar(nlayer), stat=ier(16))
!=============
! NGP: Aug 29, 2017
!=============
  allocate(kappaZ(nlayer), stat=ier(33))
  allocate(dkappadZ(nlayer), stat=ier(34))
  

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
       (ier(8) .gt. 0) .or.    &
       (ier(9) .gt. 0) .or.    &
       (ier(10) .gt. 0) .or.    &
       (ier(11) .gt. 0) .or.    &
       (ier(12) .gt. 0) .or.    &
       (ier(13) .gt. 0) .or.    &
       (ier(14) .gt. 0) .or.    &
       (ier(15) .gt. 0) .or.    &
       (ier(16) .gt. 0) )   then
     iret = -1
     return
  end if

!---------------------
! Initialise variables
!---------------------
  t = 0.0d0
  r = 0.0d0
  u = 0.0d0
  v = 0.0d0
  w = 0.0d0
  bx = 0.0d0
  by = 0.0d0
  bz = 0.0d0
  pol = 0.0d0
  tor = 0.0d0
  jsq = 0.0d0
  ex = 0.0d0
  ey = 0.0d0
  ez = 0.0d0
  bxbar = 0.0d0
  bybar = 0.0d0
  bzbar = 0.0d0
!=============
! NGP: Aug 29, 2017
!=============
  kappaZ=0.0d0
  dkappadZ=0.0d0

!---------------------
! Allocate derivatives
!---------------------
  allocate(tdot(nx,ny,nlayer,3), stat=ier(17))
  allocate(rdot(nx,ny,nlayer,3), stat=ier(18))
  allocate(udot(nx,ny,nlayer,3), stat=ier(19))
  allocate(vdot(nx,ny,nlayer,3), stat=ier(20))
  allocate(wdot(nx,ny,nlayer,3), stat=ier(21))
  allocate(poldot(nx,ny,nlayer,3), stat=ier(22))
  allocate(tordot(nx,ny,nlayer,3), stat=ier(23))
  allocate(bxbardot(nlayer,3), stat=ier(24))
  allocate(bybardot(nlayer,3), stat=ier(25))

!-------------
! Check result
!-------------
  if ( (ier(17) .gt. 0) .or.    &
       (ier(18) .gt. 0) .or.    &
       (ier(19) .gt. 0) .or.    &
       (ier(20) .gt. 0) .or.    &
       (ier(21) .gt. 0) .or.    &
       (ier(22) .gt. 0) .or.    &
       (ier(23) .gt. 0) .or.    &
       (ier(24) .gt. 0) .or.    &
       (ier(33) .gt. 0) .or.    &
       (ier(34) .gt. 0) .or.    &
       (ier(25) .gt. 0) )   then
     iret = -2
     return
  end if

!-----------------------
! Initialise derivatives
!-----------------------
  tdot  = 0.0d0
  rdot  = 0.0d0
  udot  = 0.0d0
  vdot  = 0.0d0
  wdot  = 0.0d0
  poldot = 0.0d0
  tordot = 0.0d0
  bxbardot = 0.0d0
  bybardot = 0.0d0

!-----------------------------------------
! Allocate work arrays
!
! Note that the arrays on "rank 0" need to 
! be as big as those on any of the other
! processors. This is done through nlymax
!
! Also do real work array here
!-----------------------------------------
  iret = 0
  if ( myrank .ne. 0 ) then
    allocate(wk1(nx,ny,nlayer), stat=ier(26))
    allocate(wk2(nx,ny,nlayer), stat=ier(27))
    allocate(wk3(nx,ny,nlayer), stat=ier(28))
    allocate(wk4(nx,ny,nlayer), stat=ier(29))
    allocate(rwk(nx), stat=ier(30))
  else 
    nlymax = nlayer
    do i=1, nproc-1
      ist = ((i*nz)/nproc) + 1
      ifn = ((i+1)*nz)/nproc
      nlymax = max(nlymax, ifn-ist+1)
    end do
    if ( nlymax .ne. nlayer ) then
      print*,'NOTE: increased size of work array 3 dimension ',nlayer, nlymax
    end if
    allocate(wk1(nx,ny,nlymax), stat=ier(26))
    allocate(wk2(nx,ny,nlymax), stat=ier(27))
    allocate(wk3(nx,ny,nlymax), stat=ier(28))
    allocate(wk4(nx,ny,nlymax), stat=ier(29))
    allocate(rwk(nx), stat=ier(30))
  end if

!-------------
! Check result
!-------------
  if ( (ier(26) .gt. 0) .or.    &
       (ier(27) .gt. 0) .or.    &
       (ier(28) .gt. 0) .or.    &
       (ier(29) .gt. 0) .or.    &
       (ier(30) .gt. 0) )   then
     iret = -3
     return
  end if

!-----------------------
! Initialise work arrays
!-----------------------
  wk1 = 0.0d0
  wk2 = 0.0d0
  wk3 = 0.0d0
  wk4 = 0.0d0
  rwk = 0.0

!---------------------
! Allocate derivatives
!---------------------
  allocate(wk1_1D(nlayer), stat=ier(31))
  allocate(wk2_1D(nlayer), stat=ier(32))

!-------------
! Check result
!-------------
  if ( (ier(31) .gt. 0) .or.    &
       (ier(32) .gt. 0) )   then
     iret = -4
     return
  end if

!--------------------------
! Initialise 1D work arrays
!--------------------------
  wk1_1D  = 0.0d0
  wk2_1D  = 0.0d0

END SUBROUTINE alloc_arrays


!======================
! Deallocate the memory
!======================
SUBROUTINE close_arrays
  implicit none

!-------------------------------
! Dellocate memory for variables
!-------------------------------
  deallocate(t)
  deallocate(r)
  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(bx)
  deallocate(by)
  deallocate(bz)
  deallocate(pol)
  deallocate(tor)
  deallocate(jsq)
  deallocate(ex)
  deallocate(ey)
  deallocate(ez)
  deallocate(bxbar)
  deallocate(bybar)
!=============
! NGP: Aug 29, 2017
!=============
  deallocate(kappaZ)
  deallocate(dkappadZ)

!----------------------------------
! Deallocate memory for derivatives
!----------------------------------
  deallocate(tdot)
  deallocate(rdot)
  deallocate(udot)
  deallocate(vdot)
  deallocate(wdot)
  deallocate(poldot)
  deallocate(tordot)
  deallocate(bxbardot)
  deallocate(bybardot)

!----------------------------------
! Deallocate memory for work arrays
!----------------------------------
  deallocate(wk1)
  deallocate(wk2)
  deallocate(wk3)
  deallocate(wk4)

!-------------------------------------
! Deallocate memory for 1D work arrays
!-------------------------------------
  deallocate(wk1_1D)
  deallocate(wk2_1D)

END SUBROUTINE close_arrays

END MODULE arrays
