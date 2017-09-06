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
  real(kind=dp), allocatable :: unew(:,:,:)
  real(kind=dp), allocatable :: vnew(:,:,:)
  real(kind=dp), allocatable :: wnew(:,:,:)
  real(kind=dp), allocatable :: tnew(:,:,:)
  real(kind=dp), allocatable :: rnew(:,:,:)
  real(kind=dp), allocatable :: bxnew(:,:,:)
  real(kind=dp), allocatable :: bynew(:,:,:)
  real(kind=dp), allocatable :: bznew(:,:,:)

!--------------------
! Declare work arrays
!--------------------
  real(kind=dp), allocatable :: wk1(:,:,:)
  real(kind=dp), allocatable :: wk2(:,:,:)
  real(kind=dp), allocatable :: wk3(:,:,:)
  real(kind=dp), allocatable :: wk4(:,:,:)
  real(kind=dp), allocatable :: wk1new(:,:,:)
  real(kind=dp), allocatable :: wk2new(:,:,:)
  real(kind=dp), allocatable :: wk3new(:,:,:)
  real(kind=dp), allocatable :: wk4new(:,:,:)
  real, allocatable :: rwk(:)

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
  integer               :: ier(25)
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
  allocate(tnew(2*nx,2*ny,nlayer), stat=ier(9))
  allocate(rnew(2*nx,2*ny,nlayer), stat=ier(10))
  allocate(unew(2*nx,2*ny,nlayer), stat=ier(11))
  allocate(vnew(2*nx,2*ny,nlayer), stat=ier(12))
  allocate(wnew(2*nx,2*ny,nlayer), stat=ier(13))
  allocate(bxnew(2*nx,2*ny,nlayer), stat=ier(14))
  allocate(bynew(2*nx,2*ny,nlayer), stat=ier(15))
  allocate(bznew(2*nx,2*ny,nlayer), stat=ier(16))
  
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
  tnew = 0.0d0
  rnew = 0.0d0
  unew = 0.0d0
  vnew = 0.0d0
  wnew = 0.0d0
  bxnew = 0.0d0
  bynew = 0.0d0
  bznew = 0.0d0

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
    allocate(wk1(nx,ny,nlayer), stat=ier(17))
    allocate(wk2(nx,ny,nlayer), stat=ier(18))
    allocate(wk3(nx,ny,nlayer), stat=ier(19))
    allocate(wk4(nx,ny,nlayer), stat=ier(20))
    allocate(wk1new(2*nx,2*ny,nlayer), stat=ier(21))
    allocate(wk2new(2*nx,2*ny,nlayer), stat=ier(22))
    allocate(wk3new(2*nx,2*ny,nlayer), stat=ier(23))
    allocate(wk4new(2*nx,2*ny,nlayer), stat=ier(24))
    allocate(rwk(nx), stat=ier(25))
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
    allocate(wk1(nx,ny,nlymax), stat=ier(17))
    allocate(wk2(nx,ny,nlymax), stat=ier(18))
    allocate(wk3(nx,ny,nlymax), stat=ier(19))
    allocate(wk4(nx,ny,nlymax), stat=ier(20))
    allocate(wk1new(2*nx,2*ny,nlymax), stat=ier(21))
    allocate(wk2new(2*nx,2*ny,nlymax), stat=ier(22))
    allocate(wk3new(2*nx,2*ny,nlymax), stat=ier(23))
    allocate(wk4new(2*nx,2*ny,nlymax), stat=ier(24))
    allocate(rwk(nx), stat=ier(25))
  end if

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
       (ier(25) .gt. 0) )   then
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
  wk1new = 0.0d0
  wk2new = 0.0d0
  wk3new = 0.0d0
  wk4new = 0.0d0
  rwk = 0.0

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
  deallocate(tnew)
  deallocate(rnew)
  deallocate(unew)
  deallocate(vnew)
  deallocate(wnew)
  deallocate(bxnew)
  deallocate(bynew)
  deallocate(bznew)

!----------------------------------
! Deallocate memory for work arrays
!----------------------------------
  deallocate(wk1)
  deallocate(wk2)
  deallocate(wk3)
  deallocate(wk4)
  deallocate(wk1new)
  deallocate(wk2new)
  deallocate(wk3new)
  deallocate(wk4new)
  deallocate(rwk)

END SUBROUTINE close_arrays

END MODULE arrays
