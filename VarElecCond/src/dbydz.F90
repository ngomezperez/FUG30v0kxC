!========================================================
! This module deals with Z-derivatives and communications
!========================================================
MODULE dbydz
  use types
  use dimens
  use comms
  use params
  use trans_params
  use constants
  use control
  use arrays
  use timing

  implicit none

!--------------------
! Global declarations
!--------------------
  real(kind=dp)              :: ho12
  real(kind=dp)              :: sc1

!-------------------------
! Communications variables
!-------------------------
#ifdef MPI
  integer                    :: rbufm_counter
  integer                    :: rbufp_counter
  integer                    :: rbufn_counter
  integer, allocatable       :: sendup_tag(:)
  integer, allocatable       :: senddown_tag(:)
  integer, allocatable       :: busendup_tag(:)
  integer, allocatable       :: busenddown_tag(:)
  integer, allocatable       :: bvsendup_tag(:)
  integer, allocatable       :: bwsendup_tag(:)
  integer, allocatable       :: bvsenddown_tag(:)
  integer, allocatable       :: bwsenddown_tag(:)
  real(kind=dp), allocatable :: rbufp(:,:,:)
  real(kind=dp), allocatable :: rbufm(:,:,:)
  real(kind=dp), allocatable :: rbufn(:,:,:)
  integer, allocatable       :: d1_distsend(:)
  integer, allocatable       :: d2_distsend(:)
  integer, allocatable       :: b_distsend(:)
  integer, allocatable       :: d1up_distsend(:)
  integer, allocatable       :: d1_distrecv(:)
  integer, allocatable       :: d2_distrecv(:)
  integer, allocatable       :: b_distrecv(:)
  integer, allocatable       :: d1up_distrecv(:)
  integer, allocatable       :: layer_size(:)
  integer, allocatable       :: layer_start(:)
  integer, allocatable       :: layer_finish(:)
#endif

CONTAINS

!===================================================
! Set up memory for the transmit and receive buffers
!===================================================
subroutine init_dbydz(iret)
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer, intent(out) :: iret

!-------------------
! Local declarations
!-------------------
#ifdef MPI
  integer              :: ier(11)
  integer              :: k
#endif

!----------------------
! Default value of iret
!----------------------
  iret = 0

#ifdef MPI
!---------------------------------------------------
! Allocate enough memory to send/receive five layers
!---------------------------------------------------
  allocate(rbufp(nx,ny,5), stat=ier(1))
  allocate(rbufm(nx,ny,5), stat=ier(2))
  allocate(rbufn(nx,ny,5), stat=ier(3))

!-------------
! Check result
!-------------
  if (  ier(1) .gt. 0 ) then
    iret =  iret - 1
  end if
  if (  ier(2) .gt. 0 ) then
    iret =  iret - 2
  end if
  if (  ier(3) .gt. 0 ) then
    iret =  iret - 3
  end if

!------------
! Zero arrays
!------------  
  rbufp = 0.0d0
  rbufm = 0.0d0
  rbufn = 0.0d0

!----------------------------
! Allocate communication tags
!----------------------------
  allocate(sendup_tag(nproc),   stat=ier(1))
  allocate(senddown_tag(nproc), stat=ier(2))
  allocate(busendup_tag(nproc),   stat=ier(3))
  allocate(busenddown_tag(nproc), stat=ier(4))
  allocate(bvsendup_tag(nproc),   stat=ier(5))
  allocate(bvsenddown_tag(nproc), stat=ier(6))
  allocate(bwsendup_tag(nproc),   stat=ier(7))
  allocate(bwsenddown_tag(nproc), stat=ier(8))

!-------------
! Check result
!-------------
  if (  ier(1) .gt. 0 )  iret =  -4;
  if (  ier(2) .gt. 0 )  iret =  -5;
  if (  ier(3) .gt. 0 )  iret =  -6;
  if (  ier(4) .gt. 0 )  iret =  -7;
  if (  ier(5) .gt. 0 )  iret =  -8;
  if (  ier(6) .gt. 0 )  iret =  -9;
  if (  ier(7) .gt. 0 )  iret =  -10;
  if (  ier(8) .gt. 0 )  iret =  -11;

!------------
! Define tags
!------------
  do k=1,nproc 
    sendup_tag(k) = 10000+k
    senddown_tag(k) = 20000+k
    busendup_tag(k) = 30000+k
    bvsendup_tag(k) = 40000+k
    busenddown_tag(k) = 50000+k
    bvsenddown_tag(k) = 60000+k
    bwsendup_tag(k) = 70000+k
    bwsenddown_tag(k) = 80000+k
  end do  
  
!-----------------------------------
! Allocate layer distribution arrays
!-----------------------------------
  allocate(layer_size(nproc),   stat=ier(1))
  allocate(layer_start(nproc),  stat=ier(2))
  allocate(layer_finish(nproc), stat=ier(3))
  allocate(d1_distsend(nproc),  stat=ier(4))
  allocate(d2_distsend(nproc),  stat=ier(5))
  allocate(b_distsend(nproc),   stat=ier(6))
  allocate(d1up_distsend(nproc),stat=ier(7))
  allocate(d1_distrecv(nproc),  stat=ier(8))
  allocate(d2_distrecv(nproc),  stat=ier(9))
  allocate(b_distrecv(nproc),   stat=ier(10))
  allocate(d1up_distrecv(nproc), stat=ier(11))

!-------------
! Check result
!-------------
  if (  ier(1) .gt. 0 )  iret =  -12;
  if (  ier(2) .gt. 0 )  iret =  -13;
  if (  ier(3) .gt. 0 )  iret =  -14;
  if (  ier(4) .gt. 0 )  iret =  -15;
  if (  ier(5) .gt. 0 )  iret =  -16;
  if (  ier(6) .gt. 0 )  iret =  -17;
  if (  ier(7) .gt. 0 )  iret =  -18;
  if (  ier(8) .gt. 0 )  iret =  -19;
  if (  ier(9) .gt. 0 )  iret =  -20;
  if (  ier(10) .gt. 0 )  iret =  -21;
  if (  ier(11) .gt. 0 )  iret =  -22;

!---------------------------
! Define layer distributions
!---------------------------
  do k=1, nproc
    layer_start(k) = ((k-1)*nz)/nproc + 1
    layer_finish(k) = ((k*nz)/nproc) 
    layer_size(k) = layer_finish(k) - layer_start(k) + 1
  end do

!-------------------------
! Zero distribution arrays
!-------------------------
  do k=1, nproc
    d1_distsend(k)=0
    d2_distsend(k)=0
    b_distsend(k)=0
    d1up_distsend(k)=0
    d1_distrecv(k)=0
    d2_distrecv(k)=0
    b_distrecv(k)=0
    d1up_distrecv(k)=0
  end do  

!-------------------------
! Zero other MPI variables
!-------------------------
  rbufm_counter = 0
  rbufp_counter = 0 
  rbufn_counter = 0
  ncount        = 0 
  msg           = 0 

!-------------------------------------
! Now allocate send/recv distributions
!
! Nz>6 is assumed!
!-------------------------------------
! Only if nproc > 1
!-------------------------------------
  if (nproc.gt.1) then
    do k=1, nproc

!-------------------
! Do top layer first
!-------------------
    if (layer_start(k) .eq. 1) then
      if (layer_size(k) .eq. 1) then
        call getbelow_init(3, k, d1_distsend, d1_distrecv)
        call getbelow_init(5, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
        call getbelow_init(3, k, b_distsend, b_distrecv)
      else if (layer_size(k) .eq. 2) then
        call getbelow_init(2, k, d1_distsend, d1_distrecv)
        call getbelow_init(4, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
        call getbelow_init(2, k, b_distsend, b_distrecv)
      else if (layer_size(k) .eq. 3) then
        call getbelow_init(2, k, d1_distsend, d1_distrecv)
        call getbelow_init(3, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
        call getbelow_init(1, k, b_distsend, b_distrecv)
      else if (layer_size(k) .ge. 4) then
        call getbelow_init(2, k, d1_distsend, d1_distrecv)
        call getbelow_init(2, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
      end if

!------------------------------
! If layer starts at nz=2 then:
!------------------------------
    else if (layer_start(k) .eq. 2) then
      call getabove_init(1, k, d1_distsend, d1_distrecv)
      call getabove_init(1, k, d2_distsend, d2_distrecv)
      call getabove_init(1, k, d1up_distsend, d1up_distrecv)
      if (layer_size(k) .eq. 1) then
        call getbelow_init(2, k, d1_distsend, d1_distrecv)
        call getbelow_init(4, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
      else if (layer_size(k) .eq. 2) then
        call getbelow_init(2, k, d1_distsend, d1_distrecv)
        call getbelow_init(3, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
      else if (layer_size(k) .ge. 3) then
        call getbelow_init(2, k, d1_distsend, d1_distrecv)
        call getbelow_init(2, k, d2_distsend, d2_distrecv)
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
      end if

!------------------------------
! If layer starts at nz=3 then:
!------------------------------
    else if (layer_start(k) .eq. 3) then
      call getabove_init(2, k, d1_distsend, d1_distrecv)
      call getabove_init(2, k, d2_distsend, d2_distrecv)
      call getabove_init(2, k, d1up_distsend, d1up_distrecv)
      call getbelow_init(2, k, d1_distsend, d1_distrecv)
      call getbelow_init(2, k, d2_distsend, d2_distrecv)
      if (layer_finish(k) .eq. nz-2) then
        call getbelow_init(2, k, d1up_distsend, d1up_distrecv) 
      else
        call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
      end if

!--------------------------------
! If layer finishes at nz-2 then:
!--------------------------------
    else if (layer_finish(k) .eq. nz-2) then
      call getbelow_init(2, k, d1_distsend, d1_distrecv)
      call getbelow_init(2, k, d2_distsend, d2_distrecv)
      call getbelow_init(2, k, d1up_distsend, d1up_distrecv)
      call getabove_init(2, k, d1_distsend, d1_distrecv)
      call getabove_init(2, k, d2_distsend, d2_distrecv)
      call getabove_init(3, k, d1up_distsend, d1up_distrecv)

!--------------------------------
! If layer finishes at nz-1 then:
!--------------------------------
    else if (layer_finish(k) .eq. nz-1) then
      call getbelow_init(1, k, d1_distsend, d1_distrecv)
      call getbelow_init(1, k, d2_distsend, d2_distrecv)
      call getbelow_init(1, k, d1up_distsend, d1up_distrecv)
      if (layer_size(k) .eq. 1) then
        call getabove_init(2, k, d1_distsend, d1_distrecv)
        call getabove_init(4, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
      else if (layer_size(k) .eq. 2) then
        call getabove_init(2, k, d1_distsend, d1_distrecv)
        call getabove_init(3, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
      else if (layer_size(k) .ge. 3) then
        call getabove_init(2, k, d1_distsend, d1_distrecv)
        call getabove_init(2, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
      end if

!-------------------------------------------
! If layer finishes at the bottom of the box
!-------------------------------------------
    else if (layer_finish(k) .eq. nz) then
      if (layer_size(k) .eq. 1) then
        call getabove_init(3, k, d1_distsend, d1_distrecv)
        call getabove_init(5, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
        call getabove_init(3, k, b_distsend, b_distrecv)
      else if (layer_size(k) .eq. 2) then
        call getabove_init(2, k, d1_distsend, d1_distrecv)
        call getabove_init(4, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
        call getabove_init(2, k, b_distsend, b_distrecv)
      else if (layer_size(k) .eq. 3) then
        call getabove_init(2, k, d1_distsend, d1_distrecv)
        call getabove_init(3, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
        call getabove_init(1, k, b_distsend, b_distrecv)
      else if (layer_size(k) .ge. 4) then
        call getabove_init(2, k, d1_distsend, d1_distrecv)
        call getabove_init(2, k, d2_distsend, d2_distrecv)
        call getabove_init(3, k, d1up_distsend, d1up_distrecv)
      end if

!---------------
! Interior layer
!---------------
    else
      call getabove_init(2, k, d1_distsend, d1_distrecv)
      call getabove_init(2, k, d2_distsend, d2_distrecv)
      call getabove_init(3, k, d1up_distsend, d1up_distrecv)
      call getbelow_init(2, k, d1_distsend, d1_distrecv)
      call getbelow_init(2, k, d2_distsend, d2_distrecv)
      call getbelow_init(3, k, d1up_distsend, d1up_distrecv)
    end if
    end do
  end if

#endif

!---------------------------------------------
! Define constants for derivative calculations
!---------------------------------------------
  if ( nz .le. 1 ) then
    ho12 = 1.0d0
    sc1  = 1.0d0
    iret = iret - 4
  else 
    ho12 = dble(nz-1)/12.0d0
    sc1  = dble((nz-1)*(nz-1))/12.0d0
  end if

end subroutine init_dbydz

!=====================================
! Interior points for first derivative
!-------------------------------------
! Note: Error is O(h^4)
!=====================================
subroutine interior( fout,  fm2, fm1, fp1, fp2 )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!--------------------------------
! Loop over all values of x and y
!--------------------------------
  do j=1,ny
    do i=1,nx
      fout(i,j) = ho12*(         fm2(i,j)                  &
                         - 8.0d0*fm1(i,j)                  &
                         + 8.0d0*fp1(i,j)                  &
                         -       fp2(i,j)  )
    end do
  end do

end subroutine interior

!========================================
! Interior points for 1D first derivative
!----------------------------------------
! Note: Error is O(h^4)
!========================================
subroutine interior1D( fout,  fm2, fm1, fp1, fp2 )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  real(kind=dp), intent(in)  :: fm2
  real(kind=dp), intent(in)  :: fm1
  real(kind=dp), intent(in)  :: fp1
  real(kind=dp), intent(in)  :: fp2

  fout = ho12*( fm2-8.0d0*fm1+8.0d0*fp1-fp2  )

end subroutine interior1D

!=======================================
! Interior points for second derivatives
!---------------------------------------
! Note: Error is O(h^4)
!=======================================
subroutine d2interior( fout, fm2, fm1, f0, fp1, fp2 )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!--------------------------------
! Loop over all values of x and y
!--------------------------------
  do j=1,ny
    do i=1,nx
      fout(i,j) = sc1*(  -        fm2(i,j)                  &
                         + 16.0d0*fm1(i,j)                  &
                         - 30.0d0*f0(i,j)                   &
                         + 16.0d0*fp1(i,j)                  &
                         -        fp2(i,j)  )
    end do
  end do

end subroutine d2interior

!==========================================
! Interior points for 1D second derivatives
!------------------------------------------
! Note: Error is O(h^4)
!==========================================
subroutine d2interior1D( fout, fm2, fm1, f0, fp1, fp2 )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  real(kind=dp), intent(in)  :: fm2
  real(kind=dp), intent(in)  :: fm1
  real(kind=dp), intent(in)  :: f0
  real(kind=dp), intent(in)  :: fp1
  real(kind=dp), intent(in)  :: fp2

  fout = sc1*(  -fm2+16.0d0*fm1-30.0d0*f0+16.0d0*fp1-fp2  )

end subroutine d2interior1D
!======================================================
! Interior points for upwind derivatives (fourth order)
!======================================================
subroutine dupinterior(fout, fm3, fm2, fm1, f0, fp1, fp2, fp3, zpos)
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  real(kind=dp), intent(in)  :: fm3(nx, ny )
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )
  real(kind=dp), intent(in)  :: fp3(nx, ny )
  integer, intent(in)        :: zpos

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!-------------------------------------------------------
! Calculate one-sided derivatives depending on sign of w
!-------------------------------------------------------
  do j=1,ny
    do i=1,nx
      if (w(i,j,zpos) .ge. 0.0d0) then
        fout(i,j) = ho12*( -        fm3(i,j)               &
                           +  6.0d0*fm2(i,j)               &
                           - 18.0d0*fm1(i,j)               &
                           + 10.0d0*f0(i,j)                &
                           +  3.0d0*fp1(i,j) )
      else
        fout(i,j) = ho12*(          fp3(i,j)               &
                           -  6.0d0*fp2(i,j)               &
                           + 18.0d0*fp1(i,j)               &
                           - 10.0d0*f0(i,j)                &
                           -  3.0d0*fm1(i,j) )
      end if
    end do
  end do

end subroutine dupinterior

!==============================================================
! Boundary points for upwind derivatives (third point from top)
!==============================================================
subroutine dup_up(fout, fm2, fm1, f0, fp1, fp2, fp3, zpos)
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )
  real(kind=dp), intent(in)  :: fp3(nx, ny )
  integer, intent(in)        :: zpos

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!-------------------------------------------------------
! Calculate one-sided derivatives depending on sign of w
!
! Note: Lower accuracy near boundary 
!-------------------------------------------------------
  do j=1,ny
    do i=1,nx
      if (w(i,j,zpos) .ge. 0.0d0) then
        fout(i,j) = ho12*(    2.0d0*fm2(i,j)               &
                           - 12.0d0*fm1(i,j)               &
                           +  6.0d0*f0(i,j)                 &
                           +  4.0d0*fp1(i,j) )
      else
        fout(i,j) = ho12*(          fp3(i,j)               &
                           -  6.0d0*fp2(i,j)               &
                           + 18.0d0*fp1(i,j)               &
                           - 10.0d0*f0(i,j)                 &
                           -  3.0d0*fm1(i,j) )
      end if
    end do
  end do

end subroutine dup_up

!=================================================================
! Boundary points for upwind derivatives (third point from bottom)
!-----------------------------------------------------------------
subroutine dup_low(fout, fm3, fm2, fm1, f0, fp1, fp2, zpos)
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  real(kind=dp), intent(in)  :: fm3(nx, ny )
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )
  integer, intent(in)        :: zpos

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!-------------------------------------------------------
! Calculate one-sided derivatives depending on sign of w
!
! Note: Lower accuracy near boundary 
!-------------------------------------------------------
  do j=1,ny
    do i=1,nx
      if (w(i,j,zpos) .ge. 0.0d0) then
        fout(i,j) = ho12*( -        fm3(i,j)               &
                           +  6.0d0*fm2(i,j)               &
                           - 18.0d0*fm1(i,j)               &
                           + 10.0d0*f0(i,j)                &
                           +  3.0d0*fp1(i,j) )
      else
        fout(i,j) = ho12*( -  2.0d0*fp2(i,j)               &
                           + 12.0d0*fp1(i,j)               &
                           -  6.0d0*f0(i,j)                &
                           -  4.0d0*fm1(i,j) )
      end if
    end do
  end do

end subroutine dup_low

!=========================================================
! Bottom boundary layer for first derivative (lower point)
!=========================================================
subroutine lastdz(fout, f0, fm1, fm2, fm3, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm3(nx, ny )
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

  if ( iflag .eq. 0 ) then
!----------------------------------------------------------------
! Only the interior ponts are considered and error term is O(h^3)
! rather than O(h^4) - apparently this improves stability to 
! sound waves coming up through boundary as the error term is 
! diffusive rather than dispersive 
!----------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j) = 0.0d0
      end do
    end do

  else if ( iflag .eq. 1 ) then
!--------------------------------------------------------------------
! Conserves mass - The integral of z derivative over the whole of z
! is evaluated for objects that vanish at top and bottom (eg vertical
! velocity and momentum). This integral equals zero. Simpson's rule
! calculation gives coefficients. Interior points are
! taken care of by "interior" subroutine formula
! 
! Equivalently - this formula can be derived using a cubic 
! interpolation from the (zero) endpoints. This doesnt agree on the 
! endpoint value, but it is zero anyway (pjb 9/2/04).
! Interpolation formula is : c1 z + c2 z^2 + c3 z^3
!
! works out point above to third order as above (for stability)
!--------------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = ho12*(-  14.0d0*f0(i,j)                    &
                           -  36.0d0*fm1(i,j)                   &
                           +  18.0d0*fm2(i,j)                   &
                           -   4.0d0*fm3(i,j)       )
      end do
    end do

  else if ( iflag .eq. 2 ) then
!------------------------------------------------------------
! Third order interior point, outer point has zero derivative
!------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = 0.0d0
      end do
    end do

  else if ( iflag .eq. 5 ) then
!----------------
! All third order
!----------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = ho12*(   22.0d0*f0(i,j)                    &
                           -  36.0d0*fm1(i,j)                   &
                           +  18.0d0*fm2(i,j)                   &
                           -   4.0d0*fm3(i,j)       ) 
      end do
    end do

  else 
    print*,'UNKNOWN FLAG for boundary boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine lastdz


!============================================================
! Bottom boundary layer for 1D first derivative (lower point)
!============================================================
subroutine lastdz1D(fout, f0, fm1, fm2, fm3, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm3
  real(kind=dp), intent(in)  :: fm2
  real(kind=dp), intent(in)  :: fm1
  real(kind=dp), intent(in)  :: f0

  if ( iflag .eq. 0 ) then
!----------------------------------------------------------------
! Only the interior ponts are considered and error term is O(h^3)
! rather than O(h^4) - apparently this improves stability to 
! sound waves coming up through boundary as the error term is 
! diffusive rather than dispersive 
!----------------------------------------------------------------
    fout = 0.0d0

  else if ( iflag .eq. 1 ) then
!--------------------------------------------------------------------
! Conserves mass - The integral of z derivative over the whole of z
! is evaluated for objects that vanish at top and bottom (eg vertical
! velocity and momentum). This integral equals zero. Simpson's rule
! calculation gives coefficients. Interior points are
! taken care of by "interior" subroutine formula
! 
! Equivalently - this formula can be derived using a cubic 
! interpolation from the (zero) endpoints. This doesnt agree on the 
! endpoint value, but it is zero anyway (pjb 9/2/04).
! Interpolation formula is : c1 z + c2 z^2 + c3 z^3
!
! works out point above to third order as above (for stability)
!--------------------------------------------------------------------
    fout = ho12*( -14.0d0*f0-36.0d0*fm1+18.0d0*fm2-4.0d0*fm3 )

  else if ( iflag .eq. 2 ) then
!------------------------------------------------------------
! Third order interior point, outer point has zero derivative
!------------------------------------------------------------
    fout = 0.0d0

  else if ( iflag .eq. 5 ) then
!----------------
! All third order
!----------------
    fout = ho12*( 22.0d0*f0-36.0d0*fm1+18.0d0*fm2-4.0d0*fm3 ) 

  else 
    print*,'UNKNOWN FLAG for boundary boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine lastdz1D

!=========================================================
! Bottom boundary layer for first derivative (upper point)
!=========================================================
subroutine lastdz_up(fout, f0, fm1, fm2, fm3, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm3(nx, ny )
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

  if ( iflag .eq. 0 ) then
!----------------------------------------------------------------
! Only the interior ponts are considered and error term is O(h^3)
! rather than O(h^4) - apparently this improves stability to 
! sound waves coming up through boundary as the error term is 
! diffusive rather than dispersive 
!----------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = ho12*(    4.0d0*f0(i,j)                    &
                           +   6.0d0*fm1(i,j)                   &
                           -  12.0d0*fm2(i,j)                   &
                           +   2.0d0*fm3(i,j)       )
      end do
    end do
  else if ( iflag .eq. 1 ) then
!--------------------------------------------------------------------
! Conserves mass - The integral of z derivative over the whole of z
! is evaluated for objects that vanish at top and bottom (eg vertical
! velocity and momentum). This integral equals zero. Simpson's rule
! calculation gives coefficients. Interior points are
! taken care of by "interior" subroutine formula
! 
! Equivalently - this formula can be derived using a cubic 
! interpolation from the (zero) endpoints. This doesnt agree on the 
! endpoint value, but it is zero anyway (pjb 9/2/04).
! Interpolation formula is : c1 z + c2 z^2 + c3 z^3
!
! works out point above to third order as above (for stability)
!--------------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = ho12*(    4.0d0*f0(i,j)                    &
                           +   6.0d0*fm1(i,j)                   &
                           -  12.0d0*fm2(i,j)                   &
                           +   2.0d0*fm3(i,j)       )
      end do
    end do
  else if ( iflag .eq. 2 ) then
!------------------------------------------------------------
! Third order interior point, outer point has zero derivative
!------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = ho12*(    4.0d0*f0(i,j)                    &
                           +   6.0d0*fm1(i,j)                   &
                           -  12.0d0*fm2(i,j)                   &
                           +   2.0d0*fm3(i,j)       )
      end do
    end do
  else if ( iflag .eq. 5 ) then
!----------------
! All third order
!----------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = ho12*(    4.0d0*f0(i,j)                    &
                           +   6.0d0*fm1(i,j)                   &
                           -  12.0d0*fm2(i,j)                   &
                           +   2.0d0*fm3(i,j)       )
      end do
    end do
  else 
    print*,'UNKNOWN FLAG for boundary boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine lastdz_up

!============================================================
! Bottom boundary layer for 1D first derivative (upper point)
!============================================================
subroutine lastdz_up1D(fout, f0, fm1, fm2, fm3, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm3
  real(kind=dp), intent(in)  :: fm2
  real(kind=dp), intent(in)  :: fm1
  real(kind=dp), intent(in)  :: f0

  if ( iflag .eq. 0 ) then
!----------------------------------------------------------------
! Only the interior ponts are considered and error term is O(h^3)
! rather than O(h^4) - apparently this improves stability to 
! sound waves coming up through boundary as the error term is 
! diffusive rather than dispersive 
!----------------------------------------------------------------
    fout  = ho12*( 4.0d0*f0+6.0d0*fm1-12.0d0*fm2+2.0d0*fm3 )

  else if ( iflag .eq. 1 ) then
!--------------------------------------------------------------------
! Conserves mass - The integral of z derivative over the whole of z
! is evaluated for objects that vanish at top and bottom (eg vertical
! velocity and momentum). This integral equals zero. Simpson's rule
! calculation gives coefficients. Interior points are
! taken care of by "interior" subroutine formula
! 
! Equivalently - this formula can be derived using a cubic 
! interpolation from the (zero) endpoints. This doesnt agree on the 
! endpoint value, but it is zero anyway (pjb 9/2/04).
! Interpolation formula is : c1 z + c2 z^2 + c3 z^3
!
! works out point above to third order as above (for stability)
!--------------------------------------------------------------------
    fout = ho12*( 4.0d0*f0+6.0d0*fm1-12.0d0*fm2+2.0d0*fm3 )

  else if ( iflag .eq. 2 ) then
!------------------------------------------------------------
! Third order interior point, outer point has zero derivative
!------------------------------------------------------------
    fout = ho12*( 4.0d0*f0+6.0d0*fm1-12.0d0*fm2+2.0d0*fm3 )

  else if ( iflag .eq. 5 ) then
!----------------
! All third order
!----------------
    fout = ho12*( 4.0d0*f0+6.0d0*fm1-12.0d0*fm2+2.0d0*fm3 )

  else 
    print*,'UNKNOWN FLAG for boundary boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine lastdz_up1D

!======================================================
! Top boundary layer for first derivative (upper point)
!======================================================
subroutine topdz( fout, f1, f2, f3, f4, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f1(nx, ny )
  real(kind=dp), intent(in)  :: f2(nx, ny )
  real(kind=dp), intent(in)  :: f3(nx, ny )
  real(kind=dp), intent(in)  :: f4(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!----------------------------------------------------------
! fout is the layer derivative
! f1, f2, f3, f4 are fin(top), fin(p1), fin(p2) and fin(p3) 
!----------------------------------------------------------
  if ( iflag .eq. 0 ) then
!-------------
! top deriv =0
!-------------
    do j=1,ny
      do i=1,nx
        fout(i,j) =   0.0d0
      end do
    end do

  else if ( iflag .eq. 1 ) then
!---------------------------------------------------------
! Conserving mass again at boundary, third order otherwise
! Or equivalently cubic interpolation - see lastdz
!---------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j) =   ho12*(   14.0d0*f1(i,j)          &
                             + 36.0d0*f2(i,j)          &
                             - 18.0d0*f3(i,j)          &
                             +  4.0d0*f4(i,j)    )
      end do
    end do

  else if ( iflag .eq. 2 ) then
    do j=1,ny
      do i=1,nx
        fout(i,j) =  0.0d0
      end do
    end do

  else if ( iflag .eq. 5 ) then
!--------------------
! one-sided 3rd order 
!--------------------
    do j=1,ny
      do i=1,nx
        fout(i,j) =   ho12*( - 22.0d0*f1(i,j)          &
                             + 36.0d0*f2(i,j)          &
                             - 18.0d0*f3(i,j)          &
                             +  4.0d0*f4(i,j)  )
      end do
    end do
  else 

    print*,'UNKNOWN FLAG for top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine topdz

!=========================================================
! Top boundary layer for 1D first derivative (upper point)
!=========================================================
subroutine topdz1D( fout, f1, f2, f3, f4, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f1
  real(kind=dp), intent(in)  :: f2
  real(kind=dp), intent(in)  :: f3
  real(kind=dp), intent(in)  :: f4

!----------------------------------------------------------
! fout is the layer derivative
! f1, f2, f3, f4 are fin(top), fin(p1), fin(p2) and fin(p3) 
!----------------------------------------------------------
  if ( iflag .eq. 0 ) then
!-------------
! top deriv =0
!-------------
    fout = 0.0d0

  else if ( iflag .eq. 1 ) then
!---------------------------------------------------------
! Conserving mass again at boundary, third order otherwise
! Or equivalently cubic interpolation - see lastdz
!---------------------------------------------------------
    fout = ho12*( 14.0d0*f1+36.0d0*f2-18.0d0*f3+4.0d0*f4 )

  else if ( iflag .eq. 2 ) then
    fout = 0.0d0

  else if ( iflag .eq. 5 ) then
!--------------------
! one-sided 3rd order 
!--------------------
    fout = ho12*( -22.0d0*f1+36.0d0*f2-18.0d0*f3+4.0d0*f4 )

  else 

    print*,'UNKNOWN FLAG for top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine topdz1D

!======================================================
! Top boundary layer for first derivative (lower point)
!======================================================
subroutine topdz_low( fout, f1, f2, f3, f4, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f1(nx, ny )
  real(kind=dp), intent(in)  :: f2(nx, ny )
  real(kind=dp), intent(in)  :: f3(nx, ny )
  real(kind=dp), intent(in)  :: f4(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

!----------------------------------------------------------
! fout is the layer derivative
! f1, f2, f3, f4 are fin(top), fin(p1), fin(p2) and fin(p3) 
!----------------------------------------------------------
  if ( iflag .eq. 0 ) then
!---------------------------------
! derivative is set to third order
!---------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j) =   ho12*( -  4.0d0*f1(i,j)          &
                             -  6.0d0*f2(i,j)          &
                             + 12.0d0*f3(i,j)          &
                             -  2.0d0*f4(i,j)  )
      end do
    end do
  else if ( iflag .eq. 1 ) then
!---------------------------------------------------------
! Conserving mass again at boundary, third order otherwise
! Or equivalently cubic interpolation - see lastdz
!---------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j) =   ho12*( -  4.0d0*f1(i,j)          &
                             -  6.0d0*f2(i,j)          &
                             + 12.0d0*f3(i,j)          &
                             -  2.0d0*f4(i,j)  )
      end do
    end do
  else if ( iflag .eq. 2 ) then
    do j=1,ny
      do i=1,nx
        fout(i,j) =   ho12*( -  4.0d0*f1(i,j)          &
                             -  6.0d0*f2(i,j)          &
                             + 12.0d0*f3(i,j)          &
                             -  2.0d0*f4(i,j)  )
      end do
    end do
  else if ( iflag .eq. 5 ) then
!--------------------
! one-sided 3rd order 
!--------------------
    do j=1,ny
      do i=1,nx
        fout(i,j) =   ho12*( -  4.0d0*f1(i,j)          &
                             -  6.0d0*f2(i,j)          &
                             + 12.0d0*f3(i,j)          &
                             -  2.0d0*f4(i,j)  )
      end do
    end do
  else 
    print*,'UNKNOWN FLAG for top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine topdz_low

!=========================================================
! Top boundary layer for 1D first derivative (lower point)
!=========================================================
subroutine topdz_low1D( fout, f1, f2, f3, f4, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f1
  real(kind=dp), intent(in)  :: f2
  real(kind=dp), intent(in)  :: f3
  real(kind=dp), intent(in)  :: f4

!----------------------------------------------------------
! fout is the layer derivative
! f1, f2, f3, f4 are fin(top), fin(p1), fin(p2) and fin(p3) 
!----------------------------------------------------------
  if ( iflag .eq. 0 ) then
!---------------------------------
! derivative is set to third order
!---------------------------------
    fout = ho12*( -4.0d0*f1-6.0d0*f2+12.0d0*f3-2.0d0*f4 )

  else if ( iflag .eq. 1 ) then
!---------------------------------------------------------
! Conserving mass again at boundary, third order otherwise
! Or equivalently cubic interpolation - see lastdz
!---------------------------------------------------------
    fout = ho12*( -4.0d0*f1-6.0d0*f2+12.0d0*f3-2.0d0*f4 )

  else if ( iflag .eq. 2 ) then
    fout = ho12*( -4.0d0*f1-6.0d0*f2+12.0d0*f3-2.0d0*f4 )

  else if ( iflag .eq. 5 ) then
!--------------------
! one-sided 3rd order 
!--------------------
    fout = ho12*( -4.0d0*f1-6.0d0*f2+12.0d0*f3-2.0d0*f4 )

  else 
    print*,'UNKNOWN FLAG for top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine topdz_low1D

!=========================================================
! Second derivatives of the bottom boundary layers (lower)
!=========================================================
subroutine d2lastdz(fout, f0, fm1, fm2, fm3, fm4, fm5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm5(nx, ny )
  real(kind=dp), intent(in)  :: fm4(nx, ny )
  real(kind=dp), intent(in)  :: fm3(nx, ny )
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

  if ( iflag .eq. 0 ) then
!-----------------------------------------------------------------
! Set outer value to zero and fourth order formula for inner value
!-----------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = 0.0d0
      end do
    end do
  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = sc1*(     45.0d0*f0(i,j)                     &
                           -  154.0d0*fm1(i,j)                   &
                           +  214.0d0*fm2(i,j)                   &
                           -  156.0d0*fm3(i,j)                   &
                           +   61.0d0*fm4(i,j)                   &
                           -   10.0d0*fm5(i,j)      )
      end do
    end do
  else 
    print*,'UNKNOWN FLAG for 2nd derivative on bottom boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2lastdz

!============================================================
! 1D Second derivatives of the bottom boundary layers (lower)
!============================================================
subroutine d2lastdz1D(fout, f0, fm1, fm2, fm3, fm4, fm5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm5
  real(kind=dp), intent(in)  :: fm4
  real(kind=dp), intent(in)  :: fm3
  real(kind=dp), intent(in)  :: fm2
  real(kind=dp), intent(in)  :: fm1
  real(kind=dp), intent(in)  :: f0

  if ( iflag .eq. 0 ) then
!-----------------------------------------------------------------
! Set outer value to zero and fourth order formula for inner value
!-----------------------------------------------------------------
    fout = 0.0d0

  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    fout = sc1*( 45.0d0*f0-154.0d0*fm1+214.0d0*fm2                &
                -156.0d0*fm3+61.0d0*fm4-10.0d0*fm5 )

  else 
    print*,'UNKNOWN FLAG for 2nd derivative on bottom boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2lastdz1D

!=========================================================
! Second derivatives of the bottom boundary layers (upper)
!=========================================================
subroutine d2lastdz_up(fout, f0, fm1, fm2, fm3, fm4, fm5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm5(nx, ny )
  real(kind=dp), intent(in)  :: fm4(nx, ny )
  real(kind=dp), intent(in)  :: fm3(nx, ny )
  real(kind=dp), intent(in)  :: fm2(nx, ny )
  real(kind=dp), intent(in)  :: fm1(nx, ny )
  real(kind=dp), intent(in)  :: f0(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

  if ( iflag .eq. 0 ) then
!-----------------------------------------------------------------
! Set outer value to zero and fourth order formula for inner value
!-----------------------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = sc1*(    10.0d0*f0(i,j)                      &
                           -  15.0d0*fm1(i,j)                    &
                           -   4.0d0*fm2(i,j)                    &
                           +  14.0d0*fm3(i,j)                    &
                           -   6.0d0*fm4(i,j)                    &
                           +         fm5(i,j)       )
      end do
    end do
  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = sc1*(    10.0d0*f0(i,j)                      &
                           -  15.0d0*fm1(i,j)                    &
                           -   4.0d0*fm2(i,j)                    &
                           +  14.0d0*fm3(i,j)                    &
                           -   6.0d0*fm4(i,j)                    &
                           +         fm5(i,j)       )
      end do
    end do
  else 
    print*,'UNKNOWN FLAG for 2nd derivative on bottom boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2lastdz_up

!============================================================
! 1D Second derivatives of the bottom boundary layers (upper)
!============================================================
subroutine d2lastdz_up1D(fout, f0, fm1, fm2, fm3, fm4, fm5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: fm5
  real(kind=dp), intent(in)  :: fm4
  real(kind=dp), intent(in)  :: fm3
  real(kind=dp), intent(in)  :: fm2
  real(kind=dp), intent(in)  :: fm1
  real(kind=dp), intent(in)  :: f0

  if ( iflag .eq. 0 ) then
!-----------------------------------------------------------------
! Set outer value to zero and fourth order formula for inner value
!-----------------------------------------------------------------
    fout = sc1*( 10.0d0*f0-15.0d0*fm1-4.0d0*fm2+14.0d0*fm3       &
                -6.0d0*fm4+fm5 )

  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    fout = sc1*( 10.0d0*f0-15.0d0*fm1-4.0d0*fm2+14.0d0*fm3       &
                -6.0d0*fm4+fm5 )

  else 
    print*,'UNKNOWN FLAG for 2nd derivative on bottom boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2lastdz_up1D

!======================================================
! Second derivatives of the top boundary layers (upper)
!======================================================
subroutine d2topdz(fout, f0, fp1, fp2, fp3, fp4, fp5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f0(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )
  real(kind=dp), intent(in)  :: fp3(nx, ny )
  real(kind=dp), intent(in)  :: fp4(nx, ny )
  real(kind=dp), intent(in)  :: fp5(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

  if ( iflag .eq. 0 ) then
!--------------------------------------------------
! set top to 0.0 and other is fourth order accurate
!
! Only the interior points are considered 
!--------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = 0.0d0
      end do
    end do
  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = sc1*(     45.0d0*f0(i,j)                     &
                           -  154.0d0*fp1(i,j)                   &
                           +  214.0d0*fp2(i,j)                   &
                           -  156.0d0*fp3(i,j)                   &
                           +   61.0d0*fp4(i,j)                   &
                           -   10.0d0*fp5(i,j)       )
      end do
    end do
  else 
    print*,'UNKNOWN FLAG for 2nd derivative on top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2topdz
!=========================================================
! 1D Second derivatives of the top boundary layers (upper)
!=========================================================
subroutine d2topdz1D(fout, f0, fp1, fp2, fp3, fp4, fp5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f0
  real(kind=dp), intent(in)  :: fp1
  real(kind=dp), intent(in)  :: fp2
  real(kind=dp), intent(in)  :: fp3
  real(kind=dp), intent(in)  :: fp4
  real(kind=dp), intent(in)  :: fp5

  if ( iflag .eq. 0 ) then
!--------------------------------------------------
! set top to 0.0 and other is fourth order accurate
!
! Only the interior points are considered 
!--------------------------------------------------
    fout = 0.0d0

  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    fout = sc1*( 45.0d0*f0-154.0d0*fp1+214.0d0*fp2               &
                -156.0d0*fp3+61.0d0*fp4-10.0d0*fp5 )

  else 
    print*,'UNKNOWN FLAG for 2nd derivative on top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2topdz1D

!======================================================
! Second derivatives of the top boundary layers (lower)
!======================================================
subroutine d2topdz_low(fout, f0, fp1, fp2, fp3, fp4, fp5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout(nx, ny)
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f0(nx, ny )
  real(kind=dp), intent(in)  :: fp1(nx, ny )
  real(kind=dp), intent(in)  :: fp2(nx, ny )
  real(kind=dp), intent(in)  :: fp3(nx, ny )
  real(kind=dp), intent(in)  :: fp4(nx, ny )
  real(kind=dp), intent(in)  :: fp5(nx, ny )

!-------------------
! Local declarations
!-------------------
  integer                    :: i, j

  if ( iflag .eq. 0 ) then
!--------------------------------------------------
! set top to 0.0 and other is fourth order accurate
!
! Only the interior points are considered 
!--------------------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = sc1*(    10.0d0*f0(i,j)                      &
                           -  15.0d0*fp1(i,j)                    &
                           -   4.0d0*fp2(i,j)                    &
                           +  14.0d0*fp3(i,j)                    &
                           -   6.0d0*fp4(i,j)                    &
                           +         fp5(i,j)       )
      end do
    end do
  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    do j=1,ny
      do i=1,nx
        fout(i,j)  = sc1*(    10.0d0*f0(i,j)                      &
                           -  15.0d0*fp1(i,j)                    &
                           -   4.0d0*fp2(i,j)                    &
                           +  14.0d0*fp3(i,j)                    &
                           -   6.0d0*fp4(i,j)                    &
                           +         fp5(i,j)       )
      end do
    end do
  else 
    print*,'UNKNOWN FLAG for 2nd derivative on top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2topdz_low

!=========================================================
! 1D Second derivatives of the top boundary layers (lower)
!=========================================================
subroutine d2topdz_low1D(fout, f0, fp1, fp2, fp3, fp4, fp5, iflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=dp), intent(out) :: fout
  integer,       intent(in)  :: iflag
  real(kind=dp), intent(in)  :: f0
  real(kind=dp), intent(in)  :: fp1
  real(kind=dp), intent(in)  :: fp2
  real(kind=dp), intent(in)  :: fp3
  real(kind=dp), intent(in)  :: fp4
  real(kind=dp), intent(in)  :: fp5

  if ( iflag .eq. 0 ) then
!--------------------------------------------------
! set top to 0.0 and other is fourth order accurate
!
! Only the interior points are considered 
!--------------------------------------------------
    fout = sc1*( 10.0d0*f0-15.0d0*fp1-4.0d0*fp2                   &
                +14.0d0*fp3-6.0d0*fp4+fp5 )

  else if ( iflag .eq. 3 ) then
!------------------------------------
! fourth order one-sided formulations
!------------------------------------
    fout = sc1*( 10.0d0*f0-15.0d0*fp1-4.0d0*fp2                   &
                +14.0d0*fp3-6.0d0*fp4+fp5 )

  else 
    print*,'UNKNOWN FLAG for 2nd derivative on top boundary ', iflag
    call abort_comms(7)
  end if 

end subroutine d2topdz_low1D

!====================================================================
! <<<< d1bydz >>>>
!
! Sets Fout to the first derivative in the Z direction of fin
!
! If DFLAG=0 only interior points are considered.
! If DFLAG=1 then the derivative is evaluated at z=0,1 also,
!            using conservation of mass. (for d RW /dz dnd dW/dz)
!            For an quantity which is 0 on boundaries, integral of
!            deriv. = 0 - equivalently cubic interpolation
! If DFLAG=2 the boundary derivatives are set to zero for (dUdz)
!            quartic interpolation used for inner point
! If DFLAG=5 3rd order sideways differencing is used at the boundary.
!
! Note: DFLAG=3 and DFLAG=4 used to be used for magnetic boundary
! conditions - now obsolete.
!--------------------------------------------------------------------
! As the data is layed out in slabs the neighbouring layers have
! to be communicated to this layer.
!
!
!        myrank-1      rubfm
!                    ---------
!        myrank        layer
!                    ---------
!        myrank+1      rubfp
!
!===================================================================
Subroutine d1bydz( fin, fout, dflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=8), intent(in)      :: fin(nx,ny,nlayer)
  integer,      intent(in)      :: dflag
  real(kind=8), intent(out)     :: fout(nx,ny,nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: k
#ifdef MPI
  integer                       :: ier
#endif

!----------------------------
! Starts timing off if needed
!----------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

#ifdef MPI
!------------------------------
! Do communication if necessary
!------------------------------
  if (nproc .gt. 1) then

!-----------------------
! Initialise all buffers
!-----------------------
    rbufm=0.0d0
    rbufp=0.0d0
    msg = 0
    mstat = 0
    ireq = 0 
    rbufp_counter = 0 
    rbufm_counter = 0

!----------------------------
! Post receive commands first
!----------------------------
    do k=1, myrank
      if (d1_distrecv(k) .gt. 0) then
        ncount = nx*ny*d1_distrecv(k)
        call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,         &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        senddown_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufm_counter = rbufm_counter + d1_distrecv(k)
      end if
    end do

    do k=myrank+2,nproc
      if (d1_distrecv(k) .gt. 0) then
        ncount = nx*ny*d1_distrecv(k)
        call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,         &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        sendup_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufp_counter = rbufp_counter + d1_distrecv(k)
      end if
    end do    

!----------------------------
! Carry out MPI send commands
!----------------------------
    do k=1, myrank
      if (d1_distsend(k) .gt. 0) then
        ncount = nx*ny*d1_distsend(k)
        call MPI_Isend( fin(1,1,1), ncount,                         &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        sendup_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

    do k=myrank+2,nproc
      if (d1_distsend(k) .gt. 0) then
        ncount = nx*ny*d1_distsend(k)
        call MPI_Isend( fin(1,1,nlayer-d1_distsend(k)+1), ncount,   &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        senddown_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

  if (nlayer .ge. 4) then

!--------------------------------------------
! First do the simplest case, for nlayer >= 4
!--------------------------------------------
! Do the top two layers on a node
!--------------------------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               fin(1,1,3), fin(1,1,4), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               fin(1,1,3), fin(1,1,4), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                               fin(1,1,2), fin(1,1,3), dflag )
      call interior( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),        &
                               fin(1,1,3), fin(1,1,4) )
#endif
    else 
#ifdef MPI
      call interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                               fin(1,1,2),   fin(1,1,3 ) )
      call interior( fout(1,1,2), rbufm(1,1,2), fin(1,1,1),        &
                               fin(1,1,3), fin(1,1,4) ) 
#endif
    end if

!--------------------------------
! Do the middle layers on a node
!--------------------------------
    if (nlayer .ge. 5) then
      do k=3, nlayer-2
        call interior( fout(1,1,k), fin(1,1,k-2), fin(1,1,k-1),    &
                               fin(1,1,k+1), fin(1,1,k+2) )
      end do
    end if

!-----------------------------------
! Do the bottom two layers on a node
!-----------------------------------
    if ( izfinish .eq. nz ) then
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), fin(1,1,nlayer-2),        &
                      fin(1,1,nlayer-3), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                   fin(1,1,nlayer-3), dflag )
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call interior( fout(1,1,nlayer-1), fin(1,1,nlayer-3),        & 
                     fin(1,1,nlayer-2), fin(1,1,nlayer),           &
                     rbufp(1,1,1) )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      fin(1,1,nlayer-2), dflag )
#endif
    else 
#ifdef MPI
     call interior( fout(1,1,nlayer-1),fin(1,1,nlayer-3),          &
                    fin(1,1,nlayer-2), fin(1,1,nlayer),            &
                    rbufp(1,1,1) )
     call interior( fout(1,1,nlayer),  fin(1,1,nlayer-2),          &
                    fin(1,1,nlayer-1), rbufp(1,1,1),               & 
                    rbufp(1,1,2) )
#endif
    end if

#ifdef MPI

  else if ( nlayer .eq. 3) then

!------------------------------------
! Now do the nlayer=3 case
!------------------------------------
! Start with top two layers on a node
!------------------------------------
    if ( izstart .eq. 1 ) then
     call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                              fin(1,1,3), rbufp(1,1,1), dflag )
     call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                              fin(1,1,3), rbufp(1,1,1), dflag )
     call interior( fout(1,1,3), fin(1,1,1), fin(1,1,2),          &
                              rbufp(1,1,1), rbufp(1,1,2) )
   else if (izstart .eq. 2) then
#ifdef MPI
     call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                              fin(1,1,2), fin(1,1,3), dflag )
     call interior( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),        &
                              fin(1,1,3), rbufp(1,1,1) )
     call interior( fout(1,1,3), fin(1,1,1), fin(1,1,2),          &
                              rbufp(1,1,1), rbufp(1,1,2) )
#endif
   else if (izfinish .eq. nz) then
     call interior( fout(1,1,nlayer-2),rbufm(1,1,1),rbufm(1,1,2), &
                    fin(1,1,nlayer-1), fin(1,1,nlayer) )
     call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                  fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                  rbufm(1,1,2), dflag )
     call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                  fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                  rbufm(1,1,2), dflag )
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call interior( fout(1,1,nlayer-2),rbufm(1,1,1),rbufm(1,1,2), &
                     fin(1,1,nlayer-1), fin(1,1,nlayer) )
      call interior( fout(1,1,nlayer-1),rbufm(1,1,2),              &
                     fin(1,1,nlayer-2),fin(1,1,nlayer),rbufp(1,1,1) )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      fin(1,1,nlayer-2), dflag )
#endif
    else 
#ifdef MPI
      call interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                                  fin(1,1,2),   fin(1,1,3 ) )
      call interior( fout(1,1,2), rbufm(1,1,2), fin(1,1,1),        &
                                  fin(1,1,3), rbufp(1,1,1) ) 
      call interior( fout(1,1,3), fin(1,1,1), fin(1,1,2),          &
                                  rbufp(1,1,1), rbufp(1,1,2) )
#endif
    end if

  else if ( nlayer .eq. 2) then

!-------------------------
! Now do the nlayer=2 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               rbufp(1,1,1), rbufp(1,1,2), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               rbufp(1,1,1), rbufp(1,1,2), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                              fin(1,1,2), rbufp(1,1,1), dflag )

      call interior( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),        &
                              rbufp(1,1,1), rbufp(1,1,2) )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call interior( fout(1,1,nlayer-1),rbufm(1,1,1),              &
                     rbufm(1,1,2),fin(1,1,nlayer),rbufp(1,1,1) )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      rbufm(1,1,2), dflag )
#endif
    else if (izfinish .eq. nz) then
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), rbufm(1,1,2),             &
                      rbufm(1,1,1), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), rbufm(1,1,2),                &
                   rbufm(1,1,1), dflag )

    else
      call interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                                  fin(1,1,2),   rbufp(1,1,1) )
      call interior( fout(1,1,2), rbufm(1,1,2), fin(1,1,1),        &
                                  rbufp(1,1,1), rbufp(1,1,2) ) 
    end if

  else if ( nlayer .eq. 1) then

!-------------------------
! Now do the nlayer=1 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), rbufp(1,1,1),           &
                               rbufp(1,1,2), rbufp(1,1,3), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                               rbufp(1,1,1), rbufp(1,1,2), dflag )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), rbufm(1,1,2),               &
                      rbufm(1,1,1), dflag )
#endif
    else if (izfinish .eq. nz) then
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   rbufm(1,1,3), rbufm(1,1,2),                     &
                   rbufm(1,1,1), dflag )

    else
      call interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                                  rbufp(1,1,1), rbufp(1,1,2) )
    end if
#endif
  else
    print*,' not set up for case nlayer = ',nlayer
    call abort_comms(2)
  end if

!-------------------------------------------------------
! This measures the time taken to get to this stage from
! the start of the subroutine and stores as time_dz
!-------------------------------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_dz = time_dz + (time_tmp2-time_tmp)
  end if

end subroutine d1bydz

!====================================================================
! <<<< d1bydz1D >>>>
!
! Sets Fout to the 1D first derivative in the Z direction of fin
!
! If DFLAG=0 only interior points are considered.
! If DFLAG=1 then the derivative is evaluated at z=0,1 also,
!            using conservation of mass. (for d RW /dz dnd dW/dz)
!            For an quantity which is 0 on boundaries, integral of
!            deriv. = 0 - equivalently cubic interpolation
! If DFLAG=2 the boundary derivatives are set to zero for (dUdz)
!            quartic interpolation used for inner point
! If DFLAG=5 3rd order sideways differencing is used at the boundary.
!
! Note: DFLAG=3 and DFLAG=4 used to be used for magnetic boundary
! conditions - now obsolete.
!--------------------------------------------------------------------
! As the data is layed out in slabs the neighbouring layers have
! to be communicated to this layer.
!
!
!        myrank-1      rubfm
!                    ---------
!        myrank        layer
!                    ---------
!        myrank+1      rubfp
!
!===================================================================
Subroutine d1bydz1D( fin, fout, dflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=8), intent(in)      :: fin(nlayer)
  integer,      intent(in)      :: dflag
  real(kind=8), intent(out)     :: fout(nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: k
#ifdef MPI
  integer                       :: ier
#endif

!----------------------------
! Starts timing off if needed
!----------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

#ifdef MPI
!------------------------------
! Do communication if necessary
!------------------------------
  if (nproc .gt. 1) then

!-----------------------
! Initialise all buffers
!-----------------------
    recvbuf=0.0d0
    sendbuf=0.0d0
    msg = 0
    mstat = 0
    ireq = 0 
    rbufp_counter = 0 
    rbufm_counter = 0

!----------------------------
! Post receive commands first
!----------------------------
    do k=1, myrank
      if (d1_distrecv(k) .gt. 0) then
        ncount = d1_distrecv(k)
        call MPI_Irecv( recvbuf(1+rbufm_counter), ncount,         &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        senddown_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufm_counter = rbufm_counter + d1_distrecv(k)
      end if
    end do

    do k=myrank+2,nproc
      if (d1_distrecv(k) .gt. 0) then
        ncount = d1_distrecv(k)
        call MPI_Irecv( sendbuf(1+rbufp_counter), ncount,         &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        sendup_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufp_counter = rbufp_counter + d1_distrecv(k)
      end if
    end do    

!----------------------------
! Carry out MPI send commands
!----------------------------
    do k=1, myrank
      if (d1_distsend(k) .gt. 0) then
        ncount = d1_distsend(k)
        call MPI_Isend( fin(1), ncount,                         &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        sendup_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

    do k=myrank+2,nproc
      if (d1_distsend(k) .gt. 0) then
        ncount = d1_distsend(k)
        call MPI_Isend( fin(nlayer-d1_distsend(k)+1), ncount,   &
                        MPI_DOUBLE_PRECISION, k-1,                  &
                        senddown_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

  if (nlayer .ge. 4) then

!--------------------------------------------
! First do the simplest case, for nlayer >= 4
!--------------------------------------------
! Do the top two layers on a node
!--------------------------------------------
    if ( izstart .eq. 1 ) then
      call topdz1D( fout(1), fin(1), fin(2),             &
                               fin(3), fin(4), dflag )
      call topdz_low1D( fout(2), fin(1), fin(2),         &
                               fin(3), fin(4), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call topdz_low1D( fout(1), recvbuf(1), fin(1),     &
                               fin(2), fin(3), dflag )
      call interior1D( fout(2), recvbuf(1), fin(1),      &
                               fin(3), fin(4) )
#endif
    else 
#ifdef MPI
      call interior1D( fout(1), recvbuf(1), recvbuf(2),  &
                               fin(2),   fin(3) )
      call interior1D( fout(2), recvbuf(2), fin(1),      &
                               fin(3), fin(4) ) 
#endif
    end if

!--------------------------------
! Do the middle layers on a node
!--------------------------------
    if (nlayer .ge. 5) then
      do k=3, nlayer-2
        call interior1D( fout(k), fin(k-2), fin(k-1),    &
                               fin(k+1), fin(k+2) )
      end do
    end if

!-----------------------------------
! Do the bottom two layers on a node
!-----------------------------------
    if ( izfinish .eq. nz ) then
      call lastdz_up1D( fout(nlayer-1), fin(nlayer),     & 
                      fin(nlayer-1), fin(nlayer-2),    &
                      fin(nlayer-3), dflag )
      call lastdz1D( fout(nlayer), fin(nlayer),          & 
                   fin(nlayer-1), fin(nlayer-2),       &
                   fin(nlayer-3), dflag )
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call interior1D( fout(nlayer-1), fin(nlayer-3),    & 
                     fin(nlayer-2), fin(nlayer),       &
                     sendbuf(1) )
      call lastdz_up1D( fout(nlayer), sendbuf(1),        & 
                      fin(nlayer), fin(nlayer-1),      &
                      fin(nlayer-2), dflag )
#endif
    else 
#ifdef MPI
     call interior1D( fout(nlayer-1),fin(nlayer-3),      &
                    fin(nlayer-2), fin(nlayer),        &
                    sendbuf(1) )
     call interior1D( fout(nlayer),  fin(nlayer-2),      &
                    fin(nlayer-1), sendbuf(1),         & 
                    sendbuf(2) )
#endif
    end if

#ifdef MPI

  else if ( nlayer .eq. 3) then

!------------------------------------
! Now do the nlayer=3 case
!------------------------------------
! Start with top two layers on a node
!------------------------------------
    if ( izstart .eq. 1 ) then
     call topdz1D( fout(1), fin(1), fin(2),              &
                              fin(3), sendbuf(1), dflag )
     call topdz_low1D( fout(2), fin(1), fin(2),          &
                              fin(3), sendbuf(1), dflag )
     call interior1D( fout(3), fin(1), fin(2),           &
                              sendbuf(1), sendbuf(2) )
   else if (izstart .eq. 2) then
#ifdef MPI
     call topdz_low1D( fout(1), recvbuf(1), fin(1),      &
                              fin(2), fin(3), dflag )
     call interior1D( fout(2), recvbuf(1), fin(1),       &
                              fin(3), sendbuf(1) )
     call interior1D( fout(3), fin(1), fin(2),           &
                              sendbuf(1), sendbuf(2) )
#endif
   else if (izfinish .eq. nz) then
     call interior1D( fout(nlayer-2),recvbuf(1),recvbuf(2), &
                    fin(nlayer-1), fin(nlayer) )
     call lastdz_up1D( fout(nlayer-1), fin(nlayer),         & 
                  fin(nlayer-1), fin(nlayer-2),           &
                  recvbuf(2), dflag )
     call lastdz1D( fout(nlayer), fin(nlayer),              & 
                  fin(nlayer-1), fin(nlayer-2),           &
                  recvbuf(2), dflag )
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call interior1D( fout(nlayer-2),recvbuf(1),recvbuf(2),   &
                     fin(nlayer-1), fin(nlayer) )
      call interior1D( fout(nlayer-1),recvbuf(2),              &
                     fin(nlayer-2),fin(nlayer),sendbuf(1) )
      call lastdz_up1D( fout(nlayer), sendbuf(1),              & 
                      fin(nlayer), fin(nlayer-1),            &
                      fin(nlayer-2), dflag )
#endif
    else 
#ifdef MPI
      call interior1D( fout(1), recvbuf(1), recvbuf(2),        &
                                  fin(2),   fin(3 ) )
      call interior1D( fout(2), recvbuf(2), fin(1),            &
                                  fin(3), sendbuf(1) ) 
      call interior1D( fout(3), fin(1), fin(2),                &
                                  sendbuf(1), sendbuf(2) )
#endif
    end if

  else if ( nlayer .eq. 2) then

!-------------------------
! Now do the nlayer=2 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz1D( fout(1), fin(1), fin(2),                   &
                               sendbuf(1), sendbuf(2), dflag )
      call topdz_low1D( fout(2), fin(1), fin(2),               &
                               sendbuf(1), sendbuf(2), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call topdz_low1D( fout(1), recvbuf(1), fin(1),           &
                              fin(2), sendbuf(1), dflag )

      call interior1D( fout(2), recvbuf(1), fin(1),            &
                              sendbuf(1), sendbuf(2) )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call interior1D( fout(nlayer-1),recvbuf(1),              &
                     recvbuf(2),fin(nlayer),sendbuf(1) )
      call lastdz_up1D( fout(nlayer), sendbuf(1),              & 
                      fin(nlayer), fin(nlayer-1),            &
                      recvbuf(2), dflag )
#endif
    else if (izfinish .eq. nz) then
      call lastdz_up1D( fout(nlayer-1), fin(nlayer),           & 
                      fin(nlayer-1), recvbuf(2),             &
                      recvbuf(1), dflag )
      call lastdz1D( fout(nlayer), fin(nlayer),                & 
                   fin(nlayer-1), recvbuf(2),                &
                   recvbuf(1), dflag )

    else
      call interior1D( fout(1), recvbuf(1), recvbuf(2),        &
                                  fin(2),   sendbuf(1) )
      call interior1D( fout(2), recvbuf(2), fin(1),            &
                                  sendbuf(1), sendbuf(2) ) 
    end if

  else if ( nlayer .eq. 1) then

!-------------------------
! Now do the nlayer=1 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz1D( fout(1), fin(1), sendbuf(1),               &
                               sendbuf(2), sendbuf(3), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call topdz_low1D( fout(1), recvbuf(1), fin(1),           &
                               sendbuf(1), sendbuf(2), dflag )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call lastdz_up1D( fout(nlayer), sendbuf(1),              & 
                      fin(nlayer), recvbuf(2),               &
                      recvbuf(1), dflag )
#endif
    else if (izfinish .eq. nz) then
      call lastdz1D( fout(nlayer), fin(nlayer),                & 
                   recvbuf(3), recvbuf(2),                   &
                   recvbuf(1), dflag )

    else
      call interior1D( fout(1), recvbuf(1), recvbuf(2),        &
                                  sendbuf(1), sendbuf(2) )
    end if
#endif
  else
    print*,' not set up for case nlayer = ',nlayer
    call abort_comms(2)
  end if

!-------------------------------------------------------
! This measures the time taken to get to this stage from
! the start of the subroutine and stores as time_dz
!-------------------------------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_dz = time_dz + (time_tmp2-time_tmp)
  end if

end subroutine d1bydz1D

!====================================================================
! <<<< d1upbydz >>>>
!
! Sets Fout to the first upwind derivative in the Z direction of fin
!
! If DFLAG=0 only interior points are considered.
! If DFLAG=1 then the derivative is evaluated at z=0,1 also,
!            using conservation of mass. (for d RW /dz dnd dW/dz)
!            For an quantity which is 0 on boundaries, integral of
!            deriv. = 0 - equivalently cubic interpolation
! If DFLAG=2 the boundary derivatives are set to zero for (dUdz)
!            quartic interpolation used for inner point
! If DFLAG=5 3rd order sideways differencing is used at the boundary.
!
! Note: DFLAG=3, DFLAG=4 used to be used for magnetic boundary
! conditions - now obsolete.
!--------------------------------------------------------------------
! As the data is layed out in slabs the neighbouring layers have
! to be communicated to this layer.
!
!
!        myrank-1      rubfm
!                    ---------
!        myrank        layer
!                    ---------
!        myrank+1      rubfp
!
!====================================================================
Subroutine d1upbydz( fin, fout, dflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=8), intent(in)      :: fin(nx,ny,nlayer)
  integer,      intent(in)      :: dflag
  real(kind=8), intent(out)     :: fout(nx,ny,nlayer)

!-------------------
! Local declarations
!-------------------
  integer                       :: k
#ifdef MPI 
  integer                       :: ier
#endif

!----------------------------
! Starts timing off if needed
!----------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

#ifdef MPI
!------------------------------
! Do communication if necessary
!------------------------------
  if (nproc .gt. 1) then

!-------------------
! Initialise buffers
!-------------------
    rbufm=0.0d0
    rbufp=0.0d0
    msg = 0
    mstat = 0
    ireq = 0 
    rbufp_counter = 0 
    rbufm_counter = 0

!----------------------------
! Post receive commands first
!----------------------------
    do k=1, myrank
      if (d1up_distrecv(k) .gt. 0) then
        ncount = nx*ny*d1up_distrecv(k)
        call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        senddown_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufm_counter = rbufm_counter + d1up_distrecv(k)
      end if
    end do

    do k=myrank+2,nproc
      if (d1up_distrecv(k) .gt. 0) then
        ncount = nx*ny*d1up_distrecv(k)
        call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        sendup_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufp_counter = rbufp_counter + d1up_distrecv(k)
      end if
    end do

!----------------------------
! Carry out MPI send commands
!----------------------------
    do k=1, myrank
      if (d1up_distsend(k) .gt. 0) then
        ncount = nx*ny*d1up_distsend(k)
        call MPI_Isend( fin(1,1,1), ncount,                           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        sendup_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

    do k=myrank+2,nproc    
      if (d1up_distsend(k) .gt. 0) then
        ncount = nx*ny*d1up_distsend(k)
        call MPI_Isend( fin(1,1,nlayer-d1up_distsend(k)+1), ncount,   &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        senddown_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

  if (nlayer .ge. 6) then

!--------------------------------------------
! First do the simplest case, for nlayer >= 6
!--------------------------------------------
! Start with the top three layers on a node
!--------------------------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               fin(1,1,3), fin(1,1,4), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               fin(1,1,3), fin(1,1,4), dflag )
      call dup_up( fout(1,1,3), fin(1,1,1), fin(1,1,2),            &
                 fin(1,1,3), fin(1,1,4), fin(1,1,5), fin(1,1,6), 3 )
#ifdef MPI
    else if (izstart .eq. 2) then
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                               fin(1,1,2), fin(1,1,3), dflag )
      call dup_up( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),          &
                   fin(1,1,2), fin(1,1,3), fin(1,1,4), fin(1,1,5), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,1), fin(1,1,1),     & 
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        fin(1,1,5), fin(1,1,6), 3 ) 
    else if (izstart .eq. 3) then
      call dup_up( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),        &
                   fin(1,1,1), fin(1,1,2), fin(1,1,3), fin(1,1,4), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,1), rbufm(1,1,2),   & 
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        fin(1,1,4), fin(1,1,5), 2 ) 
      call dupinterior( fout(1,1,3), rbufm(1,1,2), fin(1,1,1),     & 
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        fin(1,1,5), fin(1,1,6), 3 ) 
#endif
    else 
#ifdef MPI
      call dupinterior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),   & 
                        rbufm(1,1,3), fin(1,1,1), fin(1,1,2),      &
                        fin(1,1,3), fin(1,1,4), 1 ) 
      call dupinterior( fout(1,1,2), rbufm(1,1,2), rbufm(1,1,3),   & 
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        fin(1,1,4), fin(1,1,5), 2 ) 
      call dupinterior( fout(1,1,3), rbufm(1,1,3), fin(1,1,1),     & 
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        fin(1,1,5), fin(1,1,6), 3 ) 
#endif
    end if

!--------------------------------
! Do the middle layers on a node
!--------------------------------
    if (nlayer .ge. 7) then
      do k=4, nlayer-3
        call dupinterior( fout(1,1,k), fin(1,1,k-3), fin(1,1,k-2),&
                         fin(1,1,k-1), fin(1,1,k), fin(1,1,k+1),  &
                         fin(1,1,k+2), fin(1,1,k+3), k )
      end do
    end if

!-----------------------------------
! Do the bottom three layers on a node
!-----------------------------------
    if ( izfinish .eq. nz ) then
      call dup_low( fout(1,1,nlayer-2), fin(1,1,nlayer-5),         &
                    fin(1,1,nlayer-4), fin(1,1,nlayer-3),          &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), nlayer-2 )
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), fin(1,1,nlayer-2),        &
                      fin(1,1,nlayer-3), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                   fin(1,1,nlayer-3), dflag )
#ifdef MPI
    else if (izfinish .eq. nz-1) then
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      fin(1,1,nlayer-2), dflag )
      call dup_low( fout(1,1,nlayer-1), fin(1,1,nlayer-4),         &
                    fin(1,1,nlayer-3), fin(1,1,nlayer-2),          &
                    fin(1,1,nlayer-1), fin(1,1,nlayer),            &
                    rbufp(1,1,1), nlayer-1 ) 
      call dupinterior( fout(1,1,nlayer-2), fin(1,1,nlayer-5),     & 
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
    else if (izfinish .eq. nz-2) then
      call dup_low( fout(1,1,nlayer), fin(1,1,nlayer-3),           &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), rbufp(1,1,1), rbufp(1,1,2), nlayer)   
      call dupinterior( fout(1,1,nlayer-1), fin(1,1,nlayer-4),     &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer),        &
                        rbufp(1,1,1), rbufp(1,1,2), nlayer-1 ) 
      call dupinterior( fout(1,1,nlayer-2), fin(1,1,nlayer-5),     &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 ) 
#endif
    else 
#ifdef MPI
      call dupinterior( fout(1,1,nlayer-2),fin(1,1,nlayer-5),      &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 ) 
      call dupinterior( fout(1,1,nlayer-1),fin(1,1,nlayer-4),      &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer),        &
                        rbufp(1,1,1), rbufp(1,1,2), nlayer-1 ) 
      call dupinterior( fout(1,1,nlayer),fin(1,1,nlayer-3),        &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1),             &
                        rbufp(1,1,2), rbufp(1,1,3), nlayer ) 
#endif
    end if

#ifdef MPI

  else if ( nlayer .eq. 5) then

!-------------------------
! Now do the nlayer=5 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               fin(1,1,3), fin(1,1,4), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               fin(1,1,3), fin(1,1,4), dflag )
      call dup_up( fout(1,1,3), fin(1,1,1), fin(1,1,2),            &
                   fin(1,1,3), fin(1,1,4), fin(1,1,5), rbufp(1,1,1), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), fin(1,1,5),        &
                        rbufp(1,1,1), rbufp(1,1,2), 4 )
      call dupinterior( fout(1,1,5), fin(1,1,2), fin(1,1,3),       &
                        fin(1,1,4), fin(1,1,5), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 5 )
    else if ( izstart .eq. 2) then
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                      fin(1,1,2), fin(1,1,3), dflag )
      call dup_up( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),          &
                   fin(1,1,2), fin(1,1,3), fin(1,1,4), fin(1,1,5), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,1), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        fin(1,1,5), rbufp(1,1,1), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), fin(1,1,5),        &
                        rbufp(1,1,1), rbufp(1,1,2), 4 )
      call dupinterior( fout(1,1,5), fin(1,1,2), fin(1,1,3),       &
                        fin(1,1,4), fin(1,1,5), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 5 )
    else if ( izstart .eq. 3) then
      call dup_up( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),        &
                   fin(1,1,1), fin(1,1,2), fin(1,1,3), fin(1,1,4), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,1), rbufm(1,1,2),   &
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        fin(1,1,4), fin(1,1,5), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,2), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        fin(1,1,5), rbufp(1,1,1), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), fin(1,1,5),        &
                        rbufp(1,1,1), rbufp(1,1,2), 4 )
      call dupinterior( fout(1,1,5), fin(1,1,2), fin(1,1,3),       &
                        fin(1,1,4), fin(1,1,5), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 5 )
    else if ( izfinish .eq. nz-2) then
      call dupinterior( fout(1,1,nlayer-4), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3)     ,           &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1), nlayer-4 )
      call dupinterior( fout(1,1,nlayer-3), rbufm(1,1,2),          &
                        rbufm(1,1,3), fin(1,1,nlayer-4),           &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer), nlayer-3 )
      call dupinterior( fout(1,1,nlayer-2), rbufm(1,1,3),          &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
      call dupinterior( fout(1,1,nlayer-1), fin(1,1,nlayer-4),     &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer),        &
                        rbufp(1,1,1), rbufp(1,1,2), nlayer-1 )
      call dup_low( fout(1,1,nlayer), fin(1,1,nlayer-3),           &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), rbufp(1,1,1), rbufp(1,1,2), nlayer)
    else if ( izfinish .eq. nz-1) then
      call dupinterior( fout(1,1,nlayer-4), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3)     ,           &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1), nlayer-4 )
      call dupinterior( fout(1,1,nlayer-3), rbufm(1,1,2),          &
                        rbufm(1,1,3), fin(1,1,nlayer-4),           &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer), nlayer-3 )
      call dupinterior( fout(1,1,nlayer-2), rbufm(1,1,3),          &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
      call dup_low( fout(1,1,nlayer-1), fin(1,1,nlayer-4),         &
                    fin(1,1,nlayer-3), fin(1,1,nlayer-2),          &
                    fin(1,1,nlayer-1), fin(1,1,nlayer),            &
                    rbufp(1,1,1), nlayer-1 )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      fin(1,1,nlayer-2), dflag )
    else if ( izfinish .eq. nz) then
      call dupinterior( fout(1,1,nlayer-4), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3)     ,           &
                        fin(1,1,nlayer-4), fin(1,1,nlayer-3),      &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1), nlayer-4 )
      call dupinterior( fout(1,1,nlayer-3), rbufm(1,1,2),          &
                        rbufm(1,1,3), fin(1,1,nlayer-4),           &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer), nlayer-3 )
      call dup_low( fout(1,1,nlayer-2), rbufm(1,1,3),              &
                    fin(1,1,nlayer-4), fin(1,1,nlayer-3),          &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), nlayer-2 )
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), fin(1,1,nlayer-2),        &
                      fin(1,1,nlayer-3), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                   fin(1,1,nlayer-3), dflag )
    else
      call dupinterior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),   &
                        rbufm(1,1,3), fin(1,1,1), fin(1,1,2),      &
                        fin(1,1,3), fin(1,1,4), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,2), rbufm(1,1,3),   &
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        fin(1,1,4), fin(1,1,5), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,3), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        fin(1,1,5), rbufp(1,1,1), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), fin(1,1,5),        &
                        rbufp(1,1,1), rbufp(1,1,2), 4 )
      call dupinterior( fout(1,1,5), fin(1,1,2), fin(1,1,3),       &
                        fin(1,1,4), fin(1,1,5), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 5 )
    end if

  else if ( nlayer .eq. 4) then

!-------------------------
! Now do the nlayer=4 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               fin(1,1,3), fin(1,1,4), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               fin(1,1,3), fin(1,1,4), dflag )
      call dup_up( fout(1,1,3), fin(1,1,1), fin(1,1,2),            &
                   fin(1,1,3), fin(1,1,4), rbufp(1,1,1),rbufp(1,1,2), 3)
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 4 )
    else if ( izstart .eq. 2) then
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                      fin(1,1,2), fin(1,1,3), dflag )
      call dup_up( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),          &
                   fin(1,1,2), fin(1,1,3), fin(1,1,4), rbufp(1,1,1), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,1), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        rbufp(1,1,1), rbufp(1,1,2), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 4 )
    else if ( izstart .eq. 3) then
      call dup_up( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),        &
                   fin(1,1,1), fin(1,1,2), fin(1,1,3), fin(1,1,4), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,1), rbufm(1,1,2),   &
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        fin(1,1,4), rbufp(1,1,1), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,2), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        rbufp(1,1,1), rbufp(1,1,2), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 4 )
    else if ( izfinish .eq. nz-2) then
      call dupinterior( fout(1,1,nlayer-3), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3),                &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer), nlayer-3 )
      call dupinterior( fout(1,1,nlayer-2), rbufm(1,1,2),          &
                        rbufm(1,1,3), fin(1,1,nlayer-3),           &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
      call dupinterior( fout(1,1,nlayer-1), rbufm(1,1,3),          &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer),        &
                        rbufp(1,1,1), rbufp(1,1,2), nlayer-1 )
      call dup_low( fout(1,1,nlayer), fin(1,1,nlayer-3),           &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), rbufp(1,1,1), rbufp(1,1,2), nlayer)
    else if ( izfinish .eq. nz-1) then
      call dupinterior( fout(1,1,nlayer-3), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3),                &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer), nlayer-3 )
      call dupinterior( fout(1,1,nlayer-2), rbufm(1,1,2),          &
                        rbufm(1,1,3), fin(1,1,nlayer-3),           &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
      call dup_low( fout(1,1,nlayer-1), rbufm(1,1,3),              &
                    fin(1,1,nlayer-3), fin(1,1,nlayer-2),          &
                    fin(1,1,nlayer-1), fin(1,1,nlayer),            &
                    rbufp(1,1,1), nlayer-1 )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                   fin(1,1,nlayer), fin(1,1,nlayer-1),             &
                   fin(1,1,nlayer-2), dflag )
    else if ( izfinish .eq. nz) then
      call dupinterior( fout(1,1,nlayer-3), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3),                &
                        fin(1,1,nlayer-3), fin(1,1,nlayer-2),      &
                        fin(1,1,nlayer-1), fin(1,1,nlayer), nlayer-3 )
      call dup_low( fout(1,1,nlayer-2), rbufm(1,1,2),              &
                    rbufm(1,1,3), fin(1,1,nlayer-3),               &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), nlayer-2 )
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), fin(1,1,nlayer-2),        &
                      fin(1,1,nlayer-3), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                   fin(1,1,nlayer-3), dflag )
    else
      call dupinterior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),   &
                        rbufm(1,1,3), fin(1,1,1), fin(1,1,2),      &
                        fin(1,1,3), fin(1,1,4), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,2), rbufm(1,1,3),   &
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        fin(1,1,4), rbufp(1,1,1), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,3), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                        rbufp(1,1,1), rbufp(1,1,2), 3 )
      call dupinterior( fout(1,1,4), fin(1,1,1), fin(1,1,2),       &
                        fin(1,1,3), fin(1,1,4), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 4 )
    end if

  else if ( nlayer .eq. 3) then

!-------------------------
! Now do the nlayer=3 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               fin(1,1,3), rbufp(1,1,1), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               fin(1,1,3), rbufp(1,1,1), dflag )
      call dup_up( fout(1,1,3), fin(1,1,1), fin(1,1,2),            &
                   fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2),rbufp(1,1,3), 3)
    else if ( izstart .eq. 2) then
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                      fin(1,1,2), fin(1,1,3), dflag )
      call dup_up( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),          &
                   fin(1,1,2), fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,1), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 3 )
    else if ( izstart .eq. 3) then
      call dup_up( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),        &
                   fin(1,1,1), fin(1,1,2), fin(1,1,3), rbufp(1,1,1), 1)
      call dupinterior( fout(1,1,2), rbufm(1,1,1), rbufm(1,1,2),   &
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        rbufp(1,1,1), rbufp(1,1,2), 2 )
      if (izfinish .eq. nz-2) then
        call dup_low( fout(1,1,3), rbufm(1,1,2), fin(1,1,1),       &
                      fin(1,1,2), fin(1,1,3), rbufp(1,1,1),        &
                      rbufp(1,1,2), 3 )
      else
        call dupinterior( fout(1,1,3), rbufm(1,1,2), fin(1,1,1),   &
                          fin(1,1,2), fin(1,1,3), rbufp(1,1,1),    &
                          rbufp(1,1,2), rbufp(1,1,3), 3 )
      end if 
    else if ( izfinish .eq. nz-2) then
      call dupinterior( fout(1,1,nlayer-2), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3),                &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
      call dupinterior( fout(1,1,nlayer-1), rbufm(1,1,2),          &
                        rbufm(1,1,3), fin(1,1,nlayer-2),           &
                        fin(1,1,nlayer-1), fin(1,1,nlayer),        &
                        rbufp(1,1,1), rbufp(1,1,2), nlayer-1 )
      call dup_low( fout(1,1,nlayer), rbufm(1,1,3),                &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), rbufp(1,1,1), rbufp(1,1,2), nlayer)
    else if ( izfinish .eq. nz-1) then
      call dupinterior( fout(1,1,nlayer-2), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3),                &
                        fin(1,1,nlayer-2), fin(1,1,nlayer-1),      &
                        fin(1,1,nlayer), rbufp(1,1,1), nlayer-2 )
      call dup_low( fout(1,1,nlayer-1), rbufm(1,1,2),              &
                    rbufm(1,1,3), fin(1,1,nlayer-2),               &
                    fin(1,1,nlayer-1), fin(1,1,nlayer),            &
                    rbufp(1,1,1), nlayer-1 )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      fin(1,1,nlayer-2), dflag )
    else if ( izfinish .eq. nz) then
      call dup_low( fout(1,1,nlayer-2), rbufm(1,1,1),              &
                    rbufm(1,1,2), rbufm(1,1,3),                    &
                    fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                    fin(1,1,nlayer), nlayer-2 )
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), fin(1,1,nlayer-2),        &
                      rbufm(1,1,3), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                   rbufm(1,1,3), dflag )

    else
      call dupinterior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),   &
                        rbufm(1,1,3), fin(1,1,1), fin(1,1,2),      &
                        fin(1,1,3), rbufp(1,1,1), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,2), rbufm(1,1,3),   &
                        fin(1,1,1), fin(1,1,2), fin(1,1,3),        &
                        rbufp(1,1,1), rbufp(1,1,2), 2 )
      call dupinterior( fout(1,1,3), rbufm(1,1,3), fin(1,1,1),     &
                        fin(1,1,2), fin(1,1,3), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 3 )
    end if

  else if ( nlayer .eq. 2) then

!-------------------------
! Now do the nlayer=2 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                               rbufp(1,1,1), rbufp(1,1,2), dflag )
      call topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                               rbufp(1,1,1), rbufp(1,1,2), dflag )
    else if ( izstart .eq. 2) then
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                      fin(1,1,2), rbufp(1,1,1), dflag )
      call dup_up( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),          &
                   fin(1,1,2), rbufp(1,1,1), rbufp(1,1,2),         &
                   rbufp(1,1,3), 2 )
    else if ( izstart .eq. 3) then
      call dup_up( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),        &
                   fin(1,1,1), fin(1,1,2), rbufp(1,1,1),           &
                   rbufp(1,1,2), 1)
      if (izfinish .eq. nz-2) then
        call dup_low(fout(1,1,2), rbufm(1,1,1), rbufm(1,1,2),      &
                     fin(1,1,1), fin(1,1,2), rbufp(1,1,1),         &
                     rbufp(1,1,2), 2 )
      else
        call dupinterior( fout(1,1,2), rbufm(1,1,1), rbufm(1,1,2), &
                          fin(1,1,1), fin(1,1,2), rbufp(1,1,1),    &
                          rbufp(1,1,2), rbufp(1,1,3), 2 )
      end if
    else if ( izfinish .eq. nz-2) then
      call dupinterior( fout(1,1,nlayer-1), rbufm(1,1,1),          &
                        rbufm(1,1,2), rbufm(1,1,3),                &
                        fin(1,1,nlayer-1), fin(1,1,nlayer),        &
                        rbufp(1,1,1), rbufp(1,1,2), nlayer-1 )
      call dup_low( fout(1,1,nlayer), rbufm(1,1,2),                &
                    rbufm(1,1,3), fin(1,1,nlayer-1),               &
                    fin(1,1,nlayer), rbufp(1,1,1), rbufp(1,1,2), nlayer)
    else if ( izfinish .eq. nz-1) then
      call dup_low( fout(1,1,nlayer-1), rbufm(1,1,1),              &
                    rbufm(1,1,2), rbufm(1,1,3),                    &
                    fin(1,1,nlayer-1), fin(1,1,nlayer),            &
                    rbufp(1,1,1), nlayer-1 )
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                      rbufm(1,1,3), dflag )
    else if ( izfinish .eq. nz) then
      call lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                      fin(1,1,nlayer-1), rbufm(1,1,3),             &
                      rbufm(1,1,2), dflag )
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   fin(1,1,nlayer-1), rbufm(1,1,3),                &
                   rbufm(1,1,2), dflag )
    else
      call dupinterior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),   &
                        rbufm(1,1,3), fin(1,1,1), fin(1,1,2),      &
                        rbufp(1,1,1), rbufp(1,1,2), 1 )
      call dupinterior( fout(1,1,2), rbufm(1,1,2), rbufm(1,1,3),   &
                        fin(1,1,1), fin(1,1,2), rbufp(1,1,1),      &
                        rbufp(1,1,2), rbufp(1,1,3), 2 )
    end if

  else if ( nlayer .eq. 1) then

!-------------------------
! Now do the nlayer=1 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call topdz( fout(1,1,1), fin(1,1,1), rbufp(1,1,1),           &
                               rbufp(1,1,2), rbufp(1,1,3), dflag )
    else if ( izstart .eq. 2) then
      call topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                      rbufp(1,1,1), rbufp(1,1,2), dflag )
    else if ( izstart .eq. 3) then
      call dup_up( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),        &
                   fin(1,1,1), rbufp(1,1,1), rbufp(1,1,2), rbufp(1,1,3), 1)
    else if ( izfinish .eq. nz-2) then
      call dup_low( fout(1,1,nlayer), rbufm(1,1,1),                &
                    rbufm(1,1,2), rbufm(1,1,3),                    &
                    fin(1,1,nlayer), rbufp(1,1,1), rbufp(1,1,2), nlayer)
    else if ( izfinish .eq. nz-1) then
      call lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                      fin(1,1,nlayer), rbufm(1,1,3),               &
                      rbufm(1,1,2), dflag )
    else if ( izfinish .eq. nz) then
      call lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                   rbufm(1,1,3), rbufm(1,1,2),                     &
                   rbufm(1,1,1), dflag )

    else
      call dupinterior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),   &
                        rbufm(1,1,3), fin(1,1,1), rbufp(1,1,1),    &
                        rbufp(1,1,2), rbufp(1,1,3), 1 )
    end if
#endif
  else
    print*,' not set up for case nlayer = ',nlayer
    call abort_comms(2)
  end if

!-------------------------------------------------------
! This measures the time taken to get to this stage from
! the start of the subroutine and stores as time_dzup
!-------------------------------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_dzup = time_dzup + (time_tmp2-time_tmp)
  end if

end subroutine d1upbydz

!====================================================================
! <<< d2bydz >>>
!
! Sets fout to the second derivative of fin in the Z direction
! If DFLAG=0 the endpoint values are set to zero.
! If DFLAG=3 then use fourth order one-sided formulation
!
! IF DFLAG = 0 then the six neighbouring layers are need for the 
!               boundary point.
!
! Note: DFLAG=1, 2 used to be used for magnetic boundary points
! These flags are now obsolete (pjb 9/2/04)
!--------------------------------------------------------------------
! As the data is layed out in slabs the neighbouring layers have
! to be communicated to this layer.
!
!
!        myrank-1      rubfm
!                    ---------
!        myrank        layer
!                    ---------
!        myrank+1      rubfp
!
!--------------------------------------------------------------------
!====================================================================
Subroutine d2bydz( fin, fout, dflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=8), intent(in)      :: fin(nx,ny,nlayer)
  real(kind=8), intent(out)     :: fout(nx,ny,nlayer)
  integer,      intent(in)      :: dflag

!-------------------
! Local declarations
!-------------------
  integer                       :: k
#ifdef MPI 
  integer                       :: ier
#endif

!----------------------------
! Starts timing off if needed
!---------------------------- 
  if (timer) then
    call findtime(time_tmp)
  end if

#ifdef MPI

!------------------------------
! Do communication if necessary
!------------------------------
  if (nproc .gt. 1) then

!--------------------------------
! Initialise communcation buffers
!--------------------------------
    rbufm=0.0d0
    rbufp=0.0d0
    msg = 0
    mstat = 0
    ireq = 0 
    rbufp_counter = 0 
    rbufm_counter = 0

!----------------------------
! Post receive commands first
!----------------------------
    do k=1, myrank
      if (d2_distrecv(k) .gt. 0) then
        ncount = nx*ny*d2_distrecv(k)
        call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        senddown_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufm_counter = rbufm_counter + d2_distrecv(k)
      end if
    end do

    do k=myrank+2,nproc
      if (d2_distrecv(k) .gt. 0) then
        ncount = nx*ny*d2_distrecv(k)
        call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        sendup_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufp_counter = rbufp_counter + d2_distrecv(k)
      end if
    end do 

    do k=1, myrank
      if (d2_distsend(k) .gt. 0) then
        ncount = nx*ny*d2_distsend(k)
        call MPI_Isend( fin(1,1,1), ncount,                           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        sendup_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

    do k=myrank+2,nproc
      if (d2_distsend(k) .gt. 0) then
        ncount = nx*ny*d2_distsend(k)
        call MPI_Isend( fin(1,1,nlayer-d2_distsend(k)+1), ncount,     &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        senddown_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

  if (nlayer .ge. 4) then
!--------------------------------------------
! First do the simplest case, for nlayer >= 4
!--------------------------------------------
! Do the top two layers on a node
!--------------------------------------------
    if ( izstart .eq. 1 ) then
      if (nlayer .ge. 6) then
        call d2topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),           &
                      fin(1,1,3), fin(1,1,4), fin(1,1,5),            &
                      fin(1,1,6), dflag )
        call d2topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),       &
                          fin(1,1,3), fin(1,1,4), fin(1,1,5),        &
                          fin(1,1,6), dflag )
#ifdef MPI
      else if (nlayer .eq. 5) then
        call d2topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),           &
                      fin(1,1,3), fin(1,1,4), fin(1,1,5),            &
                      rbufp(1,1,1), dflag )
        call d2topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),       &
                          fin(1,1,3), fin(1,1,4), fin(1,1,5),        &
                          rbufp(1,1,1), dflag )
      else
        call d2topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),           &
                      fin(1,1,3), fin(1,1,4), rbufp(1,1,1),          &
                      rbufp(1,1,2), dflag )
        call d2topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),       &
                          fin(1,1,3), fin(1,1,4), rbufp(1,1,1),      &
                          rbufp(1,1,2), dflag )
#endif
      end if
#ifdef MPI
    else if (izstart .eq. 2) then
      if (nlayer .ge. 5) then
        call d2topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),     &
                          fin(1,1,2), fin(1,1,3), fin(1,1,4),        &
                          fin(1,1,5), dflag )
      else
        call d2topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),     &
                          fin(1,1,2), fin(1,1,3), fin(1,1,4),        & 
                          rbufp(1,1,1), dflag )
      end if
      call d2interior( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),        &
                       fin(1,1,2), fin(1,1,3), fin(1,1,4) )
#endif
    else 
#ifdef MPI
      call d2interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                       fin(1,1,1), fin(1,1,2),   fin(1,1,3 ) )
      call d2interior( fout(1,1,2), rbufm(1,1,2), fin(1,1,1),        &
                       fin(1,1,2), fin(1,1,3), fin(1,1,4) ) 
#endif
    end if

!--------------------------------
! Do the middle layers on a node
!--------------------------------
    if (nlayer .ge. 5) then
      do k=3, nlayer-2
        call d2interior( fout(1,1,k), fin(1,1,k-2), fin(1,1,k-1),    &
                         fin(1,1,k), fin(1,1,k+1), fin(1,1,k+2) )
      end do
    end if

!-----------------------------------
! Do the bottom two layers on a node
!-----------------------------------
    if ( izfinish .eq. nz ) then
      if (nlayer .ge. 6) then
        call d2lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),       & 
                          fin(1,1,nlayer-1), fin(1,1,nlayer-2),      &
                          fin(1,1,nlayer-3), fin(1,1,nlayer-4),      &
                          fin(1,1,nlayer-5), dflag )
        call d2lastdz( fout(1,1,nlayer), fin(1,1,nlayer),            & 
                       fin(1,1,nlayer-1), fin(1,1,nlayer-2),         &
                       fin(1,1,nlayer-3), fin(1,1,nlayer-4),         &
                       fin(1,1,nlayer-5), dflag )
#ifdef MPI
      else if (nlayer .eq. 5) then
        call d2lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),       & 
                          fin(1,1,nlayer-1), fin(1,1,nlayer-2),      &
                          fin(1,1,nlayer-3), fin(1,1,nlayer-4),      &
                          rbufm(1,1,2), dflag )
        call d2lastdz( fout(1,1,nlayer), fin(1,1,nlayer),            & 
                       fin(1,1,nlayer-1), fin(1,1,nlayer-2),         &
                       fin(1,1,nlayer-3), fin(1,1,nlayer-4),         &
                       rbufm(1,1,2), dflag )
      else
        call d2lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),       & 
                          fin(1,1,nlayer-1), fin(1,1,nlayer-2),      &
                          fin(1,1,nlayer-3), rbufm(1,1,2),           &
                          rbufm(1,1,1), dflag )
        call d2lastdz( fout(1,1,nlayer), fin(1,1,nlayer),            & 
                       fin(1,1,nlayer-1), fin(1,1,nlayer-2),         &
                       fin(1,1,nlayer-3), rbufm(1,1,2),              &
                       rbufm(1,1,1), dflag )
#endif
     end if
#ifdef MPI
   else if (izfinish .eq. nz-1) then
     if (nlayer .ge. 5) then
       call d2lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),             & 
                         fin(1,1,nlayer), fin(1,1,nlayer-1),         &
                         fin(1,1,nlayer-2), fin(1,1,nlayer-3),       &
                         fin(1,1,nlayer-4), dflag )
     else
       call d2lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),             & 
                         fin(1,1,nlayer), fin(1,1,nlayer-1),         &
                         fin(1,1,nlayer-2), fin(1,1,nlayer-3),       &
                         rbufm(1,1,2), dflag )
     end if
     call d2interior( fout(1,1,nlayer-1), fin(1,1,nlayer-3),         &
                      fin(1,1,nlayer-2), fin(1,1,nlayer-1),          &
                      fin(1,1,nlayer), rbufp(1,1,1) )
#endif
   else 
#ifdef MPI
     call d2interior( fout(1,1,nlayer-1),fin(1,1,nlayer-3),          &
                      fin(1,1,nlayer-2), fin(1,1,nlayer-1),          & 
                      fin(1,1,nlayer), rbufp(1,1,1) )
     call d2interior( fout(1,1,nlayer),  fin(1,1,nlayer-2),          &
                      fin(1,1,nlayer-1), fin(1,1,nlayer),            &
                      rbufp(1,1,1), rbufp(1,1,2) )
#endif
    end if

#ifdef MPI
  else if ( nlayer .eq. 3) then

!--------------------------------
! Now do the nlayer=3 case
!--------------------------------
! Do the top two layers on a node
!--------------------------------
    if ( izstart .eq. 1 ) then
      call d2topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                    fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2),          &
                    rbufp(1,1,3), dflag )
      call d2topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                        fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2),      &
                        rbufp(1,1,3), dflag )
      call d2interior( fout(1,1,3), fin(1,1,1), fin(1,1,2),          &
                       fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2) )
    else if (izstart .eq. 2) then
#ifdef MPI
      call d2topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                        fin(1,1,2), fin(1,1,3), rbufp(1,1,1),        &
                        rbufp(1,1,2), dflag )
      call d2interior( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),        &
                       fin(1,1,2), fin(1,1,3), rbufp(1,1,1) )
      call d2interior( fout(1,1,3), fin(1,1,1), fin(1,1,2),          &
                       fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2) )
#endif
    else if (izfinish .eq. nz) then
      call d2interior( fout(1,1,nlayer-2),rbufm(1,1,2),rbufm(1,1,3), &
                       fin(1,1,nlayer-2), fin(1,1,nlayer-1),         &
                       fin(1,1,nlayer) )
      call d2lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                        fin(1,1,nlayer-1), fin(1,1,nlayer-2),        &
                        rbufm(1,1,3), rbufm(1,1,2), rbufm(1,1,1), dflag )
      call d2lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                     fin(1,1,nlayer-1), fin(1,1,nlayer-2),           &
                     rbufm(1,1,3), rbufm(1,1,2), rbufm(1,1,1), dflag )
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call d2interior( fout(1,1,nlayer-2),rbufm(1,1,1),rbufm(1,1,2), &
                       fin(1,1,nlayer-2), fin(1,1,nlayer-1),         &
                       fin(1,1,nlayer) )
      call d2interior( fout(1,1,nlayer-1),rbufm(1,1,2),              &
                       fin(1,1,nlayer-2), fin(1,1,nlayer-1),         &
                       fin(1,1,nlayer),rbufp(1,1,1) )
      call d2lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                        fin(1,1,nlayer), fin(1,1,nlayer-1),          &
               fin(1,1,nlayer-2), rbufm(1,1,2),rbufm(1,1,1), dflag )
#endif
    else 
#ifdef MPI
      call d2interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                       fin(1,1,1),  fin(1,1,2),   fin(1,1,3 ) )
      call d2interior( fout(1,1,2), rbufm(1,1,2), fin(1,1,1),        &
                       fin(1,1,2),  fin(1,1,3), rbufp(1,1,1) ) 
      call d2interior( fout(1,1,3), fin(1,1,1), fin(1,1,2),          &
                       fin(1,1,3), rbufp(1,1,1), rbufp(1,1,2) )
#endif
    end if

  else if ( nlayer .eq. 2) then

!-------------------------
! Now do the nlayer=2 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call d2topdz( fout(1,1,1), fin(1,1,1), fin(1,1,2),             &
                    rbufp(1,1,1), rbufp(1,1,2), rbufp(1,1,3),        &
                    rbufp(1,1,4) ,dflag )
      call d2topdz_low( fout(1,1,2), fin(1,1,1), fin(1,1,2),         &
                        rbufp(1,1,1), rbufp(1,1,2), rbufp(1,1,3),    &
                        rbufp(1,1,4) ,dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call d2topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                        fin(1,1,2), rbufp(1,1,1), rbufp(1,1,2),      &
                        rbufp(1,1,3), dflag )
      call d2interior( fout(1,1,2), rbufm(1,1,1), fin(1,1,1),        &
                       fin(1,1,2), rbufp(1,1,1), rbufp(1,1,2) )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call d2interior( fout(1,1,nlayer-1),rbufm(1,1,2),              &
                       rbufm(1,1,3),fin(1,1,nlayer-1),               &
                       fin(1,1,nlayer),rbufp(1,1,1) )
      call d2lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                        fin(1,1,nlayer), fin(1,1,nlayer-1),          &
                        rbufm(1,1,3), rbufm(1,1,2), rbufm(1,1,1), dflag )
#endif
    else if (izfinish .eq. nz) then
      call d2lastdz_up( fout(1,1,nlayer-1), fin(1,1,nlayer),         & 
                        fin(1,1,nlayer-1), rbufm(1,1,4),             &
                        rbufm(1,1,3), rbufm(1,1,2), rbufm(1,1,1), dflag )
      call d2lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                     fin(1,1,nlayer-1), rbufm(1,1,4),                &
                     rbufm(1,1,3), rbufm(1,1,2), rbufm(1,1,1), dflag )
    else
      call d2interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                       fin(1,1,1), fin(1,1,2),   rbufp(1,1,1) )
      call d2interior( fout(1,1,2), rbufm(1,1,2), fin(1,1,1),        &
                       fin(1,1,2), rbufp(1,1,1), rbufp(1,1,2) ) 
    end if
  else if ( nlayer .eq. 1) then

!-------------------------
! Now do the nlayer=1 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call d2topdz( fout(1,1,1), fin(1,1,1), rbufp(1,1,1),           &
                    rbufp(1,1,2), rbufp(1,1,3), rbufp(1,1,4),        &
                    rbufp(1,1,5), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call d2topdz_low( fout(1,1,1), rbufm(1,1,1), fin(1,1,1),       &
                        rbufp(1,1,1), rbufp(1,1,2), rbufp(1,1,3),    &
                        rbufp(1,1,4), dflag )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call d2lastdz_up( fout(1,1,nlayer), rbufp(1,1,1),              & 
                        fin(1,1,nlayer), rbufm(1,1,4),               &
                        rbufm(1,1,3), rbufm(1,1,2), rbufm(1,1,1), dflag )
#endif
    else if (izfinish .eq. nz) then
      call d2lastdz( fout(1,1,nlayer), fin(1,1,nlayer),              & 
                     rbufm(1,1,5), rbufm(1,1,4),rbufm(1,1,3),        &
                     rbufm(1,1,2), rbufm(1,1,1), dflag )
    else
      call d2interior( fout(1,1,1), rbufm(1,1,1), rbufm(1,1,2),      &
                       fin(1,1,1), rbufp(1,1,1), rbufp(1,1,2) )
    end if
#endif
  else
    print*,' not set up for case nlayer = ',nlayer
    call abort_comms(2)
  end if

!-------------------------------------------------------
! This measures the time taken to get to this stage from
! the start of the subroutine and stores as time_d2z
!-------------------------------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_d2z = time_d2z + (time_tmp2-time_tmp)
  end if

end subroutine d2bydz

!====================================================================
! <<< d2bydz1D >>>
!
! Sets fout to the 1D second derivative of fin in the Z direction
! If DFLAG=0 the endpoint values are set to zero.
! If DFLAG=3 then use fourth order one-sided formulation
!
! IF DFLAG = 0 then the six neighbouring layers are need for the 
!               boundary point.
!
! Note: DFLAG=1, 2 used to be used for magnetic boundary points
! These flags are now obsolete (pjb 9/2/04)
!--------------------------------------------------------------------
! As the data is layed out in slabs the neighbouring layers have
! to be communicated to this layer.
!
!
!        myrank-1      rubfm
!                    ---------
!        myrank        layer
!                    ---------
!        myrank+1      rubfp
!
!--------------------------------------------------------------------
!====================================================================
Subroutine d2bydz1D( fin, fout, dflag )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  real(kind=8), intent(in)      :: fin(nlayer)
  real(kind=8), intent(out)     :: fout(nlayer)
  integer,      intent(in)      :: dflag

!-------------------
! Local declarations
!-------------------
  integer                       :: k
#ifdef MPI 
  integer                       :: ier
#endif

!----------------------------
! Starts timing off if needed
!---------------------------- 
  if (timer) then
    call findtime(time_tmp)
  end if

#ifdef MPI

!------------------------------
! Do communication if necessary
!------------------------------
  if (nproc .gt. 1) then

!--------------------------------
! Initialise communcation buffers
!--------------------------------
    recvbuf=0.0d0
    sendbuf=0.0d0
    msg = 0
    mstat = 0
    ireq = 0 
    rbufp_counter = 0 
    rbufm_counter = 0

!----------------------------
! Post receive commands first
!----------------------------
    do k=1, myrank
      if (d2_distrecv(k) .gt. 0) then
        ncount = d2_distrecv(k)
        call MPI_Irecv( recvbuf(1+rbufm_counter), ncount,             &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        senddown_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufm_counter = rbufm_counter + d2_distrecv(k)
      end if
    end do

    do k=myrank+2,nproc
      if (d2_distrecv(k) .gt. 0) then
        ncount = d2_distrecv(k)
        call MPI_Irecv( sendbuf(1+rbufp_counter), ncount,             &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        sendup_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufp_counter = rbufp_counter + d2_distrecv(k)
      end if
    end do 

    do k=1, myrank
      if (d2_distsend(k) .gt. 0) then
        ncount = d2_distsend(k)
        call MPI_Isend( fin(1), ncount,                               &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        sendup_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

    do k=myrank+2,nproc
      if (d2_distsend(k) .gt. 0) then
        ncount = d2_distsend(k)
        call MPI_Isend( fin(nlayer-d2_distsend(k)+1), ncount,         &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        senddown_tag(myrank+1),comm,ireq(msg+1),ier)
        msg = msg + 1
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

  if (nlayer .ge. 4) then
!--------------------------------------------
! First do the simplest case, for nlayer >= 4
!--------------------------------------------
! Do the top two layers on a node
!--------------------------------------------
    if ( izstart .eq. 1 ) then
      if (nlayer .ge. 6) then
        call d2topdz1D( fout(1), fin(1), fin(2),           &
                      fin(3), fin(4), fin(5),            &
                      fin(6), dflag )
        call d2topdz_low1D( fout(2), fin(1), fin(2),       &
                          fin(3), fin(4), fin(5),        &
                          fin(6), dflag )
#ifdef MPI
      else if (nlayer .eq. 5) then
        call d2topdz1D( fout(1), fin(1), fin(2),           &
                      fin(3), fin(4), fin(5),            &
                      sendbuf(1), dflag )
        call d2topdz_low1D( fout(2), fin(1), fin(2),       &
                          fin(3), fin(4), fin(5),        &
                          sendbuf(1), dflag )
      else
        call d2topdz1D( fout(1), fin(1), fin(2),           &
                      fin(3), fin(4), sendbuf(1),        &
                      sendbuf(2), dflag )
        call d2topdz_low1D( fout(2), fin(1), fin(2),       &
                          fin(3), fin(4), sendbuf(1),    &
                          sendbuf(2), dflag )
#endif
      end if
#ifdef MPI
    else if (izstart .eq. 2) then
      if (nlayer .ge. 5) then
        call d2topdz_low1D( fout(1), recvbuf(1), fin(1),   &
                          fin(2), fin(3), fin(4),        &
                          fin(5), dflag )
      else
        call d2topdz_low1D( fout(1), recvbuf(1), fin(1),   &
                          fin(2), fin(3), fin(4),        & 
                          sendbuf(1), dflag )
      end if
      call d2interior1D( fout(2), recvbuf(1), fin(1),      &
                       fin(2), fin(3), fin(4) )
#endif
    else 
#ifdef MPI
      call d2interior1D( fout(1), recvbuf(1), recvbuf(2),  &
                       fin(1), fin(2),   fin(3) )
      call d2interior1D( fout(2), recvbuf(2), fin(1),      &
                       fin(2), fin(3), fin(4) ) 
#endif
    end if

!--------------------------------
! Do the middle layers on a node
!--------------------------------
    if (nlayer .ge. 5) then
      do k=3, nlayer-2
        call d2interior1D( fout(k), fin(k-2), fin(k-1),    &
                         fin(k), fin(k+1), fin(k+2) )
      end do
    end if

!-----------------------------------
! Do the bottom two layers on a node
!-----------------------------------
    if ( izfinish .eq. nz ) then
      if (nlayer .ge. 6) then
        call d2lastdz_up1D( fout(nlayer-1), fin(nlayer),       & 
                          fin(nlayer-1), fin(nlayer-2),      &
                          fin(nlayer-3), fin(nlayer-4),      &
                          fin(nlayer-5), dflag )
        call d2lastdz1D( fout(nlayer), fin(nlayer),            & 
                       fin(nlayer-1), fin(nlayer-2),         &
                       fin(nlayer-3), fin(nlayer-4),         &
                       fin(nlayer-5), dflag )
#ifdef MPI
      else if (nlayer .eq. 5) then
        call d2lastdz_up1D( fout(nlayer-1), fin(nlayer),       & 
                          fin(nlayer-1), fin(nlayer-2),      &
                          fin(nlayer-3), fin(nlayer-4),      &
                          recvbuf(2), dflag )
        call d2lastdz1D( fout(nlayer), fin(nlayer),            & 
                       fin(nlayer-1), fin(nlayer-2),         &
                       fin(nlayer-3), fin(nlayer-4),         &
                       recvbuf(2), dflag )
      else
        call d2lastdz_up1D( fout(nlayer-1), fin(nlayer),       & 
                          fin(nlayer-1), fin(nlayer-2),      &
                          fin(nlayer-3), recvbuf(2),         &
                          recvbuf(1), dflag )
        call d2lastdz1D( fout(nlayer), fin(nlayer),            & 
                       fin(nlayer-1), fin(nlayer-2),         &
                       fin(nlayer-3), recvbuf(2),            &
                       recvbuf(1), dflag )
#endif
     end if
#ifdef MPI
   else if (izfinish .eq. nz-1) then
     if (nlayer .ge. 5) then
       call d2lastdz_up1D( fout(nlayer), sendbuf(1),           & 
                         fin(nlayer), fin(nlayer-1),         &
                         fin(nlayer-2), fin(nlayer-3),       &
                         fin(nlayer-4), dflag )
     else
       call d2lastdz_up1D( fout(nlayer), sendbuf(1),           & 
                         fin(nlayer), fin(nlayer-1),         &
                         fin(nlayer-2), fin(nlayer-3),       &
                         recvbuf(2), dflag )
     end if
     call d2interior1D( fout(nlayer-1), fin(nlayer-3),         &
                      fin(nlayer-2), fin(nlayer-1),          &
                      fin(nlayer), sendbuf(1) )
#endif
   else 
#ifdef MPI
     call d2interior1D( fout(nlayer-1),fin(nlayer-3),          &
                      fin(nlayer-2), fin(nlayer-1),          & 
                      fin(nlayer), sendbuf(1) )
     call d2interior1D( fout(nlayer),  fin(nlayer-2),          &
                      fin(nlayer-1), fin(nlayer),            &
                      sendbuf(1), sendbuf(2) )
#endif
    end if

#ifdef MPI
  else if ( nlayer .eq. 3) then

!--------------------------------
! Now do the nlayer=3 case
!--------------------------------
! Do the top two layers on a node
!--------------------------------
    if ( izstart .eq. 1 ) then
      call d2topdz1D( fout(1), fin(1), fin(2),                 &
                    fin(3), sendbuf(1), sendbuf(2),          &
                    sendbuf(3), dflag )
      call d2topdz_low1D( fout(2), fin(1), fin(2),             &
                        fin(3), sendbuf(1), sendbuf(2),      &
                        sendbuf(3), dflag )
      call d2interior1D( fout(3), fin(1), fin(2),              &
                       fin(3), sendbuf(1), sendbuf(2) )
    else if (izstart .eq. 2) then
#ifdef MPI
      call d2topdz_low1D( fout(1), recvbuf(1), fin(1),         &
                        fin(2), fin(3), sendbuf(1),          &
                        sendbuf(2), dflag )
      call d2interior1D( fout(2), recvbuf(1), fin(1),          &
                       fin(2), fin(3), sendbuf(1) )
      call d2interior1D( fout(3), fin(1), fin(2),              &
                       fin(3), sendbuf(1), sendbuf(2) )
#endif
    else if (izfinish .eq. nz) then
      call d2interior1D( fout(nlayer-2),recvbuf(2),recvbuf(3), &
                       fin(nlayer-2), fin(nlayer-1),         &
                       fin(nlayer) )
      call d2lastdz_up1D( fout(nlayer-1), fin(nlayer),         & 
                        fin(nlayer-1), fin(nlayer-2),        &
                        recvbuf(3), recvbuf(2), recvbuf(1), dflag )
      call d2lastdz1D( fout(nlayer), fin(nlayer),              & 
                     fin(nlayer-1), fin(nlayer-2),           &
                     recvbuf(3), recvbuf(2), recvbuf(1), dflag )
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call d2interior1D( fout(nlayer-2),recvbuf(1),recvbuf(2), &
                       fin(nlayer-2), fin(nlayer-1),         &
                       fin(nlayer) )
      call d2interior1D( fout(nlayer-1),recvbuf(2),            &
                       fin(nlayer-2), fin(nlayer-1),         &
                       fin(nlayer),sendbuf(1) )
      call d2lastdz_up1D( fout(nlayer), sendbuf(1),            & 
                        fin(nlayer), fin(nlayer-1),          &
               fin(nlayer-2), recvbuf(2),recvbuf(1), dflag )
#endif
    else 
#ifdef MPI
      call d2interior1D( fout(1), recvbuf(1), recvbuf(2),      &
                       fin(1),  fin(2),   fin(3 ) )
      call d2interior1D( fout(2), recvbuf(2), fin(1),          &
                       fin(2),  fin(3), sendbuf(1) ) 
      call d2interior1D( fout(3), fin(1), fin(2),              &
                       fin(3), sendbuf(1), sendbuf(2) )
#endif
    end if

  else if ( nlayer .eq. 2) then

!-------------------------
! Now do the nlayer=2 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call d2topdz1D( fout(1), fin(1), fin(2),                 &
                    sendbuf(1), sendbuf(2), sendbuf(3),      &
                    sendbuf(4) ,dflag )
      call d2topdz_low1D( fout(2), fin(1), fin(2),             &
                        sendbuf(1), sendbuf(2), sendbuf(3),  &
                        sendbuf(4) ,dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call d2topdz_low1D( fout(1), recvbuf(1), fin(1),         &
                        fin(2), sendbuf(1), sendbuf(2),      &
                        sendbuf(3), dflag )
      call d2interior1D( fout(2), recvbuf(1), fin(1),          &
                       fin(2), sendbuf(1), sendbuf(2) )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call d2interior1D( fout(nlayer-1),recvbuf(2),            &
                       recvbuf(3),fin(nlayer-1),             &
                       fin(nlayer),sendbuf(1) )
      call d2lastdz_up1D( fout(nlayer), sendbuf(1),            & 
                        fin(nlayer), fin(nlayer-1),          &
                        recvbuf(3), recvbuf(2), recvbuf(1), dflag )
#endif
    else if (izfinish .eq. nz) then
      call d2lastdz_up1D( fout(nlayer-1), fin(nlayer),         & 
                        fin(nlayer-1), recvbuf(4),           &
                        recvbuf(3), recvbuf(2), recvbuf(1), dflag )
      call d2lastdz1D( fout(nlayer), fin(nlayer),              & 
                     fin(nlayer-1), recvbuf(4),              &
                     recvbuf(3), recvbuf(2), recvbuf(1), dflag )
    else
      call d2interior1D( fout(1), recvbuf(1), recvbuf(2),      &
                       fin(1), fin(2),   sendbuf(1) )
      call d2interior1D( fout(2), recvbuf(2), fin(1),          &
                       fin(2), sendbuf(1), sendbuf(2) ) 
    end if
  else if ( nlayer .eq. 1) then

!-------------------------
! Now do the nlayer=1 case
!-------------------------
    if ( izstart .eq. 1 ) then
      call d2topdz1D( fout(1), fin(1), sendbuf(1),             &
                    sendbuf(2), sendbuf(3), sendbuf(4),      &
                    sendbuf(5), dflag )
    else if (izstart .eq. 2) then
#ifdef MPI
      call d2topdz_low1D( fout(1), recvbuf(1), fin(1),         &
                        sendbuf(1), sendbuf(2), sendbuf(3),  &
                        sendbuf(4), dflag )
#endif
    else if (izfinish .eq. nz-1) then
#ifdef MPI
      call d2lastdz_up1D( fout(nlayer), sendbuf(1),            & 
                        fin(nlayer), recvbuf(4),             &
                        recvbuf(3), recvbuf(2), recvbuf(1), dflag )
#endif
    else if (izfinish .eq. nz) then
      call d2lastdz1D( fout(nlayer), fin(nlayer),              & 
                     recvbuf(5), recvbuf(4),recvbuf(3),      &
                     recvbuf(2), recvbuf(1), dflag )
    else
      call d2interior1D( fout(1), recvbuf(1), recvbuf(2),      &
                       fin(1), sendbuf(1), sendbuf(2) )
    end if
#endif
  else
    print*,' not set up for case nlayer = ',nlayer
    call abort_comms(2)
  end if

!-------------------------------------------------------
! This measures the time taken to get to this stage from
! the start of the subroutine and stores as time_d2z
!-------------------------------------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_d2z = time_d2z + (time_tmp2-time_tmp)
  end if

end subroutine d2bydz1D

!================================================================
! Boundary:
!----------------------------------------------------------------
! Applies the boundary conditions:
! W = 0 at Z = 0 and z = 1;
! T =1 at Z = 0.
!----------------------------------------------------------------
! Vertical derivatives of U, V are zero at z=0 and 1.0
! Either T=1 + theta or dt/dz = theta at z=1.
! Use 3rd order finite difference formula to determine U1, V1 etc
!================================================================
Subroutine boundary( )
  implicit none 

!-------------------
! Local declarations
!-------------------
  integer          :: i, j
  integer          :: inloop
  real(kind=dp)    :: t1old
  real(kind=dp)    :: ft1, dft1
#ifdef MPI 
  integer          :: k
  integer          :: ier
#endif

!----------------------------
! Starts timing off if needed
!----------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!------------------------------------------
! Do communication for u and v if necessary
!------------------------------------------
#ifdef MPI
  if (nproc .gt. 1) then

!-------------------
! Initialise buffers
!-------------------
    rbufp = 0.0d0
    rbufm = 0.0d0
    rbufp_counter = 0
    rbufm_counter = 0
    msg = 0
    mstat = 0
    ireq  = 0 

!--------------------------------
! Post all receive commands first
!--------------------------------
    if (myrank .eq. 0) then
      rbufp_counter = 0
      rbufm_counter = 0
      do k=2, nproc, 1
        if (b_distrecv(k) .gt. 0) then
        ncount = nx*ny*b_distrecv(k)
        call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        bvsendup_tag(k), comm, ireq(msg+1), ier )
        call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        busendup_tag(k), comm, ireq(msg+2), ier ) 
        msg = msg + 2
        rbufp_counter = rbufp_counter + b_distrecv(k)
        rbufm_counter = rbufm_counter + b_distrecv(k)
        end if
      end do
    end if

    if (myrank .eq. nproc-1) then
      rbufp_counter = 0
      rbufm_counter = 0
      do k=1, nproc-1
        if (b_distrecv(k) .gt. 0) then
        ncount = nx*ny*b_distrecv(k)
        call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        bvsenddown_tag(k), comm, ireq(msg+1), ier )
        call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        busenddown_tag(k), comm, ireq(msg+2), ier ) 
        msg = msg + 2
        rbufp_counter = rbufp_counter + b_distrecv(k)
        rbufm_counter = rbufm_counter + b_distrecv(k)
        end if
      end do
    end if 

!-------------------
! Send commands here
!-------------------
    do k=2, nproc, 1
      if (myrank .eq. k-1) then
        if (b_distsend(1) .gt. 0) then
          ncount = nx*ny*b_distsend(1)
          call MPI_Isend( u(1,1,1), ncount, MPI_DOUBLE_PRECISION, 0,    &
                          busendup_tag(k), comm, ireq(msg+1), ier )
          call MPI_Isend( v(1,1,1), ncount, MPI_DOUBLE_PRECISION, 0,    &
                          bvsendup_tag(k), comm, ireq(msg+2), ier )
          msg = msg + 2
        end if
      end if
    end do  

    do k=1, nproc-1
      if (myrank .eq. k-1) then
        if (b_distsend(nproc) .gt. 0) then
          ncount = nx*ny*b_distsend(nproc)
          call MPI_Isend( u(1,1,nlayer-b_distsend(nproc)+1), ncount,    &
                          MPI_DOUBLE_PRECISION, nproc-1,                &
                          busenddown_tag(k), comm, ireq(msg+1), ier )
          call MPI_Isend( v(1,1,nlayer-b_distsend(nproc)+1), ncount,    & 
                          MPI_DOUBLE_PRECISION, nproc-1,                &
                          bvsenddown_tag(k), comm, ireq(msg+2), ier )
          msg = msg + 2
        end if
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

!-------------------------------------------------------
! Do the top layer boundary conditions if myrank is zero
!-------------------------------------------------------
  if ( myrank .eq. 0 ) then

!-------------------------------------
! Sets vertical velocity equal to zero
!-------------------------------------
    do j=1,ny
      do i=1,nx
         w(i,j,1) = 0.0d0
      end do
    end do

!----------------------------------
! Stress free conditions on u and v
!----------------------------------
    if (nlayer .ge. 4) then
      do j=1,ny
        do i=1,nx
          u(i,j,1) = elf*(  18.0d0*u(i,j,2)      &
                          - 9.0d0*u(i,j,3)       &
                          + 2.0d0*u(i,j,4)       )
          v(i,j,1) = elf*(  18.0d0*v(i,j,2)      &
                          - 9.0d0*v(i,j,3)       &
                          + 2.0d0*v(i,j,4)       )
        end do
      end do

#ifdef MPI
    else if (nlayer .eq. 3) then
      do j=1,ny
        do i=1,nx
          u(i,j,1) = elf*(  18.0d0*u(i,j,2)      &
                          - 9.0d0*u(i,j,3)       &
                          + 2.0d0*rbufm(i,j,1)   )
          v(i,j,1) = elf*(  18.0d0*v(i,j,2)      &
                          - 9.0d0*v(i,j,3)       &
                          + 2.0d0*rbufp(i,j,1)   )
        end do
      end do
    else if (nlayer .eq. 2) then
      do j=1,ny
        do i=1,nx
          u(i,j,1) = elf*(  18.0d0*u(i,j,2)      &
                          - 9.0d0*rbufm(i,j,1)   &
                          + 2.0d0*rbufm(i,j,2)   )
          v(i,j,1) = elf*(  18.0d0*v(i,j,2)      &
                          - 9.0d0*rbufp(i,j,1)   &
                          + 2.0d0*rbufp(i,j,2)   )
        end do
      end do
    else if (nlayer .eq. 1) then
      do j=1,ny
        do i=1,nx
          u(i,j,1) = elf*(  18.0d0*rbufm(i,j,1)  &
                          - 9.0d0*rbufm(i,j,2)   &
                          + 2.0d0*rbufm(i,j,3)   )
          v(i,j,1) = elf*(  18.0d0*rbufp(i,j,1)  &
                          - 9.0d0*rbufp(i,j,2)   &
                          + 2.0d0*rbufp(i,j,3)   )
        end do
      end do
#endif
    else
      print*,'Not set up for case nlayer = ',nlayer
      call abort_comms(0)
    end if
  end if

!---------------------------------------------
! Do the bottom layer if myrank + 1 equals the
! number of processors 
!---------------------------------------------
  if ( myrank + 1 .eq. nproc )  then

!-------------------------
! Boundary condition for W
!-------------------------
    do j=1,ny
      do i=1,nx
         w(i,j,nlayer) = 0.0d0
      end do
    end do

!---------------------------------------
! Bottom boundary conditions for U and V
!---------------------------------------
    if (nlayer .ge. 4) then
      do j=1,ny
        do i=1,nx
          u(i,j,nlayer) = elf*(  18.0d0*u(i,j,nlayer-1)     &
                               - 9.0d0*u(i,j,nlayer-2)      &
                               + 2.0d0*u(i,j,nlayer-3)      )
          v(i,j,nlayer) = elf*(  18.0d0*v(i,j,nlayer-1)     &
                               - 9.0d0*v(i,j,nlayer-2)      &
                               + 2.0d0*v(i,j,nlayer-3)      )
       end do
     end do

#ifdef MPI
    else if (nlayer .eq. 3) then
      do j=1,ny
        do i=1,nx
          u(i,j,nlayer) = elf*(  18.0d0*u(i,j,nlayer-1)     &
                               - 9.0d0*u(i,j,nlayer-2)      &
                               + 2.0d0*rbufm(i,j,1)         )
          v(i,j,nlayer) = elf*(  18.0d0*v(i,j,nlayer-1)     &
                               - 9.0d0*v(i,j,nlayer-2)      &
                               + 2.0d0*rbufp(i,j,1)         )
        end do
      end do
    else if (nlayer .eq. 2) then
      do j=1,ny
        do i=1,nx
          u(i,j,nlayer) = elf*(  18.0d0*u(i,j,nlayer-1)     &
                               - 9.0d0*rbufm(i,j,2)         &
                               + 2.0d0*rbufm(i,j,1)         )
          v(i,j,nlayer) = elf*(  18.0d0*v(i,j,nlayer-1)     &
                               - 9.0d0*rbufp(i,j,2)         &
                               + 2.0d0*rbufp(i,j,1)         )
        end do
      end do
    else if (nlayer .eq. 1) then
      do j=1,ny
        do i=1,nx
          u(i,j,nlayer) = elf*(  18.0d0*rbufm(i,j,3)        &
                               - 9.0d0*rbufm(i,j,2)         &
                               + 2.0d0*rbufm(i,j,1)         )
          v(i,j,nlayer) = elf*(  18.0d0*rbufp(i,j,3)        &
                               - 9.0d0*rbufp(i,j,2)         &
                               + 2.0d0*rbufp(i,j,1)         )
        end do
      end do
#endif
    else
      print*,'Not set up for case nlayer = ',nlayer
      call abort_comms(0)
    end if
  end if

!-------------------------------------------
! Finally lets do the temperature conditions
!-------------------------------------------
#ifdef MPI

!-------------------
! Initialise buffers
!-------------------
  msg = 0
  mstat = 0
  ireq = 0 
  rbufm = 0.0d0
  rbufp = 0.0d0
  rbufm_counter = 0
 
  if (nproc .gt. 1) then
!---------------------------------------------------
! If we're using radiative or constant flux boundary 
! conditions then more communication is needed.
! Use u buffers:
!---------------------------------------------------
    if ( radbc ) then
      if (myrank .eq. 0) then
        rbufm_counter = 0
        do k=2, nproc, 1
          if (b_distrecv(k) .gt. 0) then
            ncount = nx*ny*b_distrecv(k)
            call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,         &
                            MPI_DOUBLE_PRECISION, k-1,                  &
                            busendup_tag(k), comm, ireq(msg+1), ier) 
            msg = msg + 1
            rbufm_counter = rbufm_counter + b_distrecv(k)
          end if
        end do
      end if
      
      do k=2, nproc, 1
        if (myrank .eq. k-1) then
          if (b_distsend(1) .gt. 0) then
            ncount = nx*ny*b_distsend(1)
            call MPI_Isend( t(1,1,1), ncount, MPI_DOUBLE_PRECISION, 0,  &
                            busendup_tag(k), comm, ireq(msg+1), ier)
            msg = msg + 1
          end if
        end if
      end do  
    end if

!-----------------------------------------------------
! If we're using a constant flux boundary condition at 
! the lower surface then more communication is needed.
! Use u buffers:
!-----------------------------------------------------
    if ( constflux ) then
      if (myrank .eq. nproc-1) then
        rbufm_counter = 0
        do k=1, nproc-1
          if (b_distrecv(k) .gt. 0) then
            ncount = nx*ny*b_distrecv(k)
            call MPI_Irecv( rbufm(1,1,1+rbufm_counter), ncount,         &
                            MPI_DOUBLE_PRECISION, k-1,                  &
                            busenddown_tag(k),comm,ireq(msg+1),ier) 
            msg = msg + 1
            rbufm_counter = rbufm_counter + b_distrecv(k)
          end if
        end do
      end if

      do k=1, nproc-1
        if (myrank .eq. k-1) then
          if (b_distsend(nproc) .gt. 0) then
            ncount = nx*ny*b_distsend(nproc)
            call MPI_Isend( t(1,1,nlayer-b_distsend(nproc)+1), ncount,  &
                            MPI_DOUBLE_PRECISION, nproc-1,              &
                            busenddown_tag(k),comm,ireq(msg+1),ier)
            msg = msg + 1
          end if
        end if
      end do
    end if

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif

!------------------------------------
! Now do the upper boundary condition
!------------------------------------
  if (myrank .eq. 0) then
    if ( radbc ) then
!-------------------------------------------------------------
! This uses a radiative boundary condition which matches the 
! heat flux to a blackbody term. Initial conditions imply that 
! dT/dz= theta T^4. Solved by Newton-Raphson
!-------------------------------------------------------------
      if ( nlayer .ge. 4 ) then
        do j=1,ny
          do i=1,nx
            t1old = t(i,j,1)
            do inloop =1,10
              ft1 = 6.0d0*theta*dz*t1old*t1old*t1old*t1old+11.0d0*t1old &
                   -18.0d0*t(i,j,2)+9.0d0*t(i,j,3)-2.0d0*t(i,j,4)    
              dft1 = 24.0d0*theta*dz*t1old*t1old*t1old + 11.0d0
              t1old = t1old - (ft1/dft1)
            end do
            if ( abs(ft1/dft1) .gt. 1.0d-6 ) then 
              print*, 'ERROR: radbc not converged',i,j,(ft1/dft1)
            end if
            t(i,j,1) = t1old
          end do
        end do
#ifdef MPI
      else if ( nlayer .eq. 3 ) then
        do j=1,ny
          do i=1,nx
            t1old = t(i,j,1)
            do inloop =1,10
              ft1 = 6.0d0*theta*dz*t1old*t1old*t1old*t1old+11.0d0*t1old &
                   -18.0d0*t(i,j,2)+9.0d0*t(i,j,3)-2.0d0*rbufm(i,j,1)  
              dft1 = 24.0d0*theta*dz*t1old*t1old*t1old + 11.0d0
              t1old = t1old - (ft1/dft1)
            end do
            if ( abs(ft1/dft1) .gt. 1.0d-6 ) then 
              print*, 'ERROR: radbc not converged',i,j,(ft1/dft1)
            end if
            t(i,j,1) = t1old
          end do
        end do
      else if ( nlayer .eq. 2 ) then
        do j=1,ny
          do i=1,nx
            t1old = t(i,j,1)
            do inloop =1,10
              ft1 = 6.0d0*theta*dz*t1old*t1old*t1old*t1old+11.0d0*t1old &
                   -18.0d0*t(i,j,2)+9.0d0*rbufm(i,j,1)-2.0d0*rbufm(i,j,2)  
              dft1 = 24.0d0*theta*dz*t1old*t1old*t1old + 11.0d0
              t1old = t1old - (ft1/dft1)
            end do
            if ( abs(ft1/dft1) .gt. 1.0d-6 ) then 
              print*, 'ERROR: radbc not converged',i,j,(ft1/dft1)
            end if
            t(i,j,1) = t1old
          end do
        end do
      else if ( nlayer .eq. 1 ) then
        do j=1,ny
          do i=1,nx
            t1old = t(i,j,1)
            do inloop =1,10
              ft1 = 6.0d0*theta*dz*t1old*t1old*t1old*t1old+11.0d0*t1old &
                   -18.0d0*rbufm(i,j,1)+9.0d0*rbufm(i,j,2)-2.0d0*rbufm(i,j,3) 
              dft1 = 24.0d0*theta*dz*t1old*t1old*t1old + 11.0d0
              t1old = t1old - (ft1/dft1)
            end do
            if ( abs(ft1/dft1) .gt. 1.0d-6 ) then 
              print*, 'ERROR: radbc not converged',i,j,(ft1/dft1)
            end if
            t(i,j,1) = t1old
          end do
        end do
#endif
      else
        print*,'Not set up for case nlayer = ',nlayer
        call abort_comms(0)
      end if
    else
      do j=1,ny
        do i=1,nx
!--------
! fixed T
!--------
          t(i,j,1) = 1.0d0
        end do
      end do
    end if
  end if

!-------------------------------------------------
! Now do the lower boundary condition
!-------------------------------------------------
! two cases:
!  1)  fixed temperature gradient (fixed at theta)
!  2)  fixed temperature
!-------------------------------------------------
  if (myrank .eq. nproc-1) then
    if ( constflux ) then

!---------
! Case (1)
!---------
      if ( nlayer .ge. 4 ) then
        do j=1,ny
          do i=1,nx
            t(i,j,nlayer) = elf*(18.0d0*t(i,j,nlayer-1)               & 
                               - 9.0d0*t(i,j,nlayer-2)                &
                               + 2.0d0*t(i,j,nlayer-3) + 6.0d0*theta*dz)
          end do
        end do
#ifdef MPI
      else if ( nlayer .eq. 3 ) then
        do j=1,ny
          do i=1,nx
            t(i,j,nlayer) = elf*(18.0d0*t(i,j,nlayer-1)               & 
                               - 9.0d0*t(i,j,nlayer-2)                &
                               + 2.0d0*rbufm(i,j,1) + 6.0d0*theta*dz)
          end do
        end do
      else if ( nlayer .eq. 2 ) then
        do j=1,ny
          do i=1,nx
            t(i,j,nlayer) = elf*(18.0d0*t(i,j,nlayer-1)               & 
                               - 9.0d0*rbufm(i,j,2)                   &
                               + 2.0d0*rbufm(i,j,1) + 6.0d0*theta*dz)
          end do
        end do
      else if ( nlayer .eq. 1 ) then
        do j=1,ny
          do i=1,nx
            t(i,j,nlayer) = elf*(18.0d0*rbufm(i,j,3)                  & 
                               - 9.0d0*rbufm(i,j,2)                   &
                               + 2.0d0*rbufm(i,j,1) + 6.0d0*theta*dz)
          end do
        end do
#endif
      else
        print*,'Not set up for case nlayer = ',nlayer
        call abort_comms(0)
      end if
    else

!---------
! Case (2)
!---------
      do j=1,ny
        do i=1,nx
          t(i,j,nlayer) = 1.0d0 + theta
        end do
      end do
    endif
  end if
 
!--------------------------
! Increment timer if needed
!--------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_bound = time_bound + (time_tmp2-time_tmp)
  end if

end subroutine boundary
!================================================================
! Magnetic Boundary Conditions
!----------------------------------------------------------------
! Applies the boundary conditions to the magnetic field
! Perfect conductor: dBx/dz=dBy/dz=0, Bz=0
! Vertical field: Bx=By=0, dBz/dz=0 
!================================================================
Subroutine mag_boundary( )
  implicit none 

!-------------------
! Local declarations
!-------------------
  integer          :: i, j
#ifdef MPI 
  integer          :: k
  integer          :: ier
#endif

!----------------------------
! Starts timing off if needed
!----------------------------
  if (timer) then
    call findtime(time_tmp)
  end if

!--------------------------------
! Do Perfect conductor case first
!--------------------------------
  if (perfect) then

!--------------------------------------------
! Do communication for Bx and By if necessary
!--------------------------------------------
#ifdef MPI
    if (nproc .gt. 1) then

!-------------------------
! Start with mean profiles
!-------------------------

!-------------------
! Initialise buffers
!-------------------
      recvbuf = 0.0d0
      sendbuf = 0.0d0
      rbufp_counter = 0
      rbufm_counter = 0
      msg = 0
      mstat = 0
      ireq  = 0

!--------------------------------
! Post all receive commands first
!--------------------------------
      if (myrank .eq. 0) then
        rbufp_counter = 0
        rbufm_counter = 0
        do k=2, nproc, 1
          if (b_distrecv(k) .gt. 0) then
          ncount = b_distrecv(k)
          call MPI_Irecv( recvbuf(1+rbufp_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                  &
                          bvsendup_tag(k), comm, ireq(msg+1), ier )
          call MPI_Irecv( sendbuf(1+rbufm_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                  &
                          busendup_tag(k), comm, ireq(msg+2), ier ) 
          msg = msg + 2
          rbufp_counter = rbufp_counter + b_distrecv(k)
          rbufm_counter = rbufm_counter + b_distrecv(k)
          end if
        end do
      end if

      if (myrank .eq. nproc-1) then
        rbufp_counter = 0
        rbufm_counter = 0
        do k=1, nproc-1
          if (b_distrecv(k) .gt. 0) then
          ncount = b_distrecv(k)
          call MPI_Irecv( recvbuf(1+rbufp_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                  &
                          bvsenddown_tag(k), comm, ireq(msg+1), ier )
          call MPI_Irecv( sendbuf(1+rbufm_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                  &
                          busenddown_tag(k), comm, ireq(msg+2), ier ) 
          msg = msg + 2
          rbufp_counter = rbufp_counter + b_distrecv(k)
          rbufm_counter = rbufm_counter + b_distrecv(k)
          end if
        end do
      end if 

!-------------------
! Send commands here
!-------------------
      do k=2, nproc, 1
        if (myrank .eq. k-1) then
          if (b_distsend(1) .gt. 0) then
            ncount = b_distsend(1)
            call MPI_Isend( bxbar(1), ncount, MPI_DOUBLE_PRECISION, 0,    &
                            busendup_tag(k), comm, ireq(msg+1), ier )
            call MPI_Isend( bybar(1), ncount, MPI_DOUBLE_PRECISION, 0,    &
                            bvsendup_tag(k), comm, ireq(msg+2), ier )
            msg = msg + 2
          end if
        end if
      end do  

      do k=1, nproc-1
        if (myrank .eq. k-1) then
          if (b_distsend(nproc) .gt. 0) then
            ncount = b_distsend(nproc)
            call MPI_Isend( bxbar(nlayer-b_distsend(nproc)+1), ncount,    &
                            MPI_DOUBLE_PRECISION, nproc-1,                &
                            busenddown_tag(k), comm, ireq(msg+1), ier )
            call MPI_Isend( bybar(nlayer-b_distsend(nproc)+1), ncount,    & 
                            MPI_DOUBLE_PRECISION, nproc-1,                &
                            bvsenddown_tag(k), comm, ireq(msg+2), ier )
            msg = msg + 2
          end if
        end if
      end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
      call MPI_Waitall( msg, ireq, mstat, ier )

!--------------------------
! Now do standard variables
!--------------------------

!-------------------
! Initialise buffers
!-------------------
      rbufn = 0.0d0
      rbufn_counter = 0
      msg = 0 
      mstat = 0
      ireq = 0

!----------------------------
! Post receive commands first
!----------------------------
      if (myrank .eq. 0) then
        rbufn_counter = 0
        do k=2, nproc, 1
          if (b_distrecv(k) .gt. 0) then
          ncount = nx*ny*b_distrecv(k)
          call MPI_Irecv( rbufn(1,1,1+rbufn_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                    &
                          bwsendup_tag(k), comm, ireq(msg+1), ier ) 
          msg = msg + 1
          rbufn_counter = rbufn_counter + b_distrecv(k)
          end if
        end do
      end if
      
      if (myrank .eq. nproc-1) then
        rbufn_counter = 0
        do k=1, nproc-1
          if (b_distrecv(k) .gt. 0) then
          ncount = nx*ny*b_distrecv(k)
          call MPI_Irecv( rbufn(1,1,1+rbufn_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                    &
                          bwsenddown_tag(k), comm, ireq(msg+1), ier ) 
          msg = msg + 1
          rbufn_counter = rbufn_counter + b_distrecv(k)
          end if
        end do
      end if 

!-------------------
! Send commands here
!-------------------
      do k=2, nproc, 1
        if (myrank .eq. k-1) then
          if (b_distsend(1) .gt. 0) then
            ncount = nx*ny*b_distsend(1)
            call MPI_Isend( tor(1,1,1), ncount, MPI_DOUBLE_PRECISION, 0,   &
                            bwsendup_tag(k), comm, ireq(msg+1), ier )
            msg = msg + 1
          end if
        end if
      end do  

      do k=1, nproc-1
        if (myrank .eq. k-1) then
          if (b_distsend(nproc) .gt. 0) then
            ncount = nx*ny*b_distsend(nproc)
            call MPI_Isend( tor(1,1,nlayer-b_distsend(nproc)+1), ncount,   & 
                            MPI_DOUBLE_PRECISION, nproc-1,                 &
                            bwsenddown_tag(k), comm, ireq(msg+1), ier )
            msg = msg + 1
          end if
        end if
      end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
      call MPI_Waitall( msg, ireq, mstat, ier )
    end if
#endif

!-------------------------------------------------------
! Do the top layer boundary conditions if myrank is zero
!-------------------------------------------------------
    if ( myrank .eq. 0 ) then

!----------------------
! Sets bz equal to zero
!----------------------
      do j=1,ny
        do i=1,nx
          pol(i,j,1) = 0.0d0
        end do
      end do

!----------------------------------------
! zero derivative conditions on bx and by
!----------------------------------------
      if (nlayer .ge. 4) then
        bxbar(1) = elf*(  18.0d0*bxbar(2)      &
                        - 9.0d0*bxbar(3)       &
                        + 2.0d0*bxbar(4)       )
        bybar(1) = elf*(  18.0d0*bybar(2)      &
                        - 9.0d0*bybar(3)       &
                        + 2.0d0*bybar(4)       )
        do j=1,ny
          do i=1,nx
            tor(i,j,1) = elf*(  18.0d0*tor(i,j,2)    &
                              - 9.0d0*tor(i,j,3)     &
                              + 2.0d0*tor(i,j,4)     )
          end do
        end do

#ifdef MPI
      else if (nlayer .eq. 3) then
        bxbar(1) = elf*(  18.0d0*bxbar(2)        &
                        - 9.0d0*bxbar(3)         &
                        + 2.0d0*sendbuf(1)       )
        bybar(1) = elf*(  18.0d0*bybar(2)        &
                        - 9.0d0*bybar(3)         &
                        + 2.0d0*recvbuf(1)       )
        do j=1,ny
          do i=1,nx
            tor(i,j,1) = elf*(  18.0d0*tor(i,j,2)    &
                             - 9.0d0*tor(i,j,3)      &
                             + 2.0d0*rbufn(i,j,1)    )
          end do
        end do
      else if (nlayer .eq. 2) then
        bxbar(1) = elf*(  18.0d0*bxbar(2)        &
                        - 9.0d0*sendbuf(1)       &
                        + 2.0d0*sendbuf(2)       )
        bybar(1) = elf*(  18.0d0*bybar(2)        &
                        - 9.0d0*recvbuf(1)       &
                        + 2.0d0*recvbuf(2)       )
        do j=1,ny
          do i=1,nx
            tor(i,j,1) = elf*(  18.0d0*tor(i,j,2)    &
                            - 9.0d0*rbufn(i,j,1)     &
                            + 2.0d0*rbufn(i,j,2)     )
          end do
        end do
      else if (nlayer .eq. 1) then
        bxbar(1) = elf*(  18.0d0*sendbuf(1)      &
                        - 9.0d0*sendbuf(2)       &
                        + 2.0d0*sendbuf(3)       )
        bybar(1) = elf*(  18.0d0*recvbuf(1)      &
                        - 9.0d0*recvbuf(2)       &
                        + 2.0d0*recvbuf(3)       )
        do j=1,ny
          do i=1,nx
            tor(i,j,1) = elf*(  18.0d0*rbufn(i,j,1) &
                            - 9.0d0*rbufn(i,j,2)    &
                            + 2.0d0*rbufn(i,j,3)    )
          end do
        end do
#endif
      else
        print*,'Not set up for case nlayer = ',nlayer
        call abort_comms(0)
      end if

    end if

!---------------------------------------------
! Do the bottom layer if myrank + 1 equals the
! number of processors 
!---------------------------------------------
    if ( myrank + 1 .eq. nproc )  then

!--------------------------
! Boundary condition for bz
!--------------------------
      do j=1,ny
        do i=1,nx
          pol(i,j,nlayer) = 0.0d0
        end do
      end do

!-----------------------------------------
! Bottom boundary conditions for bx and by
!-----------------------------------------
      if (nlayer .ge. 4) then
        bxbar(nlayer) = elf*(  18.0d0*bxbar(nlayer-1)     &
                        - 9.0d0*bxbar(nlayer-2)      &
                        + 2.0d0*bxbar(nlayer-3)      )
        bybar(nlayer) = elf*(  18.0d0*bybar(nlayer-1)     &
                        - 9.0d0*bybar(nlayer-2)      &
                        + 2.0d0*bybar(nlayer-3)      )
        do j=1,ny
          do i=1,nx
            tor(i,j,nlayer) = elf*(  18.0d0*tor(i,j,nlayer-1)   &
                                  - 9.0d0*tor(i,j,nlayer-2)     &
                                  + 2.0d0*tor(i,j,nlayer-3)     )
          end do
        end do

#ifdef MPI
      else if (nlayer .eq. 3) then
        bxbar(nlayer) = elf*(  18.0d0*bxbar(nlayer-1)     &
                        - 9.0d0*bxbar(nlayer-2)      &
                        + 2.0d0*sendbuf(1)           )
        bybar(nlayer) = elf*(  18.0d0*bybar(nlayer-1)     &
                        - 9.0d0*bybar(nlayer-2)      &
                        + 2.0d0*recvbuf(1)           )
        do j=1,ny
          do i=1,nx
            tor(i,j,nlayer) = elf*(  18.0d0*tor(i,j,nlayer-1)   &
                                  - 9.0d0*tor(i,j,nlayer-2)     &
                                  + 2.0d0*rbufn(i,j,1)          )
          end do
        end do
      else if (nlayer .eq. 2) then
        bxbar(nlayer) = elf*(  18.0d0*bxbar(nlayer-1)     &
                        - 9.0d0*sendbuf(2)           &
                        + 2.0d0*sendbuf(1)           )
        bybar(nlayer) = elf*(  18.0d0*bybar(nlayer-1)     &
                        - 9.0d0*recvbuf(2)           &
                        + 2.0d0*recvbuf(1)           )
        do j=1,ny
          do i=1,nx
            tor(i,j,nlayer) = elf*(  18.0d0*tor(i,j,nlayer-1)   &
                                  - 9.0d0*rbufn(i,j,2)          &
                                  + 2.0d0*rbufn(i,j,1)          )
          end do
        end do
      else if (nlayer .eq. 1) then
        bxbar(nlayer) = elf*(  18.0d0*sendbuf(3)          &
                        - 9.0d0*sendbuf(2)           &
                        + 2.0d0*sendbuf(1)           )
        bybar(nlayer) = elf*(  18.0d0*recvbuf(3)          &
                        - 9.0d0*recvbuf(2)           &
                        + 2.0d0*recvbuf(1)           )
        do j=1,ny
          do i=1,nx
            tor(i,j,nlayer) = elf*(  18.0d0*rbufn(i,j,3)       &
                                  - 9.0d0*rbufn(i,j,2)         &
                                  + 2.0d0*rbufn(i,j,1)         )
          end do
        end do
#endif
      else
        print*,'Not set up for case nlayer = ',nlayer
        call abort_comms(0)
      end if
    end if

!--------------------------------
! Now do vertical field condition
!--------------------------------
  else 

!-------------------------------------
! Do communication for Bz if necessary
!-------------------------------------
#ifdef MPI
    if (nproc .gt. 1) then

!-------------------
! Initialise buffers
!-------------------
      recvbuf = 0.0d0
      rbufp = 0.0d0
      rbufp_counter = 0
      msg = 0
      mstat = 0
      ireq  = 0

!----------------------
! Mean bz not evolved!
! Do standard variables
!----------------------

!----------------------------
! Post receive commands first
!----------------------------
      if (myrank .eq. 0) then
        rbufp_counter = 0
        do k=2, nproc, 1
          if (b_distrecv(k) .gt. 0) then
          ncount = nx*ny*b_distrecv(k)
          call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                    &
                          bvsendup_tag(k), comm, ireq(msg+1), ier )
          msg = msg + 1
          rbufp_counter = rbufp_counter + b_distrecv(k)
          end if
        end do
      end if
      
      if (myrank .eq. nproc-1) then
        rbufp_counter = 0
        do k=1, nproc-1
          if (b_distrecv(k) .gt. 0) then
          ncount = nx*ny*b_distrecv(k)
          call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                          MPI_DOUBLE_PRECISION, k-1,                    &
                          bvsenddown_tag(k), comm, ireq(msg+1), ier )
          msg = msg + 1
          rbufp_counter = rbufp_counter + b_distrecv(k)
          end if
        end do
      end if 

!-------------------
! Send commands here
!-------------------
      do k=2, nproc, 1
        if (myrank .eq. k-1) then
          if (b_distsend(1) .gt. 0) then
            ncount = nx*ny*b_distsend(1)
            call MPI_Isend( pol(1,1,1), ncount, MPI_DOUBLE_PRECISION, 0,    &
                            bvsendup_tag(k), comm, ireq(msg+1), ier )
            msg = msg + 1
          end if
        end if
      end do  

      do k=1, nproc-1
        if (myrank .eq. k-1) then
          if (b_distsend(nproc) .gt. 0) then
            ncount = nx*ny*b_distsend(nproc)
            call MPI_Isend( pol(1,1,nlayer-b_distsend(nproc)+1), ncount,    & 
                            MPI_DOUBLE_PRECISION, nproc-1,                 &
                            bvsenddown_tag(k), comm, ireq(msg+1), ier )
            msg = msg + 1
          end if
        end if
      end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
      call MPI_Waitall( msg, ireq, mstat, ier )
    end if
#endif

!-------------------------------------------------------
! Do the top layer boundary conditions if myrank is zero
!-------------------------------------------------------
    if ( myrank .eq. 0 ) then

!---------------------------
! Sets bx=by=0 equal to zero
!---------------------------
      do j=1,ny
        do i=1,nx
          bxbar(1) = 0.0d0
          bybar(1) = 0.0d0
          tor(i,j,1) = 0.0d0
        end do
      end do

!---------------------------------
! zero derivative conditions on bz
!---------------------------------
      if (nlayer .ge. 4) then
        do j=1,ny
          do i=1,nx
            pol(i,j,1) = elf*(  18.0d0*pol(i,j,2)    &
                             - 9.0d0*pol(i,j,3)      &
                             + 2.0d0*pol(i,j,4)      )
          end do
        end do

#ifdef MPI
      else if (nlayer .eq. 3) then
        do j=1,ny
          do i=1,nx
            pol(i,j,1) = elf*(  18.0d0*pol(i,j,2)    &
                             - 9.0d0*pol(i,j,3)      &
                             + 2.0d0*rbufp(i,j,1)    )
          end do
        end do
      else if (nlayer .eq. 2) then
        do j=1,ny
          do i=1,nx
            pol(i,j,1) = elf*(  18.0d0*pol(i,j,2)    &
                            - 9.0d0*rbufp(i,j,1)     &
                            + 2.0d0*rbufp(i,j,2)     )
          end do
        end do
      else if (nlayer .eq. 1) then
        do j=1,ny
          do i=1,nx
            pol(i,j,1) = elf*(  18.0d0*rbufp(i,j,1) &
                            - 9.0d0*rbufp(i,j,2)    &
                            + 2.0d0*rbufp(i,j,3)   )
          end do
        end do
#endif
      else
        print*,'Not set up for case nlayer = ',nlayer
        call abort_comms(0)
      end if

    end if

!---------------------------------------------
! Do the bottom layer if myrank + 1 equals the
! number of processors 
!---------------------------------------------
    if ( myrank + 1 .eq. nproc )  then

!------------------------------
! Boundary condition for bx, by
!------------------------------
      do j=1,ny
        do i=1,nx
          bxbar(nlayer) = 0.0d0
          bybar(nlayer) = 0.0d0
          tor(i,j,nlayer) = 0.0d0
        end do
      end do

!-----------------------------------------
! Bottom boundary conditions for bx and by
!-----------------------------------------
      if (nlayer .ge. 4) then
        do j=1,ny
          do i=1,nx
            pol(i,j,nlayer) = elf*(  18.0d0*pol(i,j,nlayer-1)   &
                                  - 9.0d0*pol(i,j,nlayer-2)     &
                                  + 2.0d0*pol(i,j,nlayer-3)     )
          end do
        end do

#ifdef MPI
      else if (nlayer .eq. 3) then
        do j=1,ny
          do i=1,nx
            pol(i,j,nlayer) = elf*(  18.0d0*pol(i,j,nlayer-1)   &
                                  - 9.0d0*pol(i,j,nlayer-2)     &
                                  + 2.0d0*rbufp(i,j,1)          )
          end do
        end do
      else if (nlayer .eq. 2) then
        do j=1,ny
          do i=1,nx
            pol(i,j,nlayer) = elf*(  18.0d0*pol(i,j,nlayer-1)   &
                                  - 9.0d0*rbufp(i,j,2)          &
                                  + 2.0d0*rbufp(i,j,1)          )
          end do
        end do
      else if (nlayer .eq. 1) then
        do j=1,ny
          do i=1,nx
            pol(i,j,nlayer) = elf*(  18.0d0*rbufp(i,j,3)       &
                                  - 9.0d0*rbufp(i,j,2)         &
                                  + 2.0d0*rbufp(i,j,1)         )
          end do
        end do
#endif
      else
        print*,'Not set up for case nlayer = ',nlayer
        call abort_comms(0)
      end if
    end if
  end if

!--------------------------
! Increment timer if needed
!--------------------------
  if (timer) then
    call findtime(time_tmp2)
    time_bound = time_bound + (time_tmp2-time_tmp)
  end if

end subroutine mag_boundary

!=============================================
! This subroutine works out the Nusselt number
!=============================================
Subroutine nusselt(Nuss)
  implicit none

!-------------------------
! Declare output variables
!-------------------------
  real(kind=dp)   :: Nuss

!-------------------
! Local declarations
!-------------------
  real(kind=dp)   :: tzbar
  integer         :: i, j
#ifdef MPI
  integer         :: k
  integer         :: ier
#endif

!----------------------------
! Initialise output and tzbar
!----------------------------
  Nuss = 0.0d0
  tzbar = 0.0d0

!------------------------------------
! Do communication for T if necessary
!------------------------------------
#ifdef MPI
  if (nproc .gt. 1) then

!-------------------
! Initialise buffers
!-------------------
    rbufp = 0.0d0
    rbufp_counter = 0
    msg = 0
    mstat = 0
    ireq  = 0 

!----------------------------
! Post receive commands first
!----------------------------
    if (myrank .eq. nproc-1) then
      rbufp_counter = 0
      do k=1, nproc-1
        if (b_distrecv(k) .gt. 0) then
        ncount = nx*ny*b_distrecv(k)
        call MPI_Irecv( rbufp(1,1,1+rbufp_counter), ncount,           &
                        MPI_DOUBLE_PRECISION, k-1,                    &
                        bvsenddown_tag(k), comm, ireq(msg+1), ier )
        msg = msg + 1
        rbufp_counter = rbufp_counter + b_distrecv(k)
        end if
      end do
    end if 

!------------------
! Send command here
!------------------
    do k=1, nproc-1
      if (myrank .eq. k-1) then
        if (b_distsend(nproc) .gt. 0) then
          ncount = nx*ny*b_distsend(nproc)
          call MPI_Isend( t(1,1,nlayer-b_distsend(nproc)+1), ncount,    &
                          MPI_DOUBLE_PRECISION, nproc-1,                &
                          bvsenddown_tag(k), comm, ireq(msg+1), ier )
          msg = msg + 1
        end if
      end if
    end do

!--------------------------------------
! Check that communication has finished
!--------------------------------------
    call MPI_Waitall( msg, ireq, mstat, ier )
  end if
#endif 

!---------------------------------------
! Calculate mean value of dtdz at bottom
!---------------------------------------
  if ( myrank + 1 .eq. nproc )  then
    if (nlayer .ge. 4) then
      do j=1,ny
        do i=1,nx
          tzbar = tzbar + ( ( 1.0d0/( 6.0d0*dz ) )*         &
           ( 11.0d0*t(i,j,nlayer) - 18.0d0*t(i,j,nlayer-1)  &
           + 9.0d0*t(i,j,nlayer-2) - 2.0d0*t(i,j,nlayer-3)) )
        end do
      end do                    
#ifdef MPI
    else if (nlayer .eq. 3) then
      do j=1,ny
        do i=1,nx
          tzbar = tzbar + ( ( 1.0d0/( 6.0d0*dz ) )*         &
           ( 11.0d0*t(i,j,nlayer) - 18.0d0*t(i,j,nlayer-1)  &
           + 9.0d0*t(i,j,nlayer-2) - 2.0d0*rbufp(i,j,1)) )
        end do
      end do  
    else if (nlayer .eq. 2) then
      do j=1,ny
        do i=1,nx
          tzbar = tzbar + ( ( 1.0d0/( 6.0d0*dz ) )*         &
           ( 11.0d0*t(i,j,nlayer) - 18.0d0*t(i,j,nlayer-1)  &
           + 9.0d0*rbufp(i,j,2) - 2.0d0*rbufp(i,j,1)) )
        end do
      end do  
    else if (nlayer .eq. 1) then
      do j=1,ny
        do i=1,nx
          tzbar = tzbar + ( ( 1.0d0/( 6.0d0*dz ) )*         &
           ( 11.0d0*t(i,j,nlayer) - 18.0d0*rbufp(i,j,3)  &
           + 9.0d0*rbufp(i,j,2) - 2.0d0*rbufp(i,j,1)) )
        end do
      end do  
#endif
    else
      print*,'Not set up for case nlayer = ',nlayer
      call abort_comms(0)
    end if
    tzbar = tzbar * xn

!-------------------------
! Calculate Nusselt number
!-------------------------
    Nuss = (tzbar - nusstmp)/(theta - nusstmp)   
  end if


End Subroutine nusselt
#ifdef MPI
!===========================================================
! Initialisation routine for communication - layers below
! Works out layer distribution and stores it in given arrays
!===========================================================
Subroutine getbelow_init(reqd_layers, rankcounter, dsend, drecv )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  integer, intent(in)      :: reqd_layers
  integer, intent(in)      :: rankcounter
  integer, intent(inout)   :: dsend(nproc)
  integer, intent(inout)   :: drecv(nproc)

!-------------------
! Local declarations
!-------------------
  integer                  :: j
  integer                  :: templayers

!---------------------------------------
! Loop over processors below rankcounter
!---------------------------------------
  templayers=reqd_layers

!-------------
! Sanity check
!-------------
  if (rankcounter.eq.nproc) then
    return
  else
    j=rankcounter+1
    do while (templayers .gt. 0)
!-------------
! Sanity check
!-------------
      if (j.gt.nproc) then
        return
      end if
      if (layer_size(j) .ge. templayers) then
        if (myrank.eq.rankcounter-1) then
          drecv(j) = templayers
        else if (myrank.eq.j-1) then
          dsend(rankcounter) = templayers
        end if
        templayers = -1
      else 
        if (myrank.eq.rankcounter-1) then
          drecv(j) = layer_size(j)
        else if (myrank.eq.j-1) then
          dsend(rankcounter) = layer_size(j)
        end if
        templayers = templayers - layer_size(j)
      end if
      j=j+1
    end do
  end if

end subroutine getbelow_init

!===========================================================
! Initialisation routine for communication - layers above
! Works out layer distribution and stores it in given arrays
!===========================================================
Subroutine getabove_init(reqd_layers, rankcounter, dsend, drecv )
  implicit none

!-----------------------------------
! Declare input and output variables
!-----------------------------------
  integer, intent(in)      :: reqd_layers
  integer, intent(in)      :: rankcounter
  integer, intent(inout)   :: dsend(nproc)
  integer, intent(inout)   :: drecv(nproc)

!-------------------
! Local declarations
!-------------------
  integer                  :: j
  integer                  :: templayers

!---------------------------------------
! Loop over processors above rankcounter
!---------------------------------------
  templayers=reqd_layers

!-------------
! Sanity check
!-------------
  if (rankcounter.eq.1) then
    return
  else
    j=rankcounter-1
    do while (templayers .gt. 0)

!-------------
! Sanity check
!-------------
      if (j.lt.1) then
        return
      end if
      if (layer_size(j) .ge. templayers) then
        if (myrank.eq.rankcounter-1) then
          drecv(j) = templayers
        else if (myrank.eq.j-1) then
          dsend(rankcounter) = templayers
        end if
        templayers = -1
      else 
        if (myrank.eq.rankcounter-1) then
          drecv(j) = layer_size(j)
        else if (myrank.eq.j-1) then
          dsend(rankcounter) = layer_size(j)
        end if
        templayers = templayers - layer_size(j)
      end if
      j=j-1
    end do
  end if

end subroutine getabove_init
#endif

END MODULE dbydz
