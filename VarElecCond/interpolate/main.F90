!-----------------------------------------------------------------------
! F90 version or p3d
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 3-Dimension compressible magnetoconvection program:
! 
! Solves the equations of compressible magnetohydrodynamics in 3D
! cartesian geometry. Box is periodic in x and y; z=0 corresponds to
! the top of the box and z=1 corresponds to the bottom. 
! 
! Horizontal derivatives are evaluated in Fourier space (FFTW)
! Vertical derivatives are 4th Order finite-difference
!
! Time-stepping is done via 3rd order Adams-Bashforth
!
! Rotation added - magnitude determined by mid-layer Taylor number
! ("TAYL") - orientation determined by "PSI". 
! Unit rotation vector = ( 0 , cos (PSI) , -sin (PSI) )
!
! Static state is a polytrope with all the usual parameters
! (See Matthews, Proctor and Weiss (1995) JFM, 305, 281)
! INPUT contains all the necessary input parameters
!
! "SF" is the time-stepping safety factor (typically SF=0.2 is needed)
! "TOTAL" is the number of time-steps to be taken
!
! Data is dumped to terminal every "NDU06" steps
! Full data files dumped in NetCDF format every "NDU16" steps
!
! Every "NDURESTART" data dumps, the files are stored as
! double precision for clean restart
! 
! If "NFRAME" equals 0 then static start
! 
! If "NFRAME" equals 1 then reads initial condition from 'restart.nc'
!
! If "NFRAME" equals 2 then reads intial condition from 'restart.nc',
! other parameters from INPUT
!  
! Boundary conditions: (top and bottom)
! =====================================
! Velocity: Impermeable and stress-free
! Magnetic: Vertical field (default)
!            or perfect conductor (UNTESTED!)
! Temperature: Fixed top and bottom (default)
!              or radiative ("radbc": top) 
!              or constant flux ("constflux": bottom)
!
! ---------------------------------------------------------------------- 
! Paralellisation information (MPI):
!
! Distribution over nodes is done via:
! isum = 0
! do iproc=1,nproc
!   istart = ((iproc-1)*nz)/nproc + 1
!   ifinish = (iproc*nz)/nproc
!   print*,'iproc ', iproc, 'istart ', istart, ' ifinish ', &
!    ifinish, ifinish - istart+1
!   isum = isum + (ifinish -istart+1)
! end do
!
!
!----------------------------------------------------------------------- 
! NOTE: Needs pre-processing in order to sort out compiler flags
!-----------------------------------------------------------------------
PROGRAM pd3
  use types
  use comms
  use params
  use dimens
  use control
  use constants
  use fft
  use arrays
  use timing
  use initial
  use utils 
  use netcdf

  implicit none

  integer       :: ier
  integer       :: i, j, k
  integer       :: nit
  integer       :: iter
  integer       :: iproc
  integer       :: ist, ifn
  integer       :: isum

  real(kind=dp) :: dt, dt1, dt2
  real(kind=dp) :: c0, c1, c2
  real(kind=dp) :: ttime, z1
  real(kind=dp) :: ubar, vbar, kenergy, menergy, Nuss
  real(kind=dp) :: divbmax, divbmin, divbtest
  real(kind=dp) :: rmin, tmin, umin, vmin, wmin
  real(kind=dp) :: rmax, tmax, umax, vmax, wmax 
  real(kind=dp) :: bxmin, bymin, bzmin
  real(kind=dp) :: bxmax, bymax, bzmax
  real(kind=dp) :: machmin, machmax
  real(kind=dp) :: alfvenmin, alfvenmax
  real(kind=dp) :: rmsmin, rmsmax

!--------------------------
! Set the control variables
!--------------------------
  dalias = .true.
  kinematic  = .false.
  first   = .true.
  second  = .false.
  radbc = .false.
  constflux = .true.
  perfect = .false.
  timer = .false.
  dynamo = .false.

!---------------
! Version number
!---------------
  version = 3

!--------------------------
! This sets up MPI routines
!--------------------------
  call init_comms(ier)
  if ( ier .lt. 0 ) then
    print*,'TROUBLE allocating comms arrays', ier
    call abort_comms(ier)
  end if

!------------------------------------------------------
! Greetings - prints out info regarding parallelisation
!------------------------------------------------------
  if ( myrank .eq. 0 ) then
    print*,'---------------------------------------------------'
    print*,'Simulations of compressible 3-D magnetoconvection  '  
    print*,'---------------------------------------------------'

    if ( nproc .eq. 1 ) then
      print*,'Serial version of the code'
    else 
      print*,'Parallel version of the code running on ' , nproc, ' nodes'
    end if
    print *,' '
    print*, '3rd-order Adams-Bashforth time-stepping scheme'
    print*, 'Momentum equation in non-conservative form'
    if ( dalias) then 
      print*, 'De-aliasing is switched on'
    else
      print*, 'No de-aliasing'
    end if 
    if ( constflux ) then
      print *, 'Fixed-flux temperature boundary condition'
    else
      print *, 'Fixed temperature lower boundary condition'
    end if
    if ( radbc ) then
      print *, 'Radiative upper boundary condition'
    else
      print *, 'Fixed temperature upper boundary condition'
    end if
    print *,''
    print *,'FFTW libraries used for horizontal derivatives'
  end if 

!---------------------
! Start timing routine
!---------------------
  call init_timing()

!-----------------------------------------------------
! Read parameters from input file or from restart file
!-----------------------------------------------------
  call readparam( iter, ttime, dt1, ier )
  if ( ier .lt. 0 ) then
    print*,'ERROR: Unable to read parameters or restart file ', ier
    call abort_comms(1)
  end if

!-----------------
! Print parameters
!-----------------
  if ( myrank .eq. 0 ) then
    print*,' --------------------------------------'
    print*,' Array dimensions for this calculation:'
    print*,'    nx = ', nx
    print*,'    ny = ', ny
    print*,'    nz = ', nz
    print*,' --------------------------------------'
    print*, ' GAMMA     = ', gamma
    print*, ' SIGMA     = ', sigma
    print*, ' THETA     = ', theta
    print*, ' CK        = ', ck
    print*, ' RAYL      = ', rayl
    print*, ' CM        = ', cm
    print*, ' F         = ', f
    print*, ' CHAND     = ', chand
    print*, ' TAYL      = ', tayl
    print*, ' PSI       = ', psi
    print*, ' ZETA      = ', zeta
    print*, ' XMAX      = ', xmax
    print*, ' YMAX      = ', ymax
    print*,' --------------------------------------'
  end if

!--------------------------------------------------
! Set distribution of layers amongst the processors 
!--------------------------------------------------
  izstart  = ((myrank)*nz)/nproc + 1
  izfinish = ((myrank+1)*nz)/nproc
  nlayer   = izfinish - izstart + 1

  if (( nproc .eq. 1 ).and.( nlayer .lt. 6 )) then
    print*,'ERROR: need at least six layers for d2dybdz '
    call abort_comms(1)
  end if

  if (( nproc .eq. 2 ).and.( nlayer .lt. 3 )) then
    print*,'ERROR: need at least six layers for d2dybdz '
    call abort_comms(1)
  end if

!----------------------------------------
! Confirms the dimensioning of the arrays
!----------------------------------------
! Main arrays are of length: nx*ny*nlayer 
!----------------------------------------
  if ( myrank .eq. 0 ) then
    print*,'DIMENSION of arrays ', nx*ny*nlayer
  end if

!------------------
! Set the constants
!------------------
   call set_constants()

!---------------
! Initialize FFT
!---------------
  call init_fft(ier)
  if ( ier .lt. 0 ) then
    print*,'TROUBLE allocating FFT arrays', ier
    call abort_comms(ier)
  end if

!-------------------------------------
! Allocate and initialise large arrays
!-------------------------------------
  call alloc_arrays(ier)
  if ( ier .lt. 0 ) then
    print*,'TROUBLE allocating large arrays', ier
    call abort_comms(ier)
  end if

!---------------------------
! Read in initial conditions
!---------------------------
  call ReadNetCDF(ier)
  if ( ier .ne. 0 ) then
    print*,'Trouble Reading NetCDF files '
    call abort_comms(ier)
  end if

!-------------------------------------------
! Print the inital greetings and information
!-------------------------------------------
  if ( myrank .eq. 0 ) then
    print*,' -------------------------------'
    print*,' nx ', nx, ' ny ', ny, ' nz ', nz
    print*,' -------------------------------'
    print*,' '
    print*,'Distribution over processors is by layer.'
    print*,'Processor   Initial   Final layer   Total '
    print*,'---------   -------   -----------   ----- '
    isum = 0
    do iproc=1,nproc
      ist = ((iproc-1)*nz)/nproc + 1
      ifn = (iproc*nz)/nproc
      write(*,'(2x,i4,7x,i4,8x,i4,10x,i4)') iproc,ist,ifn,ifn - ist + 1
      isum = isum + ( ifn - ist + 1 )
    end do
    write(*,'(35x,i4)') isum
    print*,' ' 
    print*,' --------------------------------------------------------- '
    print*,' ' 
  end if

!------------------------------------------------------------------
! Calculate timestep DT: find velocity + sound speed + Alfven speed
!------------------------------------------------------------------
  call calcdt( dt, 0 )
  if ( myrank .eq. 0 ) then
    print*,'Initial Timestep, dt = ', dt
    print*,''
  end if

!----------------------
! Do interpolation here
!----------------------
  call expand(r,rnew)
  call expand(t,tnew)
  call expand(u,unew)
  call expand(v,vnew)
  call expand(w,wnew)
  call expand(bx,bxnew)
  call expand(by,bynew)
  call expand(bz,bznew)
  
!------------------
! Print some output
!------------------
  if ( myrank .eq. 0 ) then
    call fdate(adate)
    adate = adate(12:19)
    write(*,*)
    write(*,*) adate,' TTIME ', ttime
  endif

  call GlobalAverages( ubar, vbar, kenergy, menergy)
  if ( myrank .eq. 0 ) then
    print*,' ubar= ', ubar, ' vbar= ', vbar
  end if
      
  call sminmax(rnew, rmin, rmax)
  call sminmax(tnew, tmin, tmax)
  call sminmax(unew, umin, umax)
  call sminmax(vnew, vmin, vmax)
  call sminmax(wnew, wmin, wmax)
  call sminmax(bxnew, bxmin, bxmax )
  call sminmax(bynew, bymin, bymax )
  call sminmax(bznew, bzmin, bzmax )
  wk1new = 0.0d0
  wk2new = 0.0d0
  wk3new = 0.0d0
  do k=1,nlayer
    do j=1,2*ny
      do i=1,2*nx
        wk1new(i,j,k)= ( unew(i,j,k)*unew(i,j,k) + &
                         vnew(i,j,k)*vnew(i,j,k) + &
                         wnew(i,j,k)*wnew(i,j,k) ) / tnew(i,j,k)
        wk2new(i,j,k)= ( bxnew(i,j,k)*bxnew(i,j,k) + &
                         bynew(i,j,k)*bynew(i,j,k) + &
                         bznew(i,j,k)*bznew(i,j,k) ) / rnew(i,j,k)
        wk3new(i,j,k)= ( unew(i,j,k)*unew(i,j,k) + &
                         vnew(i,j,k)*vnew(i,j,k) + &
                         wnew(i,j,k)*wnew(i,j,k) )
      end do
    end do
  end do
  call sminmax(wk1new, machmin, machmax)
  call sminmax(wk2new, alfvenmin, alfvenmax)
  call sminmax(wk3new, rmsmin, rmsmax)
  if (myrank .eq. 0) then
    print*,' alfvenmin= ', sqrt( f*alfvenmin ), &
           ' alfvenmax= ', sqrt( f*alfvenmax )
    print*,' machmin= ', sqrt(machmin/gamma), &
           ' machmax= ', sqrt(machmax/gamma)
    print*,' rmsmin= ', sqrt(rmsmin), &
           ' rmsmax= ', sqrt(rmsmax)
    print*,' rmin= ', rmin, ' rmax= ', rmax
    print*,' tmin= ', tmin, ' tmax= ', tmax
    print*,' umin= ', umin, ' umax= ', umax
    print*,' vmin= ', vmin, ' vmax= ', vmax
    print*,' wmin= ', wmin, ' wmax= ', wmax
    print*,' bxmin= ', bxmin, ' bxmax= ', bxmax
    print*,' bymin= ', bymin, ' bymax= ', bymax
    print*,' bzmin= ', bzmin, ' bzmax= ', bzmax
    print*,''
    print*,' magenergy= ', menergy
    print*,' kinenergy= ', kenergy
    print*,' ' 
  end if

!------------------------------------------------------------------------
! Setting itype to 1 dimensions the big arrays as (nx,ny,nz) 
!------------------------------------------------------------------------
  itype = 1

!------------------------------------------------------------------------
! Set the default output precision to 8 byte reals
!------------------------------------------------------------------------
  iprec = 8

!------------------
! Generate filename
!------------------
  name = 'output.nc'

!-------------------------------
! Write output file
!-------------------------------
  call CreateNetCDF(ier)
  if ( ier .ne. 0 ) then
    call abort_comms(13)
  end if

  call WriteNetCDF(iter, ttime, dt1, ier)
  if ( ier .ne. 0 ) then
    call abort_comms(14)
  end if

  call CloseNetCDF(ier)
  if ( ier .ne. 0 ) then
    call abort_comms(15)
  end if
   
!----------------------------------------------------------------------
! Tidy up and finish
!----------------------------------------------------------------------

   if ( myrank .eq. 0 ) then
     call print_timing()
   end if

!---------------------------------------------------
! These routines shut the code down
!---------------------------------------------------

   call close_arrays()
   call close_fft()
   call exit_comms()


 END PROGRAM pd3
