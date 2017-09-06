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
  use trans_params
  use dimens
  use control
  use constants
  use fft
  use dbydz
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
  real(kind=dp) :: ttime
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
  constflux = .false.
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
    print*, ' SF        = ', sf
    print*, ' NFRAME    = ', nframe
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

!-----------------
! Initialise dbydz
!-----------------
  call init_dbydz(ier) 
  if ( ier .lt. 0 ) then
    print*,'TROUBLE allocating dybdz arrays', ier
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

  if ( myrank .eq. 0 ) then
    print*,'NFRAME = ', nframe
    print*,'STARTING VALUE FOR ITER = ', iter
  end if

!---------------------------
! Read in initial conditions
!---------------------------
  if ( nframe .eq. 0 ) then
    call static()
!========================
! NGP Aug 30 2017
! initialize kappaZ and dkappadZ
    call    ini_trans_params()
!========================
  else 
    call ReadNetCDF(ier)
    if ( ier .ne. 0 ) then
      print*,'Trouble Reading NetCDF files '
      call abort_comms(ier)
    end if
  end if

!-------------------------------------------
! Print the inital greetings and information
!-------------------------------------------
  if ( myrank .eq. 0 ) then
    print*,' -------------------------------'
    print*,' nx ', nx, ' ny ', ny, ' nz ', nz
    print*,' -------------------------------'
    print*,'This calculation requires ', ntotal, ' iterations.'
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

!-----------------------------------------
! Set up poloidal toroidal components of B
!-----------------------------------------
  if (f .eq. 0) then
    if (kinematic) then
      call mag_init()
    end if
  else 
    call mag_init()
  end if

!---------------------------------------
! This is the beginning of the main loop
!---------------------------------------
  nit = 0

1000 nit = nit + 1

  call findtime(start_iteration)

  iter  = iter + 1 
  ttime = ttime + dt

!------------------------------
! Evolve hydrodynamic variables
!------------------------------
  call hydro_deriv()

! ----------------------------------------------------------------
! Check to see whether or not we need to evolve any magnetic terms
! ----------------------------------------------------------------
  if (f .eq. 0) then

!-----------------------------------
! Evolve magnetic terms if kinematic
!-----------------------------------
    if (kinematic) then
      call mag_deriv()  
    else
      do k=1,nlayer
        do j=1,ny
          do i=1,nx
            poldot(i,j,k,3) = 0.0d0
            tordot(i,j,k,3) = 0.0d0
          end do
        end do
        bxbardot(k,3) = 0.0d0
        bybardot(k,3) = 0.0d0
      end do
    end if
  else     
    call mag_deriv()
  end if

!---------------------------------
! Calculate time-step coefficients
!---------------------------------
  if ( first ) then
    first = .false.
    second = .true.
    dt1 = 0.0d0
    dt2 = 0.0d0
    c0 = dt
    c1  = 0.0d0
    c2 = 0.0d0
  else if ( second ) then 
    first = .false.
    second = .false.
    dt2 = 0.0d0
    c0 = dt*( 1.0d0 + 0.5*(dt/dt1))
    c1 = -0.5d0*dt*(dt/dt1)
    c2 = 0.0d0
  else
    c0 = dt + ( ( (dt*dt)/(dt1*dt2) )*( (dt/3.0d0)       &
         + 0.5d0*(dt1+dt2) ) )
    c1 = ( (dt*dt)/( dt1*(dt1-dt2) ) )*( (dt/3.0d0)      &
         + 0.5d0*dt2 )
    c2 = -( (dt*dt)/( dt2*(dt1-dt2) ) )*( (dt/3.0d0)     &
         + 0.5d0*dt1 )
  end if

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!-------------------------------
! The time-stepping is done here
!-------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx

!--------------------
! Update 3D variables
!--------------------
        u(i,j,k) = u(i,j,k) + c0*udot(i,j,k,3)                &
                   + c1*udot(i,j,k,2) + c2*udot(i,j,k,1)
        v(i,j,k) = v(i,j,k) + c0*vdot(i,j,k,3)                &
                   + c1*vdot(i,j,k,2) + c2*vdot(i,j,k,1)
        w(i,j,k) = w(i,j,k) + c0*wdot(i,j,k,3)                &
                   + c1*wdot(i,j,k,2) + c2*wdot(i,j,k,1)
        r(i,j,k) = r(i,j,k) + c0*rdot(i,j,k,3)                &
                   + c1*rdot(i,j,k,2) + c2*rdot(i,j,k,1)
        t(i,j,k) = t(i,j,k) + c0*tdot(i,j,k,3)                &
                   + c1*tdot(i,j,k,2) + c2*tdot(i,j,k,1)
        pol(i,j,k) = pol(i,j,k) + c0*poldot(i,j,k,3)          &
                     + c1*poldot(i,j,k,2) + c2*poldot(i,j,k,1)
        tor(i,j,k) = tor(i,j,k) + c0*tordot(i,j,k,3)          &
                     + c1*tordot(i,j,k,2) + c2*tordot(i,j,k,1)

!--------------------------------------------
! Make sure the density doesn't get too small
!--------------------------------------------
        if (r(i,j,k) .lt. 0.05d0) then
          r(i,j,k)=0.05d0
          rdot(i,j,k,1)=0.0d0
          rdot(i,j,k,2)=0.0d0
	  rdot(i,j,k,3)=0.0d0
        end if

!----------------------
! Update 3D derivatives
!----------------------
        udot(i,j,k,1) = udot(i,j,k,2)
        udot(i,j,k,2) = udot(i,j,k,3)
        vdot(i,j,k,1) = vdot(i,j,k,2)
        vdot(i,j,k,2) = vdot(i,j,k,3)
        wdot(i,j,k,1) = wdot(i,j,k,2)
        wdot(i,j,k,2) = wdot(i,j,k,3)
        rdot(i,j,k,1) = rdot(i,j,k,2)
        rdot(i,j,k,2) = rdot(i,j,k,3)
        tdot(i,j,k,1) = tdot(i,j,k,2)
        tdot(i,j,k,2) = tdot(i,j,k,3)
        poldot(i,j,k,1) = poldot(i,j,k,2)
        poldot(i,j,k,2) = poldot(i,j,k,3)
        tordot(i,j,k,1) = tordot(i,j,k,2)
        tordot(i,j,k,2) = tordot(i,j,k,3)
      end do
    end do

!-------------------------------------------------------
! Update 1D (mean) fields and mean-field derivatives
! Note: bzbar is constant with these boundary conditions
!-------------------------------------------------------
    bxbar(k) = bxbar(k) + c0*bxbardot(k,3)               &
                    + c1*bxbardot(k,2) + c2*bxbardot(k,1)
    bybar(k) = bybar(k) + c0*bybardot(k,3)               &
                    + c1*bybardot(k,2) + c2*bybardot(k,1)
    bxbardot(k,1) = bxbardot(k,2)
    bxbardot(k,2) = bxbardot(k,3)
    bybardot(k,1) = bybardot(k,2)
    bybardot(k,2) = bybardot(k,3)
  end do

!-------------------
! Increment t15_time
!-------------------
  call findtime(time_tmp2)
  t15_time = t15_time + (time_tmp2 - time_tmp)

!-----------------------------------
! Work out new values of dt1 and dt2
!-----------------------------------
   dt2 = dt + dt1
   dt1 = dt

!--------------------------------------------
! Dealias arrays every 10 timesteps if needed
! Do magnetic dealiasing every timestep
!--------------------------------------------
   if ( dalias ) then
     if ( mod(nit,10) .eq. 0 ) then 
       call dealias(r)
       call dealias(t)
       call dealias(u)
       call dealias(v)
       call dealias(w)
     end if
     if (f .eq. 0) then
       if (kinematic) then
         call magdealias()
       end if
     else
       call magdealias()
     end if
   end if

!---------------------
! Reconstitute B and J
!---------------------
   if (f .eq. 0) then
     if (kinematic) then
       call bandj()
     end if   
   else  
     call bandj()
   end if

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!------------------------------
! Now apply boundary conditions
!------------------------------
   call boundary()
   if (f .eq. 0) then
     if (kinematic) then
       call mag_boundary()
     end if
   else 
     call mag_boundary()
   end if

!-------------------
! Increment t17_time
!-------------------
  call findtime(time_tmp2)
  t17_time = t17_time + (time_tmp2 - time_tmp)

!------------------------------------------------------
! Check whether it is time for some output
!----------------------------------------------------------------------------
! Summary data goes to standard output
! Unit 16: Dump of all 8 arrays for restarts and movies
!----------------------------------------------------------------------------
    if ( mod(nit, ndu06) .eq. 0 ) then
!----------------------------------------------------------------------------

      if ( myrank .eq. 0 ) then
          call fdate(adate)
          adate = adate(12:19)
          write(*,*)
          write(*,*) adate,' NIT= ', nit, ' TTIME ', ttime
      endif

      call GlobalAverages( ubar, vbar, kenergy, menergy)
      if ( myrank .eq. 0 ) then
        print*,' ubar= ', ubar, ' vbar= ', vbar
      end if

      if (kinematic) then
        menergy = menergy
      else
        menergy = f*menergy
      end if
      
      call nusselt(nuss)
#ifdef MPI
       CALL MPI_Bcast( nuss, 1, MPI_DOUBLE_PRECISION, nproc-1, comm, ier )
#endif

      call checkdivb( divbmax, divbmin )
      if ( myrank .eq. 0 ) then
        print*,' divbmin= ', divbmin, ' divbmax= ', divbmax
      end if

      divbtest = max( abs(divbmax), abs(divbmin) )
      if ( divbtest .gt. 0.01d0 ) then
        if ( myrank .eq. 0 ) then
          print*,' ERROR: DIV B is not zero '
          print*,' DIV B = ', divbtest
        end if
        call abort_comms(12)
      end if


      call sminmax(r, rmin, rmax)
      call sminmax(t, tmin, tmax)
      call sminmax(u, umin, umax)
      call sminmax(v, vmin, vmax)
      call sminmax(w, wmin, wmax)
      call sminmax(bx, bxmin, bxmax )
      call sminmax(by, bymin, bymax )
      call sminmax(bz, bzmin, bzmax )
      wk1 = 0.0d0
      wk2 = 0.0d0
      wk3 = 0.0d0
      do k=1,nlayer
       do j=1,ny
         do i=1,nx
           wk1(i,j,k)= ( u(i,j,k)*u(i,j,k) + &
                         v(i,j,k)*v(i,j,k) + &
                         w(i,j,k)*w(i,j,k) ) / t(i,j,k)
           wk2(i,j,k)= ( bx(i,j,k)*bx(i,j,k) + &
                         by(i,j,k)*by(i,j,k) + &
                         bz(i,j,k)*bz(i,j,k) ) / r(i,j,k)
           wk3(i,j,k)= ( u(i,j,k)*u(i,j,k) + &
                         v(i,j,k)*v(i,j,k) + &
                         w(i,j,k)*w(i,j,k) )
        end do
       end do
      end do
      call sminmax(wk1, machmin, machmax)
      call sminmax(wk2, alfvenmin, alfvenmax)
      call sminmax(wk3, rmsmin, rmsmax)
      if ( myrank .eq. 0 ) then
        if (kinematic) then
          print*,' alfvenmin= ', sqrt( alfvenmin ), &
                 ' alfvenmax= ', sqrt( alfvenmax )
        else   
          print*,' alfvenmin= ', sqrt( f*alfvenmin ), &
                 ' alfvenmax= ', sqrt( f*alfvenmax )
        end if
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
        print*,' Nusselt= ', Nuss
        print*,' ' 
      end if
     
!----------------------------------------------------------------------------
    end if
!----------------------------------------------------------------------------
! Check if it is time to dump to unit 16, the restart/movie file
! A different one is written each time
!
! If iprec is set to 8 the output is written as doubles and 
! as reals if iprec=4.
!
! By default the results will be saved as "single precision" however
! when nit .ge. ntotal they will be saved as doubles. Ideally the
! derivatives and previous time step should also be saved so that
! the job can be restarted as though it was not stopped.
!
! Every "ndurestart" storages, double precision data is dumped to 
! file alongside double precision derivatives - can now do a clean restart
! 
!----------------------------------------------------------------------------
! First this is a good place to calculate the new time-step.
! Currently re-evaluated every ten steps and prints out whenever 
! nit = ndu06
!--------------------------------------------------------------------------

   if ( mod(nit,10) .eq. 0 ) then

     call calcdt( dt, mod(nit,ndu06) )
     if ( mod(nit,ndu06) .eq. 0 ) then
 
       if ( myrank .eq. 0 ) then
         print*,' New value of dt ', dt
       end if

     end if

   end if

!------------------------------------------------------------------------
! Setting itype to 1 dimensions the big arrays as (nx,ny,nz) 
!------------------------------------------------------------------------

    itype = 1

!------------------------------------------------------------------------
! Set the default output precision to 8 byte reals
!------------------------------------------------------------------------

    iprec = 4
    if ( mod(nit, ndu16) .eq. 0 ) then

        if ( nit .ge. ntotal ) then
          iprec = 8
        end if

        if ( mod(nit, ndurestart*ndu16) .eq. 0 ) then
          iprec = 8
        end if

        !------------------------------
        ! Generate the name of the file
        !------------------------------


        i = mod(iter/ndu16,1000)
        hundreds = char(i/100 + 48 )
        tens     = char(mod(i,100)/10 + 48 )
        units    = char(mod(i,10) + 48 )
        name = 'T16_'//HUNDREDS//TENS//UNITS//'.nc'
        call fdate(adate)
        call hostnm(hname)
        history = 'Created '//adate//' on '//hname//' by Paul Bushby'

!-------------------------------
! Write T16 file
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

    end if
   
!-------------------------------------------------------------
! Finish time-step and go back to start of loop if necessary
!-------------------------------------------------------------
   call findtime(finish_iteration)

   if  ( nit .lt. ntotal ) go to 1000
   
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
