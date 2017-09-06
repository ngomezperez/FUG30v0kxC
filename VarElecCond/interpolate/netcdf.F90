!======================================
! This module holds the NetCDF routines
!======================================
MODULE netcdf
  use types
  use dimens
  use params
  use arrays
  use comms

  implicit none
  include 'netcdf.inc'

  character(1)  :: hundreds, tens, units 
  character(4)  :: varname
  character(10) :: restart_file
  character(80) :: name 
  character(80) :: history 
  character(24) :: adate 
  character(25) :: hname 

  integer       :: ncid 
  integer       :: itype 
  integer       :: iprec 
  integer       :: ncdfstatus
  integer       :: nxid, nyid, nzid
  integer       :: varid
  integer       :: version  
  integer       :: icount(3)
  integer       :: istarting(3)

CONTAINS

!================================================================
! This subroutine creates a NetCDF file
! It also creates the array dimensions and puts it into data mode
!================================================================
SUBROUTINE CreateNetCDF( ier ) 
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer,       intent(out) ::  ier

!-------------------
! Local declarations
!-------------------
  integer :: nxdim, nydim, nzdim
  integer :: onedim
  integer :: dims(3)
  integer :: lenhis

!------------------------
! Sets dummy value of ier
!------------------------
  ier = -1

!------------------------------------
! Checks to make sure that myrank = 0
!------------------------------------
  if ( myrank .ne. 0 ) then
    ier = 0
    return
  end if

!---------------------------------------------------------------------
! Create NetCDF file and put into define mode
! ncdfstatus command automatically puts into define mode
! name is file name, NF_CLOBBER overwrites existing files of same name
! ncid is NetCDF id. Returns NF_NOERR if all goes well.
!---------------------------------------------------------------------
  ncdfstatus = nf_create( name, IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!------------------------------------------------------------
! Now Create attributes
! NF_GLOBAL indicates a global attribute
! NF_INT implies that it is of integer type, "1" gives length
!------------------------------------------------------------
  ncdfstatus = nf_put_att_int( ncid, NF_GLOBAL, 'version',      &
                                   NF_INT, 1, version )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  lenhis = 79
  ncdfstatus = nf_put_att_text( ncid, NF_GLOBAL, 'history',     &
                                   lenhis, history )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'gamma',     &
                                   NF_DOUBLE, 1, gamma )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'sigma',     &
                                   NF_DOUBLE, 1, sigma )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'theta',     &
                                   NF_DOUBLE, 1, theta )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'ck',        &
                                   NF_DOUBLE, 1, ck )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'cm',        &
                                   NF_DOUBLE, 1, cm )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'f',         &
                                   NF_DOUBLE, 1, f )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'tayl',      &
                                   NF_DOUBLE, 1, tayl )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'psi',       &
                                   NF_DOUBLE, 1, psi )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'zeta',      &
                                   NF_DOUBLE, 1, zeta )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'xmax',      &
                                   NF_DOUBLE, 1, xmax )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_double( ncid, NF_GLOBAL, 'ymax',      &
                                   NF_DOUBLE, 1, ymax )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!---------------------------------------------------
! Create dimensions
! nf_def_dim returns an integer (NF_NOERR if all OK)
! ncid is NetCDF ID, "name" of dimension
! nx is length of dimension
! Final argument (integer) is dimension ID
!---------------------------------------------------
  ncdfstatus = nf_def_dim( ncid, 'nx', 2*nx, nxdim )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_def_dim( ncid, 'ny', 2*ny, nydim )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_def_dim( ncid, 'nz', nz, nzdim )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_def_dim( ncid, 'one', 1, onedim )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  
!----------------------
! Define the dimensions
!----------------------
  dims(1)  = onedim  
  ncdfstatus = nf_def_var( ncid, 'time', NF_DOUBLE, 1, dims, varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_put_att_text( ncid, varid, 'units', 7,             &
                                'UNKNOWN' )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  dims(1)  = onedim  
  ncdfstatus = nf_def_var( ncid, 'iteration', NF_INT, 1, dims, varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  dims(1)  = onedim  
  ncdfstatus = nf_def_var( ncid, 'timestep', NF_DOUBLE, 1, dims, varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  if ( itype .eq. 1 ) then
    dims(1)  = nxdim  
    dims(2)  = nydim  
    dims(3)  = nzdim  
  else 
    dims(1)  = nzdim  
    dims(2)  = nydim  
    dims(3)  = nxdim  
  end if

!---------------------
! Define the variables
!---------------------
  if ( iprec .eq. 8 ) then
    ncdfstatus = nf_def_var( ncid, 'R', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'T', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'U', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'V', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'W', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'Bx', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'By', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'Bz', NF_DOUBLE, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  else if ( iprec .eq. 4 ) then
    ncdfstatus = nf_def_var( ncid, 'R', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'T', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'U', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'V', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'W', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'Bx', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'By', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_def_var( ncid, 'Bz', NF_REAL, 3, dims, varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  else 
    print*,'UNKNOWN type ', iprec ,' Type must be 4 or 8'
    ier = -1
    return
  end if

!-------------------------------------------------------
! Change mode so that data can be written to NetCDF file
!-------------------------------------------------------
  ncdfstatus = nf_enddef(ncid)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!------------------------------------------------------------------
! nf_sync synchronises disk writes with memory buffers (immediately 
! available after writing)
!------------------------------------------------------------------
  ncdfstatus = nf_sync(ncid)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
 
!-------------------------
! Sets return value of ier
!-------------------------
  ier = 0

  return

END SUBROUTINE CreateNetCDF

!========================================
! This subroutine writes to a NetCDF file
!========================================
SUBROUTINE WriteNetCDF( iter, ttime, dt1, ier )
  implicit none

!-----------------------------------
! Declare output and input variables
!-----------------------------------
  integer,       intent(out) :: ier
  integer,       intent(in)  :: iter
  real(kind=dp), intent(in)  :: ttime               
  real(kind=dp), intent(in)  :: dt1

!-------------------
! Local declarations
!-------------------
  integer            :: i, j, k
#ifdef MPI
  integer, parameter :: rtag  = 11
  integer, parameter :: ttag  = 12
  integer, parameter :: utag  = 13
  integer, parameter :: vtag  = 14
  integer, parameter :: wtag  = 15
  integer, parameter :: bxtag = 16
  integer, parameter :: bytag = 17
  integer, parameter :: bztag = 18  
#endif

!--------------------
! Dummy value for ier
!--------------------
  ier = -1
 
!---------------------------------------------
! Set default values of ireq, mstat and ncount
!---------------------------------------------
#ifdef MPI
  ireq = 0 
  mstat = 0
  ncount = 0
#endif

!----------------------------------------
! Save the simple variables on myrank = 0
!----------------------------------------
  if (myrank .eq. 0) then
    ncdfstatus = nf_inq_varid( ncid, 'time', varid )
    if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    ncdfstatus = nf_put_var_double( ncid, varid, ttime )
    if ( ncdfstatus .ne. NF_NOERR) call handle_error()

    ncdfstatus = nf_inq_varid( ncid, 'iteration', varid )
    if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    ncdfstatus = nf_put_var_int( ncid, varid, iter )
    if ( ncdfstatus .ne. NF_NOERR) call handle_error()

    ncdfstatus = nf_inq_varid( ncid, 'timestep', varid )
    if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    ncdfstatus = nf_put_var_double( ncid, varid, dt1 )
    if ( ncdfstatus .ne. NF_NOERR) call handle_error()
  end if

#ifdef MPI

!-----------------------
! Do parallel case first
! Note: work on myrank=0
!-----------------------
  if (myrank .eq. 0) then
    call MPI_Barrier(comm,ier)
    varname = 'R'
    call GatherWrite( rtag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'T'
    call GatherWrite( ttag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'U'
    call GatherWrite( utag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'V'
    call GatherWrite( vtag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'W'
    call GatherWrite( wtag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'Bx'
    call GatherWrite( bxtag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'By'
    call GatherWrite( bytag, ier )

    call MPI_Barrier(comm,ier)
    varname = 'Bz'
    call GatherWrite( bztag, ier )

    print*,'Finished writing to NetCDF file'
  else
    ncount = nlayer*4*nx*ny

    call MPI_Barrier(comm,ier)
    call MPI_Send( rnew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, rtag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( tnew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, ttag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( unew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, utag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( vnew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, vtag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( wnew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, wtag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( bxnew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, bxtag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( bynew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, bytag, comm, ier )

    call MPI_Barrier(comm,ier)
    call MPI_Send( bznew, ncount, MPI_DOUBLE_PRECISION,              &
                       0, bztag, comm, ier )    

  end if

#else

!----------------------
! Now do serial version
!----------------------
  print*,'Writing R'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'R', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, rnew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing T'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'T', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, tnew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing U'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'U', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                      icount, unew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing V'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'V', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, vnew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing W'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'W', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, wnew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing Bx'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'Bx', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, bxnew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing By'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'By', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, bynew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Writing Bz'
  do k=1,nlayer
    do j=1,2*ny
      istarting(1) = 1
      istarting(2) = j
      istarting(3) = k
      icount(1) = 2*nx
      icount(2) = 1
      icount(3) = 1
      ncdfstatus = nf_inq_varid( ncid, 'Bz', varid )
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
      ncdfstatus = nf_put_vara_double(ncid,varid,istarting,   &
                                        icount, bznew(1,j,k))
      if ( ncdfstatus .ne. NF_NOERR) call handle_error()
    end do
  end do

  print*,'Finished writing to NetCDF file'

#endif

!---------------------
! Return value for ier
!---------------------
  ier = 0 

  return

END SUBROUTINE WriteNetCDF

!==============================================
! This subroutine opens a netcdf file
! It also reads attributes and array dimensions
!==============================================
SUBROUTINE OpenNetCDF( iter,ttime,dt1,ier)
  implicit none

!-------------------------
! Declare output variables
!-------------------------
  integer,       intent(out) ::  iter               
  integer,       intent(out) ::  ier
  real(kind=dp), intent(out) ::  ttime               
  real(kind=dp), intent(out) ::  dt1

!----------------------------
! Set dummy values for output
!----------------------------
  ncid = -1
  ier = -1
  ttime = 0.0d0
  ier = 0.0d0

!-------------------------------------------------------------
! NF_OPEN reads named file, NF_WRITE implies file is writeable
! ncid is returned NetCDF file ID.
!-------------------------------------------------------------
  ncdfstatus =  NF_OPEN( restart_file, NF_WRITE, ncid )
  if ( ncdfstatus .ne. NF_NOERR ) then
    print*,'ERROR: Unable to open file ', restart_file
    call  handle_error()
    return
  end if 

!----------------------------------------------------------------------
! Now read the attributes
!----------------------------------------------------------------------
  ncdfstatus = nf_get_att_int( ncid,NF_GLOBAL,'version',version)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_text( ncid,NF_GLOBAL,'history',history)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'gamma', gamma)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'sigma', sigma)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL,'theta', theta)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'ck', ck)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'cm', cm )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'f', f )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'tayl', tayl )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'psi', psi )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL, 'zeta', zeta )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL,'xmax', xmax )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_att_double( ncid, NF_GLOBAL,'ymax', ymax )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!---------------------------------
! Read the time and store in ttime
!---------------------------------
  ncdfstatus = nf_inq_varid( ncid, 'time', varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_var_double( ncid, varid, ttime)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!--------------------------------------
! Read the iterations and store in iter
!--------------------------------------
  ncdfstatus = nf_inq_varid( ncid, 'iteration', varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_var_int( ncid, varid, iter)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!-----------------------------------
! Read the timestep and store in dt1
!-----------------------------------
  ncdfstatus = nf_inq_varid( ncid, 'timestep', varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_get_var_double( ncid, varid, dt1)
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!----------------------------------------------------------------
! Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid
!---------------------------------------------------------------- 
  ncdfstatus = nf_inq_dimid( ncid, 'nx', nxid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_inq_dimid( ncid, 'ny', nyid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_inq_dimid( ncid, 'nz', nzid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!------------------------------------------
! Read the actual dimensions 
! nf_inq_dimlen returns length of dimension
!------------------------------------------
  ncdfstatus = nf_inq_dimlen( ncid, nxid, nx )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_inq_dimlen( ncid, nyid, ny )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_inq_dimlen( ncid, nzid, nz )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!------------------------
! Set return value of ier
!------------------------
  ier = 0

  return

END SUBROUTINE OpenNetCDF

!============================================
! This subroutine reads in the data variables
!============================================
SUBROUTINE ReadNetCDF( ier )
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer,       intent(out) ::  ier

!-------------------
! Local declarations
!-------------------
  integer            :: i, j, k
  integer            :: lnx, lny, lnz
#ifdef MPI
  integer, parameter :: rtag  = 1
  integer, parameter :: ttag  = 2
  integer, parameter :: utag  = 3
  integer, parameter :: vtag  = 4
  integer, parameter :: wtag  = 5
  integer, parameter :: bxtag = 6
  integer, parameter :: bytag = 7
  integer, parameter :: bztag = 8
#endif  

!-----------------------
! Set dummy value of ier
!-----------------------
  ier = -1

!---------------------------------------------
! Set default values of ireq, mstat and ncount
!---------------------------------------------
#ifdef MPI
  ireq = 0 
  mstat = 0
  ncount = 0
#endif
!-----------------------------------------------------
! First get the dimensions of the arrays on myrank = 0
!-----------------------------------------------------
  if (myrank .eq. 0) then
    ncdfstatus = nf_inq_dimid( ncid, 'nx', nxid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_inq_dimid( ncid, 'ny', nyid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_inq_dimid( ncid, 'nz', nzid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

    ncdfstatus = nf_inq_dimlen( ncid, nxid, lnx )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_inq_dimlen( ncid, nyid, lny )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    ncdfstatus = nf_inq_dimlen( ncid, nzid, lnz )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

    if (( lnx .ne. nx ).or.( lny.ne.ny).or.(lnz.ne.nz)) then
      write(*,*)'ERROR: Dimensions disagree ', lnx, nx,lny,ny,lnz,nz
      ier = -1
      return
    end if
  end if

!-----------------------------
! Now read data and distribute
!-----------------------------
! First the MPI routines
!-----------------------------
#ifdef MPI
  call MPI_BARRIER(comm, ier )
  if (myrank .eq. 0) then

!------------------
! Check data format
!------------------
    call GetType()
    if ( itype .ne. 1 ) then
      write(*,*)'ERROR: Arrays in restart file are not dimensioned as '
      write(*,*)'(nx,ny,nz). --- Please convert to this form '
      ier = -1
      return
    end if

!---------------------------------------------
! Broadcast data from myrank = 0 to myrank > 0
!---------------------------------------------
    varname = 'R'
    call ReadBcast( rtag, ier )
    varname = 'T'
    call ReadBcast( ttag, ier )
    varname = 'U'
    call ReadBcast( utag, ier )
    varname = 'V'
    call ReadBcast( vtag, ier )
    varname = 'W'
    call ReadBcast( wtag, ier )
    varname = 'Bx'
    call ReadBcast( bxtag, ier )
    varname = 'By'
    call ReadBcast( bytag, ier )
    varname = 'Bz'
    call ReadBcast( bztag, ier )
    print *,'Finished Reading NetCDF file'    
  else

!----------------------------
! Receive data from ReadBcast
!----------------------------
    ncount = nx*ny*nlayer
    call MPI_Irecv( r(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, rtag, comm, ireq(1), ier )
    call MPI_Irecv( t(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, ttag, comm, ireq(2), ier )
    call MPI_Irecv( u(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, utag, comm, ireq(3), ier )
    call MPI_Irecv( v(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, vtag, comm, ireq(4), ier )
    call MPI_Irecv( w(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, wtag, comm, ireq(5), ier )
    call MPI_Irecv( bx(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, bxtag, comm, ireq(6), ier )
    call MPI_Irecv( by(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, bytag, comm, ireq(7), ier )
    call MPI_Irecv( bz(1,1,1), ncount, MPI_DOUBLE_PRECISION,   &
                         0, bztag, comm, ireq(8), ier )
    call MPI_Waitall( 8, ireq, mstat, ier )
  end if

!-------------------
! Now do serial case
!-------------------
#else

!----------------
! Check data type
!----------------
  call GetType()
  if ( itype .ne. 1 ) then
    write(*,*)'ERROR: Arrays in restart file are not dimensioned as '
    write(*,*)'(nx,ny,nz). --- Please convert to this form '
    ier = -1
    return
  end if
 
!-------------------------------
! Do double precision case first
!-------------------------------
  if (iprec .eq. 8) then

    write(*,*) 'Reading R'
    ncdfstatus = nf_inq_varid( ncid, 'R', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, r(1,j,k))
      end do
    end do

    write(*,*) 'Reading T'
    ncdfstatus = nf_inq_varid( ncid, 'T', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, t(1,j,k))
      end do
    end do

    write(*,*) 'Reading U'
    ncdfstatus = nf_inq_varid( ncid, 'U', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, u(1,j,k))
      end do
    end do

    write(*,*) 'Reading V'
    ncdfstatus = nf_inq_varid( ncid, 'V', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, v(1,j,k))
      end do
    end do

    write(*,*) 'Reading W'
    ncdfstatus = nf_inq_varid( ncid, 'W', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, w(1,j,k))
      end do
    end do

    write(*,*) 'Reading Bx'
    ncdfstatus = nf_inq_varid( ncid, 'Bx', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, bx(1,j,k))
      end do
    end do

    write(*,*) 'Reading By'
    ncdfstatus = nf_inq_varid( ncid, 'By', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, by(1,j,k))
      end do
    end do

    write(*,*) 'Reading Bz'
    ncdfstatus = nf_inq_varid( ncid, 'Bz', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_double(ncid, varid, istarting,    &
                                        icount, bz(1,j,k))
      end do
    end do

  else if (iprec .eq. 4) then

!----------------------------------------
! Read as real data and convert to double
!----------------------------------------
    write(*,*) 'Reading R'
    ncdfstatus = nf_inq_varid( ncid, 'R', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,     &
                                      icount, rwk)
        do i=1,nx
          r(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading T'
    ncdfstatus = nf_inq_varid( ncid, 'T', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          t(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading U'
    ncdfstatus = nf_inq_varid( ncid, 'U', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          u(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading V'
    ncdfstatus = nf_inq_varid( ncid, 'V', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          v(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading W'
    ncdfstatus = nf_inq_varid( ncid, 'W', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          w(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading Bx'
    ncdfstatus = nf_inq_varid( ncid, 'Bx', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          bx(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading By'
    ncdfstatus = nf_inq_varid( ncid, 'By', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          by(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

    write(*,*) 'Reading Bz'
    ncdfstatus = nf_inq_varid( ncid, 'Bz', varid )
    if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
    do k=1,nlayer
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
        ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                      icount, rwk)
        do i=1,nx
          bz(i,j,k) = dble(rwk(i))
        end do
      end do
    end do

  else
    ier = -2 
    write(*,*)'ERROR: WHAT TYPE of READ '
  end if

  print *,'Finished Reading NetCDF file'
#endif

!------------------------
! Set return value of ier
!------------------------
  ier = 0

  return

END SUBROUTINE ReadNetCDF

!=============================================================
! Closes off NetCDF file if on the master processor (myrank=0)
!=============================================================
SUBROUTINE CloseNetCDF(ier)
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer,       intent(out) ::  ier

!-----------------------
! Set dummy value of ier
!-----------------------
  ier = -1

!---------------------------
! Check to see if myrank = 0
!---------------------------
  if ( myrank .ne. 0 ) then
    ier = 0
    return  
  end if

!--------------------
! Synchronise buffers
!--------------------
  ncdfstatus = nf_sync(ncid)
  if ( ncdfstatus .ne. NF_NOERR ) call handle_error()

!-----------------------
! Closes off NetCDF file
!-----------------------
  ncdfstatus = nf_close(ncid)
  if ( ncdfstatus .ne. NF_NOERR ) call handle_error()

!------------------------
! Set return value of ier
!------------------------
  ier = 0

  return

END SUBROUTINE CloseNetCDF

!=========================================================
! Returns type = 1 if the arrays are dimension (nx,ny,nz),
!         type = -1 if they are dimensioned (nz,ny,nx) and
!         type = 0 otherwise.
! Returns iprec = 8 or iprec = 4 depending upon data types
!=========================================================
SUBROUTINE GetType
  implicit none

!-------------------
! Local declarations
!-------------------
  integer        :: dimids(3)
  integer        :: xtype

!---------------------------------------
! Set default values for iprec and itype
!---------------------------------------
  itype = 0
  iprec = 0

!-----------------------------------------
! Work out array ordering using "R"
! nf_inq_varid returns variable ID for "R"
!-----------------------------------------
  ncdfstatus = nf_inq_varid( ncid, 'R', varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!-------------------------------------------------------
! Get the variable ids for "R"
! nf_inq_vardimid returns dim Ids for varid in dimids(3)
!-------------------------------------------------------
  ncdfstatus = nf_inq_vardimid( ncid, varid, dimids )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!------------------------------------------
! Read the id's of the variables nx, ny, nz
!------------------------------------------
  ncdfstatus = nf_inq_dimid( ncid, 'nx', nxid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_inq_dimid( ncid, 'ny', nyid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()
  ncdfstatus = nf_inq_dimid( ncid, 'nz', nzid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!---------------------------------------------------
! Case 1 corresponds to nx,ny,nz; Case 2 to nz,ny,nx
!--------------------------------------------------- 
  if      (( dimids(1) .eq. nxid ).and.   &
           ( dimids(2) .eq. nyid ).and.   &
           ( dimids(3) .eq. nzid ) ) then
    itype = 1
  else if (( dimids(1) .eq. nzid ).and.   &
          ( dimids(2) .eq. nyid ).and.   &
          ( dimids(3) .eq. nxid ) ) then
    itype =-1
  end if

!-------------------------------------
! Find out the precision of "R"
! nf_inq_vartype returns type of varid
!-------------------------------------
  ncdfstatus = nf_inq_vartype( ncid, varid, xtype )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!--------------------------------
! Allocate correct value to xtype 
!--------------------------------
  if ( xtype .eq. NF_DOUBLE ) then
    iprec = 8
  else if ( xtype .eq. NF_REAL ) then
    iprec = 4
  end if

END SUBROUTINE GetType

!=================================================================
! This subroutine uses the global variable ncdfstatus to determine 
! whether or not an error has occurred. NF_NOERR is an integer 
! which implies success. NF_STRERROR returns a string indicating 
! that there has been an error
!=================================================================
SUBROUTINE handle_error 
  implicit none

  if ( ncdfstatus .ne. NF_NOERR ) then
    print*, NF_STRERROR(ncdfstatus)
  end if
END SUBROUTINE handle_error

!-------------------------------------------------------
! Now we have two routines that are only needed with MPI
!-------------------------------------------------------
#ifdef MPI

!-----------------------
! Read Broadcast routine
!-----------------------
SUBROUTINE ReadBcast( dummytag, ier )
  implicit none

!----------------------
! Argument declarations
!----------------------
  integer, intent(in)   :: dummytag
  integer, intent(out)  :: ier

!-------------------
! Local declarations
!-------------------
  integer :: i, j, k
  integer :: iproc 
  integer :: kp
  integer :: istt, ifin

  print*,'Reading ', varname 

!------------------------
! Set dummy value for ier
!------------------------
  ier = -1

!-------------------------------
! Find varid for current varname
!-------------------------------
  ncdfstatus = nf_inq_varid( ncid, varname, varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

!---------------------------------------------
! Work out what to broadcast to each processor
!---------------------------------------------
  do iproc=0,nproc-1
    istt = ((iproc)*nz)/nproc +1
    ifin = ((iproc+1)*nz)/nproc
    ncount = nx*ny*(ifin-istt+1)
    do k=istt,ifin
      kp = k - istt  + 1
      do j=1,ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k
        icount(1) = nx
        icount(2) = 1
        icount(3) = 1
 
        if ( iprec .eq. 8 ) then
          ncdfstatus = nf_get_vara_double(ncid, varid, istarting,  &
                                          icount, wk1(1,j,kp) )
        else 
          ncdfstatus = nf_get_vara_real(ncid, varid, istarting,    &
                                          icount, rwk )
          do i=1,nx
            wk1(i,j,kp) = dble(rwk(i))
          end do
        end if
      end do
    end do

!-----------------------------------------------------
! Carry out broadcast (or read directly if myrank = 0)
!-----------------------------------------------------
    if ( iproc .eq. 0 ) then
      do k=istt,ifin
        kp = k - istt  + 1
        do j=1,ny
          do i=1,nx
            if (varname .eq. 'R') r(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'T') t(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'U') u(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'V') v(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'W') w(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'Bx') bx(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'By') by(i,j,kp) = wk1(i,j,kp)
            if (varname .eq. 'Bz') bz(i,j,kp) = wk1(i,j,kp)
          end do
        end do
      end do
    else
      call MPI_Send( wk1, ncount, MPI_DOUBLE_PRECISION, &
                        iproc, dummytag, comm, ier )
    end if
  end do
    
!-------------------------
! Set return value for ier
!------------------------- 
  ier = 0 
 
  return

END SUBROUTINE ReadBcast

!=======================================
! Gatherwrite routine for NetCDF writing
!=======================================
SUBROUTINE GatherWrite( dummytag, ier )
  implicit none

!----------------------
! Argument declarations
!----------------------
  integer, intent(in)   :: dummytag
  integer, intent(out)  :: ier

!-------------------
! Local declarations
!-------------------
  integer :: i, j, k
  integer :: iproc 
  integer :: istt, ifin, nlayer_size

  print*,'Writing ', varname 

!------------------------
! Set dummy value for ier
!------------------------
  ier = -1

!-----------------------------
! Gets Variable ID for varname
!-----------------------------
  ncdfstatus = nf_inq_varid( ncid, varname, varid )
  if ( ncdfstatus .ne. NF_NOERR ) call  handle_error()

  do iproc=0,nproc-1

!---------------------------
! Work out data distribution
!---------------------------
    istt = ((iproc)*nz)/nproc +1
    ifin = ((iproc+1)*nz)/nproc
    nlayer_size = ifin - istt + 1
    ncount = 4*nx*ny*nlayer_size

!-----------------------------------------------------
! If on leading processor, copy variable to work array
!-----------------------------------------------------
    if (iproc .eq. 0) then
      do k=1, nlayer_size
        do j=1,2*ny
          do i=1,2*nx
            if (varname .eq. 'R') wk1new(i,j,k) = rnew(i,j,k)
            if (varname .eq. 'T') wk1new(i,j,k) = tnew(i,j,k)
            if (varname .eq. 'U') wk1new(i,j,k) = unew(i,j,k)
            if (varname .eq. 'V') wk1new(i,j,k) = vnew(i,j,k)
            if (varname .eq. 'W') wk1new(i,j,k) = wnew(i,j,k)
            if (varname .eq. 'Bx') wk1new(i,j,k) = bxnew(i,j,k)
            if (varname .eq. 'By') wk1new(i,j,k) = bynew(i,j,k)
            if (varname .eq. 'Bz') wk1new(i,j,k) = bznew(i,j,k)
          end do
        end do
      end do
    else

!------------------------------------------------------
! Otherwise carry out a blocking receive to fill wk1new      
!------------------------------------------------------
      call MPI_Recv( wk1new, ncount, MPI_DOUBLE_PRECISION, &
                        iproc, dummytag, comm, status, ier )
    end if

!------------------
! Now write to file
!------------------
    do k=1,nlayer_size
      do j=1,2*ny
        istarting(1) = 1
        istarting(2) = j
        istarting(3) = k+istt-1
        icount(1) = 2*nx
        icount(2) = 1
        icount(3) = 1

        ncdfstatus = nf_put_vara_double(ncid, varid, istarting,  &
                                          icount, wk1new(1,j,k) )
	if ( ncdfstatus .ne. NF_NOERR ) call handle_error()
      end do
    end do
  end do

!-------------------------
! Set return value for ier
!------------------------- 
  ier = 0 
 
  return

END SUBROUTINE GatherWrite
#endif

END MODULE netcdf
