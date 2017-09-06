!================================
! Module to handle communications 
!================================
MODULE comms
  use types

  implicit none

!-----------------------------------
! Put in MPI Include file
! 
! NOTE: Compile path needs to be set
!-----------------------------------
#ifdef MPI
  INCLUDE 'mpif.h'
#endif

!--------------------
! Global declarations
!--------------------
  integer  :: myrank
  integer  :: nproc
  integer  :: comm

#ifdef MPI
!-----------------------------
! MPI specific status variable
!-----------------------------
  integer              :: status(mpi_status_size)
  integer              :: ncount
  integer              :: msg
  integer, allocatable :: ireq(:)
  integer, allocatable :: mstat(:,:)

!--------------------------------
! Global send and receive buffers
!--------------------------------
  integer       :: isendbuf(11)
  real(kind=dp) :: sendbuf(20)
  real(kind=dp) :: recvbuf(5)
#endif

CONTAINS

!===========================================
! This subroutine initialises communications
!===========================================
SUBROUTINE init_comms(iret)
  implicit none

!------------------------
! Declare output variable
!------------------------
  integer, intent(out) :: iret

!------------------
! Local declaration
!------------------
#ifdef MPI
  integer  :: ier(6) 
#endif

!----------------------
! Default value of iret
!----------------------
  iret = 0

!-----------------------------------------------------
! MPI initialisation
!
! NOTE: dummy values set for variables if MPI not used
!-----------------------------------------------------
#ifdef MPI
  call mpi_init(ier(1))
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ier(2))
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ier(3))
  call mpi_comm_dup( MPI_COMM_WORLD, comm, ier(4))
#else
  nproc = 1
  myrank = 0
  comm = -1
#endif

#ifdef MPI
!------------------------------------
! Allocate and define wait parameters
!------------------------------------
  allocate(mstat(MPI_STATUS_SIZE,40), stat=ier(5))
  allocate(ireq(40),                  stat=ier(6))

!-------------
! Check result
!-------------
  if ( (ier(5) .gt. 0) .or.    &
       (ier(6) .gt. 0) )   then
    iret = -1
    return
  end if

!---------------------
! Initialise variables
!---------------------
  msg = 0 
  ncount = 0
  mstat = 0
  ireq = 0
  status = 0
  isendbuf = 0
  sendbuf = 0.0d0
  recvbuf = 0.0d0
#endif

  return

END SUBROUTINE init_comms

!=====================
! Clean exit from code
!=====================
SUBROUTINE exit_comms()
  implicit none

!-----------------------
! Stop all MPI processes
!-----------------------
#ifdef MPI
  integer  ::   ier
  call mpi_finalize(ier)
#endif

  return

END SUBROUTINE exit_comms

!=====================
! Panic exit from code
!=====================
SUBROUTINE abort_comms(error)
  implicit none

!------------------------
! Declare input variable
!------------------------
  integer, intent(in) ::   error

!--------------------
! Abort communication
!--------------------
#ifdef MPI
  integer  :: ier
  call mpi_abort(comm, error, ier )
#else
  print*,'ERROR = ',error
  call exit(1)
#endif

  return

END SUBROUTINE abort_comms

END MODULE comms





