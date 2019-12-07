! $Id: depth_check.F,v 1.2.2.4 2008/04/03 16:23:28 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.2.2.4 $
!==============================================================================|
!     ENSURE DEPTH IS GREATER THAN MIN_DEPTH                                   |
!     IF THIS CONDITION IS VIOLATED, HALT THE PROGRAM AND WRITE A WARNING      |
!     MESSAGE TO THE SCREEN                                                    |
!==============================================================================|

SUBROUTINE DEPTH_CHECK

  !==============================================================================|
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_PAR  
  IMPLICIT NONE
  INTEGER, DIMENSION(NPROCS) :: SBUF,RBUF
  INTEGER  :: I,II,MLOC,IERR, RECV
  REAL(SP) :: DMIN
  !==============================================================================|

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: depth_check"
!--Calculate Minimum Depth and Set Global Node Number if Min Depth < MIN_DEPTH
   SBUF = 0
   MLOC = 0
   IERR = 0
   DMIN = MINVAL(D(1:M))
   ! NGID NOW EXISTS FOR BOTH SERIAL AND PARALLEL
   IF(DMIN < MIN_DEPTH) MLOC = NGID(MINLOC(D(1:M),DIM=1))

!--Reduce in Master Processor Array and Dump To Screen
   SBUF(MYID) = MLOC
   RBUF = SBUF
#  if defined (MULTIPROCESSOR)
   RECV=MSRID-1
   IF(PAR)CALL MPI_REDUCE(SBUF,RBUF,NPROCS,MPI_INTEGER,MPI_MAX,RECV,MPI_FVCOM_GROUP,IERR)
   IF(IERR/=0) CALL FATAL_ERROR("COULD NOT MPI_REDUCE IN DEPTH CHECK")
#  endif


!--If Depth Condition is Violated Write Warning and Halt User
   IF(MSR)THEN
     IF(SUM(RBUF) /= 0) THEN
        DO I=1,NPROCS
           II = RBUF(I)
           IF(II /= 0)THEN
              !         WRITE(*,*)'DEPTH IN NODE: ',II,' AT ',XG(II)+VXMIN,YG(II)+VYMIN, &
              !                        ' IS LESS THAN MIN_DEPTH'
              
              WRITE(*,*)'DEPTH IN NODE: ',II,'; IS LESS THAN MIN_DEPTH'
              WRITE(*,*)'ADJUST BATHYMETRY AT THIS (THESE) LOCATION(S) OR'
              WRITE(*,*)'RECOMPILE FVCOM WITH FLOODING/DRYING FORMULATION'                 
              WRITE(*,*)'STOPPING....'
              CALL PSTOP
           END IF
        END DO
     END IF
  END IF

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: depth_check"

   END SUBROUTINE DEPTH_CHECK
!==============================================================================|
