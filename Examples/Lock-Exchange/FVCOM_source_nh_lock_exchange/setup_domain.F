! $Id: setup_domain.F,v 1.1.2.10 2008/06/05 14:57:29 jqi Exp $
! $Name: New_Input $
! $Revision: 1.1.2.10 $
!==============================================================
SUBROUTINE SETUP_DOMAIN
!==============================================================
! IF(PAR): 
!          Uses METIS to determine local domains
!          Create global to local and local to global maps
!==============================================================
  USE MOD_PAR
  USE ALL_VARS
  USE MOD_NCDIO
  USE MOD_SETUP
  USE MOD_NESTING
  
  !==============================================================================|
  !  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
  !  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
  !  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
  !  TRANSFER GLOBAL COORDINATES TO LOCAL DOMAINS
  !==============================================================================|

  IMPLICIT NONE
  INTEGER :: SENDER,STATUS
  INTEGER :: nvals, IERR, I, EGL
  

     if (dbg_set(dbg_log)) &
          & write(IPT,*) "SETTING UP FVCOM PROCESSOR ENVIRONMENT"
  
  ! SETUP OUTPUT MODE OR TURN IT OFF
  CALL SETUP_MPI_IO_MODE(USE_MPI_IO_MODE,MPI_FVCOM_GROUP)
  
  IF(SERIAL) THEN
     
     ! DO NOT ALLOCATE EL_PID FOR NON-MULTIPROCESSOR CASES
          
     ! SHIFT GRID DIMENSIONS TO LOCAL FOR ANY SERIAL CASE
     MT=MGL
     NT=NGL
     
     M = MGL
     N = NGL
     
     ! TRANSFER CONNECTIVTY TO LOCAL DOMAIN NAME VARIABLES
     ALLOCATE(NV(0:NT,4))
     ALLOCATE(NVGL(0:NT,3))
     NV   = NVG
     ! DUPLICATE FOR NC OUTPUT - (IF PAR HAS GLOBAL NODE NUMBERING)
     NVGL = NVG(:,1:3)

     NV(0,:) = 0
     NVGL(0,:) = 0

     DEALLOCATE(NVG)
     

     ! MUST ALLOCATE EL_PID FOR OUTPUT OF DATA
     ALLOCATE(EL_PID(NGL))  ;  EL_PID = 1


     ! SETUP NGID
     ALLOCATE(NGID(0:MGL))  ; 
     DO i=0,MGL
        NGID(i)=i
     END DO

     NGID_X => NGID

     ! SETUP NLID
     ALLOCATE(NLID(0:MGL))
     DO i=0,MGL
        NLID(i) = i
     END DO

     NLID_X => NLID

     ! SETUP EGID
     ALLOCATE(EGID(0:NGL))  ; 
     DO i=0,NGL
        EGID(i)=i
     END DO

     EGID_X => EGID

     ! SETUP ELID
     ALLOCATE(ELID(0:NGL))
     DO i=0,NGL
        ELID(i) = i
     END DO

     ELID_X => ELID

     ALLOCATE(NDE_ID(MGL)); NDE_ID=0

  ELSE ! IF PAR

#  if defined (MULTIPROCESSOR) 
     
     !
     !  DECOMPOSE DOMAIN BY ELEMENTS USING METIS
     !
     ALLOCATE(EL_PID(NGL))  ;  EL_PID = 1
     CALL DOMDEC(NVG,NPROCS,EL_PID,MSR)
     
     ! GENERATE MAP TO EXCHANGE DATA BETWEEN GLOBAL AND LOCAL DOMAINS
     CALL GENMAP

     IF(.NOT. IOPROC) THEN
        ! TRANSFER CONNECTIVTY TO LOCAL DOMAIN NAME VARIABLES
        ALLOCATE(NV(0:NT,4))
        ! DUPLICATE FOR NC OUTPUT - HAS GLOBAL NODE NUMBERING
        ALLOCATE(NVGL(0:NT,3))
        
        DO I=1,N
           EGL  = EGID(I)
           NVGL(I,1:3) = NVG(EGL,1:3)! NVGL is used to create the data files
           !(MUST HAVE GLOBAL NODE NUMBERING)
           
           NV(I,1:4) = NLID(NVG(EGL,1:4))! NV is the local connectivity
        END DO
        
        DO I=1,NHE
           EGL  = HE_LST(I)
           NVGL(I+N,1:3) = NVG(EGL,1:3)! NVGL is used to create the data files
           !(MUST HAVE GLOBAL NODE NUMBERING)
           
           NV(I+N,1:4) = NLID_X(NVG(EGL,1:4))! NV is the local connectivity
        END DO

        NV(0,:) = 0
        NVGL(0,:) = 0

     END IF

# endif
  END IF
  
  ! PASS OBC TO LOCAL VARIABLES AND GENERATE MAP FROM GLOBAL

  IF(OBC_ON) CALL GENMAP_OBC
  IF(OBC_LONGSHORE_FLOW_ON) CALL GENMAP_LSF
  IF(NESTING_ON) CALL GENMAP_NEST
  
  IF (USE_MPI_IO_MODE) then
# if defined (MULTIPROCESSOR)
     CALL MPI_IO_LOOP
# else
     CALL FATAL_ERROR("'USE_MPI_IO_MODE' SHOULD NEVER BE TRUE WHEN C&
          &OMPILED WITHOUT MULTIPROCESSOR")
# endif
  END IF
  
  if(DBG_SET(dbg_log))  write(IPT,*) "Finished SETUP_LOCAL_DOMAINS"
  

END SUBROUTINE SETUP_DOMAIN
