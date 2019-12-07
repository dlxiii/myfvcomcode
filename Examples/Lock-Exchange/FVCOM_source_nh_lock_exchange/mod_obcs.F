! $Id: mod_obcs.F,v 1.7.2.11 2008/05/09 19:49:22 jqi Exp $
! $Name: New_Input $
! $Revision: 1.7.2.11 $
!==============================================================================!
! TYPE_OBC = 1: Surface Elevation Specified (Tidal Forcing) (ASL)              !
! TYPE_OBC = 2: AS TYPE_OBC=1 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 3: Zero Surface Elevation Boundary Condition (ASL_CLP)            !
! TYPE_OBC = 4: AS TYPE_OBC=3 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 5: GRAVITY-WAVE RADIATION IMPLICIT OPEN BOUNDARY CONDITION (GWI)  !
! TYPE_OBC = 6: AS TYPE_OBC=5 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 7: BLUMBERG AND KHANTA IMPLICIT OPEN BOUNDARY CONDITION (BKI)     !
! TYPE_OBC = 8: AS TYPE_OBC=7 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
! TYPE_OBC = 9: ORLANSKI RADIATION EXPLICIT OPEN BOUNDARY CONDITION (ORE)      !
! TYPE_OBC =10: AS TYPE_OBC=9 AND NON-LINEAR FLUX FOR CURRENT AT OPEN BOUNDARY !
!                                                                              !
! TYPE_TSOBC = 1: THE PERTURBATION OF TEMPERATURE AND SALINITY AT OPEN BOUNDARY!
!                 ARE EQUAL TO THAT AT NEXT_OBC                                !
! TYPE_TSOBC = 2: GRAVITY-WAVE RADIATION IMPLICIT OPEN BOUNDARY CONDITION FOR  !
!                 THE PERTURBATION OF TEMPERATURE AND SALINITY                 !
! TYPE_TSOBC = 3: BLUMBERG AND KHANTA IMPLICIT OPEN BOUNDARY CONDITION FOR     !
!                 THE PERTURBATION OF TEMPERATURE AND SALINITY                 !
! TYPE_TSOBC = 4: ORLANSKI RADIATION EXPLICIT OPEN BOUNDARY CONDITION FOR      !
!                 THE PERTURBATION OF TEMPERATURE AND SALINITY                 !
!==============================================================================!

MODULE MOD_OBCS
  USE MOD_PREC
  USE MOD_UTILS
  USE CONTROL,  type_tsobc => obc_ts_type
  USE MOD_TIME
  USE MOD_FORCE
  USE BCS
  IMPLICIT NONE
  SAVE

  !--Open Boundary Types, Lists, Pointers
# if defined (TIDE_OUTPUT)
  INTEGER               :: IOBCNODE_GL      !! = ntidenode_GL
  INTEGER,  ALLOCATABLE :: I_OBCNODE_GL(:)  !! = I_TIDENODE_GL
  INTEGER               :: IOBCELL_GL       !! = ntidecell_GL
  INTEGER,  ALLOCATABLE :: I_OBCELL_GL(:)   !! = I_TIDECELL_GL
# endif

  !   INTEGER               :: IOBCN_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY NODES
  !   INTEGER               :: IOBCN            !!LOCAL NUMBER OF OPEN BOUNDARY NODES
  !   INTEGER,  ALLOCATABLE :: I_OBC_GL(:)      !!GLOBAL ID OF OPEN BOUNDARY NODES
  !   INTEGER,  ALLOCATABLE :: I_OBC_N(:)       !!OPEN BOUNDARY NODE LIST
  INTEGER,  ALLOCATABLE :: NEXT_OBC(:)      !!INTERIOR NEIGHBOR OF OPEN BOUNDARY NODE
  INTEGER,  ALLOCATABLE :: NEXT_OBC2(:)     !!INTERIOR NEIGHBOR OF NEXT_OBC
  !   INTEGER,  ALLOCATABLE :: TYPE_OBC(:)      !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
  !   INTEGER,  ALLOCATABLE :: TYPE_OBC_GL(:)   !!OUTER BOUNDARY NODE TYPE (FOR SURFACE ELEVATION)
  INTEGER               :: IBCN(5)          !!NUMBER OF EACH TYPE OF OBN IN LOCAL  DOM
  INTEGER               :: IBCN_GL(5)       !!NUMBER OF EACH TYPE OF OBN IN GLOBAL DOM
  INTEGER,  ALLOCATABLE :: OBC_LST(:,:)     !!MAPPING OF OPEN BOUNDARY ARRAYS TO EACH TYPE
  INTEGER,  ALLOCATABLE :: NADJN_OBC(:)     !!NUMBER OF ADJACENT OPEN BOUNDARY NODES TO OBN
  INTEGER,  ALLOCATABLE :: ADJN_OBC(:,:)    !!ADJACENT OBN's of OBN
  INTEGER,  ALLOCATABLE :: NADJC_OBC(:)     !!NUMBER OF ADJACENT OPEN BOUNDARY CELLS TO OBN
  INTEGER,  ALLOCATABLE :: ADJC_OBC(:,:)    !!ADJACENT OPEN BOUNDARY CELLS

  !-- Open Boundary TEMPERATURE AND SALINITY FOR NUDGING
  REAL(SP), ALLOCATABLE :: TEMP_OBC(:,:)      !!OPEN BOUNDARY TEMPERATURE (see mod_force.F)     
  REAL(SP), ALLOCATABLE :: SALT_OBC(:,:)      !!OPEN BOUNDARY SALT (see mod_force.F)   

  !--Open Boundary Metrics
  INTEGER,  ALLOCATABLE :: NFLUXF_OBC(:)    !!NUMBER OF FLUX SEGMENTS TO OBN
  REAL(SP), ALLOCATABLE :: FLUXF_OBC(:,:)   !!FLUX FRACTION ON EACH SIDE OF OBN
  REAL(SP), ALLOCATABLE :: NXOBC(:)         !!INWARD POINTING X-NORMAL OF OBN
  REAL(SP), ALLOCATABLE :: NYOBC(:)         !!INWARD POINTING Y-NORMAL OF OBN
  REAL(SP), ALLOCATABLE :: DLTN_OBC(:)      !!DISTANCE BETWEEN NEXT_OBC AND OBN NORMAL TO BOUNDARY

  !--Previous Time Level Free Surface Fields
  REAL(SP), ALLOCATABLE :: ELM1(:)          !!SURFACE ELEV FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: ELM2(:)          !!SURFACE ELEV FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: T1M1(:,:)          !!TEMPERATURE FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: T1M2(:,:)          !!TEMPERATURE FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: S1M1(:,:)          !!SALINITY FROM PREVIOUS TIME LEVEL (ORLANSKI COND)
  REAL(SP), ALLOCATABLE :: S1M2(:,:)          !!SALINITY FROM PREVIOUS TIME LEVEL (ORLANSKI COND)

  !--Nonlinear Velocity Open Boundary Condition Arrays
  REAL(SP), ALLOCATABLE :: FLUXOBN(:)
  REAL(SP), ALLOCATABLE :: IUCP(:)
  REAL(SP), ALLOCATABLE :: XFLUX_OBCN(:)
  REAL(SP), ALLOCATABLE :: UARD_OBCN(:)
  REAL(SP), ALLOCATABLE :: XFLUX_OBC(:,:)


CONTAINS
  !==========================================================================|
  SUBROUTINE ALLOC_OBC_DATA

    !--------------------------------------------------------------------------|
    !  ALLOCATE AND INITIALIZE SURFACE ELEVATION ARRAYS FOR                    |
    !  TIME STEPS N-1 AND N-2                                                  |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    ALLOCATE(ELM1(0:MT))             ;ELM1       = ZERO 
    ALLOCATE(ELM2(0:MT))             ;ELM2       = ZERO 
    ALLOCATE(NEXT_OBC(IOBCN))        ;NEXT_OBC   = 0
    ALLOCATE(NEXT_OBC2(IOBCN))       ;NEXT_OBC2  = 0
    ALLOCATE(FLUXOBN(1:NT))          ;FLUXOBN    = ZERO
    ALLOCATE(IUCP(0:NT))             ;IUCP       = 1
    ALLOCATE(XFLUX_OBCN(IOBCN+1))    ;XFLUX_OBCN = ZERO
    ALLOCATE(UARD_OBCN(IOBCN+1))     ;UARD_OBCN  = ZERO
    ALLOCATE(XFLUX_OBC(IOBCN+1,KBM1));XFLUX_OBC  = ZERO

    ALLOCATE(TEMP_OBC(IOBCN,KBM1))   ;TEMP_OBC   = ZERO
    ALLOCATE(SALT_OBC(IOBCN,KBM1))   ;SALT_OBC   = ZERO

    IF(TYPE_TSOBC == 4)THEN
       ALLOCATE(T1M1(0:MT,KBM1))        ;T1M1       = ZERO 
       ALLOCATE(T1M2(0:MT,KBM1))        ;T1M2       = ZERO 
       ALLOCATE(S1M1(0:MT,KBM1))        ;S1M1       = ZERO 
       ALLOCATE(S1M2(0:MT,KBM1))        ;S1M2       = ZERO 
    END IF
    RETURN
  END SUBROUTINE ALLOC_OBC_DATA
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE ASSIGN_ELM1_TO_ELM2

    !--------------------------------------------------------------------------|
    !  Assign ELM1 to ELM2 and EL to ELM1                                      |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    ELM2 = ELM1
    ELM1 = EL

    IF(TYPE_TSOBC == 4)THEN
       T1M2 = T1M1
       T1M1 = T1

       S1M2 = S1M1
       S1M1 = S1
    END IF

    RETURN
  END SUBROUTINE ASSIGN_ELM1_TO_ELM2
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE ELEVATION_EQUI ! ONLY USED if defined (EQUI_TIDE)
    !--------------------------------------------------------------------------|
    !  Surface Elevation of EQUILIBRIUM TIDE                                   |
    !--------------------------------------------------------------------------|
    USE ALL_VARS
    IMPLICIT NONE

    TYPE(TIME):: TIME_ELAPSED
    INTEGER :: I,J
!    REAL(SP), ALLOCATABLE :: FORCE(:)
    REAL(DP):: TIME_FLT,PHAI_IJ


# if !defined(SPHERICAL)
    IF(.not.USE_PROJ) CALL FATAL_ERROR &
         & ("ELEVATION_EQI: ILLEGAL ATTEMPT TO USE EQUILIBRIUM TIDE IN CARTESIAN MODE",&
         & "YOU MUST BE USING PROJ4 TO DO THIS",&
         & "THE POLICE HAVE BEEN NOTIFIED AND WILL BE THERE SHORTLY")
# endif

    ! THIS CODE IS CALLED DURING EACH RK STAGE
    ! DECIDE WHICH TIME YOU WANT TO USE!
    TIME_ELAPSED = RKTime - SpecTime

!! IN MY JUDGEMENT, A TIMESERIES FORCING DOES NOT MAKE SENSE FOR ATMO/EQUITIDE


!    SELECT CASE(TIDE_FORCING_TYPE)
!    CASE(TIDE_FORCING_TIMESERIES)
       !
       !-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
       !
!       CALL FATAL_ERROR("ELEVATION_EQUI: Not set up for Julian Timeseries forcing yet?")

!    CASE(TIDE_FORCING_SPECTRAL)
       !
       !-Non-Julian: Set Elevation of Equilibrium Tide -----------------------------|
       !

       ELF_EQI=0.0_SP

       DO J = 1,nTideComps

          TIME_FLT = SECONDS(TIME_ELAPSED * PI2/PERIOD(J))             
          SELECT CASE(TIDE_TYPE(J))
          CASE(SEMIDIURNAL)

             DO I = 1, MT
                PHAI_IJ = LON(I)*PI2/360.0_SP
                ELF_EQI(I) = ELF_EQI(I) + BETA_EQI(J)*APT_EQI(J)*COS(LAT(I)*PI2/360.0_SP)**2         &
                     *DCOS(TIME_FLT+2.0_DP*PHAI_IJ)
             END DO
          CASE(DIURNAL)
             
             DO I = 1, MT   
                PHAI_IJ = LON(I)*PI2/360.0_SP
                ELF_EQI(I) = ELF_EQI(I) + BETA_EQI(J)*APT_EQI(J)*SIN(2.0_SP*LAT(I)*PI2/360.0_SP)     &
                     *DCOS(TIME_FLT+PHAI_IJ)
             END DO
             
          CASE DEFAULT
             CALL FATAL_ERROR("ELEVATION_EQUI: INVALID TIDE_TYPE:",&
                  &"Choices:"//TRIM(DIURNAL)//" or "//TRIM(SEMIDIURNAL) )
          END SELECT
       END DO
       
       ELF_EQI = ELF_EQI * RAMP
       
!    CASE DEFAULT
!       CALL FATAL_ERROR("ELEVATION_EQUI: INVALID TIDE FORCING TYPE ?")
!    END SELECT

    RETURN
  END SUBROUTINE ELEVATION_EQUI
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE ELEVATION_ATMO ! ONLY USED if defined (ATMO_TIDE)
    !--------------------------------------------------------------------------|
    !  Surface Elevation of ATMOSPHERIC TIDE                                   |
    !--------------------------------------------------------------------------|
    USE ALL_VARS
    IMPLICIT NONE

    REAL(SP),PARAMETER :: APT_ATMO = 0.0113_SP       
    REAL(DP),PARAMETER :: ALFA_ATMO = 112.0_SP*PI2/360.0_SP
    REAL(SP),PARAMETER :: S2PERIOD = 43200.0_SP

    INTEGER :: I,J
    TYPE(TIME):: TIME_ELAPSED
    REAL(DP):: TIME_FLT,PHAI_IJ
    REAL(SP):: FORCE

# if !defined(SPHERICAL)
    IF(.not.USE_PROJ) CALL FATAL_ERROR &
         & ("ELEVATION_ATMO: ILLEGAL ATTEMPT TO USE ATMOSPHERIC TIDE IN CARTESIAN MODE",&
         & "YOU MUST BE USING PROJ4 TO DO THIS",&
         & "THE POLICE HAVE BEEN NOTIFIED AND WILL BE THERE SHORTLY")
# endif

    ! THIS CODE IS CALLED DURING EACH RK STAGE
    ! DECIDE WHICH TIME YOU WANT TO USE!
    TIME_ELAPSED = RKTime - SpecTime

!! IN MY JUDGEMENT, A TIMESERIES FORCING DOES NOT MAKE SENSE FOR ATMO/EQUITIDE

!    SELECT CASE(TIDE_FORCING_TYPE)
!    CASE(TIDE_FORCING_TIMESERIES)
       !
       !-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
       !
!       CALL FATAL_ERROR("ELEVATION_EQUI: Not set up for Julian Timeseries forcing yet?")

!    CASE(TIDE_FORCING_SPECTRAL)
       !
       !-Non-Julian: Set Elevation of Equilibrium Tide -----------------------------|
       !

       TIME_FLT = SECONDS(TIME_ELAPSED * PI2/S2PERIOD)

       DO I = 1, MT
          PHAI_IJ = LON(I)*PI2/360.0_SP
          FORCE = APT_ATMO*DCOS(TIME_FLT+2.0_DP*PHAI_IJ-ALFA_ATMO)
          ELF_ATMO(I) = FORCE * RAMP
       END DO

!    CASE DEFAULT
!       CALL FATAL_ERROR("ELEVATION_ATMO: INVALID TIDE FORCING TYPE ?")
!    END SELECT



    RETURN
  END SUBROUTINE ELEVATION_ATMO
  !==========================================================================|

!==============================================================================|
!     SET BOUNDARY CONDITIONS FOR CALCULATING atmospher pressure external model|
!                                                                              |
!==============================================================================|

   SUBROUTINE BCOND_PA_AIR
#  if defined (AIR_PRESSURE)
!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_PAR

#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE

   REAL(SP) :: FACT,UFACT,TX,TY
!    for ice model
   REAL(SP) :: pa_air2
   INTEGER  I,J,L1,L2,IERR

!==============================================================================|
!   SURFACE AIR PRESSURE FORCING FOR INTERNAL MODE                                          !
!==============================================================================|
   ELF_AIR = -(PA_AIR-SLP0)/(1025._SP*GRAV_N)

# if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,ELF_AIR)
# endif
 
   RETURN
#  endif
   END SUBROUTINE BCOND_PA_AIR
   
   
!==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_ASL

    !--------------------------------------------------------------------------|
    !  Surface Elevation Boundary Condition (Tidal Forcing)                    |
    !--------------------------------------------------------------------------|

    USE ALL_VARS, only : ELF
    IMPLICIT NONE

    INTEGER :: I,II,J,JN
    TYPE(TIME):: TIME_ELAPSED
    REAL(DP):: TIME_FLT,PHAI_IJ
    REAL(SP):: FORCE
    REAL(SP), ALLOCATABLE :: BND_ELV(:)


    ALLOCATE(BND_ELV(IOBCN))
    !
    !-Julian: Set Elevation Based on Linear Interpolation Between Two Data Times-|
    !
    SELECT CASE(TIDE_FORCING_TYPE)
    CASE(TIDE_FORCING_TIMESERIES)

       ! THIS CODE IS CALLED DURING EACH RK STAGE
       ! DECIDE WHICH TIME YOU WANT TO USE!
       CALL UPDATE_TIDE(RKTime,BND_ELV)
       !      CALL UPDATE_TIDE(ExtTime,BND_ELV)
       DO J = 1, IBCN(1)
          JN = OBC_LST(1,J)
          II = I_OBC_N(JN)

          ELF(II) = BND_ELV(J)*RAMP
       END DO

       !
       !-Non-Julian: Set Elevation Based on Input Amplitude and Phase of Tidal Comps-|
       !

    CASE(TIDE_FORCING_SPECTRAL)

       ! THIS CODE IS CALLED DURING EACH RK STAGE
       ! DECIDE WHICH TIME YOU WANT TO USE!
       TIME_ELAPSED = RKTime - SpecTime

       !      write(ipt,*) "-----------------------------------------"
       DO I = 1, IBCN(1)
          JN = OBC_LST(1,I)
          II = I_OBC_N(JN)
          FORCE = 0.0_SP
          DO J = 1,nTideComps
             PHAI_IJ = PHAI(I,J)*PI2/360.0_SP
             TIME_FLT = SECONDS(TIME_ELAPSED * PI2/PERIOD(J))
             ! Take cosine of a double precision number to preserve accuracy!
             FORCE = APT(I,J)*DCOS(TIME_FLT -PHAI_IJ) + FORCE
          END DO
          FORCE = FORCE + EMEAN(I)
          ELF(II) = FORCE * RAMP
          !         write(ipt,*) "ELF=",ELF(II)

       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("BCOND_ASL: INVALID TIDE FORCING TYPE ?")
    END SELECT

    DEALLOCATE(BND_ELV)

    RETURN
  END SUBROUTINE BCOND_ASL
!!$!==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_ASL_CLP

    !--------------------------------------------------------------------------|
    !  Zero Surface Elevation Boundary Condition                               |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I,II,JN

    DO I = 1, IBCN(2)
       JN = OBC_LST(2,I)
       II = I_OBC_N(JN)  
       ELF(II) = 0.0_SP
    END DO

    RETURN
  END SUBROUTINE BCOND_ASL_CLP
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_GWI

    !--------------------------------------------------------------------------|
    !  GRAVITY-WAVE RADIATION IMPLICIT OPEN BOUNDARY CONDITION (GWI)           |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I1,I2,J,JN
    REAL(SP):: CC,CP

    DO J = 1,IBCN(3)
       JN = OBC_LST(3,J)
       I1 =  I_OBC_N(JN) 
       I2 = NEXT_OBC(JN) 
       CC = SQRT(GRAV_N(I1)*H(I1))*DTE/DLTN_OBC(JN)
       CP = CC + 1.0_SP
       ELF(I1) = (CC*ELF(I2) + EL(I1))/CP
    END DO

    RETURN
  END SUBROUTINE BCOND_GWI
  !==========================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==========================================================================|
  SUBROUTINE BCOND_BKI

    !--------------------------------------------------------------------------|
    !  BLUMBERG AND KHANTA IMPLICIT OPEN BOUNDARY CONDITION                    |
    !--------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I1,I2,J,JN
    REAL(SP):: CC,CP

    DO J = 1,IBCN(4)
       JN = OBC_LST(4,J)
       I1 = I_OBC_N(JN)
       I2 = NEXT_OBC(JN)
       CC = SQRT(GRAV_N(I1)*H(I1))*DTE/DLTN_OBC(JN)
       CP = CC + 1.0_SP

       ELF(I1) = (CC*ELF(I2) + EL(I1)*(1.0_SP-DTE/10800.0_SP))/CP
    END DO

    RETURN
  END SUBROUTINE BCOND_BKI
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE BCOND_ORE

    !------------------------------------------------------------------------------|
    !  ORLANSKI RADIATION EXPLICIT OPEN BOUNDARY CONDITION (ORE)                   |
    !------------------------------------------------------------------------------|

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER  :: I1,I2,I3,J,JN
    REAL(SP) :: CL, MU

    DO J = 1, IBCN(5)
       JN = OBC_LST(5,J)
       I1 = I_OBC_N(JN)
       I2 = NEXT_OBC(JN)
       I3 = NEXT_OBC2(JN)

       CL = (ELM2(I2)-EL(I2))/(EL(I2)+ELM2(I2)-2.0*ELM1(I3))
       IF(CL >= 1.0)THEN
          MU = 1.0
       ELSE IF(CL > 0.0 .AND. CL < 1.0)THEN
          MU = CL
       ELSE
          MU = 0.0
       END IF

       ELF(I1)=(ELM1(I1)*(1.0-MU)+2.0*MU*EL(I2))/(1.0+MU)
    END DO

    RETURN
  END SUBROUTINE BCOND_ORE
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !========================================================================

  SUBROUTINE SETUP_OBCTYPES
    USE LIMS
    IMPLICIT NONE

    INTEGER              :: I,I1,I2,NCNT,IERR,J
    INTEGER, ALLOCATABLE :: TEMP1(:),TEMP2(:),TEMP3(:),TEMP4(:),TEMP5(:),TEMP6(:),TEMP7(:),ITEMP(:)

    ! CREATE THE BC TYPE INTEGERS USED IN THE MODEL
    CALL SEPARATE_OBC ! (MOD_OBCS.F)

    !==============================================================================|
    !   REPORT AND CHECK OBC RESULTS                                               |
    !==============================================================================|
    ALLOCATE(TEMP1(NPROCS),TEMP2(NPROCS),TEMP3(NPROCS),TEMP4(NPROCS))
    ALLOCATE(TEMP5(NPROCS),TEMP6(NPROCS),TEMP7(NPROCS))
    TEMP1(1)  = IOBCN
    TEMP2(1) = IBCN(1)
    TEMP3(1) = IBCN(2)
    TEMP4(1) = IBCN(3)
    TEMP5(1) = IBCN(4)
    TEMP6(1) = IBCN(5)
    !   TEMP7(1) = NOBCGEO

# if defined (MULTIPROCESSOR)
    CALL MPI_ALLGATHER(IOBCN,  1,MPI_INTEGER,TEMP1,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(1),1,MPI_INTEGER,TEMP2,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(2),1,MPI_INTEGER,TEMP3,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(3),1,MPI_INTEGER,TEMP4,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(4),1,MPI_INTEGER,TEMP5,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    CALL MPI_ALLGATHER(IBCN(5),1,MPI_INTEGER,TEMP6,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)
    !   CALL MPI_GATHER(NOBCGEO,1,MPI_INTEGER,TEMP7,1,MPI_INTEGER,0,MPI_FVCOM_GROUP,IERR)
# endif


    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*)  '! BCOND TYPE             : TOTAL BC NODES; BC NODES IN EACH DOMAIN'
       WRITE(IPT,100)'!  IOBCN                 :',IOBCN_GL,   (TEMP1(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(1)               :',IBCN_GL(1), (TEMP2(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(2)               :',IBCN_GL(2), (TEMP3(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(3)               :',IBCN_GL(3), (TEMP4(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(4)               :',IBCN_GL(4), (TEMP5(I),I=1,NPROCS)
       WRITE(IPT,100)'!  IBCN(5)               :',IBCN_GL(5), (TEMP6(I),I=1,NPROCS)
       !   IF(MSR)WRITE(IPT,100)'!  NOBCGEO               :'
       !,NOBCGEO_GL, (TEMP7(I),I=1,NPROCS)
    END IF

    DEALLOCATE(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7)

    RETURN
100 FORMAT(1X,A26,I6," =>",2X,4(I5,1H,))

  END SUBROUTINE SETUP_OBCTYPES
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE SEPARATE_OBC

    !------------------------------------------------------------------------------|
    ! Accumulate separately the amounts of nodes for 11 types of open boundaries   |
    !------------------------------------------------------------------------------|
    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I,I1,I2,I3,I4,I5,II,J

    IBCN = 0
    IBCN_GL = 0

    DO I = 1, IOBCN_GL
       IF(TYPE_OBC_GL(I) == 1 .OR. TYPE_OBC_GL(I) == 2)  IBCN_GL(1) = IBCN_GL(1) + 1
       IF(TYPE_OBC_GL(I) == 3 .OR. TYPE_OBC_GL(I) == 4)  IBCN_GL(2) = IBCN_GL(2) + 1
       IF(TYPE_OBC_GL(I) == 5 .OR. TYPE_OBC_GL(I) == 6)  IBCN_GL(3) = IBCN_GL(3) + 1
       IF(TYPE_OBC_GL(I) == 7 .OR. TYPE_OBC_GL(I) == 8)  IBCN_GL(4) = IBCN_GL(4) + 1
       IF(TYPE_OBC_GL(I) == 9 .OR. TYPE_OBC_GL(I) == 10) IBCN_GL(5) = IBCN_GL(5) + 1
    END DO

    DO I = 1, IOBCN
       IF(TYPE_OBC(I) == 1 .OR. TYPE_OBC(I) == 2)  IBCN(1) = IBCN(1) + 1
       IF(TYPE_OBC(I) == 3 .OR. TYPE_OBC(I) == 4)  IBCN(2) = IBCN(2) + 1
       IF(TYPE_OBC(I) == 5 .OR. TYPE_OBC(I) == 6)  IBCN(3) = IBCN(3) + 1
       IF(TYPE_OBC(I) == 7 .OR. TYPE_OBC(I) == 8)  IBCN(4) = IBCN(4) + 1
       IF(TYPE_OBC(I) == 9 .OR. TYPE_OBC(I) == 10) IBCN(5) = IBCN(5) + 1
    END DO

    I1 = 0
    I2 = 0
    I3 = 0
    I4 = 0
    I5 = 0


    ALLOCATE(OBC_LST(5,MAXVAL(IBCN))) ; OBC_LST = 0

    DO I=1,IOBCN
       IF(TYPE_OBC(I) == 1 .OR. TYPE_OBC(I) == 2)THEN
          I1 = I1 + 1
          OBC_LST(1,I1) = I
       ELSE IF(TYPE_OBC(I) == 3 .OR. TYPE_OBC(I) == 4)THEN
          I2 = I2 + 1
          OBC_LST(2,I2) = I
       ELSE IF(TYPE_OBC(I) == 5 .OR. TYPE_OBC(I) == 6)THEN
          I3 = I3 + 1
          OBC_LST(3,I3) = I
       ELSE IF(TYPE_OBC(I) == 7 .OR. TYPE_OBC(I) == 8)THEN
          I4 = I4 + 1
          OBC_LST(4,I4) = I
       ELSE IF(TYPE_OBC(I) == 9 .OR. TYPE_OBC(I) == 10)THEN
          I5 = I5 + 1
          OBC_LST(5,I5) = I
       END IF
    END DO

    RETURN
  END SUBROUTINE SEPARATE_OBC
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE SETUP_OBC

    !------------------------------------------------------------------------------!
    USE ALL_VARS
#  if defined (SPHERICAL)
    USE MOD_SPHERICAL
#  endif
# if defined (MULTIPROCESSOR)
    USE MOD_PAR  
# endif
    IMPLICIT NONE

    REAL(SP) :: DXC,DYC,DXN,DYN,CROSS,E1,E2,DOTMAX,DOT,DX,DY,DXN_TMP,DYN_TMP
    INTEGER  :: I,J,JJ,INODE,JNODE,I1,I2,IC,N1,N2,N3
    LOGICAL  :: DEBUG

    REAL(SP), ALLOCATABLE :: NXOBC_TMP(:),NYOBC_TMP(:)

# if defined (SPHERICAL)
    REAL(DP) x1_dp,y1_dp,x2_dp,y2_dp,side
# endif
    !------------------------------------------------------------------------------!


    IF(DBG_SET(DBG_LOG))WRITE(IPT,*) "!"
    IF(DBG_SET(DBG_LOG))WRITE(IPT,*) "!           SETTING UP OPEN BOUNDARY CONDITIONS"
    IF(DBG_SET(DBG_LOG))WRITE(IPT,*) "!"

    !--Determine Adjacent Open Boundary Points-------------------------------------!

    ALLOCATE(NADJN_OBC(IOBCN))  ; NADJN_OBC = 0
    ALLOCATE(ADJN_OBC(IOBCN,2)) ; ADJN_OBC = 0

    DO I=1,IOBCN
       INODE = I_OBC_N(I)
       DO J=1,NTSN(INODE)-1
          JNODE = NBSN(INODE,J)
          IF(ISONB(JNODE) == 2 .AND. INODE /= JNODE)THEN
             NADJN_OBC(I) = NADJN_OBC(I) + 1
             ADJN_OBC(I,NADJN_OBC(I)) = JNODE
          END IF
       END DO
    END DO


    DO I=1,IOBCN
       IF(NADJN_OBC(I) == 0)THEN
          WRITE(*,*)'NO ADJACENT NODE FOUND FOR BOUNDARY NODE',I
          WRITE(*,*)'IN PROCESSOR',MYID
          CALL PSTOP
       END IF
    END DO

    !--Determine Adjacent Cells-(Nonlinear Only)-----------------------------------!
    !--Simultaneously Determine INWARD Pointing Normal NXOBC,NYOBC                 !

    ALLOCATE(NADJC_OBC(IOBCN))  ; NADJC_OBC = 0
    ALLOCATE(ADJC_OBC(IOBCN,2)) ; ADJC_OBC = 0
    ALLOCATE(NXOBC(IOBCN)) ; NXOBC = 0
    ALLOCATE(NYOBC(IOBCN)) ; NYOBC = 0
    ALLOCATE(NXOBC_TMP(IOBCN)) ; NXOBC_TMP = 0
    ALLOCATE(NYOBC_TMP(IOBCN)) ; NYOBC_TMP = 0

    DO I=1,IOBCN
       I1 = I_OBC_N(I)

       !!Mark First Cell on Boundary Edge Adjacent to Node I
       I2 = ADJN_OBC(I,1)
       DO J = 1, NTVE(I1)
          IC = NBVE(I1,J)
          N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
          IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
# if defined (SPHERICAL)
             x1_dp=VX(I1); y1_dp=VY(I1)
             x2_dp=VX(I2); y2_dp=VY(I2)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXN = side; DYN = (VY(I2)-VY(I1))*TPI

             x2_dp=XC(IC); y2_dp=YC(IC)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXC = side; DYC = (YC(IC)-VY(I1))*TPI
# else
             DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
             DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)
# endif
             CROSS = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
             NXOBC_TMP(I) =  CROSS*DYN/SQRT(DXN**2 +DYN**2)
             NYOBC_TMP(I) = -CROSS*DXN/SQRT(DXN**2 +DYN**2)
             NXOBC(I) = NXOBC_TMP(I)
             NYOBC(I) = NYOBC_TMP(I)
             NADJC_OBC(I) = NADJC_OBC(I) + 1
             ADJC_OBC(I,NADJC_OBC(I)) = IC
             IF(MOD(TYPE_OBC(I),2) == 1)THEN
                !!Node is Linear, Mark Cell as Linear for Flux Update
                IUCP(IC) = 0
             END IF
          END IF
       END DO

       IF(NADJN_OBC(I) > 1)THEN
          I2 = ADJN_OBC(I,2)
          DO J = 1, NTVE(I1)
             IC = NBVE(I1,J)
             N1 = NV(IC,1) ; N2 = NV(IC,2) ; N3 = NV(IC,3)
             IF( N1-I2 == 0 .OR. N2-I2 == 0 .OR. N3-I2 == 0)THEN
# if defined (SPHERICAL)
                x1_dp=VX(I1); y1_dp=VY(I1)
                x2_dp=VX(I2); y2_dp=VY(I2)
                CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
                DXN = side; DYN = (VY(I2)-VY(I1))*TPI

                x2_dp=XC(IC); y2_dp=YC(IC)
                CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
                DXC = side; DYC = (YC(IC)-VY(I1))*TPI
# else
                DXN = VX(I2)-VX(I1) ; DYN = VY(I2)-VY(I1)
                DXC = XC(IC)-VX(I1) ; DYC = YC(IC)-VY(I1)
# endif
                CROSS = SIGN(1.0_SP,DXC*DYN - DYC*DXN)
                NXOBC_TMP(I) = NXOBC_TMP(I) + CROSS*DYN/SQRT(DXN**2 +DYN**2)
                NYOBC_TMP(I) = NYOBC_TMP(I) - CROSS*DXN/SQRT(DXN**2 +DYN**2)
                NXOBC(I) = NXOBC_TMP(I)/SQRT(NXOBC_TMP(I)**2 + NYOBC_TMP(I)**2)
                NYOBC(I) = NYOBC_TMP(I)/SQRT(NXOBC_TMP(I)**2 + NYOBC_TMP(I)**2)

                NADJC_OBC(I) = NADJC_OBC(I) + 1
                ADJC_OBC(I,NADJC_OBC(I)) = IC
                IF(MOD(TYPE_OBC(I),2) == 1)THEN
                   !!Node is Linear, Mark Cell as Linear for Flux Update
                   IUCP(IC) = 0
                END IF
             END IF
          END DO
       END IF
    END DO

    DEALLOCATE(NXOBC_TMP,NYOBC_TMP)

    !--Determine Adjacent FluxFractions--------------------------------------------!

    ALLOCATE(NFLUXF_OBC(IOBCN))  ; NFLUXF_OBC = 0
    ALLOCATE(FLUXF_OBC(IOBCN,2)) ; FLUXF_OBC  = 0
    DO I=1,IOBCN
       IF(NADJN_OBC(I) == 1)THEN
          NFLUXF_OBC(I) = 1
          FLUXF_OBC(I,1) = 1.
          FLUXF_OBC(I,2) = 0.
       ELSE
          NFLUXF_OBC(I) = 2
          N1   = I_OBC_N(I)
          N2   = ADJN_OBC(I,1)
          N3   = ADJN_OBC(I,2)
# if defined (SPHERICAL)
          x1_dp=VX(N2)
          y1_dp=VY(N2)
          x2_dp=VX(N1)
          y2_dp=VY(N1)
          CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
          E1 = SQRT( side**2 + ((VY(N1)-VY(N2))*TPI)**2)

          x1_dp=VX(N3)
          y1_dp=VY(N3)
          CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
          E2 = SQRT( side**2 + ((VY(N1)-VY(N3))*TPI)**2)
# else
          E1 = SQRT( (VX(N1)-VX(N2))**2 + (VY(N1)-VY(N2))**2)
          E2 = SQRT( (VX(N1)-VX(N3))**2 + (VY(N1)-VY(N3))**2)
# endif
          FLUXF_OBC(I,1) =  E1/(E1+E2)
          FLUXF_OBC(I,2) =  E2/(E1+E2)
       END IF
    END DO
    !--Determine 1st Layer Neighbor for Open Boundary Points-----------------------!
    !--Node Chosen is Node That is Connected to OBC Node and is Oriented           !
    !--Most Normal to the Boundary.  It is not Necessarily the Closest Node        !
    !--Determine also DLTN_OBC, the normal component of the distance between       !
    !--Next_obc and the open boundary node                                         !

    DO I=1,IOBCN
       DOTMAX =  -1.0
       INODE = I_OBC_N(I)
       DO J=1,NTSN(INODE)-1
          JNODE = NBSN(INODE,J)
          IF(ISONB(JNODE) /= 2 .AND. INODE /= JNODE)THEN
# if defined (SPHERICAL)
             x1_dp=VX(INODE)
             y1_dp=VY(INODE)
             x2_dp=VX(JNODE)
             y2_dp=VY(JNODE)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXN_TMP = side
             DYN_TMP = (VY(JNODE)-VY(INODE))*TPI
# else
             DXN_TMP = VX(JNODE)-VX(INODE)
             DYN_TMP = VY(JNODE)-VY(INODE)
# endif
             DXN = DXN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
             DYN = DYN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
             DOT = DXN*NXOBC(I) + DYN*NYOBC(I)
             IF(DOT > DOTMAX)THEN
                DOTMAX = DOT
                NEXT_OBC(I) = JNODE
             END IF
          END IF
       END DO
    END DO


    !--Determine 2nd Layer Neighbor for Open Boundary Points-----------------------!

    DO I=1,IOBCN
       DOTMAX =  -1.0
       INODE = NEXT_OBC(I)
       DO J=1,NTSN(INODE)-1
          JNODE = NBSN(INODE,J)
          IF(ISONB(JNODE) /= 2)THEN
# if defined (SPHERICAL)
             x1_dp=VX(INODE)
             y1_dp=VY(INODE)
             x2_dp=VX(JNODE)
             y2_dp=VY(JNODE)
             CALL ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
             DXN_TMP = side
             DYN_TMP = (VY(JNODE)-VY(INODE))*TPI
# else
             DXN_TMP = VX(JNODE)-VX(INODE)
             DYN_TMP = VY(JNODE)-VY(INODE)
# endif
             DXN = DXN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
             DYN = DYN_TMP/SQRT(DXN_TMP**2 + DYN_TMP**2)
             DOT = DXN*NXOBC(I) + DYN*NYOBC(I)
             IF(DOT > DOTMAX)THEN
                DOTMAX = DOT
                NEXT_OBC2(I) = JNODE
             END IF
          END IF
       END DO
    END DO

    !--Determine DLTN_OBC----------------------------------------------------------!
    ALLOCATE(DLTN_OBC(IOBCN))
    DO I=1,IOBCN
       I1 = I_OBC_N(I)
       I2 = NEXT_OBC(I)

# if defined (SPHERICAL)
       x1_dp=VX(I1)
       y1_dp=VY(I1)
       x2_dp=VX(I2)
       y2_dp=VY(I2)
       call ARCX(x1_dp,y1_dp,x2_dp,y2_dp,side)
       DX = side
       DY = (VY(I2)-VY(I1))*TPI
# else
       DX = VX(I2)-VX(I1)
       DY = VY(I2)-VY(I1)
# endif
       DLTN_OBC(I) = ABS(DX*NXOBC(I) + DY*NYOBC(I))
    END DO

!!$   RETURN
!!$!--Dump Information to Matlab Files for Checking-------------------------------!
!!$
!!$   OPEN(UNIT=81,FILE='mesh.scatter',FORM='formatted')
!!$   DO I=1,M
!!$     WRITE(81,*)vx(i),vy(i)
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='nextobc.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = NEXT_OBC(I)
!!$     WRITE(81,*)VX(I1),VY(I1)
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='nextobc2.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = NEXT_OBC2(I)
!!$     WRITE(81,*)VX(I1),VY(I1)
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='iobcn.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = I_OBC_N(I)
!!$     WRITE(81,*)VX(I1),VY(I1) 
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='obcnorm.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     I1 = I_OBC_N(I)
!!$     WRITE(81,*)NXOBC(I1),NYOBC(I1) 
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='nonlinear.scatter',FORM='formatted')
!!$   DO I=1,IOBCN
!!$     IF(NADJC_OBC(I) > 0) WRITE(81,*)XC(ADJC_OBC(I,1)),YC(ADJC_OBC(I,1)) 
!!$     IF(NADJC_OBC(I) > 1) WRITE(81,*)XC(ADJC_OBC(I,2)),YC(ADJC_OBC(I,2)) 
!!$   END DO
!!$   CLOSE(81)
!!$   OPEN(UNIT=81,FILE='linear.scatter',FORM='formatted')
!!$   DO I=1,N
!!$     IF(IUCP(I)==0)THEN
!!$       WRITE(81,*)XC(I),YC(I) 
!!$     END IF
!!$   END DO
!!$   CLOSE(81)

    RETURN
  END SUBROUTINE SETUP_OBC
  !==============================================================================|

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
# if !defined (SEMI_IMPLICIT)
  SUBROUTINE FLUX_OBN(K)

    USE ALL_VARS
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: K
    INTEGER              :: I,J,C1,C2
    REAL(SP)             :: FF,FLUX 

    FLUXOBN = 0.0_SP

    DO I = 1, IOBCN

       J = I_OBC_N(I)
       !Compute Boundary Flux From Continuity Flux Defect
       FLUX = -(ELF(J)-ELRK(J))*ART1(J)/(ALPHA_RK(K)*DTE)-XFLUX_OBCN(I)

       !Set Flux In Adjacent Nonlinear BC Element 1 (If Exists)
       IF(NADJC_OBC(I) > 0) THEN
          C1 = ADJC_OBC(I,1) 
          FF = FLUXF_OBC(I,1)
          FLUXOBN(C1) = FLUXOBN(C1) + FF*FLUX 
       END IF

       !Set Flux In Adjacent Nonlinear BC Element 2 (If Exists)
       IF(NADJC_OBC(I) > 1) THEN
          C2 =  ADJC_OBC(I,2) 
          FF = FLUXF_OBC(I,2)
          FLUXOBN(C2) = FLUXOBN(C2) + FF*FLUX 
       END IF

    END DO

    RETURN
  END SUBROUTINE FLUX_OBN
# endif
  !========================================================================

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE BCOND_T_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)
    !==============================================================================|
    ! Calculate the OBC for temperature perturbation                               |
    !==============================================================================|
    USE ALL_VARS
    IMPLICIT NONE

    !   INTEGER :: I1,I2,J,JN
    INTEGER :: I,J,J1,J2,K
    REAL(SP):: CC,CP,MU,CL
    REAL(SP):: PERT_NEXT,PERT,T2D_NEXT,T2D
    REAL(SP):: T2D_NEXT1,TM12D_NEXT2,TM12D_NEXT1,TM22D_NEXT1
    REAL(SP):: TTMP(IOBCN,KBM1)

    SELECT CASE(TYPE_TSOBC)

    CASE(1)
       DO K=1,KBM1
          TTMP(I,K) = TF1(J1,K) - T2D_NEXT
       END DO
    CASE(2)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = TF1(J1,K) - T2D_NEXT
          PERT      = T1(J,K) - T2D
          TTMP(I,K) = (CC*PERT_NEXT + PERT)/CP
       END DO
    CASE(3)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = TF1(J1,K) - T2D_NEXT
          PERT      = T1(J,K) - T2D
          TTMP(I,K) = (CC*PERT_NEXT + PERT*(1.0_SP - DTI/10800.0_SP))/CP
       END DO
    CASE(4)
       J2 = NEXT_OBC2(I)
       T2D_NEXT1  =0.0_SP
       TM12D_NEXT2=0.0_SP
       TM12D_NEXT1=0.0_SP
       TM22D_NEXT1=0.0_SP
       DO K=1,KBM1
          T2D_NEXT1  =T2D_NEXT1  +T1(J1,K)*DZ(J1,K)
          TM12D_NEXT2=TM12D_NEXT2+T1M1(J2,K)*DZ(J2,K)
          TM12D_NEXT1=TM12D_NEXT1+T1M1(J,K)*DZ(J,K)
          TM22D_NEXT1=TM22D_NEXT1+T1M2(J1,K)*DZ(J1,K)
       END DO

       DO K=1,KBM1
          CL = ((T1M2(J1,K)-TM22D_NEXT1)-(T1(J1,K)-T2D_NEXT1))/   &
               ((T1(J1,K)-T2D_NEXT1)+(T1M2(J1,K)-TM22D_NEXT1)     &
               -2.0*(T1M1(J2,K)-TM12D_NEXT2))
          IF(CL >= 1.0)THEN
             MU = 1.0
          ELSE IF(CL > 0.0 .AND. CL < 1.0)THEN
             MU = CL
          ELSE
             MU = 0.0
          END IF

          TTMP(I,K)=((T1M1(J,K)-TM12D_NEXT1)*(1.0-MU)     &
               +2.0*MU*(T1(J1,K)-T2D_NEXT1))/(1.0+MU)
       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("INVALID OBC_TS_TYPE IN NML_OPEN_BOUNDARY_CONTROL"&
            &, "See mod_obcs.F")

    END SELECT

    RETURN
  END SUBROUTINE BCOND_T_PERTURBATION
  !========================================================================

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

  !==============================================================================|
  SUBROUTINE BCOND_S_PERTURBATION(S2D_NEXT,S2D,STMP,I,J,J1)
    !==============================================================================|
    ! Calculate the OBC for salinity perturbation                                  |
    !==============================================================================|
    USE ALL_VARS
    IMPLICIT NONE

    INTEGER :: I,J,J1,J2,K
    REAL(SP):: CC,CP,MU,CL
    REAL(SP):: PERT_NEXT,PERT,S2D_NEXT,S2D
    REAL(SP):: S2D_NEXT1,SM12D_NEXT2,SM12D_NEXT1,SM22D_NEXT1
    REAL(SP):: STMP(IOBCN,KBM1)

    SELECT CASE(TYPE_TSOBC)

    CASE(1)
       DO K=1,KBM1
          STMP(I,K) = SF1(J1,K) - S2D_NEXT
       END DO
    CASE(2)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = SF1(J1,K) - S2D_NEXT
          PERT      = S1(J,K) - S2D
          STMP(I,K) = (CC*PERT_NEXT + PERT)/CP
       END DO
    CASE(3)
       CC = SQRT(GRAV_N(J)*H(J))*DTI/DLTN_OBC(I)
       CP = CC + 1.0_SP
       DO K=1,KBM1
          PERT_NEXT = SF1(J1,K) - S2D_NEXT
          PERT      = S1(J,K) - S2D
          STMP(I,K) = (CC*PERT_NEXT + PERT*(1.0_SP - DTI/10800.0_SP))/CP
       END DO
    CASE(4)
       J2 = NEXT_OBC2(I)
       S2D_NEXT1  =0.0_SP
       SM12D_NEXT2=0.0_SP
       SM12D_NEXT1=0.0_SP
       SM22D_NEXT1=0.0_SP
       DO K=1,KBM1
          S2D_NEXT1  =S2D_NEXT1  +S1(J1,K)*DZ(J1,K)
          SM12D_NEXT2=SM12D_NEXT2+S1M1(J2,K)*DZ(J2,K)
          SM12D_NEXT1=SM12D_NEXT1+S1M1(J,K)*DZ(J,K)
          SM22D_NEXT1=SM22D_NEXT1+S1M2(J1,K)*DZ(J1,K)
       END DO

       DO K=1,KBM1
          CL = ((S1M2(J1,K)-SM22D_NEXT1)-(S1(J1,K)-S2D_NEXT1))/   &
               ((S1(J1,K)-S2D_NEXT1)+(S1M2(J1,K)-SM22D_NEXT1)     &
               -2.0*(S1M1(J2,K)-SM12D_NEXT2))
          IF(CL >= 1.0)THEN
             MU = 1.0
          ELSE IF(CL > 0.0 .AND. CL < 1.0)THEN
             MU = CL
          ELSE
             MU = 0.0
          END IF

          STMP(I,K)=((S1M1(J,K)-SM12D_NEXT1)*(1.0-MU)     &
               +2.0*MU*(S1(J1,K)-S2D_NEXT1))/(1.0+MU)
       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("INVALID OBC_TS_TYPE IN NML_OPEN_BOUNDARY_CONTROL"&
            &, "See mod_obcs.F")

    END SELECT

    RETURN
  END SUBROUTINE BCOND_S_PERTURBATION
  !========================================================================

!!$

END MODULE MOD_OBCS
