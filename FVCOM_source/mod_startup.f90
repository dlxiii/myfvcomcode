










!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

MODULE MOD_STARTUP
  USE MOD_UTILS
  USE MOD_NCTOOLS
  USE MOD_INPUT
  USE ALL_VARS
  USE EQS_OF_STATE
  USE MOD_WD
  USE SINTER



  


  USE MOD_DYE, ONLY: DYE
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: READ_SSH
  PUBLIC :: READ_UV
  PUBLIC :: READ_TURB
  PUBLIC :: READ_TS
  PUBLIC :: READ_WETDRY

  PUBLIC :: STARTUP
!# if defined (DATA_ASSIM)
!  PUBLIC :: STARTUP_ASSIM
!# endif 
CONTAINS
!==============================================================================!
  SUBROUTINE STARTUP
    IMPLICIT NONE

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: Startup"

    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*  )'!'
       WRITE(IPT,*  )'!           SETTING INITIAL CONDITIONS '
       WRITE(IPT,*  )'!'
    END IF

    IF (ASSOCIATED(NC_START))THEN
       IF(DBG_SET(DBG_LOG)) write(ipt,*) "! NC_START FILE INFO:"
       CALL PRINT_FILE(NC_START)
    END IF
       
!    SET THE SEA SURFACE HEIGHT
    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING SEA SURFACE HEIGHT"
    SELECT CASE(STARTUP_TYPE)
! =================================================
    CASE(STARTUP_TYPE_HOTSTART)
! =================================================

       CALL READ_SSH

       IF(WETTING_DRYING_ON) CALL READ_WETDRY





! ------- New: Karsten Lettmann, June 2016 -----------------------
! ---------- end  new --------------------------------------------

!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA ICE"



!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA ICE222"




    CALL READ_DYE

! =================================================
    CASE(STARTUP_TYPE_CRASHRESTART)
! =================================================
       CALL READ_SSH
     
       IF(WETTING_DRYING_ON) CALL READ_WETDRY





! ------- New: Karsten Lettmann, June 2016 -----------------------
! ---------- end  new --------------------------------------------





    CALL READ_DYE

! =================================================
    CASE(STARTUP_TYPE_COLDSTART)
! =================================================

       CALL SET_WATER_DEPTH
       IF(WETTING_DRYING_ON) CALL SET_WD_DATA

! =================================================
    CASE DEFAULT
! =================================================
       CALL FATAL_ERROR("STARTUP: UNKNOWN STARTUP TYPE: "//TRIM(STARTUP_TYPE))
    END SELECT


    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING VELOCITY FIELDS"

! SET STARTUP VALUES FOR VELOCITY
    SELECT CASE (STARTUP_UV_TYPE)
    CASE (STARTUP_TYPE_OBSERVED)
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_LINEAR) 
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_CONSTANT)
       !CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
       CALL SET_CONSTANT_UV
    CASE(STARTUP_TYPE_DEFAULT)
       ! OKAY - it is zero
    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_UV

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_UV_TYPE")
    END SELECT

    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING TURBULENCE FIELDS"

! SET STARTUP VALUES FOR TURBULENCE
    SELECT CASE(STARTUP_TURB_TYPE)
    CASE(STARTUP_TYPE_OBSERVED) 
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_LINEAR)
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_CONSTANT)
       CALL FATAL_ERROR("I DON'T KNOW HOW TO DO THAT KIND OF STARTUP")
    CASE(STARTUP_TYPE_DEFAULT)

       CALL SET_DEFAULT_TURB

    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_TURB

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_TURB_TYPE")
    END SELECT

    if (dbg_set(dbg_log)) write(ipt,*) &
         & "! INITIALIZING TEMPERATURE AND SALINITY"

! SET STARTUP VALUES FOR TEMPERATER AND SALINITY
    SELECT CASE(STARTUP_TS_TYPE)
    CASE(STARTUP_TYPE_OBSERVED) 

       CALL SET_OBSERVED_TS

    CASE(STARTUP_TYPE_LINEAR)

       CALL SET_LINEAR_TS

    CASE(STARTUP_TYPE_CONSTANT)

       CALL SET_CONSTANT_TS

    CASE(STARTUP_TYPE_DEFAULT)
       CALL FATAL_ERROR("There is no default startup for Temperature and Salinity")
    CASE(STARTUP_TYPE_SETVALUES)

       CALL READ_TS

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN STARTUP_TS_TYPE")
    END SELECT



    ! ONCE WE HAVE STARTED THE MODEL POINT THE START FILE OBJECT TO
    ! ITS OWN OUTPUT TO RELOAD OLD TIME STATES
    IF(.not. associated(NC_START,NC_RST)) THEN
       
       ! SOME START UPS DO NOT HAVE A START FILE
       IF(Associated(NC_START)) CALL KILL_FILE(NC_START)

       NC_START => NC_RST
    END IF

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: Startup"


  END SUBROUTINE STARTUP
!==============================================================================!

!# if defined (DATA_ASSIM)
!!==============================================================================!
!  SUBROUTINE STARTUP_ASSIM
!# if defined(SEDIMENT) 
!  USE MOD_SED, ONLY: sed_hot_start
!# endif
!    IMPLICIT NONE

!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: Startup_ASSIM"
       
!!    SET THE SEA SURFACE HEIGHT
!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING SEA SURFACE HEIGHT"
!       CALL READ_SSH

!       IF(WETTING_DRYING_ON) CALL READ_WETDRY

!# if defined (NH)
!       CALL READ_NH
!# endif

!# if defined (ATMO_TIDE)
!       CALL READ_ATMO
!# endif

!# if defined (EQUI_TIDE)
!       CALL READ_EQI
!# endif

!# if defined (ICE)    
!       CALL  READ_ICE
!       CALL  AGGREGATE  !!  ggao 07312008
!# endif    

!# if defined (SEDIMENT)
!  if(sed_hot_start)then
!    CALL READ_SED
!  endif
!# endif


!       CALL READ_UV

!       CALL READ_TURB

!       CALL READ_TS
!    if (dbg_set(dbg_log)) write(ipt,*) &
!         & "! INITIALIZING TS"


!#   if defined (WATER_QUALITY)

!       CALL READ_WQM

!#   endif

!#   if defined (BioGen)

!       CALL READ_BIO

!#   endif


!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: Startup_Assim"



!  END SUBROUTINE STARTUP_ASSIM
!!==============================================================================!
!# endif

  SUBROUTINE READ_SSH
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"
    
    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD EL
    VAR => FIND_VAR(NC_START,'zeta',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zeta'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, EL)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ET
    VAR => FIND_VAR(NC_START,'et',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'et'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ET)
    CALL NC_READ_VAR(VAR,STKCNT)

    !----------------------------------------------------------------
    ! Read the most recent bathymetry if Morphodynamics is Active
    !----------------------------------------------------------------


    !----------------------------------------------------------------
    ! Given SSH and Bathy, Update the Bathymetry 
    !----------------------------------------------------------------
    D  = H + EL
    DT = H + ET

    
    CALL N2E2D(H,H1)
    CALL N2E2D(EL,EL1)
    CALL N2E2D(D,D1)
    CALL N2E2D(DT,DT1)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"

  END SUBROUTINE READ_SSH
!==============================================================================!
!==============================================================================!
  SUBROUTINE READ_WETDRY
    USE MOD_WD
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_WETDRY"

    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD ISWETN
    VAR => FIND_VAR(NC_START,'wet_nodes',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETN)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ISWETC
    VAR => FIND_VAR(NC_START,'wet_cells',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETC)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD ISWETNT
    VAR => FIND_VAR(NC_START,'wet_nodes_prev_int',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes_prev_int'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETNT)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ISWETCT
    VAR => FIND_VAR(NC_START,'wet_cells_prev_int',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_int'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, ISWETCT)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD ISWETCE
    VAR => FIND_VAR(NC_START,'wet_cells_prev_ext',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_ext'&
         & IN THE HOTSTART FILE OBJECT")
    ! ------ new: Karsten Lettmann, May 2016 -------------
    !CALL NC_CONNECT_AVAR(VAR, ISWETC)  ! original line
    CALL NC_CONNECT_AVAR(VAR, ISWETCE)
    ! ------------- end new ------------------------------
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_WETDRY"
    
  END SUBROUTINE READ_WETDRY
!==============================================================================!
!==============================================================================!
  SUBROUTINE READ_ATMO
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_ATMO" 

    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD ATMO
    VAR => FIND_VAR(NC_START,'el_atmo',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'el_atmo' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, el_atmo)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_ATMO"

  END SUBROUTINE READ_ATMO
!==============================================================================!
  SUBROUTINE READ_EQI
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_EQI" 
    
    STKCNT = NC_START%FTIME%PREV_STKCNT

    ! LOAD ATMO
    VAR => FIND_VAR(NC_START,'el_eqi',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR&
         &("COULD NOT FIND VARIABLE 'el_eqi' IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, el_eqi)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_EQI" 
   
  END SUBROUTINE READ_EQI
!==============================================================================!

! ------- New: Karsten Lettmann, June 2016 -----------------------
! ------------- end new ---------------------------------------------------------

  SUBROUTINE READ_DYE
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT, K

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_DYE"

    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD DYE
    VAR => FIND_VAR(NC_START,'DYE',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'DYE'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, DYE)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_DYE"

  END SUBROUTINE READ_DYE
!==============================================================================!
  SUBROUTINE READ_TS
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT, K
    REAL(SP), DIMENSION(0:MT,KB) :: PRESSURE
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_TS" 
   
    STKCNT = NC_START%FTIME%PREV_STKCNT


    ! LOAD TEMPERATURE
    VAR => FIND_VAR(NC_START,'temp',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'temp'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, T1)
    CALL NC_READ_VAR(VAR,STKCNT)


    ! LOAD MEAN INITIAL TEMPERATURE
    VAR => FIND_VAR(NC_START,'tmean1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tmean1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, tmean1)
    CALL NC_READ_VAR(VAR)

    ! LOAD SALINITY
    VAR => FIND_VAR(NC_START,'salinity',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'saltinity'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, S1)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD MEAN INITIAL SALINITY
    VAR => FIND_VAR(NC_START,'smean1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'smean1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, smean1)
    CALL NC_READ_VAR(VAR)

    ! LOAD DENSITY
    VAR => FIND_VAR(NC_START,'rho1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rho1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, RHO1)
    CALL NC_READ_VAR(VAR,STKCNT)

    ! LOAD MEAN DENSITY
    VAR => FIND_VAR(NC_START,'rmean1',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rmean1'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, rmean1)
    CALL NC_READ_VAR(VAR)


    ! AVERAGE FROM CELLS TO FACE CENTERS


!JQI    SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
!JQI    CASE(SW_DENS1)

!JQI       ! SET MEAN DENSITY
!JQI       DO K=1,KBM1
!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
!JQI       END DO
!JQI       CALL FOFONOFF_MILLARD(SMEAN1,TMEAN1,PRESSURE,0.0_SP,RMEAN1)
!JQI       RMEAN1(0,:)=0.0_SP
!JQI       RMEAN1(:,KB)=0.0_SP

!JQI       ! SET REAL DENSITY
!JQI       CALL DENS1 ! GENERIC CALL TO FOFONOFF_MILLARD FOR S1,T1...
       
!JQI    CASE(SW_DENS2)
!JQI       ! SET MEAN DENSITY
!JQI       CALL DENS2G(SMEAN1,TMEAN1,RMEAN1)
!JQI       RMEAN1(0,:)=0.0_SP
!JQI       RMEAN1(:,KB)=0.0_SP

!JQI       ! SET REAL DENSITY
!JQI       CALL DENS2 ! GENERIC CALL TO DENS2G FOR S1,T1...

!JQI    CASE(SW_DENS3)

!JQI       ! SET MEAN DENSITY
!JQI       DO K=1,KBM1
!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.01_SP
!JQI       END DO
!JQI       CALL JACKET_MCDOUGALL(SMEAN1,TMEAN1,PRESSURE,RMEAN1)
!JQI       RMEAN1(0,:)=0.0_SP
!JQI       RMEAN1(:,KB)=0.0_SP

!JQI       ! SET REAL DENSITY
!JQI       CALL DENS3 ! GENERIC CALL TO JACKET_MCDOUGALL FOR S1,T1.

!JQI    CASE DEFAULT
!JQI       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
!JQI            & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
!JQI    END SELECT

!J. GE
    T0=T1
    S0=S1
!J. GE    
    CALL N2E3D(T1,T)
    CALL N2E3D(S1,S)
    CALL N2E3D(RHO1,RHO)
    CALL N2E3D(Tmean1,Tmean)
    CALL N2E3D(Smean1,Smean)
    CALL N2E3D(Rmean1,Rmean)

   
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_TS" 
   
 END SUBROUTINE READ_TS
!==============================================================================!
  SUBROUTINE SET_OBSERVED_TS
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT
    INTEGER KSL                                 !!NUMBER OF STANDARD SEA LEVELS 
    REAL(SP), ALLOCATABLE, TARGET :: DPTHSL(:)          !!DEPTH AT STANDARD SEA LEVEL
    REAL(SP), ALLOCATABLE, TARGET :: TSL(:,:),SSL(:,:)  !!T/S AT STANDARD SEA LEVEL
    REAL(SP),ALLOCATABLE          :: TA(:),SA(:)

    REAL(SP),DIMENSION(KBM1)      :: TI,SI,ZI
    REAL(SP),DIMENSION(0:MT,KB) :: PRESSURE
    INTEGER :: K, I, IB, IERR
    REAL(SP) :: SCOUNT, FAC, RBUF
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_OBSERVED_TS" 

    TI       = 0.0_SP
    SI       = 0.0_SP
    ZI       = 0.0_SP
    PRESSURE = 0.0_SP

    STKCNT = NC_START%FTIME%PREV_STKCNT

    DIM => FIND_DIM(NC_START,'ksl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND DIMENSION 'ksl'&
         & IN THE STARTUP FILE OBJECT")

    KSL = DIM%DIM
    ! ALLOCATE SPACE
    ALLOCATE(TSL(MT,KSL))
    ALLOCATE(SSL(MT,KSL))
    ALLOCATE(DPTHSL(KSL))
    TSL   = 0.0_SP
    SSL   = 0.0_SP
    DPTHSL= 0.0_SP

    ! ALLOCATE SPACE FOR THE AVERAGE TEMPERATURE AT EACH DEPTH
    ALLOCATE(TA(KSL),SA(KSL))
    TA = 0.0_SP
    SA = 0.0_SP

    ! READ THE DATA
    VAR => FIND_VAR(NC_START,'tsl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tsl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, TSL)
    CALL NC_READ_VAR(VAR,STKCNT)
    CALL NC_DISCONNECT(VAR)

    VAR => FIND_VAR(NC_START,'ssl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ssl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, SSL)
    CALL NC_READ_VAR(VAR,STKCNT)
    CALL NC_DISCONNECT(VAR)


    VAR => FIND_VAR(NC_START,'zsl',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zsl'&
         & IN THE STARTUP FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, DPTHSL)
    CALL NC_READ_VAR(VAR)
    CALL NC_DISCONNECT(VAR)


!==============================================================================|
!    CALCULATE HORIZONTAL AVERAGE VALUES OF TEMPERATURE/SALINITY/DENSITY       !
!==============================================================================|
   
      where(ssl <0.0_SP) ssl=0.0
 
    ! GET THE AVERAGE T/S AT EACH DEPTH ON STANDARD LEVELS
    IF(SERIAL)THEN
       DO K = 1, KSL
          SCOUNT = 0.0_SP
          DO I = 1, MT
             IF(-H(I) <= DPTHSL(K)) THEN
                SCOUNT = SCOUNT + 1.0_SP
                TA(K) = TA(K) + TSL(I,K)
                SA(K) = SA(K) + SSL(I,K)
             END IF
          END DO
          IF(SCOUNT >= 1.0_SP)THEN
             TA(K)=TA(K)/SCOUNT
             SA(K)=SA(K)/SCOUNT
          END IF
       END DO
    END IF
    
    
    IF(PAR)THEN
       
       DO K = 1, KSL
          SCOUNT = 0.0_SP
          DO I = 1, M
             IF(-H(I) <= DPTHSL(K)) THEN
                IF(NDE_ID(I) == 0)THEN !!INTERNAL NODE
                   SCOUNT = SCOUNT + 1.0_SP
                   TA(K) = TA(K) + TSL(I,K)
                   SA(K) = SA(K) + SSL(I,K)
                ELSE  !!BOUNDARY NODE, ACCUMULATE FRACTION ONLY
                   DO IB = 1,NBN
                      IF(BN_LOC(IB) == I)THEN
                         FAC = 1.0_SP/FLOAT(BN_MLT(IB))
                         SCOUNT = SCOUNT + FAC 
                         TA(K) = TA(K) + FAC*TSL(I,K)
                         SA(K) = SA(K) + FAC*SSL(I,K)
                      END IF
                   END DO
                END IF
             END IF
          END DO
          CALL MPI_ALLREDUCE(TA(K),RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          TA(K) = RBUF
          CALL MPI_ALLREDUCE(SA(K),RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          SA(K) = RBUF
          CALL MPI_ALLREDUCE(SCOUNT,RBUF,1,MPI_F,MPI_SUM,MPI_FVCOM_GROUP,IERR)
          SCOUNT = RBUF

          IF(SCOUNT >= 1.0_SP)THEN
             TA(K)=TA(K)/SCOUNT
             SA(K)=SA(K)/SCOUNT
          END IF
          
       END DO !!LOOP OVER KSL
    END IF  !!PARALLEL
    
    ! NOW INTERPOLATE FROM STANDARD LEVELS TO THE VALUE AT EACH NODE
    DO I=1,MT
       DO K=1,KBM1
          ZI(K)=ZZ(I,K)*D(I)+EL(I)
       END DO
       
       ! LEVEL AVERAGE T AND S
       CALL SINTER_EXTRP_UP(DPTHSL,TA,ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(DPTHSL,SA,ZI,SI,KSL,KBM1)
       
       TMEAN1(I,1:KBM1) = TI(:)
       SMEAN1(I,1:KBM1) = SI(:)
              
       ! REAL T AND S
       CALL SINTER_EXTRP_UP(DPTHSL,TSL(I,:),ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(DPTHSL,SSL(I,:),ZI,SI,KSL,KBM1)
       
       T1(I,1:KBM1) = TI(:)
       S1(I,1:KBM1) = SI(:)
!J. Ge
       T0(I,1:KBM1) = TI(:)
       S0(I,1:KBM1) = SI(:)
!J. Ge
    END DO

       where(S1<0.0_SP) S1=0.0

    IF(.NOT.BAROTROPIC)THEN
     SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
     CASE(SW_DENS1)

       ! SET MEAN DENSITY
       DO K=1,KBM1
          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
       END DO
       CALL FOFONOFF_MILLARD(SMEAN1,TMEAN1,PRESSURE,0.0_SP,RMEAN1)
       RMEAN1(0,:)=0.0_SP
       RMEAN1(:,KB)=0.0_SP

       ! SET REAL DENSITY
       CALL DENS1 ! GENERIC CALL TO FOFONOFF_MILLARD FOR S1,T1...

     CASE(SW_DENS2)

       ! SET MEAN DENSITY
       CALL DENS2G(SMEAN1,TMEAN1,RMEAN1)
       RMEAN1(0,:)=0.0_SP
       RMEAN1(:,KB)=0.0_SP

       ! SET REAL DENSITY
       CALL DENS2 ! GENERIC CALL TO DENS2G FOR S1,T1...

     CASE(SW_DENS3)
     
      ! SET MEAN DENSITY
       DO K=1,KBM1
          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.01_SP
       END DO
       CALL JACKET_MCDOUGALL(SMEAN1,TMEAN1,PRESSURE,RMEAN1)
       RMEAN1(0,:)=0.0_SP
       RMEAN1(:,KB)=0.0_SP

       ! SET REAL DENSITY
       CALL DENS3 ! GENERIC CALL TO JACKET_MCDOUGALL FOR S1,T1.

     CASE DEFAULT
      CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
     END SELECT
    ELSE
     RHO1   = 2.3E-2_SP
     RHO    = 2.3E-2_SP
     RMEAN1 = 2.3E-2_SP
    END IF 
   

   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   CALL N2E3D(Tmean1,Tmean)
   CALL N2E3D(Smean1,Smean)
   CALL N2E3D(Rmean1,Rmean)
   ! DENSITY IS ALREADY INTERPOLATED TO ELEMENTS FOR RHO IN DENSX

   DEALLOCATE(TSL,SSL,DPTHSL)
   DEALLOCATE(TA,SA)


   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_OBSERVED_TS" 
    
  END SUBROUTINE SET_OBSERVED_TS
!==============================================================================!
  SUBROUTINE SET_LINEAR_TS
    IMPLICIT NONE

    ! BY DEFINITION OF LINEAR...
    INTEGER, PARAMETER :: KSL =2

    REAL(SP), DIMENSION(KBM1)   :: TI,SI,ZI
    REAL(SP), DIMENSION(KSL)   :: TA,SA,ZA
    INTEGER :: K, I

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_LINEAR_TS" 

    ZA(1) = 0.0_SP
    ZA(2) = STARTUP_DMAX
    
    TA(1) = STARTUP_T_VALS(1)
    TA(2) = STARTUP_T_VALS(2)
    
    SA(1) = STARTUP_S_VALS(1)
    SA(2) = STARTUP_S_VALS(2)
    
    ! THE HORIZONTAL AVERAGE VALUE IS ALSO THE TRUE VALUE SINCE THE
    ! LINEAR EQUATION DOES NOT VARY WITH LOCATION 

    DO I = 1, MT
       DO K=1,KBM1
          ! CALCULATE ZI relative to z=0
          ZI(K)=ZZ(I,K)*D(I)+EL(I)
       END DO
    
       CALL SINTER_EXTRP_UP(ZA,TA,ZI,TI,KSL,KBM1)
       CALL SINTER_EXTRP_UP(ZA,SA,ZI,SI,KSL,KBM1)
    
       T1(I,1:kbm1) = TI
       S1(I,1:kbm1) = SI
!J. Ge
       T0(I,1:kbm1) = TI
       S0(I,1:kbm1) = SI
!J. Ge

    END DO

    IF(.NOT.BAROTROPIC)THEN
     SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
     CASE(SW_DENS1)
       CALL DENS1 ! USE GENERIC INTERFACE TO DENSX
     CASE(SW_DENS2)
       CALL DENS2
     CASE(SW_DENS3)
       CALL DENS3
     CASE DEFAULT
       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
     END SELECT
    ELSE
     RHO1 = 2.3E-2_SP
     RHO  = 2.3E-2_SP
    END IF 
   
   ! FOR LINEAR AND CONSTANT, THE MEAN IS EQUAL TO THE INITIAL VALUE

   TMEAN1=T1
   SMEAN1=S1
   RMEAN1=RHO1
   RMEAN=RHO
   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   TMEAN=T
   SMEAN=S
   
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_LINEAR_TS" 
    
  END SUBROUTINE SET_LINEAR_TS
!==============================================================================!
  SUBROUTINE SET_CONSTANT_TS
    IMPLICIT NONE
    
    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_CONSTANT_TS" 


    T1(:,1:KBM1) = STARTUP_T_VALS(1)
    S1(:,1:KBM1) = STARTUP_S_VALS(1)
!J. Ge
    T0(:,1:KBM1) = STARTUP_T_VALS(1)
    S0(:,1:KBM1) = STARTUP_S_VALS(1)
!J. Ge

    IF(.NOT.BAROTROPIC)THEN
     SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
     CASE(SW_DENS1)
       CALL DENS1
     CASE(SW_DENS2)
       CALL DENS2
     CASE(SW_DENS3)
       CALL DENS3
     CASE DEFAULT
       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
            & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
     END SELECT
    ELSE
     RHO1 = 2.3E-2_SP
     RHO  = 2.3E-2_SP
    END IF 

   ! FOR LINEAR AND CONSTANT, THE MEAN IS EQUAL TO THE INITIAL VALUE

   TMEAN1=T1
   SMEAN1=S1
   RMEAN1=RHO1

   CALL N2E3D(T1,T)
   CALL N2E3D(S1,S)
   
   TMEAN=T
   RMEAN=RHO
   SMEAN=S

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_CONSTANT_TS" 

 END SUBROUTINE SET_CONSTANT_TS
!==============================================================================!
!==============================================================================!
  SUBROUTINE SET_CONSTANT_UV
    IMPLICIT NONE

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_CONSTANT_UV"


    U(:,1:KBM1) = STARTUP_U_VALS    
    V(:,1:KBM1) = STARTUP_V_VALS
    UA          = STARTUP_U_VALS
    VA          = STARTUP_V_VALS
    WTS         = 0.0          
    W           = 0.0

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_CONSTANT_UV"

 END SUBROUTINE SET_CONSTANT_UV
!==============================================================================!

  SUBROUTINE READ_UV
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_UV" 


    STKCNT = NC_START%FTIME%PREV_STKCNT

    VAR => FIND_VAR(NC_START,'u',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'u'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, U)
    CALL NC_READ_VAR(VAR,STKCNT)


    VAR => FIND_VAR(NC_START,'v',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'v'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, V)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'omega',FOUND)
    IF(FOUND) THEN
       CALL NC_CONNECT_AVAR(VAR, WTS)
       CALL NC_READ_VAR(VAR,STKCNT)
       
       CALL N2E3D(WTS,W)
    ELSE
       VAR => FIND_VAR(NC_START,'w',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'w' &
            & or 'omega' IN THE HOTSTART FILE OBJECT")
       CALL NC_CONNECT_AVAR(VAR, W)
       CALL NC_READ_VAR(VAR,STKCNT)
    END IF

    VAR => FIND_VAR(NC_START,'ua',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ua'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, UA)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'va',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'va'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, VA)
    CALL NC_READ_VAR(VAR,STKCNT)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_UV" 

  END SUBROUTINE READ_UV
!==============================================================================!
  SUBROUTINE READ_TURB
    IMPLICIT NONE
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCDIM),  POINTER :: DIM
    LOGICAL :: FOUND
    INTEGER :: STKCNT

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_TURB" 

    STKCNT = NC_START%FTIME%PREV_STKCNT


    VAR => FIND_VAR(NC_START,'q2',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, Q2)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'q2l',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2l'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, Q2L)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'l',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'l'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, L)
    CALL NC_READ_VAR(VAR,STKCNT)


    VAR => FIND_VAR(NC_START,'km',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'km'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, KM)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'kq',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kq'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, KQ)
    CALL NC_READ_VAR(VAR,STKCNT)

    VAR => FIND_VAR(NC_START,'kh',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kh'&
         & IN THE HOTSTART FILE OBJECT")
    CALL NC_CONNECT_AVAR(VAR, KH)
    CALL NC_READ_VAR(VAR,STKCNT)


    CALL N2E3D(KM,KM1)

    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_TURB" 

  END SUBROUTINE READ_TURB
!==============================================================================|
   SUBROUTINE SET_DEFAULT_TURB         
!==============================================================================|
!   Initialize Turbulent Kinetic Energy and Length Scale                       |
!==============================================================================|
   IMPLICIT NONE
   INTEGER :: I,K
!==============================================================================|
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_DEFAULT_TURB" 
   
   DO K=1,KBM1
        AAM(:,K) = NN_HVC(:)
   END DO

!
!------------------------BOUNDARY VALUES---------------------------------------!
!

   DO I = 1, MT
     KM(I,1)   = 0.0_SP
     KM(I,KB)  = 0.0_SP
     KH(I,1)   = 0.0_SP
     KH(I,KB)  = 0.0_SP
     KQ(I,1)   = 0.0_SP
     KQ(I,KB)  = 0.0_SP
     L(I,1)    = 0.0_SP
     L(I,KB)   = 0.0_SP
     Q2(I,1)   = 0.0_SP
     Q2(I,KB)  = 0.0_SP
     Q2L(I,1)  = 0.0_SP
     Q2L(I,KB) = 0.0_SP
   END DO


!
!------------------------INTERNAL VALUES---------------------------------------!
!
   DO  K = 2, KBM1
     DO I = 1, MT
       IF (D(I) > 0.0_SP) THEN
         Q2(I,K)  = 1.E-8
         Q2L(I,K) = 1.E-8
         L(I,K)   = 1.
         KM(I,K)  = 2.*UMOL
         KQ(I,K)  = 2.*UMOL
         KH(I,K)  = 2.*UMOL / VPRNU
       END IF
     END DO
   END DO
   CALL N2E3D(KM,KM1)

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_DEFAULT_TURB" 

 END SUBROUTINE SET_DEFAULT_TURB
!==============================================================================|

!==============================================================================|
!   READ IN STATIC WATER DEPTH AND CALCULATE RELATED QUANTITIES                |
!                                                                              |
!   INPUTS: H(NNODE) BATHYMETRIC DEPTH AT NODES				       |
!   INITIALIZES: D(NNODE) DEPTH AT NODES				       |
!   INITIALIZES: DT(NNODE) ???					               |
!   INITIALIZES: H1(NNODE) BATHYMETRIC DEPTH AT ELEMENTS		       |
!   INITIALIZES: D1(NNODE) DEPTH AT NODES                		       | 
!   INITIALIZES: DT1(NNODE) ??                                   	       |
!==============================================================================|

   SUBROUTINE SET_WATER_DEPTH         
!------------------------------------------------------------------------------|
   USE MOD_OBCS
   USE MOD_NESTING
   IMPLICIT NONE
   REAL(SP) :: TEMP
   INTEGER  :: I,K,J1,J2
!------------------------------------------------------------------------------|
   IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: SET_WATER_DEPTH" 


!
!  ADJUST STATIC HEIGHT AND CALCULATE DYNAMIC DEPTHS (D) AND (DT)
!
   H  = H + STATIC_SSH_ADJ
   D  = H + EL
   DT = H + ET

!
!  ADJUST SIGMA VALUES ON OUTER BOUNDARY 
!

   IF(IOBCN > 0 .AND. OBC_DEPTH_CONTROL_ON) THEN
     IF(.NOT. NESTING_ON)THEN
     DO I = 1,IOBCN
       J1 = I_OBC_N(I)
       J2 = NEXT_OBC(I)
       H(J1) = H(J2)
       D(J1) = D(J2)
       DT(J1) = DT(J2)
       
       DO K = 1,KB
         Z(J1,K)   = Z(J2,K)
         ZZ(J1,K)  = ZZ(J2,K)
         DZ(J1,K)  = DZ(J2,K)
         DZZ(J1,K) = DZZ(J2,K)
       END DO	 
     END DO
     END IF
   END IF

! 
!  CALCULATE FACE-CENTERED VALUES OF BATHYMETRY AND DEPTH
!
   DO I=1,NT    
     H1(I)  = (H(NV(I,1))+H(NV(I,2))+H(NV(I,3)))/3.0_SP
     D1(I)  = H1(I)+EL1(I)
     DT1(I) = H1(I)+ET1(I)
     
     DO K = 1,KB
       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0_SP
     END DO  
   END DO

!-----AFTER MODIFYING BOUNDARY SIGMA VALUES ------------------------------------!
!-----RECOMPUTE SIGMA DERIVATIVES AND INTRA SIGMA LEVELS AGAIN ON CELL----------!

   DO K=1,KB-1
     DO I=1,NT
       DZ1(I,K)  = Z1(I,K)-Z1(I,K+1)
       ZZ1(I,K)  = .5_SP*(Z1(I,K)+Z1(I,K+1))
     END DO
   END DO

   DO I=1,NT
     ZZ1(I,KB) = 2.0_SP*ZZ1(I,KB-1)-ZZ1(I,KB-2)
   END DO

   DO K=1,KBM2
     DO I=1,NT
       DZZ1(I,K) = ZZ1(I,K)-ZZ1(I,K+1)
     END DO
   END DO
   
   DZZ1(:,KBM1) = 0.0_SP
   DZ1(:,KB)    = 0.0_SP

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: SET_WATER_DEPTH" 
   RETURN
 END SUBROUTINE SET_WATER_DEPTH
!==============================================================================|








END MODULE MOD_STARTUP
