! $Id: mod_assim.F,v 1.5.2.11 2008/05/21 00:15:09 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.5.2.11 $
!!NOTES TO LIU, IS TEMP/SAL DATA SHIFTED WITH THE 744 LIKE CURRENT????

MODULE MOD_ASSIM

  USE CONTROL
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_PREC
  USE MOD_WD

#  if defined (DYE_RELEASE)
  USE MOD_DYE
#  endif 

  IMPLICIT NONE
  SAVE


  !
  !--Data Assimilation Parameters for SST Assimilation
  !
  LOGICAL  :: SST_ASSIM                    !!TRUE IF SST ASSIMILATION ACTIVE
  CHARACTER(LEN=80) SST_ASSIM_FILE         !!FILE NAME FOR SST DATA
  REAL(SP) :: SST_RADIUS                   !!SEARCH RADIUS FOR INTERPOLATION POINTS
  REAL(SP) :: SST_WEIGHT_MAX                    
  REAL(SP) :: SST_TIMESCALE
  REAL(SP) :: SST_TIME_WINDOW            !!TIME WINDOW FOR OBSERVATION ASSIMILATION
  INTEGER  :: SST_N_PER_INTERVAL           !! ASSUMING DAILY DATA:
  !! AVERAGE OVER

  NAMELIST /NML_SST_ASSIMILATION/      &
       & SST_ASSIM,                    &
       & SST_ASSIM_FILE,               &
       & SST_RADIUS,                   &
       & SST_WEIGHT_MAX,               &
       & SST_TIMESCALE,                &
       & SST_TIME_WINDOW,              &
       & SST_N_PER_INTERVAL

  !
  !--Data Assimilation Parameters for SST GRID Assimilation
  !
  LOGICAL  :: SSTGRD_ASSIM                    !!TRUE IF SST ASSIMILATION ACTIVE
  CHARACTER(LEN=80) SSTGRD_ASSIM_FILE         !!FILE NAME FOR SST DATA
  REAL(SP) :: SSTGRD_WEIGHT_MAX                    
  REAL(SP) :: SSTGRD_TIMESCALE
  REAL(SP) :: SSTGRD_TIME_WINDOW            !!TIME WINDOW FOR OBSERVATION ASSIMILATION
  INTEGER  :: SSTGRD_N_PER_INTERVAL           !! ASSUMING DAILY DATA:
  !! AVERAGE OVER

  NAMELIST /NML_SSTGRD_ASSIMILATION/      &
       & SSTGRD_ASSIM,                    &
       & SSTGRD_ASSIM_FILE,               &
       & SSTGRD_WEIGHT_MAX,               &
       & SSTGRD_TIMESCALE,                &
       & SSTGRD_TIME_WINDOW,              &
       & SSTGRD_N_PER_INTERVAL

  REAL(SP) :: SST_WM                    
  REAL(SP) :: SST_TSCALE
  REAL(SP) :: SST_TWINDOW            !!TIME WINDOW FOR OBSERVATION ASSIMILATION
  INTEGER  :: SST_SAVE_N



  !
  !--Data Assimilation Parameters for Current Assimilation
  !
  LOGICAL  :: CUR_ASSIM                 !!TRUE IF CURRENT ASSIMILATION ACTIVE
  CHARACTER(LEN=80) CUR_ASSIM_FILE         !!FILE NAME FOR CURRENT DATA
  REAL(SP) :: CUR_RADIUS                       !!SEARCH RADIUS FOR INTERPOLATION POINTS 
  REAL(SP) :: CUR_WEIGHT_MAX
  REAL(SP) :: CUR_TIMESCALE
  REAL(SP) :: CUR_TIME_WINDOW
  INTEGER  :: CUR_MAX_LAYER                 !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT

  NAMELIST /NML_CUR_ASSIMILATION/      &
       & CUR_ASSIM,                    &
       & CUR_ASSIM_FILE,               &
       & CUR_RADIUS,                   &
       & CUR_WEIGHT_MAX,               &
       & CUR_TIMESCALE,                &
       & CUR_TIME_WINDOW,              &
       & CUR_MAX_LAYER 


  !
  !--Data Assimilation Parameters for Temp/Salinity Data Assimilation
  !
  LOGICAL  :: TS_ASSIM                  !!TRUE IF TEMP/SAL ASSIMILATION ACTIVE
  CHARACTER(LEN=80) TS_ASSIM_FILE         !!FILE NAME FOR TEMP/SAL DATA
  REAL(SP) :: TS_RADIUS                       !!SEARCH RADIUS FOR INTERPOLATION POINTS
  REAL(SP) :: TS_WEIGHT_MAX  
  REAL(SP) :: TS_TIMESCALE  
  REAL(SP) :: TS_TIME_WINDOW  
  INTEGER  :: TS_MAX_LAYER                  !!MAXIMUM NUMBER OF VERTICAL DATA FROM ANY OBS POINT

  NAMELIST /NML_TS_ASSIMILATION/       &
       & TS_ASSIM,                     &
       & TS_ASSIM_FILE,                &
       & TS_RADIUS,                    &
       & TS_WEIGHT_MAX,                &
       & TS_TIMESCALE,                 &
       & TS_TIME_WINDOW,               &
       & TS_MAX_LAYER


  !
  !--Current Assimilation Object Type
  !
  TYPE ASSIM_OBJ_CUR
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS 
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ODEPTH      !!OBSERVATION DEPTHS
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UO,VO       !!OBSERVATION DATA FOR X,Y VELOCITY COMPONENTS
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC 
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING                               
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
  END TYPE ASSIM_OBJ_CUR

  !
  !--Temperature/Salinity Assimilation Object Type
  !
  TYPE ASSIM_OBJ_TS  
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR CURRENT OBSERVATION WEIGHTING
     INTEGER  :: N_LAYERS                                 !!NUMBER OF OBSERVATIONS IN THE VERTICAL
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: DEPTH                                    !!DEPTH AT OBSERVATION STATION (MOORING)
     REAL(SP) :: SITA                                     !!LOCAL ISOBATH ANGLE AT OBSERVATION
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ODEPTH      !!OBSERVATION DEPTHS
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TEMP        !!OBSERVATION DATA FOR TEMPERATURE 
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SAL         !!OBSERVATION DATA FOR SALINITY 
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC
     INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: S_INT       !!SIGMA INTERVALS SURROUNDING CURRENT MEASUREMENT
     REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: S_WEIGHT    !!SIGMA WEIGHTING
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
  END TYPE ASSIM_OBJ_TS



  !
  !--SST Assimilation Object Type
  !
  TYPE ASSIM_OBJ_SST
     INTEGER  :: N_TIMES                                  !!NUMBER OF DATA TIMES
     INTEGER  :: N_INTPTS                                 !!NUMBER OF INTERPOLATION POINTS 
     INTEGER  :: N_T_WEIGHT                               !!DATA TIME FOR SST OBSERVATION WEIGHTING
     REAL(SP) :: X,Y                                      !!X AND Y COORDINATES OF OBSERVATION
     REAL(SP) :: T_WEIGHT                                 !!TIME WEIGHT
     REAL(SP) :: T_INT_AVGD                               !!INTERPOLATE SIM DATA AVERAGED OVER OBSERVATION PERIOD
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: T_INT_HOUR  !!INTERPOLATED SIM DATA OBTAINED AT EACH HOUR DURING SIM
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TIMES       !!DATA TIMES (DAYS)
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: SST         !!OBSERVATION DATA FOR SST 
     INTEGER,  ALLOCATABLE, DIMENSION(:)   :: INTPTS      !!POINTS USED TO INTERPOLATE TO OBSERVATION LOC 
     REAL(SP), ALLOCATABLE, DIMENSION(:)   :: X_WEIGHT    !!SPATIAL WEIGHTING FOR INTERPOLATION POINTS
  END TYPE ASSIM_OBJ_SST


  !
  !--Current Data Assimilation Variables
  !
  INTEGER                             N_ASSIM_CUR !!NUMBER OF CURRENT OBSERVATIONS 
  TYPE(ASSIM_OBJ_CUR), ALLOCATABLE :: CUR_OBS(:)  !!CURRENT ASSIMILATION DATA OBJECTS
  INTEGER, ALLOCATABLE             :: DA_CUR(:)   !!FLAG IF ELEMENT IS USED FOR CURRENT DA INTERP 

  !
  !--Salinity/Temperature Data Assimilation Variables
  !
  INTEGER                                N_ASSIM_TS   !!NUMBER OF TEMPERATURE OBSERVATIONS
  TYPE(ASSIM_OBJ_TS), ALLOCATABLE     :: TS_OBS(:)    !!TEMP ASSIMILATION DATA OBJECTS
  INTEGER, ALLOCATABLE                :: DA_TS(:)     !!FLAG IF NODE IS USED FOR CURRENT DA INTERP

  !
  !--SST Data Assimilation Variables
  !

   TYPE(NCFILE), POINTER :: SST_FILE

  INTEGER                             N_ASSIM_SST !!NUMBER OF SST  OBSERVATIONS 
  INTEGER                             N_TIMES_SST !!NUMBER OF SST OBSERVATIONS TIMES

  INTEGER    :: SST_SAVE_INDEX
  TYPE(TIME) :: SST_SAVE_TIME
  INTEGER(ITIME)  :: SST_SAVE_INTERVAL
  TYPE(NCVAR), POINTER :: VAR_SST
  TYPE(TIME)  :: ASSIM_RESTART_TIME 

  INTEGER                             ASSIM_FLAG  !!TRUE IF ON ASSIMILATION SWEEP   

  REAL(SP), ALLOCATABLE, TARGET, DIMENSION(:)  :: SST_OBS   !!SST OBS ON GRID
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)  :: SST_SAVED  !!SST ON EACH HOUR DURING SIM 
  REAL(SP), ALLOCATABLE, DIMENSION(:)    :: SST_AVG   !!SIM SST AVERAGED OVER OBSERVATION PERIOD 

  ! Some IINT VARIABLES TO KEEP TRACK OF PROGRESS
  INTEGER(ITIME) :: INT_COUNT, INT_START, INT_END


  !--save all_var to the *_BUF 

  TYPE(TIME)                              ::IntTime_BUF 
  TYPE(TIME)                              ::ExtTime_BUF 
  INTEGER(ITIME)                          ::IINT_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::U_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::V_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::WTS_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::S1_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::T1_BUF
!  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::RHO_BUF

  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::KM_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::KH_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::KQ_BUF

#  if defined (EQUI_TIDE)
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: EL_EQI_BUF
#  endif

#  if defined (ATMO_TIDE)
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: EL_ATMO_BUF
#  endif

#    if defined (GOTM)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: TKE_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: TEPS_BUF
# else
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::Q2_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::Q2L_BUF  
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::L_BUF
#    endif

  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::UA_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::VA_BUF

  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::EL_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::ET_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:)     ::H_BUF

 

#  if defined (DYE_RELEASE)  
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::DYE_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::DYEF_BUF
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)   ::DYEMEAN_BUF
#  endif 


CONTAINS !------------------------------------------------------------------!
  ! SET_ASSIM_PARAM     :   READ ASSIMILATION PARAMETERS FROM INPUT  !
  ! SET_ASSIM_INTERVALS :   SET UP DATA ASSIMILATION SWEEPS          !
  ! SET_CUR_ASSIM_DATA  :   READ AND SET CURRENT ASSIMILATION DATA   ! 
  ! SET_TS_ASSIM_DATA   :   READ AND SET TEMP/SAL ASSIMILATION DATA  ! 
  ! SET_SST_ASSIM_DATA  :   READ AND SET SST ASSIMILATION DATA       ! 
  ! CURRENT_NUDGING     :   NUDGE CURRENT USING ASSIMILATION         !
  ! TEMP NUDGING        :   NUDGE TEMP USING ASSIMILATION            !
  ! SALT_NUDGING        :   NUDGE SALT USING ASSIMILATION            !
  ! SST_NUDGING         :   NUDGE SST USING ASSIMILATION             !
  ! SST_INT             :   INTERPOLATE HOURLY/AVGE SST DATA TO OBS  !
  ! HOT_START_SST       :   READ ICs FOR SST_ASSIMILATION STAGE      !
  ! ARC_SST             :   DUMP ICs FOR SST_ASSIMILATION STAGE      ! 
  ! -----------------------------------------------------------------!
  ! -----------------------------------------------------------------!
  ! SET_SST_ASSIM_DATA_GRID:   READ AND SET SST GRID ASSIM DATA      ! 
  ! SST_NUDGING_GRID       :   NUDGE SST USING ASSIMILATION          !
  ! SST_INT_GRID           :   SAVE HOURLY/AVGE SST USING SIMULATION !  
  ! -----------------------------------------------------------------!

  !==============================================================================|
  !==============================================================================|

  SUBROUTINE INITIALIZE_ASSIM_PARAM
    USE CONTROL
    IMPLICIT NONE

    !    NML_SST_ASSIMILIATION
    SST_ASSIM          = .FALSE.
    SST_ASSIM_FILE     = trim(CASENAME)//"_sst.nc"
    SST_RADIUS         = 0.0_SP
    SST_WEIGHT_MAX     = 0.0_SP
    SST_TIMESCALE      = 0.0_SP
    SST_TIME_WINDOW    = 0.0_SP
    SST_N_PER_INTERVAL = 0

    !    NML_SST_ASSIMILIATION
    SSTGRD_ASSIM          = .FALSE.
    SSTGRD_ASSIM_FILE     = trim(CASENAME)//"_sstgrd.nc"
    SSTGRD_WEIGHT_MAX     = 0.0_SP
    SSTGRD_TIMESCALE      = 0.0_SP
    SSTGRD_TIME_WINDOW    = 0.0_SP
    SSTGRD_N_PER_INTERVAL = 0


    !    NML_CUR_ASSIMILIATION
    CUR_ASSIM       = .FALSE.
    CUR_ASSIM_FILE  = trim(CASENAME)//"_cur.nc"
    CUR_RADIUS      = 0.0_SP
    CUR_WEIGHT_MAX  = 0.0_SP
    CUR_TIMESCALE   = 0.0_SP
    CUR_TIME_WINDOW = 0.0_SP
    CUR_MAX_LAYER   = 0


    !    NML_TS_ASSIMILIATION
    TS_ASSIM          = .FALSE.
    TS_ASSIM_FILE  = trim(CASENAME)//"_ts.nc"
    TS_RADIUS         = 0.0_SP
    TS_WEIGHT_MAX     = 0.0_SP
    TS_TIMESCALE      = 0.0_SP
    TS_TIME_WINDOW    = 0.0_SP
    TS_MAX_LAYER      = 0


  END SUBROUTINE INITIALIZE_ASSIM_PARAM

  SUBROUTINE SET_ASSIM_PARAM 

    !------------------------------------------------------------------------------|
    !  READ IN PARAMETERS CONTROLLING ASSIMILATION                                 |
    !------------------------------------------------------------------------------|

    USE MOD_PREC
    USE ALL_VARS

    IMPLICIT NONE
    CHARACTER(LEN=120) :: FNAME, pathnfile
    INTEGER :: IOS, charnum
    LOGICAL FEXIST


    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: SET_ASSIM_PARAM"

    charnum = index (DATA_ASSIMILATION_FILE,".nml")
    if (charnum /= len_trim(DATA_ASSIMILATION_FILE)-3)&
         & CALL WARNING("DATA ASSIMILATION FILE does not end in .nml", &
         & trim(DATA_ASSIMILATION_FILE))
    ! OPEN FILE - try both with appending input dir and without!
    pathnfile = trim(INPUT_DIR)//trim(DATA_ASSIMILATION_FILE)
    INQUIRE(FILE=PATHNFILE,EXIST=FEXIST)
    IF(FEXIST) THEN
       Call FOPEN(ASSIMUNIT,trim(pathnfile),'cfr')
    ELSE
       pathnfile = trim(DATA_ASSIMILATION_FILE)
       Call FOPEN(ASSIMUNIT,trim(pathnfile),'cfr')
    END IF


    CALL INITIALIZE_ASSIM_PARAM

    ! SAVE FILE NAME USED TO FNAME FOR ERROR CHECKING
    FNAME =  pathnfile

    ! Read SST ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_SST_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) &
            & write(IPT,*)"Can Not Read NameList NML_SST_ASSIMILATION from file: "//trim(FNAME)
    end if
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SST_ASSIMILATION)


    ! Read SST GRID ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_SSTGRD_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) &
            & write(IPT,*)"Can Not Read NameList NML_SSTGRD_ASSIMILATION from file: "//trim(FNAME)
    end if
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_SSTGRD_ASSIMILATION)


    ! Read CUR ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_CUR_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) &
            & write(IPT,*)"Can Not Read NameList NML_CUR_ASSIMILATION from file: "//trim(FNAME)
    end if
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_CUR_ASSIMILATION)


    ! Read SST ASSIMILATION Settings
    READ(UNIT=ASSIMUNIT, NML=NML_TS_ASSIMILATION,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) &
            & write(IPT,*)"Can Not Read NameList NML_TS_ASSIMILATION from file: "//trim(FNAME)
    end if
    REWIND(ASSIMUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_TS_ASSIMILATION)


    IF (DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*)''
       WRITE(IPT,*)'!        DATA ASSIMILATION PARAMETERS       '
       IF(CUR_ASSIM)THEN
          WRITE(IPT,*)'!  # CUR_ASSIM       :  ACTIVE'
          WRITE(IPT,*)'!  # CUR_ASSIM_FILE      : '//TRIM(CUR_ASSIM_FILE)
          WRITE(IPT,*)'!  # CUR_RADIUS          :',CUR_RADIUS
          WRITE(IPT,*)'!  # CUR_WEIGHT_MAX      :',CUR_WEIGHT_MAX
          WRITE(IPT,*)'!  # CUR_TIMESCALE       :',CUR_TIMESCALE   
          WRITE(IPT,*)'!  # CUR_TIME_WINDOW     :',CUR_TIME_WINDOW
          WRITE(IPT,*)'!  # CUR_MAX_LAYER       :',CUR_MAX_LAYER
       ELSE
          WRITE(IPT,*)'!  # CUR_ASSIM           :  NOT ACTIVE'
       END IF
       IF(SST_ASSIM)THEN
          WRITE(IPT,*)'!  # SST_ASSIM           :  ACTIVE'
          WRITE(IPT,*)'!  # SST_RADIUS          :',SST_RADIUS
          WRITE(IPT,*)'!  # SST_ASSIM_FILE      : '//TRIM(SST_ASSIM_FILE)
          WRITE(IPT,*)'!  # SST_TIMESCALE       :',SST_TIMESCALE
          WRITE(IPT,*)'!  # SST_WEIGHT_MAX      :',SST_WEIGHT_MAX   
          WRITE(IPT,*)'!  # SST_TIME_WINDOW     :',SST_TIME_WINDOW
          WRITE(IPT,*)'!  # SST_N_PER_INTERVAL  :',SST_N_PER_INTERVAL
       ELSE
          WRITE(IPT,*)'!  # SST_ASSIM           :  NOT ACTIVE'
       END IF
       IF(SSTGRD_ASSIM)THEN
          WRITE(IPT,*)'!  # SSTGRD_ASSIM           :  ACTIVE'
          WRITE(IPT,*)'!  # SSTGRD_ASSIM_FILE      : '//TRIM(SSTGRD_ASSIM_FILE)
          WRITE(IPT,*)'!  # SSTGRD_TIMESCALE       :',SSTGRD_TIMESCALE
          WRITE(IPT,*)'!  # SSTGRD_WEIGHT_MAX      :',SSTGRD_WEIGHT_MAX   
          WRITE(IPT,*)'!  # SSTGRD_TIME_WINDOW     :',SSTGRD_TIME_WINDOW
          WRITE(IPT,*)'!  # SSTGRD_N_PER_INTERVAL  :',SSTGRD_N_PER_INTERVAL
       ELSE
          WRITE(IPT,*)'!  # SSTGRD_ASSIM           :  NOT ACTIVE'
       END IF

       IF(TS_ASSIM)THEN
          WRITE(IPT,*)'!  # TS_ASSIM       :  ACTIVE'
          WRITE(IPT,*)'!  # TS_ASSIM_FILE      : '//TRIM(TS_ASSIM_FILE)
          WRITE(IPT,*)'!  # TS_RADIUS          :',TS_RADIUS
          WRITE(IPT,*)'!  # TS_WEIGHT_MAX      :',TS_WEIGHT_MAX
          WRITE(IPT,*)'!  # TS_TIMESCALE       :',TS_TIMESCALE   
          WRITE(IPT,*)'!  # TS_TIME_WINDOW     :',TS_TIME_WINDOW
          WRITE(IPT,*)'!  # TS_MAX_LAYER       :',TS_MAX_LAYER
       ELSE
          WRITE(IPT,*)'!  # TS_ASSIM           :  NOT ACTIVE'
       END IF
    END IF
    
    IF(SST_ASSIM .AND. SSTGRD_ASSIM) CALL FATAL_ERROR &
         & ('Using two kinds of sst assimilation does not make any sense!')

    ! SET THE FVCOM RUN MODE BASED ON WHAT TYPE OF ASSIMILATION WE
    ! ARE DOING
    IF(SST_ASSIM .OR. SSTGRD_ASSIM) THEN
       FVCOM_RUN_MODE = FVCOM_NUDGE_AVG_SST
    END IF
       
    IF(SST_ASSIM) THEN
       SST_WM = SST_WEIGHT_MAX
       SST_TSCALE = SST_TIMESCALE
       SST_TWINDOW= SST_TIME_WINDOW
       SST_SAVE_N = SST_N_PER_INTERVAL
    ELSE IF(SSTGRD_ASSIM) THEN
       SST_WM = SSTGRD_WEIGHT_MAX
       SST_TSCALE = SSTGRD_TIMESCALE
       SST_TWINDOW= SSTGRD_TIME_WINDOW
       SST_SAVE_N = SSTGRD_N_PER_INTERVAL
    END IF

    IF (SST_TSCALE * DTI >= 1.0_sp) CALL FATAL_ERROR &
         ('The efolding time scale for the sst nudging is to short for the model time step!') 


    IF (DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: SET_ASSIM_PARAM"

    RETURN
  END SUBROUTINE SET_ASSIM_PARAM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE LOAD_ASSIM_DATA
!==============================================================================!
    IMPLICIT NONE

    IF(SST_ASSIM) CALL LOAD_SST_ASSIM_DATA

    IF(SSTGRD_ASSIM) CALL LOAD_SSTGRD_ASSIM_DATA

    IF(TS_ASSIM) CALL LOAD_TS_ASSIM_DATA

    IF(CUR_ASSIM) CALL LOAD_CUR_ASSIM_DATA

    
  END SUBROUTINE LOAD_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE LOAD_SST_ASSIM_DATA
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios
    
    back = .true.
    
    ! TEST FILE NAME
    charnum = index (SST_ASSIM_FILE,".nc",back)
    if (charnum /= len_trim(SST_ASSIM_FILE)-2)&
         & CALL WARNING("SST_ASSIM_FILE name does not end in .nc", &
         & trim(SST_ASSIM_FILE))
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(SST_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
    
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if
    
  END SUBROUTINE LOAD_SST_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE LOAD_SSTGRD_ASSIM_DATA
!==============================================================================!
    IMPLICIT NONE
    
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios
    
    back = .true.
    
    ! TEST FILE NAME
    charnum = index (SSTGRD_ASSIM_FILE,".nc",back)
    if (charnum /= len_trim(SSTGRD_ASSIM_FILE)-2)&
         & CALL WARNING("SSTGRD_ASSIM_FILE name does not end in .nc", &
         & trim(SSTGRD_ASSIM_FILE))
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(SSTGRD_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
    
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if




  END SUBROUTINE LOAD_SSTGRD_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE LOAD_CUR_ASSIM_DATA
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios
    
    back = .true.
    
    ! TEST FILE NAME
    charnum = index (CUR_ASSIM_FILE,".nc",back)
    if (charnum /= len_trim(CUR_ASSIM_FILE)-2)&
         & CALL WARNING("CUR_ASSIM_FILE name does not end in .nc", &
         & trim(CUR_ASSIM_FILE))
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(CUR_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
    
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if


  END SUBROUTINE LOAD_CUR_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================!
  SUBROUTINE LOAD_TS_ASSIM_DATA
!==============================================================================!
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    integer :: ncfileind, datfileind,ios,charnum, i
    logical :: fexist,back,connected
    character(len=100) :: testchar
    character(len=160) :: pathnfile
    character(len=2) :: cios
    
    back = .true.
    
    ! TEST FILE NAME
    charnum = index (TS_ASSIM_FILE,".nc",back)
    if (charnum /= len_trim(TS_ASSIM_FILE)-2)&
         & CALL WARNING("TS_ASSIM_FILE name does not end in .nc", &
         & trim(TS_ASSIM_FILE))
    
    ! INITIALIZE TYPE TO HOLD FILE METADATA
    pathnfile= trim(INPUT_DIR)//trim(TS_ASSIM_FILE)

    CALL NC_INIT(NCF,pathnfile)
    
    ! OPEN THE FILE AND LOAD METADATA       
    If(.not. NCF%OPEN) then
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       FILEHEAD => ADD(FILEHEAD,NCF)
    end if

  END SUBROUTINE LOAD_TS_ASSIM_DATA

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!!$!==============================================================================!
!!$  SUBROUTINE SET_SST_ASSIM_INTERVALS
!!$!==============================================================================!
!!$    
!!$    !------------------------------------------------------------------------------!
!!$    !  SET UP ASSIMILATION SWEEP INTERVALS                                         |
!!$    !------------------------------------------------------------------------------!
!!$    USE ALL_VARS
!!$    IMPLICIT NONE
!!$    REAL(DP) NDAYS
!!$
!!$    ISTART_DAY = 1
!!$    IEND_DAY   = 1
!!$    ISWEEP     = 1
!!$
!!$    CALL TIME2DAYS(EndTime-StartTime,NDAYS)
!!$    
!!$    write(ipt,*) ndays
!!$
!!$    call pshutdown
!!$
    
!!$    ISTART_DAY=INT(FLOAT(ISTART-1)*DTI/86400.+0.5)
!!$    IEND_DAY  =INT(FLOAT(IEND)*DTI/86400.-1.+0.5)
!!$    PURE_SIM  = .FALSE.
!!$    ISWEEP    = 2
!!$    IF(MSR)THEN
!!$       WRITE(IPT,*)'!'
!!$       WRITE(IPT,*)'!  STARTING DAY :  ',ISTART_DAY
!!$       WRITE(IPT,*)'!  ENDING DAY   :  ',IEND_DAY
!!$       WRITE(IPT,*)'!  SST DA INT   :  ',IAV_DAY
!!$       WRITE(IPT,*)'!  #SST DA INT  :  ',(IEND_DAY-ISTART_DAY-1)/IAV_DAY
!!$    END IF
!!$    
!!$  END SUBROUTINE SET_SST_ASSIM_INTERVALS


!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$!==============================================================================!
!!$   SUBROUTINE SET_CUR_ASSIM_DATA 
!!$!==============================================================================!
!!$
!!$!------------------------------------------------------------------------------!
!!$!  SET UP ASSIMILATION DATA FOR CURRENT OBSERVATIONS                           |
!!$!------------------------------------------------------------------------------!
!!$   USE ALL_VARS
!!$   IMPLICIT NONE
!!$   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS,NLAY
!!$   CHARACTER(LEN=120) :: FNAME,ONAME
!!$   CHARACTER(LEN= 2 ) :: NAC   
!!$   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
!!$   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
!!$   REAL(SP), PARAMETER :: ALF = 0.05_SP
!!$   LOGICAL :: FEXIST
!!$   INTEGER :: MAXEL,NBD_CNT
!!$   real(sp) :: lmin 
!!$   integer jmin 
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Read Number of Current Observations and Coordinates of Each                 !
!!$!------------------------------------------------------------------------------!
!!$       
!!$   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_current.xy"
!!$!
!!$!--Make Sure Current Assimilation Data File Exists-----------------------------!
!!$!
!!$   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!!$   IF(MSR .AND. .NOT.FEXIST)THEN
!!$     WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!!$     WRITE(IPT,*)'HALTING.....'
!!$     CALL PSTOP
!!$   END IF
!!$     
!!$!
!!$!--Read Number of Current Measurement Stations---------------------------------!
!!$!
!!$   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
!!$   READ(1,*) N_ASSIM_CUR
!!$   ALLOCATE(CUR_OBS(N_ASSIM_CUR))
!!$
!!$!
!!$!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!!$!
!!$
!!$   DO I=1,N_ASSIM_CUR
!!$     READ(1,*)ITMP,CUR_OBS(I)%X,CUR_OBS(I)%Y,CUR_OBS(I)%DEPTH,NLAY,CUR_OBS(I)%SITA
!!$     CUR_OBS(I)%N_LAYERS = NLAY
!!$     ALLOCATE(CUR_OBS(I)%ODEPTH(NLAY))
!!$     DO J=1,NLAY
!!$       READ(1,*)CUR_OBS(I)%ODEPTH(J)
!!$       IF(CUR_OBS(I)%ODEPTH(J) > CUR_OBS(I)%DEPTH)THEN
!!$         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF CURRENT MOORING',I
!!$         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH'
!!$         IF(MSR)WRITE(IPT,*)'HALTING...........'
!!$         CALL PSTOP
!!$       END IF
!!$     END DO
!!$   END DO
!!$   MAX_LAYER_CUR = MAXVAL(CUR_OBS(1:N_ASSIM_CUR)%N_LAYERS)
!!$
!!$!
!!$!--Shift Coordinates-----------------------------------------------------------!
!!$!
!!$   CUR_OBS(:)%X = CUR_OBS(:)%X - VXMIN 
!!$   CUR_OBS(:)%Y = CUR_OBS(:)%Y - VYMIN 
!!$   
!!$!
!!$!--Close Current Observation Global File---------------------------------------!
!!$!
!!$   CLOSE(1)
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Open Current Observation Files for Each Observation Point and Read Data     !
!!$!------------------------------------------------------------------------------!
!!$
!!$   DO I=1,N_ASSIM_CUR
!!$
!!$!----Make Sure Current Observation File Exists---------------------------------!
!!$     WRITE(NAC,'(I2.2)')I
!!$     ONAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//'_curr'//NAC//'.dat' 
!!$     INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
!!$     IF(MSR .AND. .NOT.FEXIST)THEN
!!$       WRITE(IPT,*)'CURRENT OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
!!$       WRITE(IPT,*)'HALTING.....'
!!$       CALL PSTOP
!!$     END IF
!!$
!!$!----Open Current Observation File for Read------------------------------------!
!!$     OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)
!!$     NCNT = 0
!!$
!!$!----Count Number of Data Entries for Observation I----------------------------!
!!$     DO WHILE(.TRUE.)
!!$       READ(1,*,IOSTAT=IOS)
!!$       IF(IOS < 0)EXIT
!!$       NCNT = NCNT + 1
!!$     END DO
!!$     CUR_OBS(I)%N_TIMES = NCNT
!!$     REWIND(1)
!!$     IF(NCNT == 0)THEN
!!$       IF(MSR)WRITE(IPT,*)'NO DATA FOR CURRENT OBSERVATION',I
!!$       CALL PSTOP
!!$     END IF
!!$
!!$!----Allocate Arrays to Hold Current (UA,VA) and Time (TIME)-------------------!
!!$     ALLOCATE(CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES))
!!$     ALLOCATE(CUR_OBS(I)%UO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
!!$     ALLOCATE(CUR_OBS(I)%VO( CUR_OBS(I)%N_TIMES , CUR_OBS(I)%N_LAYERS ))
!!$
!!$!----Read in Current Data for Observation I------------------------------------!
!!$     NLAY = CUR_OBS(I)%N_LAYERS
!!$     DO J=1,CUR_OBS(I)%N_TIMES
!!$       READ(1,*)CUR_OBS(I)%TIMES(J),(CUR_OBS(I)%UO(J,K),CUR_OBS(I)%VO(J,K),K=1,NLAY)
!!$     END DO
!!$     CLOSE(1)
!!$
!!$!----Convert Time to Seconds---------------------------------------------------!
!!$!----Shift Jan 1 Based Time Data to Dec 1 Based Time Data-----CASESPECIFIC-----!
!!$     IF(trim(CASENAME) == 'gom')THEN
!!$!       CUR_OBS(I)%TIMES = ((CUR_OBS(I)%TIMES-1.0_SP)*24.0_SP+744.0_SP)*3600.0_SP
!!$!       CUR_OBS(I)%TIMES = (CUR_OBS(I)%TIMES-1.0_SP+61.0_SP)*24.0_SP*3600.0_SP !for 1995
!!$       CUR_OBS(I)%TIMES = (CUR_OBS(I)%TIMES-1.0_SP)*24.0_SP*3600.0_SP !after 1996 
!!$     ELSE   
!!$       CUR_OBS(I)%TIMES = CUR_OBS(I)%TIMES*3600.0_SP*24.0_SP
!!$     END IF
!!$ 
!!$!----Convert Current Data from cm/s to m/s-------------------------------------!
!!$     CUR_OBS(I)%UO = CUR_OBS(I)%UO * .01_SP
!!$     CUR_OBS(I)%VO = CUR_OBS(I)%VO * .01_SP
!!$  END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Count Number of Points with Bad Data (UO = 0. + V0 = 0.)         
!!$!------------------------------------------------------------------------------!
!!$  NBD_CNT = 0
!!$  DO I=1,N_ASSIM_CUR
!!$  DO J=1,CUR_OBS(I)%N_TIMES
!!$  DO K=1,CUR_OBS(I)%N_LAYERS
!!$  IF(ABS(CUR_OBS(I)%UO(J,K)) + ABS(CUR_OBS(I)%VO(J,K)) < .0001) THEN
!!$    NBD_CNT = NBD_CNT + 1
!!$  END IF
!!$  END DO
!!$  END DO
!!$  END DO
!!$ 
!!$!------------------------------------------------------------------------------!
!!$!  Compute Spatial Interpolation Weights for each Mooring Location 
!!$!------------------------------------------------------------------------------!
!!$!   DO I=1,N_ASSIM_CUR
!!$!   LMIN = 100000000.
!!$!     X0 = CUR_OBS(I)%X
!!$!     Y0 = CUR_OBS(I)%Y
!!$!     DO J=1,MGL
!!$!       DX = ABS(XG(J)-X0)
!!$!       DY = ABS(YG(J)-Y0)
!!$!       IF(SQRT(DX**2 + DY**2) < LMIN)THEN
!!$!         LMIN = SQRT(DX**2 + DY**2)
!!$!         JMIN = J
!!$!       END IF
!!$!      END DO
!!$!      CUR_OBS(I)%SITA = SITA_GD(JMIN) + 3.14159_SP/2.0_SP
!!$!    END DO
!!$
!!$
!!$   ALLOCATE(ITEMP(NGL),FTEMP(NGL),DA_CUR(NGL))     ; DA_CUR = 0
!!$   DO I=1,N_ASSIM_CUR
!!$     X0 = CUR_OBS(I)%X
!!$     Y0 = CUR_OBS(I)%Y
!!$     ISOBATH_ANGLE = CUR_OBS(I)%SITA/180.0_SP*3.1415926_SP
!!$     ECNT = 0
!!$     DO J=1,NGL
!!$       DX = ABS(XCG(J)-X0)
!!$       DY = ABS(YCG(J)-Y0)
!!$       RD = SQRT(DX**2 + DY**2)
!!$       IF(RD <= RAD_CUR)THEN
!!$         DA_CUR(J)   = 1
!!$         ECNT        = ECNT + 1      
!!$         ITEMP(ECNT) =  J
!!$         FTEMP(ECNT) = (RAD_CUR**2 - RD**2) / (RAD_CUR**2 + RD**2)
!!$         ANG_OBS_SIM = ATAN2(DY,DX)
!!$         D_ANGLE     = ANG_OBS_SIM - ISOBATH_ANGLE 
!!$         D_ANGLE     = D_ANGLE - INT(D_ANGLE/3.1415926_SP)*3.1415926_SP
!!$         D_ANGLE     = ABS(D_ANGLE)
!!$         DIR_WEIGHT  = (ABS(D_ANGLE-0.5*3.1415926_SP)+ALF*3.1415926_SP)/ &
!!$                       ((0.5_SP+ALF)*3.1415926_SP)
!!$         FTEMP(ECNT) = FTEMP(ECNT)*DIR_WEIGHT
!!$       END IF
!!$     END DO
!!$     IF(ECNT == 0)THEN
!!$       WRITE(IPT,*)'ERROR SETTING UP CURRENT DATA ASSIMILATION'
!!$       WRITE(IPT,*)'NO ELEMENTS LIE WITHIN RADIUS',RAD_CUR
!!$       WRITE(IPT,*)'OF OBSERVATION POINT',I
!!$       CALL PSTOP   
!!$     ELSE
!!$       CUR_OBS(I)%N_INTPTS = ECNT
!!$       ALLOCATE(CUR_OBS(I)%INTPTS(ECNT))
!!$       ALLOCATE(CUR_OBS(I)%X_WEIGHT(ECNT))
!!$       CUR_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
!!$       CUR_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
!!$     END IF
!!$   END DO
!!$   DEALLOCATE(FTEMP,ITEMP)
!!$
!!$
!!$     
!!$!------------------------------------------------------------------------------!
!!$!  Compute Sigma Layer Weights for Vertical Interpolation                                                                 
!!$!------------------------------------------------------------------------------!
!!$   DO I=1,N_ASSIM_CUR
!!$     NLAY = CUR_OBS(I)%N_LAYERS
!!$     ALLOCATE(CUR_OBS(I)%S_INT(NLAY,2))
!!$     ALLOCATE(CUR_OBS(I)%S_WEIGHT(NLAY,2))
!!$     DO J=1,NLAY
!!$       SIGMA_C = -CUR_OBS(I)%ODEPTH(J)/CUR_OBS(I)%DEPTH
!!$       DO K=2,KBM1
!!$         IF(ZZ(K) <= SIGMA_C .AND. ZZ(K-1) > SIGMA_C)THEN 
!!$           CUR_OBS(I)%S_INT(J,1) = K-1
!!$           CUR_OBS(I)%S_INT(J,2) = K
!!$           CUR_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ(K))/(ZZ(K-1)-ZZ(K))
!!$           CUR_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - CUR_OBS(I)%S_WEIGHT(J,1) 
!!$         END IF  
!!$       END DO
!!$       IF(ZZ(1) < SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
!!$         CUR_OBS(I)%S_INT(J,1) = 1
!!$         CUR_OBS(I)%S_INT(J,2) = 1
!!$         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
!!$         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
!!$       END IF
!!$       IF(ZZ(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
!!$         CUR_OBS(I)%S_INT(J,1) = KBM1
!!$         CUR_OBS(I)%S_INT(J,2) = KBM1
!!$         CUR_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
!!$         CUR_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
!!$       END IF
!!$
!!$     END DO
!!$   END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Report Number of Interpolation Points, Location and Number of Data 
!!$!------------------------------------------------------------------------------!
!!$   IF(.NOT. MSR)RETURN
!!$
!!$   WRITE(IPT,*)
!!$   WRITE(IPT,*)'!            CURRENT OBSERVATION DATA           '
!!$   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_EL   SITA"
!!$   DO I=1,N_ASSIM_CUR
!!$     MAXEL = MAXLOC(CUR_OBS(I)%X_WEIGHT,DIM=1)
!!$     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
!!$     I,CUR_OBS(I)%X/1000.,CUR_OBS(I)%Y/1000., &
!!$       CUR_OBS(I)%N_INTPTS,CUR_OBS(I)%N_TIMES,CUR_OBS(I)%INTPTS(MAXEL),&
!!$       CUR_OBS(I)%SITA
!!$   END DO
!!$   WRITE(IPT,*)
!!$   WRITE(IPT,*)'NUMBER OF BAD CURRENT DATA POINTS: ',NBD_CNT
!!$   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
!!$   DO I=1,N_ASSIM_CUR
!!$   WRITE(IPT,*)I,CUR_OBS(I)%TIMES(1)/(24.*3600.),&
!!$       CUR_OBS(I)%TIMES(CUR_OBS(I)%N_TIMES)/(24.*3600.)
!!$   END DO
!!$
!!$   RETURN
!!$   END SUBROUTINE SET_CUR_ASSIM_DATA
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$   SUBROUTINE CURRENT_NUDGING
!!$!==============================================================================|
!!$!  USE CURRENT OBSERVATION DATA TO ADJUST VELOCITY COMPONENTS                  |
!!$!==============================================================================|
!!$   USE MOD_PREC
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR 
!!$#  endif
!!$   IMPLICIT NONE
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UINT,VINT,UCORR,VCORR,UG,VG,TWGHT
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
!!$   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
!!$   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
!!$   INTEGER I,J,K,J1,K1,K2,NLAY,ITIME,NTIME,IERR
!!$   INTRINSIC MINLOC
!!$!==============================================================================|
!!$
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Gather U and V Fields to Master Processor                                   ! 
!!$!------------------------------------------------------------------------------!
!!$   ALLOCATE(UG(NGL,KB))
!!$   ALLOCATE(VG(NGL,KB))
!!$#  if defined (MULTIPROCESSOR)
!!$   IF(PAR)THEN
!!$     CALL GATHER(LBOUND(UF,1),  UBOUND(UF,1),  N,NGL,KB,MYID,NPROCS,EMAP,UF,UG)
!!$     CALL GATHER(LBOUND(VF,1),  UBOUND(VF,1),  N,NGL,KB,MYID,NPROCS,EMAP,VF,VG)
!!$   END IF
!!$#  endif
!!$   IF(SERIAL)THEN
!!$     UG(1:NGL,1:KBM1) = UF(1:NGL,1:KBM1)
!!$     VG(1:NGL,1:KBM1) = VF(1:NGL,1:KBM1)
!!$   END IF
!!$!------------------------------------------------------------------------------!
!!$!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   IF(MSR)THEN
!!$   CUR_OBS%T_WEIGHT = 0. 
!!$   T_THRESH         = ASTIME_WINDOW_CUR   
!!$   SIMTIME          = TIME*86400
!!$
!!$   DO I=1,N_ASSIM_CUR       
!!$     NTIME = CUR_OBS(I)%N_TIMES
!!$     ALLOCATE(FTEMP(NTIME)) 
!!$     FTEMP(1:NTIME) = ABS(SIMTIME - CUR_OBS(I)%TIMES(1:NTIME))
!!$     DT_MIN = MINVAL(FTEMP(1:NTIME))
!!$     CUR_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)
!!$
!!$     IF(DT_MIN < T_THRESH)THEN     
!!$       IF(DT_MIN < .5_SP*T_THRESH) THEN
!!$         CUR_OBS(I)%T_WEIGHT = 1.0_SP
!!$       ELSE
!!$         CUR_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
!!$       END IF
!!$     END IF
!!$
!!$     DEALLOCATE(FTEMP)
!!$   END DO
!!$   
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Interpolate Simulation Data to Local Observation Point                      ! 
!!$!------------------------------------------------------------------------------!
!!$       
!!$   ALLOCATE(UINT(N_ASSIM_CUR,MAX_LAYER_CUR)) ; UINT = 0. 
!!$   ALLOCATE(VINT(N_ASSIM_CUR,MAX_LAYER_CUR)) ; VINT = 0.
!!$
!!$   DO I=1,N_ASSIM_CUR  
!!$     DO J=1,CUR_OBS(I)%N_INTPTS
!!$       J1        = CUR_OBS(I)%INTPTS(J)
!!$       WGHT      = CUR_OBS(I)%X_WEIGHT(J)
!!$       NLAY      = CUR_OBS(I)%N_LAYERS
!!$       DO K=1,NLAY
!!$         U1 = UG(J1,CUR_OBS(I)%S_INT(K,1))
!!$         U2 = UG(J1,CUR_OBS(I)%S_INT(K,2))
!!$         V1 = VG(J1,CUR_OBS(I)%S_INT(K,1))
!!$         V2 = VG(J1,CUR_OBS(I)%S_INT(K,2))
!!$         W1 = CUR_OBS(I)%S_WEIGHT(K,1)
!!$         W2 = CUR_OBS(I)%S_WEIGHT(K,2)
!!$         UINT(I,K) = UINT(I,K) + (U1*W1 + U2*W2)*WGHT 
!!$         VINT(I,K) = VINT(I,K) + (V1*W1 + V2*W2)*WGHT 
!!$       END DO
!!$     END DO
!!$     TOT_WGHT = SUM(CUR_OBS(I)%X_WEIGHT(1:CUR_OBS(I)%N_INTPTS))
!!$     UINT(I,1:NLAY) = UINT(I,1:NLAY)/TOT_WGHT
!!$     VINT(I,1:NLAY) = VINT(I,1:NLAY)/TOT_WGHT
!!$   END DO
!!$     
!!$!------------------------------------------------------------------------------!
!!$!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   ALLOCATE(TWGHT(NGL,KBM1))   ; TWGHT = 0.
!!$   ALLOCATE(UCORR(NGL,KBM1))   ; UCORR   = 0.
!!$   ALLOCATE(VCORR(NGL,KBM1))   ; VCORR   = 0.
!!$
!!$   DO I=1,N_ASSIM_CUR
!!$     DO J=1,CUR_OBS(I)%N_INTPTS
!!$       J1     = CUR_OBS(I)%INTPTS(J)
!!$       ITIME  = CUR_OBS(I)%N_T_WEIGHT
!!$       NLAY   = CUR_OBS(I)%N_LAYERS
!!$       DO K=1,NLAY
!!$         K1           = CUR_OBS(I)%S_INT(K,1)
!!$         K2           = CUR_OBS(I)%S_INT(K,2)
!!$         W1           = CUR_OBS(I)%S_WEIGHT(K,1)
!!$         W2           = CUR_OBS(I)%S_WEIGHT(K,2)
!!$         WEIGHT1      = CUR_OBS(I)%T_WEIGHT*CUR_OBS(I)%X_WEIGHT(J)*W1
!!$         WEIGHT2      = CUR_OBS(I)%T_WEIGHT*CUR_OBS(I)%X_WEIGHT(J)*W2
!!$         TWGHT(J1,K1) = TWGHT(J1,K1) + WEIGHT1   
!!$         TWGHT(J1,K2) = TWGHT(J1,K2) + WEIGHT2   
!!$         DEFECT       = CUR_OBS(I)%UO(ITIME,K) - UINT(I,K)
!!$         CORRECTION   = GAMA_CUR*DEFECT
!!$!qxu{
!!$!         UCORR(J1,K1) = UCORR(J1,K1) + CORRECTION*WEIGHT1
!!$!         UCORR(J1,K2) = UCORR(J1,K2) + CORRECTION*WEIGHT2
!!$         UCORR(J1,K1) = UCORR(J1,K1) + CORRECTION*WEIGHT1**2
!!$         UCORR(J1,K2) = UCORR(J1,K2) + CORRECTION*WEIGHT2**2
!!$!qxu}
!!$         DEFECT       = CUR_OBS(I)%VO(ITIME,K) - VINT(I,K)
!!$         CORRECTION   = GAMA_CUR*DEFECT
!!$!QXU{
!!$!         VCORR(J1,K1) = VCORR(J1,K1) + CORRECTION*WEIGHT1
!!$!         VCORR(J1,K2) = VCORR(J1,K2) + CORRECTION*WEIGHT2
!!$
!!$         VCORR(J1,K1) = VCORR(J1,K1) + CORRECTION*WEIGHT1**2
!!$         VCORR(J1,K2) = VCORR(J1,K2) + CORRECTION*WEIGHT2**2
!!$!QXU}
!!$!        GEOFF NEW
!!$         IF(ABS(CUR_OBS(I)%UO(ITIME,K)) + ABS(CUR_OBS(I)%VO(ITIME,K)) < .0001)THEN
!!$           TWGHT(J1,K1) = 0.
!!$           TWGHT(J1,K2) = 0.
!!$         END IF
!!$       END DO
!!$     END DO
!!$   END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Nudge Simulation Data Using Local Corrections                               ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   DO I=1,NGL
!!$     DO K=1,KBM1
!!$       IF(DA_CUR(I) == 1 .AND. TWGHT(I,K) > 1.0E-08)THEN
!!$         UG(I,K) = UG(I,K) + DTI*GALPHA_CUR*UCORR(I,K)/TWGHT(I,K)
!!$         VG(I,K) = VG(I,K) + DTI*GALPHA_CUR*VCORR(I,K)/TWGHT(I,K)
!!$       END IF
!!$     END DO
!!$   END DO
!!$
!!$   DEALLOCATE(TWGHT,UCORR,VCORR,UINT,VINT)
!!$   END IF  !!MASTER
!!$!------------------------------------------------------------------------------!
!!$!  Disperse New Data Fields to Slave Processors                                ! 
!!$!------------------------------------------------------------------------------!
!!$   IF(SERIAL)THEN
!!$     UF(1:N,1:KBM1) = UG(1:N,1:KBM1)
!!$     VF(1:N,1:KBM1) = VG(1:N,1:KBM1)
!!$   END IF
!!$#  if defined (MULTIPROCESSOR) 
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$   CALL MPI_BCAST(UG,NGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
!!$   CALL MPI_BCAST(VG,NGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$
!!$   IF(PAR)THEN
!!$     DO I=1,N
!!$       UF(I,1:KBM1) = UG(EGID(I),1:KBM1)
!!$       VF(I,1:KBM1) = VG(EGID(I),1:KBM1)
!!$     END DO
!!$   END IF
!!$#  endif
!!$
!!$   DEALLOCATE(UG,VG)
!!$          
!!$
!!$   RETURN
!!$   END SUBROUTINE CURRENT_NUDGING
!!$!==============================================================================|
!!$!==============================================================================|
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$
!!$   SUBROUTINE SET_SST_ASSIM_DATA 
!!$
!!$!------------------------------------------------------------------------------!
!!$!  SET UP ASSIMILATION DATA FOR SST OBSERVATIONS                               |
!!$!------------------------------------------------------------------------------!
!!$   USE ALL_VARS
!!$   IMPLICIT NONE
!!$   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS
!!$   CHARACTER(LEN=120) :: FNAME
!!$   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
!!$   REAL(SP):: X0,Y0,DX,DY,RD,TEMPF
!!$   LOGICAL :: FEXIST
!!$       
!!$!------------------------------------------------------------------------------!
!!$     
!!$!
!!$!--Read Number SST Measurements and Data Set Size------------------------------!
!!$!
!!$
!!$
!!$
!!$   ! IT WOULD BE SOMETHING LIKE DIMENSIONS: NODES, TIME
!!$   ! VARIABLES: SST, LON, LAT
!!$
!!$
!!$
!!$   READ(1,*) N_TIMES_SST 
!!$   READ(1,*) TEMPF,N_ASSIM_SST 
!!$   ALLOCATE(SST_OBS(N_ASSIM_SST))
!!$   DO I=1,N_ASSIM_SST
!!$     SST_OBS(I)%N_TIMES = N_TIMES_SST
!!$     ALLOCATE(SST_OBS(I)%SST(N_TIMES_SST))
!!$     ALLOCATE(SST_OBS(I)%TIMES(N_TIMES_SST))
!!$     ALLOCATE(SST_OBS(I)%T_INT_HOUR(N_DA_HOURS)) ; SST_OBS(I)%T_INT_HOUR = 0.
!!$   END DO
!!$   REWIND(1)  ; READ(1,*)
!!$   
!!$!
!!$!--Read SST Data---------------------------------------------------------------!
!!$!
!!$   DO I=1,N_TIMES_SST
!!$     READ(1,*) TEMPF,N_ASSIM_SST 
!!$     DO J=1,N_ASSIM_SST
!!$       READ(1,*)SST_OBS(J)%X,SST_OBS(J)%Y,SST_OBS(J)%SST(I)
!!$       SST_OBS(J)%TIMES(I) = TEMPF
!!$     END DO
!!$   END DO
!!$
!!$!
!!$!--Shift Coordinates-of SST Observation Locations------------------------------!
!!$!
!!$   SST_OBS(:)%X = SST_OBS(:)%X - VXMIN 
!!$   SST_OBS(:)%Y = SST_OBS(:)%Y - VYMIN 
!!$   
!!$!
!!$!--Close SST Observation Data File---------------------------------------------!
!!$!
!!$   CLOSE(1)
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Compute Spatial Interpolation Weights for each Observation Location 
!!$!------------------------------------------------------------------------------!
!!$
!!$   ALLOCATE(ITEMP(MGL),FTEMP(MGL),DA_SST(MGL))     ; DA_SST = 0
!!$   DO I=1,N_ASSIM_SST
!!$     X0 = SST_OBS(I)%X
!!$     Y0 = SST_OBS(I)%Y
!!$     ECNT = 0
!!$     DO J=1,MGL
!!$       DX = ABS(XG(J)-X0)
!!$       DY = ABS(YG(J)-Y0)
!!$       RD = SQRT(DX**2 + DY**2)
!!$       IF(RD <= RAD_SST)THEN
!!$         DA_SST(J) = 1
!!$         ECNT = ECNT + 1      
!!$         ITEMP(ECNT) =  J
!!$         FTEMP(ECNT) = (RAD_SST**2 - RD**2) / (RAD_SST**2 + RD**2)
!!$       END IF
!!$     END DO
!!$     IF(ECNT == 0)THEN
!!$       SST_OBS(I)%N_INTPTS = ECNT
!!$!QXU{ 
!!$!       ALLOCATE(SST_OBS(I)%INTPTS(1))
!!$!       ALLOCATE(SST_OBS(I)%X_WEIGHT(1))
!!$       ALLOCATE(SST_OBS(I)%INTPTS(0:1))
!!$       ALLOCATE(SST_OBS(I)%X_WEIGHT(0:1))
!!$!QXU}
!!$     ELSE
!!$       SST_OBS(I)%N_INTPTS = ECNT
!!$       ALLOCATE(SST_OBS(I)%INTPTS(ECNT))
!!$       ALLOCATE(SST_OBS(I)%X_WEIGHT(ECNT))
!!$       SST_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
!!$       SST_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
!!$     END IF
!!$   END DO
!!$   DEALLOCATE(FTEMP,ITEMP)
!!$     
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Report Number of Interpolation Points, Location and Number of Data 
!!$!------------------------------------------------------------------------------!
!!$   IF(.NOT. MSR)RETURN
!!$   WRITE(IPT,*)
!!$   WRITE(IPT,*)'!                SST OBSERVATION DATA           '
!!$!   WRITE(IPT,*)"  OBS#      X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES"
!!$!   DO I=1,N_ASSIM_SST
!!$!     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X)') &
!!$!       I,SST_OBS(I)%X/1000.,SST_OBS(I)%Y/1000.,SST_OBS(I)%N_INTPTS,SST_OBS(I)%N_TIMES
!!$!   END DO
!!$   WRITE(IPT,*)'NUMBER OF NODES WITHOUT NUDGING: ',MGL-SUM(DA_SST(1:MGL))
!!$   WRITE(IPT,*)'MAXIMUM NUMBER OF INTERP POINTS: ',MAXVAL(SST_OBS(1:N_ASSIM_SST)%N_INTPTS)
!!$   WRITE(IPT,*)
!!$   RETURN
!!$   END SUBROUTINE SET_SST_ASSIM_DATA
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$   SUBROUTINE SST_NUDGING
!!$!==============================================================================|
!!$!  USE SST OBSERVATION DATA TO NUDGE SURFACE TEMPERATURE                       |
!!$!==============================================================================|
!!$   USE MOD_PREC
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR
!!$#  endif
!!$   IMPLICIT NONE
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1G
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TINT,TCORR
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TWEIGHT,FTEMP
!!$   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
!!$   REAL(SP) :: TRUTH,ttemp,TIME_24,TMP,ASTIME,tmax,TSHIFT
!!$   INTEGER I,J,J1,ITIME,NTIME,IERR,ITMP,jtmax,ibad
!!$   INTRINSIC MINLOC
!!$!==============================================================================|
!!$
!!$   IF(MSR)THEN
!!$     IF(SST_CYC > N_TIMES_SST)THEN
!!$     WRITE(IPT,*)'NUMBER OF DA SWEEPS EXCEEDS NUMBER OF SST DATA INTERVALS'
!!$     WRITE(IPT,*)'NUMBER OF SWEEPS: ',SST_CYC
!!$     WRITE(IPT,*)'NUMBER OF DATA INTERVALS: ',N_TIMES_SST
!!$     WRITE(IPT,*)'HALTING'
!!$     CALL PSTOP
!!$     END IF
!!$   END IF
!!$!------------------------------------------------------------------------------!
!!$!  Calculate Hour Number of Current Assimilation Cycle                         ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   TSHIFT = THOUR-FLOAT(ISTART_DAY*24)-FLOAT((SST_CYC-1)*IAV_DAY*24)
!!$   DA_HOUR = INT(TSHIFT+.5)
!!$   SST_OBS(:)%N_T_WEIGHT = DA_HOUR
!!$   IF(DA_HOUR == 0)DA_HOUR = 1
!!$   ASTIME = ABS(TSHIFT-DA_HOUR)*3600.
!!$
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$
!!$   IF(ASTIME < 0.5*ASTIME_WINDOW_SST) THEN
!!$     SST_OBS(:)%T_WEIGHT = 1.0
!!$   ELSE IF(ASTIME < ASTIME_WINDOW_SST) THEN
!!$     SST_OBS(:)%T_WEIGHT = (ASTIME_WINDOW_SST-ASTIME)/ASTIME_WINDOW_SST*2.0
!!$   ELSE
!!$     SST_OBS(:)%T_WEIGHT = 0.0
!!$   END IF
!!$
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Gather Surface Temperature Field to Master Processor                        ! 
!!$!------------------------------------------------------------------------------!
!!$   ALLOCATE(T1G(MGL,1:KB))
!!$#  if defined (MULTIPROCESSOR)
!!$   IF(PAR)THEN
!!$     CALL GATHER(LBOUND(T1,1),UBOUND(T1,1),M,MGL,KB,MYID,NPROCS,NMAP,T1,T1G)
!!$   END IF
!!$#  endif
!!$   IF(SERIAL)THEN
!!$     T1G(1:MGL,1) = T1(1:MGL,1)
!!$   END IF
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Interpolate Simulation Data to Local Observation Point                      ! 
!!$!------------------------------------------------------------------------------!
!!$       
!!$   IF(MSR)THEN
!!$   ALLOCATE(TINT(N_ASSIM_SST)) ; TINT = 0. 
!!$
!!$   DO I=1,N_ASSIM_SST  
!!$     IF(SST_OBS(I)%N_INTPTS > 0)THEN
!!$     DO J=1,SST_OBS(I)%N_INTPTS
!!$       J1        = SST_OBS(I)%INTPTS(J)
!!$       WGHT      = SST_OBS(I)%X_WEIGHT(J)
!!$       TINT(I)   = TINT(I) + SST_OBS(I)%X_WEIGHT(J)*T1G(J1,1) 
!!$     END DO
!!$     TOT_WGHT = SUM(SST_OBS(I)%X_WEIGHT(1:SST_OBS(I)%N_INTPTS))
!!$     TINT(I)  = TINT(I)/TOT_WGHT
!!$     END IF 
!!$   END DO
!!$     
!!$!------------------------------------------------------------------------------!
!!$!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   ALLOCATE(TWEIGHT(MGL))    ; TWEIGHT = 0.
!!$   ALLOCATE(TCORR(MGL))      ; TCORR   = 0.
!!$
!!$   DO I=1,N_ASSIM_SST
!!$     DO J=1,SST_OBS(I)%N_INTPTS
!!$       J1           = SST_OBS(I)%INTPTS(J)
!!$       WEIGHT       = SST_OBS(I)%T_WEIGHT*SST_OBS(I)%X_WEIGHT(J)
!!$       TWEIGHT(J1)  = TWEIGHT(J1) + WEIGHT
!!$       TRUTH        = &
!!$           SST_OBS(I)%T_INT_HOUR(DA_HOUR) + SST_OBS(I)%SST(SST_CYC) - SST_OBS(I)%T_INT_AVGD
!!$       DEFECT       = TRUTH - TINT(I)
!!$       CORRECTION   = GAMA_SST*DEFECT*WEIGHT**2
!!$       TCORR(J1)    = TCORR(J1) + CORRECTION
!!$!qxu{check 
!!$       if(j1.eq.24122) then
!!$!                   IINT, NO_site,SST_obs,SST_avg,SST_hour,SST_now
!!$         write(302,'(2i7,6f7.2,4f8.4)') IINT,I,&
!!$	             SST_OBS(I)%SST(SST_CYC),&
!!$                     SST_OBS(I)%T_INT_AVGD,&
!!$		     SST_OBS(I)%T_INT_HOUR(DA_HOUR),&
!!$		     TRUTH,&  
!!$		     TINT(I),&
!!$		     DEFECT,&
!!$		     WEIGHT,&
!!$		     SST_OBS(I)%X_WEIGHT(J),&
!!$		     SST_OBS(I)%T_WEIGHT,&
!!$		     CORRECTION
!!$       endif
!!$		     
!!$!qxu}
!!$     END DO
!!$   END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Nudge Simulation Data Using Local Corrections                               ! 
!!$!------------------------------------------------------------------------------!
!!$!QXU{check weight 
!!$   ibad = 24122    !bad node the temp(24122,1)=36
!!$   write(300,'(i7,3f10.3,E12.5)') iint, &
!!$                 T1G(ibad,1) + DTI*GALPHA_SST*TCORR(ibad)/TWEIGHT(ibad),&
!!$		 T1G(ibad,1),&
!!$		 TCORR(ibad),TWEIGHT(ibad)
!!$!QXU}
!!$   DO I=1,MGL
!!$     IF(DA_SST(I) == 1 .AND. TWEIGHT(I) > 1.0E-08)THEN
!!$       T1G(I,1) = T1G(I,1) + DTI*GALPHA_SST*TCORR(I)/TWEIGHT(I)
!!$     END IF
!!$   END DO
!!$
!!$   DEALLOCATE(TWEIGHT,TCORR,TINT)
!!$   END IF !!MSR
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Disperse New Data Fields to Slave Processors                                ! 
!!$!------------------------------------------------------------------------------!
!!$   IF(SERIAL)THEN
!!$     T1(1:M,1) = T1G(1:M,1)
!!$   END IF
!!$!!   goto 100  !gwc debug dont nudge
!!$#  if defined (MULTIPROCESSOR) 
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$   CALL MPI_BCAST(T1G,MGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$
!!$   IF(PAR)THEN
!!$     DO I=1,M
!!$       T1(I,1) = T1G(NGID(I),1)
!!$     END DO
!!$   END IF
!!$#  endif
!!$!!   100 continuE
!!$
!!$   DEALLOCATE(T1G)
!!$          
!!$
!!$   RETURN
!!$   END SUBROUTINE SST_NUDGING
!!$!==============================================================================|
!!$!==============================================================================|
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$   SUBROUTINE SST_INT
!!$!==============================================================================|
!!$!  INTERPOLATE HOURLY AND INTERVAL-AVERAGAED SST DATA TO OBSERVATIN POINTS     | 
!!$!==============================================================================|
!!$   USE MOD_PREC
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR
!!$#  endif
!!$   IMPLICIT NONE
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: T_TMP_L,T_TMP_G 
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: TWEIGHT,FTEMP
!!$   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
!!$   INTEGER I,J,J1,ITIME,NTIME,IERR
!!$   INTRINSIC MINLOC
!!$!==============================================================================|
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Calculate Hour Number of Current Assimilation Cycle                         ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   IF(ABS(THOUR-INT(THOUR+0.5)) < 1.E-8) THEN
!!$     DA_HOUR = MOD(INT(THOUR-FLOAT(ISTART_DAY*24)+0.5),24*IAV_DAY)
!!$     IF(DA_HOUR == 0) DA_HOUR = 24*IAV_DAY
!!$   ELSE
!!$     RETURN
!!$   END IF
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Extract Sea Surface Temp and Gather to Global Array                         ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   ALLOCATE(T_TMP_G(MGL))
!!$#  if defined (MULTIPROCESSOR)
!!$   IF(PAR)THEN
!!$     ALLOCATE(T_TMP_L(0:M))
!!$     T_TMP_L(1:M) = T1(1:M,1)
!!$!QXU{
!!$     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$     CALL GATHER(0,M,M,MGL,1,MYID,NPROCS,NMAP,T_TMP_L,T_TMP_G)
!!$   END IF
!!$#  endif
!!$   IF(SERIAL)THEN
!!$     T_TMP_G(1:MGL) = T1(1:MGL,1)
!!$   END IF
!!$
!!$
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Interpolate Hourly Data to Local Observation Point                          ! 
!!$!------------------------------------------------------------------------------!
!!$   IF(MSR)THEN
!!$   DO I=1,N_ASSIM_SST
!!$     SST_OBS(I)%T_INT_HOUR(DA_HOUR) = 0.
!!$     IF(SST_OBS(I)%N_INTPTS /= 0)THEN
!!$     DO J=1,SST_OBS(I)%N_INTPTS
!!$       J1        = SST_OBS(I)%INTPTS(J)
!!$       WGHT      = SST_ObS(I)%X_WEIGHT(J)
!!$!QXU{check no_obs=2122 why the temp = 9 ?? intep node to this obs site
!!$       IF(I.EQ.2122) THEN
!!$         write(304,'(4i7,2f10.4)') iint,DA_HOUR,J,J1,T_TMP_G(J1),WGHT
!!$       ENDIF
!!$!QXU}
!!$       SST_OBS(I)%T_INT_HOUR(DA_HOUR)   = &
!!$               SST_OBS(I)%T_INT_HOUR(DA_HOUR) + WGHT*T_TMP_G(J1)
!!$     END DO
!!$     TOT_WGHT = SUM(SST_OBS(I)%X_WEIGHT(1:SST_OBS(I)%N_INTPTS))
!!$     SST_OBS(I)%T_INT_HOUR(DA_HOUR)  = SST_OBS(I)%T_INT_HOUR(DA_HOUR)/TOT_WGHT
!!$     END IF
!!$
!!$   END DO
!!$  
!!$!   WRITE(IPT,*)'  INTERPOLATING HOURLY SST DATA',DA_HOUR,N_DA_HOURS
!!$   END IF
!!$!QXU{
!!$#  if defined (MULTIPROCESSOR)
!!$   call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$#  endif
!!$!QXU}
!!$
!!$   DEALLOCATE(T_TMP_G)
!!$   IF(PAR) DEALLOCATE(T_TMP_L)
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Average Hourly Data Over Integration Period                                 ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   IF(DA_HOUR /=  N_DA_HOURS)RETURN
!!$
!!$   IF(MSR)THEN
!!$   DO I=1,N_ASSIM_SST
!!$     SST_OBS(I)%T_INT_AVGD = SUM(SST_OBS(I)%T_INT_HOUR(1:N_DA_HOURS))/FLOAT(N_DA_HOURS)
!!$!QXU{check average
!!$     write(301,*) SST_OBS(I)%X+VXMIN,SST_OBS(I)%Y+VYMIN,SST_OBS(I)%T_INT_AVGD
!!$!QXU}
!!$   END DO
!!$!QXU{check
!!$   do i=1,N_DA_HOURS
!!$      write(303,'(2i7,2f10.3)') iint,i,SST_OBS(2122)%T_INT_HOUR(I),SST_OBS(2123)%T_INT_HOUR(I)
!!$   enddo
!!$!QXU}
!!$   
!!$   WRITE(IPT,*)'  COMPUTING AVERAGE INTERPOLATED SST DATA',DA_HOUR,N_DA_HOURS
!!$   END IF
!!$!QXU{
!!$#  if defined (MULTIPROCESSOR)
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$#  endif
!!$!QXU}   
!!$   
!!$   RETURN
!!$   END SUBROUTINE SST_INT
!!$
!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!   READ IN DATA FROM SIMULATION STAGE AND FOR ASSIMILATION STAGE STARTUP      |
!==============================================================================|

  SUBROUTINE SST_RESTORE_STATE      !=HOT_START_SST

    !------------------------------------------------------------------------------|
    USE EQS_OF_STATE
    USE ALL_VARS
    USE MOD_PAR
    IMPLICIT NONE

    if(dbg_set(dbg_sbr)) WRITE(IPT,*) "start SST_RESTORE_STATE"

    IINT   = IINT_BUF
    ExtTime=ExtTime_BUF
    IntTime=IntTime_BUF

    IF(RECALCULATE_RHO_MEAN) THEN
       RECALC_RHO_MEAN =IntTime_BUF
    END IF


    U      = U_BUF
    V      = V_BUF
    WTS    = WTS_BUF
    S1     = S1_BUF
    T1     = T1_BUF

    KM     = KM_BUF
    KH     = KH_BUF
    KQ     = KQ_BUF

    UA     = UA_BUF
    VA     = VA_BUF

    EL     = EL_BUF
    ET     = ET_BUF
    H      = H_BUF

#  if defined (EQUI_TIDE)
    EL_EQI = EL_EQI_BUF
#  endif

#  if defined (ATMO_TIDE)
    EL_ATMO= EL_ATMO_BUF
#  endif


# if defined(GOTM)
    TEPS   = TEPS_BUF
    TKE    = TKE_BUF
# else
    Q2     = Q2_BUF
    Q2L    = Q2L_BUF
    L      = L_BUF
# endif   

#  if defined (DYE_RELEASE)   
    DYE    = DYE_BUF
    DYEF   = DYEF_BUF
    DYEMEAN = DYEMEAN_BUF
#  endif    


    ! TO SAVE MEMEORY WE SAVED ONLY THE INTERNAL VALUES: CALL EXCHANGE
    ! TO SET THE HALO ELEMENTS
    IF(PAR) CALL EXCHANGE_ALL

    ! SET DEPTH AN INTERPOLATE VALUES TO CELL CENTERS
    D  = H + EL
    DT = H + ET

    CALL N2E2D(H,H1)
    CALL N2E2D(EL,EL1)
    CALL N2E2D(D,D1)
    CALL N2E2D(DT,DT1)

    ! SET THE DENSITY
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

    CALL N2E3D(T1,T)
    CALL N2E3D(S1,S)
    !DENSITY IS ALREADY INTERPOLATED TO ELEMENT IN DENSX

    ! SET THE VERTICAL SIGMA VELOCITY
    CALL N2E3D(WTS,W)

    ! SET THE TURBULENT QUANTITIES
#    if defined (GOTM)
    L = .001 
    L(1:MT,2:KBM1) = (.5544**3)*TKE(1:MT,2:KBM1)**1.5/TEPS(1:MT,2:KBM1)
# endif

    CALL N2E3D(KM,KM1)


    if(dbg_set(dbg_sbr)) WRITE(IPT,*) "END SST_RESTORE_STATE"


    RETURN
  END SUBROUTINE SST_RESTORE_STATE
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============================================================================|
!   EXCHANGE HALO'S BECAUSE WE ONLY SAVED INTERNAL VALUES TO SAVE MEMORY       |
!==============================================================================|
   SUBROUTINE EXCHANGE_ALL 

!------------------------------------------------------------------------------|

     USE ALL_VARS

!!$# if defined (WATER_QUALITY)
!!$     USE MOD_WQM
!!$# endif
     
     USE MOD_PAR

#  if defined (DYE_RELEASE)
     USE MOD_DYE
#  endif   
     
     IMPLICIT NONE

!==============================================================================|
   

#  if defined (MULTIPROCESSOR)

#  if defined (GOTM)
   CALL AEXCHANGE(NC,MYID,NPROCS,TKE,TEPS)
#  else
   CALL AEXCHANGE(NC,MYID,NPROCS,Q2,Q2L,L)
#  endif

   CALL AEXCHANGE(NC,MYID,NPROCS,KM,KQ,KH)
   CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
   CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA)

   CALL AEXCHANGE(NC,MYID,NPROCS,S1,T1,WTS)
   CALL AEXCHANGE(NC,MYID,NPROCS,EL,H,ET)

#  if defined (EQUI_TIDE)
   CALL AEXCHANGE(NC,MYID,NPROCS,EL_EQI)
#  endif

#  if defined (ATMO_TIDE)
   CALL AEXCHANGE(NC,MYID,NPROCS,EL_ATMO)
#  endif

#  if defined (WATER_QUALITY)
   CALL FATAL_ERROR("SST ASSIMILATION IS NOT SET UP FOR WATER QUALITY MODULE")
#  endif

#  if defined (DYE_RELEASE)
   CALL AEXCHANGE(NC,MYID,NPROCS,DYEMEAN)
   CALL AEXCHANGE(NC,MYID,NPROCS,DYE)
#  endif   

#  endif

   RETURN
   END SUBROUTINE EXCHANGE_ALL

!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============================================================================|
!   DUMP DATA FILE FOR ASSIMILATION RESTART                                    |
!==============================================================================|

   SUBROUTINE SST_SAVE_STATE      !=ARC_SST 

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE

   ExtTime_BUF = ExtTime
   IntTime_BUF = IntTime
   IINT_BUF   = IINT

   U_BUF      = U
   V_BUF      = V
   WTS_BUF    = WTS
   S1_BUF     = S1
   T1_BUF     = T1

   KM_BUF     = KM
   KH_BUF     = KH
   KQ_BUF     = KQ

#    if defined (GOTM)
   TKE_BUF    = TKE
   TEPS_BUF   = TEPS
# else
   Q2_BUF     = Q2
   Q2L_BUF    = Q2L
   L_BUF      = L
# endif

   UA_BUF     = UA
   VA_BUF     = VA

   EL_BUF     = EL
   ET_BUF     = ET
   H_BUF      = H
   
#  if defined (DYE_RELEASE)   
   DYE_BUF    = DYE
   DYEF_BUF   = DYEF
   DYEMEAN_BUF = DYEMEAN
#  endif

   RETURN
 END SUBROUTINE SST_SAVE_STATE
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE ALLOC_BUFFER 
!------------------------------------------------------------------------------|
   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE

   ALLOCATE(U_BUF(0:N,KB));   U_BUF = 0.0_SP
   ALLOCATE(V_BUF(0:N,KB));   V_BUF = 0.0_SP
   ALLOCATE(WTS_BUF(0:M,KB)); WTS_BUF = 0.0_SP

   ALLOCATE(S1_BUF(0:M,KB));  S1_BUF = 0.0_SP
   ALLOCATE(T1_BUF(0:M,KB));  T1_BUF = 0.0_SP

   ALLOCATE(KM_BUF(0:M,KB));  KM_BUF = 0.0_SP
   ALLOCATE(KH_BUF(0:M,KB));  KH_BUF = 0.0_SP
   ALLOCATE(KQ_BUF(0:M,KB));  KQ_BUF = 0.0_SP

   ALLOCATE(UA_BUF(0:N));     UA_BUF = 0.0_SP
   ALLOCATE(VA_BUF(0:N));     VA_BUF = 0.0_SP

   ALLOCATE(EL_BUF(0:M));     EL_BUF = 0.0_SP
   ALLOCATE(ET_BUF(0:M));     ET_BUF = 0.0_SP
   ALLOCATE(H_BUF(0:M));      H_BUF = 0.0_SP

#  if defined (EQUI_TIDE)
   ALLOCATE(EL_EQI_BUF(0:M)); EL_EQI_BUF = 0.0_SP
#  endif

#  if defined (ATMO_TIDE)
   ALLOCATE(EL_ATMO_BUF(0:M)); EL_ATMO_BUF = 0.0_SP
#  endif


#    if defined (GOTM)
   ALLOCATE(TKE_BUF(0:M,KB)); TKE_BUF = 0.0_SP
   ALLOCATE(TEPS_BUF(0:M,KB));TEPS_BUF = 0.0_SP
# else
   ALLOCATE(Q2_BUF(0:M,KB));  Q2_BUF = 0.0_SP
   ALLOCATE(Q2L_BUF(0:M,KB)); Q2L_BUF = 0.0_SP
   ALLOCATE(L_BUF(0:M,KB));   L_BUF = 0.0_SP
# endif

#  if defined (DYE_RELEASE)
   ALLOCATE(DYE_BUF(0:M,KB)); DYE_BUF = 0.0_SP
   ALLOCATE(DYEF_BUF(0:M,KB)); DYEF_BUF = 0.0_SP
   ALLOCATE(DYEMEAN_BUF(0:M,KB)); DYEMEAN_BUF = 0.0_SP
#  endif   
   
   RETURN
   END SUBROUTINE ALLOC_BUFFER 
!==============================================================================|


!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$
!!$
!!$   SUBROUTINE SET_TS_ASSIM_DATA 
!!$
!!$!------------------------------------------------------------------------------!
!!$!  SET UP ASSIMILATION DATA FOR TEMP/SAL OBSERVATIONS                          |
!!$!------------------------------------------------------------------------------!
!!$   USE ALL_VARS
!!$   IMPLICIT NONE
!!$   INTEGER I,J,K,ECNT,ITMP,NCNT,IOS,NLAY
!!$   CHARACTER(LEN=120) :: FNAME,ONAME
!!$   CHARACTER(LEN= 2 ) :: NAC   
!!$   INTEGER,  ALLOCATABLE, DIMENSION(:) :: ITEMP
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:) :: FTEMP
!!$   REAL(SP):: X0,Y0,DX,DY,RD,SIGMA_C,ISOBATH_ANGLE,D_ANGLE,ANG_OBS_SIM,DIR_WEIGHT
!!$   REAL(SP), PARAMETER :: ALF = 0.05_SP
!!$   LOGICAL :: FEXIST
!!$   INTEGER :: MAXEL,NBD_CNT
!!$   integer :: jmin 
!!$   real(sp) :: lmin
!!$
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Read Number of Scalar Observations and Coordinates of Each                  !
!!$!------------------------------------------------------------------------------!
!!$       
!!$   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_ts.xy"
!!$!
!!$!--Make Sure Temperature and Salinity Assimilation Data File Exists-----------------------------!
!!$!
!!$   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
!!$   IF(MSR .AND. .NOT.FEXIST)THEN
!!$     WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',FNAME,' DOES NOT EXIST'
!!$     WRITE(IPT,*)'HALTING.....'
!!$     CALL PSTOP
!!$   END IF
!!$     
!!$!
!!$!--Read Number of Current Measurement Stations---------------------------------!
!!$!
!!$   OPEN(1,FILE=TRIM(FNAME),STATUS='OLD')
!!$   READ(1,*) N_ASSIM_TS                  !nomber of TS Obs station
!!$   ALLOCATE(TS_OBS(N_ASSIM_TS))          !Type for TS_OBS
!!$
!!$!
!!$!--Read X,Y Coordinate of Measurement Stations---------------------------------!
!!$!
!!$
!!$   DO I=1,N_ASSIM_TS  
!!$     READ(1,*)ITMP,TS_OBS(I)%X,TS_OBS(I)%Y,TS_OBS(I)%DEPTH,NLAY,TS_OBS(I)%SITA
!!$     TS_OBS(I)%N_LAYERS = NLAY
!!$     ALLOCATE(TS_OBS(I)%ODEPTH(NLAY))
!!$     DO J=1,NLAY
!!$       READ(1,*)TS_OBS(I)%ODEPTH(J)
!!$       IF(TS_OBS(I)%ODEPTH(J) > TS_OBS(I)%DEPTH)THEN
!!$         IF(MSR)WRITE(IPT,*)'OBSERVATION DEPTH',J,'OF TEMP/SALINITY OBS',I
!!$         IF(MSR)WRITE(IPT,*)'EXCEEDS BATHYMETRIC DEPTH'
!!$         IF(MSR)WRITE(IPT,*)'HALTING...........'
!!$         CALL PSTOP
!!$       END IF
!!$     END DO
!!$   END DO
!!$   MAX_LAYER_TS = MAXVAL(TS_OBS(1:N_ASSIM_TS)%N_LAYERS)
!!$
!!$!
!!$!--Shift Coordinates-----------------------------------------------------------!
!!$!
!!$   TS_OBS(:)%X = TS_OBS(:)%X - VXMIN 
!!$   TS_OBS(:)%Y = TS_OBS(:)%Y - VYMIN 
!!$   
!!$!
!!$!--Close Current Observation Global File---------------------------------------!
!!$!
!!$   CLOSE(1)
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Open Temp/Sal Observation Files for Each Observation Point and Read Data    !
!!$!------------------------------------------------------------------------------!
!!$
!!$   DO I=1,N_ASSIM_TS
!!$
!!$!----Make Sure Temperature/Salinity Observation File Exists--------------------!
!!$     WRITE(NAC,'(I2.2)')I
!!$     ONAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//'_ts'//NAC//'.dat' 
!!$     INQUIRE(FILE=TRIM(ONAME),EXIST=FEXIST)
!!$     IF(MSR .AND. .NOT.FEXIST)THEN
!!$       WRITE(IPT,*)'TEMP/SALINITY OBSERVATION FILE: ',ONAME,' DOES NOT EXIST'
!!$       WRITE(IPT,*)'HALTING.....'
!!$       CALL PSTOP
!!$     END IF
!!$
!!$!----Open Temp/Salinity Observation File for Read------------------------------!
!!$     OPEN(1,FILE=ONAME,STATUS='old')  ; REWIND(1)
!!$     NCNT = 0
!!$
!!$!----Count Number of Data Entries for Observation I----------------------------!
!!$     DO WHILE(.TRUE.)
!!$       READ(1,*,IOSTAT=IOS)
!!$       IF(IOS < 0)EXIT
!!$       NCNT = NCNT + 1
!!$     END DO
!!$     TS_OBS(I)%N_TIMES = NCNT
!!$     REWIND(1)
!!$     IF(NCNT == 0)THEN
!!$       IF(MSR)WRITE(IPT,*)'NO DATA FOR TEMPERATURE/SALINITY OBSERVATION',I
!!$       CALL PSTOP
!!$     END IF
!!$
!!$!----Allocate Arrays to Hold Temp/Salinity (TEMP/SAL) and Time (TIME)----------!
!!$     ALLOCATE(TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES))
!!$     ALLOCATE(TS_OBS(I)%TEMP( TS_OBS(I)%N_TIMES , TS_OBS(I)%N_LAYERS ))
!!$     ALLOCATE(TS_OBS(I)%SAL(  TS_OBS(I)%N_TIMES , TS_OBS(I)%N_LAYERS ))
!!$
!!$!----Read in Current Data for Observation I------------------------------------!
!!$     NLAY = TS_OBS(I)%N_LAYERS
!!$     DO J=1,TS_OBS(I)%N_TIMES
!!$     READ(1,*)TS_OBS(I)%TIMES(J),(TS_OBS(I)%TEMP(J,K),TS_OBS(I)%SAL(J,K),K=1,NLAY)
!!$     END DO
!!$     CLOSE(1)
!!$     
!!$
!!$!----Convert Time to Seconds---------------------------------------------------!
!!$!----Shift Jan 1 Based Time Data to Dec 1 Based Time Data-----CASESPECIFIC-----!
!!$     IF(trim(CASENAME) == 'gom')THEN
!!$!      becuse the mode runing from 1th Nov last year for spinup run so add 61d
!!$!       TS_OBS(I)%TIMES = ((TS_OBS(I)%TIMES-1.0_SP )*24.0_SP + 744.0_SP)*3600.0_SP
!!$!       TS_OBS(I)%TIMES = (TS_OBS(I)%TIMES-1.0_SP + 61.0_SP)*24.0_SP*3600.0_SP !for 1995
!!$       TS_OBS(I)%TIMES = (TS_OBS(I)%TIMES-1.0_SP)*24.0_SP*3600.0_SP !after 1996
!!$     ELSE   
!!$       TS_OBS(I)%TIMES = TS_OBS(I)%TIMES*3600.0_SP*24.0_SP
!!$     END IF
!!$ 
!!$!----Convert Temperature and Salinity to PSU/Celsius-(If Necessary)------------!
!!$     TS_OBS(I)%TEMP = TS_OBS(I)%TEMP 
!!$     TS_OBS(I)%SAL  = TS_OBS(I)%SAL  
!!$  END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Count Number of Points with Bad Data (TEMP = 0. .OR. SAL = 0.)         
!!$!------------------------------------------------------------------------------!
!!$  NBD_CNT = 0
!!$  DO I=1,N_ASSIM_TS
!!$  DO J=1,TS_OBS(I)%N_TIMES
!!$  DO K=1,TS_OBS(I)%N_LAYERS
!!$  IF(TS_OBS(I)%TEMP(J,K) < -90. .OR. TS_OBS(I)%SAL(J,K) < -90.) THEN
!!$    NBD_CNT = NBD_CNT + 1
!!$  END IF
!!$  END DO
!!$  END DO
!!$  END DO
!!$ 
!!$!------------------------------------------------------------------------------!
!!$!  Compute Spatial Interpolation Weights for each Mooring Location 
!!$!------------------------------------------------------------------------------!
!!$!   LMIN = 100000000.
!!$!   DO I=1,N_ASSIM_TS
!!$!     X0 = TS_OBS(I)%X
!!$!     Y0 = TS_OBS(I)%Y
!!$!     DO J=1,MGL
!!$!       DX = ABS(XG(J)-X0)
!!$!       DY = ABS(YG(J)-Y0)
!!$!       IF(SQRT(DX**2 + DY**2) < LMIN)THEN
!!$!         LMIN = SQRT(DX**2 + DY**2)
!!$!         JMIN = J
!!$!       END IF
!!$!      END DO
!!$!      TS_OBS(I)%SITA = SITA_GD(JMIN) + 3.14159_SP/2.0_SP
!!$!    end do
!!$
!!$
!!$!JQI<NGL ==> MGL (05/11/05)
!!$!   ALLOCATE(ITEMP(NGL),FTEMP(NGL),DA_TS(NGL))     ; DA_TS = 0
!!$   ALLOCATE(ITEMP(MGL),FTEMP(MGL),DA_TS(MGL))     ; DA_TS = 0
!!$!JQI>   
!!$   DO I=1,N_ASSIM_TS
!!$     X0 = TS_OBS(I)%X
!!$     Y0 = TS_OBS(I)%Y
!!$     ISOBATH_ANGLE = TS_OBS(I)%SITA/180.0_SP*3.1415926_SP
!!$     ECNT = 0
!!$     DO J=1,MGL
!!$       DX = ABS(XG(J)-X0)
!!$       DY = ABS(YG(J)-Y0)
!!$       RD = SQRT(DX**2 + DY**2)
!!$       IF(RD <= RAD_TS)THEN
!!$         DA_TS(J)   = 1
!!$         ECNT        = ECNT + 1      
!!$         ITEMP(ECNT) =  J
!!$         FTEMP(ECNT) = (RAD_TS**2 - RD**2) / (RAD_TS**2 + RD**2)
!!$         ANG_OBS_SIM = ATAN2(DY,DX)
!!$         D_ANGLE     = ANG_OBS_SIM - ISOBATH_ANGLE
!!$         D_ANGLE     = D_ANGLE - INT(D_ANGLE/3.1415926_SP)*3.1415926_SP
!!$         D_ANGLE     = ABS(D_ANGLE)
!!$         DIR_WEIGHT  = (ABS(D_ANGLE-0.5*3.1415926_SP)+ALF*3.1415926_SP)/ &
!!$                       ((0.5_SP+ALF)*3.1415926_SP)
!!$         FTEMP(ECNT) = FTEMP(ECNT)*DIR_WEIGHT
!!$       END IF
!!$     END DO
!!$     IF(ECNT == 0)THEN
!!$       WRITE(IPT,*)'ERROR SETTING UP TEMP/SAL DATA ASSIMILATION'
!!$       WRITE(IPT,*)'NO ELEMENTS LIE WITHIN RADIUS',RAD_TS
!!$       WRITE(IPT,*)'OF OBSERVATION POINT',I
!!$       CALL PSTOP   
!!$     ELSE
!!$       TS_OBS(I)%N_INTPTS = ECNT
!!$       ALLOCATE(TS_OBS(I)%INTPTS(ECNT))
!!$       ALLOCATE(TS_OBS(I)%X_WEIGHT(ECNT))
!!$       TS_OBS(I)%INTPTS(1:ECNT)  = ITEMP(1:ECNT)
!!$       TS_OBS(I)%X_WEIGHT(1:ECNT) = FTEMP(1:ECNT)
!!$     END IF
!!$   END DO
!!$   DEALLOCATE(FTEMP,ITEMP)
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Compute Sigma Layer Weights for Vertical Interpolation            
!!$!------------------------------------------------------------------------------!
!!$   DO I=1,N_ASSIM_TS
!!$     NLAY = TS_OBS(I)%N_LAYERS
!!$     ALLOCATE(TS_OBS(I)%S_INT(NLAY,2))
!!$     ALLOCATE(TS_OBS(I)%S_WEIGHT(NLAY,2))
!!$     DO J=1,NLAY
!!$       SIGMA_C = -TS_OBS(I)%ODEPTH(J)/TS_OBS(I)%DEPTH
!!$       DO K=2,KBM1
!!$         IF(ZZ(K) <= SIGMA_C .AND. ZZ(K-1) > SIGMA_C)THEN 
!!$           TS_OBS(I)%S_INT(J,1) = K-1
!!$           TS_OBS(I)%S_INT(J,2) = K
!!$           TS_OBS(I)%S_WEIGHT(J,1) = (SIGMA_C-ZZ(K))/(ZZ(K-1)-ZZ(K))
!!$           TS_OBS(I)%S_WEIGHT(J,2) = 1.0_SP - TS_OBS(I)%S_WEIGHT(J,1) 
!!$         END IF  
!!$       END DO
!!$!JQI       IF(ZZ(1) < SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
!!$       IF(ZZ(1) <= SIGMA_C)THEN  !!OBSERVATION ABOVE CENTROID OF FIRST SIGMA LAYER
!!$         TS_OBS(I)%S_INT(J,1) = 1
!!$         TS_OBS(I)%S_INT(J,2) = 1
!!$         TS_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
!!$         TS_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
!!$       END IF
!!$       IF(ZZ(KBM1) > SIGMA_C)THEN !!OBSERVATION BELOW CENTROID OF BOTTOM SIGMA LAYER
!!$         TS_OBS(I)%S_INT(J,1) = KBM1
!!$         TS_OBS(I)%S_INT(J,2) = KBM1
!!$         TS_OBS(I)%S_WEIGHT(J,1) = 1.0_SP
!!$         TS_OBS(I)%S_WEIGHT(J,2) = 0.0_SP
!!$       END IF
!!$
!!$     END DO
!!$   END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Report Number of Interpolation Points, Location and Number of Data 
!!$!------------------------------------------------------------------------------!
!!$   IF(.NOT. MSR)RETURN
!!$
!!$   WRITE(IPT,*)
!!$   WRITE(IPT,*)'!            TEMP/SALINITY OBSERVATION DATA           '
!!$   WRITE(IPT,*)" MOORING#   X(KM)      Y(KM)  #INTERP PTS  #DATA TIMES  NEAR_NODE  SITA"
!!$   DO I=1,N_ASSIM_TS
!!$     MAXEL = MAXLOC(TS_OBS(I)%X_WEIGHT,DIM=1)
!!$     WRITE(IPT,'(2X,I5,3X,F8.1,3X,F8.1,3X,I6,5X,I6,5X,I6,5X,F8.1)') &
!!$     I,TS_OBS(I)%X/1000.,TS_OBS(I)%Y/1000., &
!!$       TS_OBS(I)%N_INTPTS,TS_OBS(I)%N_TIMES,TS_OBS(I)%INTPTS(MAXEL),&
!!$       TS_OBS(I)%SITA
!!$   END DO
!!$   WRITE(IPT,*)
!!$   WRITE(IPT,*)'NUMBER OF BAD TS DATA POINTS: ',NBD_CNT
!!$   WRITE(IPT,*)" MOORING #   BEGIN TIME  END TIME"
!!$   DO I=1,N_ASSIM_TS
!!$   WRITE(IPT,*)I,TS_OBS(I)%TIMES(1)/(24.*3600.),&
!!$       TS_OBS(I)%TIMES(TS_OBS(I)%N_TIMES)/(24.*3600.)
!!$   END DO
!!$
!!$   RETURN
!!$   END SUBROUTINE SET_TS_ASSIM_DATA
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$   SUBROUTINE TEMP_NUDGING
!!$!==============================================================================|
!!$!  USE TEMP OBSERVATION DATA TO ADJUST TEMP FIELD                              |
!!$!==============================================================================|
!!$   USE MOD_PREC
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR 
!!$#  endif
!!$   IMPLICIT NONE
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TINT,TCORR,TG
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TWGHT_T
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
!!$   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
!!$   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
!!$   INTEGER I,J,K,J1,K1,K2,NLAY,ITIME,NTIME,IERR
!!$   INTRINSIC MINLOC
!!$!==============================================================================|
!!$
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Gather T and S Fields to Master Processor                                   ! 
!!$!------------------------------------------------------------------------------!
!!$   ALLOCATE(TG(MGL,KB))
!!$#  if defined (MULTIPROCESSOR)
!!$   IF(PAR)THEN
!!$     CALL GATHER(LBOUND(TF1,1),UBOUND(TF1,1),M,MGL,KB,MYID,NPROCS,NMAP,TF1,TG)
!!$   END IF
!!$#  endif
!!$   IF(SERIAL)THEN
!!$     TG(1:MGL,1:KBM1) = TF1(1:MGL,1:KBM1)
!!$   END IF
!!$!------------------------------------------------------------------------------!
!!$!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   IF(MSR)THEN
!!$   TS_OBS%T_WEIGHT = 0. 
!!$   T_THRESH         = ASTIME_WINDOW_TS    
!!$   SIMTIME          = TIME*86400
!!$
!!$   DO I=1,N_ASSIM_TS
!!$     NTIME = TS_OBS(I)%N_TIMES
!!$     ALLOCATE(FTEMP(NTIME)) 
!!$     FTEMP(1:NTIME) = ABS(SIMTIME - TS_OBS(I)%TIMES(1:NTIME))
!!$     DT_MIN = MINVAL(FTEMP(1:NTIME))
!!$     TS_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)
!!$
!!$     IF(DT_MIN < T_THRESH)THEN     
!!$       IF(DT_MIN < .5_SP*T_THRESH) THEN
!!$         TS_OBS(I)%T_WEIGHT = 1.0_SP
!!$       ELSE
!!$         TS_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
!!$       END IF
!!$     END IF
!!$
!!$     DEALLOCATE(FTEMP)
!!$   END DO
!!$   
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Interpolate Simulation Data to Local Observation Point                      ! 
!!$!------------------------------------------------------------------------------!
!!$       
!!$   ALLOCATE(TINT(N_ASSIM_TS,MAX_LAYER_TS)) ; TINT = 0. 
!!$
!!$   DO I=1,N_ASSIM_TS   
!!$     DO J=1,TS_OBS(I)%N_INTPTS
!!$       J1        = TS_OBS(I)%INTPTS(J)
!!$       WGHT      = TS_OBS(I)%X_WEIGHT(J)
!!$       NLAY      = TS_OBS(I)%N_LAYERS
!!$       DO K=1,NLAY
!!$         U1 = TG(J1,TS_OBS(I)%S_INT(K,1))
!!$         U2 = TG(J1,TS_OBS(I)%S_INT(K,2))
!!$         W1 = TS_OBS(I)%S_WEIGHT(K,1)
!!$         W2 = TS_OBS(I)%S_WEIGHT(K,2)
!!$         TINT(I,K) = TINT(I,K) + (U1*W1 + U2*W2)*WGHT 
!!$       END DO
!!$     END DO
!!$     TOT_WGHT = SUM(TS_OBS(I)%X_WEIGHT(1:TS_OBS(I)%N_INTPTS))
!!$     TINT(I,1:NLAY) = TINT(I,1:NLAY)/TOT_WGHT
!!$   END DO
!!$     
!!$!------------------------------------------------------------------------------!
!!$!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   ALLOCATE(TWGHT_T(MGL,KBM1))   ; TWGHT_T   = 0.
!!$   ALLOCATE(TCORR(MGL,KBM1))   ; TCORR   = 0.
!!$
!!$   DO I=1,N_ASSIM_TS 
!!$     DO J=1,TS_OBS(I)%N_INTPTS
!!$       J1     = TS_OBS(I)%INTPTS(J)
!!$       ITIME  = TS_OBS(I)%N_T_WEIGHT
!!$       NLAY   = TS_OBS(I)%N_LAYERS
!!$       DO K=1,NLAY
!!$         K1             = TS_OBS(I)%S_INT(K,1)
!!$         K2             = TS_OBS(I)%S_INT(K,2)
!!$         W1             = TS_OBS(I)%S_WEIGHT(K,1)
!!$         W2             = TS_OBS(I)%S_WEIGHT(K,2)
!!$!QXU{
!!$         DEFECT         = TS_OBS(I)%TEMP(ITIME,K) - TINT(I,K)
!!$	 IF(ABS(DEFECT).LT.10.0) then
!!$!QXU}	 
!!$           WEIGHT1        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W1
!!$           WEIGHT2        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W2
!!$           TWGHT_T(J1,K1) = TWGHT_T(J1,K1) + WEIGHT1   
!!$           TWGHT_T(J1,K2) = TWGHT_T(J1,K2) + WEIGHT2   
!!$!QXU         DEFECT         = TS_OBS(I)%TEMP(ITIME,K) - TINT(I,K)
!!$           CORRECTION     = GAMA_TS*DEFECT
!!$!QXU{
!!$!           TCORR(J1,K1)   = TCORR(J1,K1) + CORRECTION*WEIGHT1
!!$!           TCORR(J1,K2)   = TCORR(J1,K2) + CORRECTION*WEIGHT2
!!$
!!$           TCORR(J1,K1)   = TCORR(J1,K1) + CORRECTION*WEIGHT1**2
!!$           TCORR(J1,K2)   = TCORR(J1,K2) + CORRECTION*WEIGHT2**2
!!$
!!$	 ENDIF
!!$!QXU}	 
!!$!        GEOFF NEW
!!$         IF(TS_OBS(I)%TEMP(ITIME,K) < -90.) then
!!$            TWGHT_T(J1,K1) = 0.0_SP
!!$            TWGHT_T(J1,K2) = 0.0_SP
!!$         ENDIF
!!$       END DO
!!$     END DO
!!$   END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Nudge Simulation Data Using Local Corrections                               ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   DO I=1,MGL
!!$     DO K=1,KBM1
!!$       IF(DA_TS(I) == 1 .AND. TWGHT_T(I,K) > 1.0E-08)THEN
!!$         TG(I,K) = TG(I,K) + DTI*GALPHA_TS*TCORR(I,K)/TWGHT_T(I,K)
!!$       END IF
!!$     END DO
!!$   END DO
!!$
!!$   DEALLOCATE(TWGHT_T,TCORR,TINT)
!!$   END IF  !!MASTER
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Disperse New Data Fields to Slave Processors                                ! 
!!$!------------------------------------------------------------------------------!
!!$   IF(SERIAL)THEN
!!$     TF1(1:M,1:KBM1) = TG(1:M,1:KBM1)
!!$   END IF
!!$#  if defined (MULTIPROCESSOR) 
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$   CALL MPI_BCAST(TG,MGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$
!!$   IF(PAR)THEN
!!$     DO I=1,M
!!$       TF1(I,1:KBM1) = TG(NGID(I),1:KBM1)
!!$     END DO
!!$   END IF
!!$#  endif
!!$
!!$   DEALLOCATE(TG)
!!$          
!!$
!!$   RETURN
!!$   END SUBROUTINE TEMP_NUDGING
!!$!==============================================================================|
!!$!==============================================================================|
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$   SUBROUTINE SALT_NUDGING()
!!$!==============================================================================|
!!$!  USE SALT OBSERVATION DATA TO ADJUST SALINITY FIELD                          |
!!$!==============================================================================|
!!$   USE MOD_PREC
!!$   USE ALL_VARS
!!$#  if defined (MULTIPROCESSOR)
!!$   USE MOD_PAR 
!!$#  endif
!!$   IMPLICIT NONE
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: SINT,SCORR,SG
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: TWGHT_S
!!$   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: FTEMP
!!$   REAL(SP) :: WEIGHT,DEFECT,CORRECTION,DT_MIN,SIMTIME,T_THRESH,WGHT,TOT_WGHT
!!$   REAL(SP) :: U1,U2,V1,V2,W1,W2,WEIGHT1,WEIGHT2
!!$   INTEGER I,J,K,J1,K1,K2,NLAY,ITIME,NTIME,IERR
!!$   INTRINSIC MINLOC
!!$!==============================================================================|
!!$
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Gather S Fields to Master Processor                                         ! 
!!$!------------------------------------------------------------------------------!
!!$   ALLOCATE(SG(MGL,KB))
!!$#  if defined (MULTIPROCESSOR)
!!$   IF(PAR)THEN
!!$     CALL GATHER(LBOUND(SF1,1),UBOUND(SF1,1),M,MGL,KB,MYID,NPROCS,NMAP,SF1,SG)
!!$   END IF
!!$#  endif
!!$   IF(SERIAL)THEN
!!$     SG(1:MGL,1:KBM1) = SF1(1:MGL,1:KBM1)
!!$   END IF
!!$!------------------------------------------------------------------------------!
!!$!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   IF(MSR)THEN
!!$   TS_OBS%T_WEIGHT = 0. 
!!$   T_THRESH         = ASTIME_WINDOW_TS    
!!$   SIMTIME          = TIME*86400
!!$
!!$   DO I=1,N_ASSIM_TS
!!$     NTIME = TS_OBS(I)%N_TIMES
!!$     ALLOCATE(FTEMP(NTIME)) 
!!$     FTEMP(1:NTIME) = ABS(SIMTIME - TS_OBS(I)%TIMES(1:NTIME))
!!$     DT_MIN = MINVAL(FTEMP(1:NTIME))
!!$     TS_OBS(I)%N_T_WEIGHT = MINLOC(FTEMP,DIM=1)
!!$
!!$     IF(DT_MIN < T_THRESH)THEN     
!!$       IF(DT_MIN < .5_SP*T_THRESH) THEN
!!$         TS_OBS(I)%T_WEIGHT = 1.0_SP
!!$       ELSE
!!$         TS_OBS(I)%T_WEIGHT = (T_THRESH-DT_MIN)/T_THRESH*2.0_SP
!!$       END IF
!!$     END IF
!!$
!!$     DEALLOCATE(FTEMP)
!!$   END DO
!!$   
!!$       
!!$!------------------------------------------------------------------------------!
!!$!  Interpolate Simulation Data to Local Observation Point                      ! 
!!$!------------------------------------------------------------------------------!
!!$       
!!$   ALLOCATE(SINT(N_ASSIM_TS,MAX_LAYER_TS)) ; SINT = 0.
!!$
!!$   DO I=1,N_ASSIM_TS   
!!$     DO J=1,TS_OBS(I)%N_INTPTS
!!$       J1        = TS_OBS(I)%INTPTS(J)
!!$       WGHT      = TS_OBS(I)%X_WEIGHT(J)
!!$       NLAY      = TS_OBS(I)%N_LAYERS
!!$       DO K=1,NLAY
!!$         V1 = SG(J1,TS_OBS(I)%S_INT(K,1))
!!$         V2 = SG(J1,TS_OBS(I)%S_INT(K,2))
!!$         W1 = TS_OBS(I)%S_WEIGHT(K,1)
!!$         W2 = TS_OBS(I)%S_WEIGHT(K,2)
!!$         SINT(I,K) = SINT(I,K) + (V1*W1 + V2*W2)*WGHT 
!!$       END DO
!!$     END DO
!!$     TOT_WGHT = SUM(TS_OBS(I)%X_WEIGHT(1:TS_OBS(I)%N_INTPTS))
!!$     SINT(I,1:NLAY) = SINT(I,1:NLAY)/TOT_WGHT
!!$   END DO
!!$     
!!$!------------------------------------------------------------------------------!
!!$!  Compute Local Correction by Interpolating Observed/Computed Defect          ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   ALLOCATE(TWGHT_S(MGL,KBM1))   ; TWGHT_S   = 0.
!!$   ALLOCATE(SCORR(MGL,KBM1))   ; SCORR   = 0.
!!$
!!$   DO I=1,N_ASSIM_TS 
!!$     DO J=1,TS_OBS(I)%N_INTPTS
!!$       J1     = TS_OBS(I)%INTPTS(J)
!!$       ITIME  = TS_OBS(I)%N_T_WEIGHT
!!$       NLAY   = TS_OBS(I)%N_LAYERS
!!$       DO K=1,NLAY
!!$         K1             = TS_OBS(I)%S_INT(K,1)
!!$         K2             = TS_OBS(I)%S_INT(K,2)
!!$         W1             = TS_OBS(I)%S_WEIGHT(K,1)
!!$         W2             = TS_OBS(I)%S_WEIGHT(K,2)
!!$!QXU{
!!$         DEFECT         = TS_OBS(I)%SAL(ITIME,K) - SINT(I,K)
!!$	 IF(ABS(DEFECT).LT.10.0) THEN
!!$!QXU}	 
!!$           WEIGHT1        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W1
!!$           WEIGHT2        = TS_OBS(I)%T_WEIGHT*TS_OBS(I)%X_WEIGHT(J)*W2
!!$           TWGHT_S(J1,K1) = TWGHT_S(J1,K1) + WEIGHT1   
!!$           TWGHT_S(J1,K2) = TWGHT_S(J1,K2) + WEIGHT2   
!!$!QXU       DEFECT         = TS_OBS(I)%SAL(ITIME,K) - SINT(I,K)
!!$           CORRECTION     = GAMA_TS*DEFECT
!!$!JQI< (05/11/05)
!!$!          SCORR(J1,K1)   = SCORR(J1,K1) + CORRECTION*WEIGHT1
!!$!          SCORR(J1,K2)   = SCORR(J1,K2) + CORRECTION*WEIGHT2
!!$           SCORR(J1,K1)   = SCORR(J1,K1) + CORRECTION*WEIGHT1**2
!!$           SCORR(J1,K2)   = SCORR(J1,K2) + CORRECTION*WEIGHT2**2
!!$!JQI>	 
!!$
!!$!QXU{
!!$         ENDIF
!!$!QXU}	 
!!$!        GEOFF NEW
!!$         IF(TS_OBS(I)%SAL(ITIME,K)  < -90.) then
!!$           TWGHT_S(J1,K1) = 0.0_SP
!!$           TWGHT_S(J1,K2) = 0.0_SP
!!$         ENDIF
!!$       END DO
!!$     END DO
!!$   END DO
!!$
!!$!------------------------------------------------------------------------------!
!!$!  Nudge Simulation Data Using Local Corrections                               ! 
!!$!------------------------------------------------------------------------------!
!!$
!!$   DO I=1,MGL
!!$     DO K=1,KBM1
!!$       IF(DA_TS(I) == 1 .AND. TWGHT_S(I,K) > 1.0E-08)THEN
!!$         SG(I,K) = SG(I,K) + DTI*GALPHA_TS*SCORR(I,K)/TWGHT_S(I,K)
!!$       END IF
!!$     END DO
!!$   END DO
!!$
!!$   DEALLOCATE(TWGHT_S,SCORR,SINT)
!!$   END IF  !!MASTER
!!$   
!!$!------------------------------------------------------------------------------!
!!$!  Disperse New Data Fields to Slave Processors                                ! 
!!$!------------------------------------------------------------------------------!
!!$   IF(SERIAL)THEN
!!$     SF1(1:M,1:KBM1) = SG(1:M,1:KBM1)
!!$   END IF
!!$#  if defined (MULTIPROCESSOR) 
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$   CALL MPI_BCAST(SG,MGL*KB,MPI_F,0,MPI_COMM_WORLD,IERR)
!!$!QXU{
!!$   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!!$!QXU}
!!$
!!$   IF(PAR)THEN
!!$     DO I=1,M
!!$#  if defined (WET_DRY)
!!$       IF(ISWETNT(I)*ISWETN(I) == 1)THEN
!!$#  endif       
!!$       SF1(I,1:KBM1) = SG(NGID(I),1:KBM1)
!!$#  if defined (WET_DRY)
!!$       END IF
!!$#  endif              
!!$     END DO
!!$   END IF
!!$#  endif
!!$
!!$   DEALLOCATE(SG)
!!$          
!!$   RETURN
!!$   END SUBROUTINE SALT_NUDGING
!!$!==============================================================================|
!!$!==============================================================================|
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!!$
!!$

   SUBROUTINE SET_SSTGRD_ASSIM_DATA

!------------------------------------------------------------------------------!
!  SET UP ASSIMILATION DATA FOR SST OBSERVATIONS                               |
!------------------------------------------------------------------------------!
   USE ALL_VARS
   USE MOD_PAR 

   IMPLICIT NONE
   INTEGER I
   REAL(DP) :: TEMP      
   ! SOME NC POINTERS

   TYPE(NCATT), POINTER :: ATT
   TYPE(NCDIM), POINTER :: DIM
   TYPE(NCVAR), POINTER :: VAR
   
   ! SOME HANDY VARIABLES TO PLAY WITH
   LOGICAL FOUND  
   TYPE(TIME) :: FSTART, FEND,TCNT
   !------------------------------------------------------------------------------!
       
   ! GET THE FILE POINTER
   SST_FILE => FIND_FILE(FILEHEAD,trim(SSTGRD_ASSIM_FILE),FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("COULD NOT FIND SSTGRD_ASSIM_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE))
   
   ! CHECK THE ATTRIBUTES
   ATT => FIND_ATT(SST_FILE,"source",FOUND)
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SST_ASSIMGRD_FILE FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
   
   IF(ATT%CHR(1)(1:5) /= "FVCOM") CALL FATAL_ERROR &
        &('IN SSTGRD_ASSIM_FILE FILE OBJECT',&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        & "THE GLOBAL ATTRIBURE 'source' SHOULD BE 'FVCOM' ")
   
   ! LOOK FOR THE DIMENSIONS 
   DIM => FIND_DIM(SST_FILE,'node',FOUND) ! NODE
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSTGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'node'")
   
   IF (DIM%DIM /= MGL) CALL FATAL_ERROR&
        & ("THE NUMBER OF NODES IN THE SST GRID ASSIM FILE DOES NOT MATCH THE GRID FILE!")
   
   DIM => FIND_DIM(SST_FILE,'time',FOUND) ! TIME
   IF(.not. FOUND) CALL FATAL_ERROR &
        & ("IN SSTGRD_ASSIM_FILE OBJECT",&
        & "FILE NAME: "//TRIM(SSTGRD_ASSIM_FILE),&
        &"COULD NOT FIND DIMENSION 'time'")
   
   N_TIMES_SST = DIM%DIM ! SET THE VARIABLE VALUE
   
   
   
   ! CHECK THE BEGIN TIME
   FSTART = GET_FILE_TIME(SST_FILE,1)
   
   IF (FSTART > STARTTIME + MuSpD/2) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FSTART,IPT,"FIRST TIME IN SST GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(STARTTIME,IPT,"MODEL START TIME")
      END IF
      CALL FATAL_ERROR("THE SST GRID ASSIM START TIME IS INCORRECT",&
           & "THE MODEL STARTS MORE THAN 12 HOURS BEFORE THE FIRST SST DATA")
   END IF
   
   ! CHECK THE END TIME
   FEND = GET_FILE_TIME(SST_FILE,N_TIMES_SST)
   
   IF (FEND < ENDTIME - MuSpD/2) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN SST GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(ENDTIME,IPT,"MODEL END TIME")
      END IF
      CALL FATAL_ERROR("THE SST GRID ASSIM END TIME IS INCORRECT",&
           & "THE MODEL ENDS MORE THAN 12 HOURS AFTER THE LAST SST DATA")
   END IF
   
   ! CHECK THE FILE INTERVAL - EXPECTING DAILY DATA
   TCNT = FSTART
   DO I =1 , (N_TIMES_SST-1)
      TCNT = TCNT + MUSPD
   END DO
   
   IF (.NOT. FEND == TCNT) THEN
      IF(DBG_SET(DBG_LOG)) THEN
         CALL PRINT_REAL_TIME(FEND,IPT,"LAST TIME IN SST GRID ASSIM FILE")
         CALL PRINT_REAL_TIME(TCNT,IPT,"EXPECTED END TIME FOR DAILY DATA")
      END IF
      CALL FATAL_ERROR("THE SST GRID ASSIM FILE TIME STEP IS NOT DAILY",&
           & "THE START TIME PLUS NUMBER OF DAYS IN THE FILE DOES N&
           &OT EQUAL THE END TIME?")
   END IF
   
   ! GET INTERVAL FOR SAVING SST STATE 
   TEMP= 1.0_DP/REAL(SST_SAVE_N,DP)
   SST_SAVE_INTERVAL = anint(TEMP * REAL(MUSPD,DP),ITIME)

   ! ALLOCATE MEMORY
   ALLOCATE(SST_OBS(M))   
   ALLOCATE(SST_AVG(M))
   ALLOCATE(SST_SAVED(M,1:SST_SAVE_N))


   VAR => FIND_VAR(SST_FILE,'sst',FOUND)
   IF (.NOT.FOUND) CALL FATAL_ERROR &
        ("THE VARIABLE 'sst' IS MISSING FROM THE SSTGRD_ASSIM FILE?")

   CALL NC_CONNECT_AVAR(VAR,SST_OBS)
   VAR_SST =>VAR
   

   ! SET THE FILE STACK COUNT
   ! LOOK FOR TIME WITHIN ONE DTI OF MID-DAY

   FSTART = StartTime + MUSPD/2

   DO I =1 ,N_TIMES_SST

      TCNT = GET_FILE_TIME(SST_FILE,I)

      ! IF TIME IS BEFORE MID POINT THE FIRST DAY
      IF (TCNT < FSTART - IMDTI) CYCLE
      
      ! IF IT IS ALSO LESS THEN THIS
      IF (TCNT < FSTART + IMDTI) THEN
         ! WE FOUND IT
         SST_FILE%FTIME%NEXT_STKCNT = I
         EXIT

      ELSE ! THRE IS NO VALID TIME IN THE FILE

         CALL FATAL_ERROR &
              & ("THERE IS NO VALID FIRST TIME IN THE SST OBSERVATION FILE")
      END IF
   END DO


   RETURN
 END SUBROUTINE SET_SSTGRD_ASSIM_DATA

 SUBROUTINE SST_OBSERVATION_UPDATE
   IMPLICIT NONE

   CALL NC_READ_VAR(VAR_SST,SST_FILE%FTIME%NEXT_STKCNT)

   SST_FILE%FTIME%NEXT_STKCNT = SST_FILE%FTIME%NEXT_STKCNT + 1


   ! RESET THE INDEX AND THE TIME FOR USE IN THE ASSIMILATION CYCLE
   SST_SAVE_INDEX = 1
   SST_SAVE_TIME = IntTime_BUF + SST_SAVE_INTERVAL
   

 END SUBROUTINE SST_OBSERVATION_UPDATE



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|


  SUBROUTINE SST_NUDGE
  
!==============================================================================|
!  USE SST OBSERVATION DATA TO NUDGE SURFACE TEMPERATURE                       |
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   REAL(SP) :: TDIFF,T_WEIGHT,TRUTH,ADJUSTMENT
   INTEGER I

!==============================================================================|


!------------------------------------------------------------------------------!
!  Calculate Temporal Weight of Measurement (I) at Time(TIME)                  ! 
!------------------------------------------------------------------------------!
   IF(IntTime > SST_SAVE_TIME + SST_SAVE_INTERVAL/2) THEN
   
      SST_SAVE_INDEX = SST_SAVE_INDEX + 1

      SST_SAVE_TIME = SST_SAVE_TIME + SST_SAVE_INTERVAL
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Setting sst state # ",SST_SAVE_INDEX
!         CALL PRINT_REAL_TIME(SST_SAVE_TIME,IPT,"Set at")
      END IF

   END IF

   ! GET TIME DIFFERENCE IN SECONDS
   TDIFF = SECONDS(IntTime - SST_SAVE_TIME)
   TDIFF = ABS(TDIFF)

   IF(TDIFF < 0.5*SST_TWINDOW) THEN
     T_WEIGHT = 1.0
   ELSE IF(TDIFF < SST_TWINDOW) THEN
     T_WEIGHT = (SST_TWINDOW-TDIFF)/SST_TWINDOW*2.0
   ELSE
     T_WEIGHT = 0.0
   END IF
  
   IF(T_WEIGHT == 0.0) RETURN

   DO I=1,M
      !TRUTH = MODEL_SAVED@INTERVAL + OBSERVED_AVG  - MODEL_AVG
      TRUTH = SST_SAVED(I,SST_SAVE_INDEX) + SST_OBS(I) - SST_AVG(I)
      ! ADJUSMENT = MAX_WEIGHT * (TRUTH-CURRENT) * (WEIGHT_BY_TIME)
      ADJUSTMENT = SST_WM *(TRUTH-T1(i,1))*T_WEIGHT
      ! RESULT = CURRENT + DTI*TIMESCALE * ADJUSTMENT 
      T1(I,1)=T1(I,1)+DTI*SST_TSCALE*ADJUSTMENT
   ENDDO    

   
#     if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,T1) !Interprocessor Exchange   !
#     endif

   
   CALL N2E3D(T1,T)
   

   RETURN
   END SUBROUTINE SST_NUDGE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

   SUBROUTINE SST_SAVE     
!==============================================================================|
!  SAVE INTERVAL-AVERAGAED SST DATA ON GRID TO USE AS OBSERVATION              | 
!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR
   IMPLICIT NONE
   INTEGER I
   
!==============================================================================|

   IF(IntTime >= SST_SAVE_TIME) THEN
   
      SST_SAVE_INDEX = SST_SAVE_INDEX + 1
      
      IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*) "Saving sst state # ",SST_SAVE_INDEX
      END IF
!------------------------------------------------------------------------------!
!  Save Hourly SST on grid in Simulation running                               ! 
!------------------------------------------------------------------------------!
      SST_SAVED(1:M,SST_SAVE_INDEX) = T1(1:M,1)

      ! INCRIMENT THE SAVE_TIME
      SST_SAVE_TIME = SST_SAVE_TIME + SST_SAVE_INTERVAL

!------------------------------------------------------------------------------!
!  Average Hourly Data Over Integration Period                                 ! 
!------------------------------------------------------------------------------!

      IF(SST_SAVE_INDEX == SST_SAVE_N )THEN
      
         DO I=1,M
            SST_AVG(I) = SUM(SST_SAVED(I,1:SST_SAVE_N))/FLOAT(SST_SAVE_N)
         END DO
         
      END IF

   END IF
   
   RETURN
 END SUBROUTINE SST_SAVE

!!$# else
!!$  IMPLICIT NONE
!!$
!!$  CONTAINS
!!$
!!$    SUBROUTINE LOAD_ASSIM_DATA
!!$      IMPLICIT NONE
!!$      CALL FATAL_ERROR("LOAD_ASSIM_DATA: THIS IS A DUMMY PROGRAM- YOU SHOULD NOT BE HERE")
!!$    END SUBROUTINE LOAD_ASSIM_DATA
!!$
!!$    SUBROUTINE SET_ASSIM_INTERVALS
!!$      IMPLICIT NONE
!!$      CALL FATAL_ERROR("SET_ASSIM_INTERVALS: THIS IS A DUMMY PROGRAM- YOU SHOULD NOT BE HERE")
!!$
!!$    END SUBROUTINE SET_ASSIM_INTERVALS
!!$
!!$#  endif
END MODULE MOD_ASSIM
