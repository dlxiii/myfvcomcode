! $Id: fvcom.F,v 1.1.2.29 2008/06/09 22:37:09 gao Exp $
! $Name: New_Input $
! $Revision: 1.1.2.29 $
!==============================================================================!
!  VERSION 3.0 Alpha 
!==============================================================================!

PROGRAM FVCOM

! CODING STATUS: Stable but many features still not activated
! Remaining Work:
! Add MPI calls for derived type TIME
! Finish lagrangian coding
! Add no limit coding
! Add non hydrostatic
! Add variable, 1/ART, 1/ART2... etc
! Fix bottom boundary condition for scalars
! Fix bottom boundary roughness
! Add semi-implicit coding
! Add sediment 
! Test 2d model
! Test 1d model
! Finish Ice model
! Add assimilation modes
! Add mean flow
! Add dye release
! Test momentum balance
! Test GOTM Interface
!
! Biological model ?




  !==============================================================================!
  !                                                                              !
  !                            USG-FVCOM                                         !
  !    The Unstructured Grid Finite Volume Coastal Ocean Model                   !
  !                                                                              !
  !    The USG-FVCOM (publically called FVCOM) was developed by Drs. Changsheng  !
  !  Chen and Hedong Liu at the Marine Ecosystem Dynamics Modeling Laboratory    !
  !  at the School of Marine Science and Technology (SMAST), University of       !
  !  Massachusetts-Dartmouth (UMASSD) and Dr. Robert C. Beardsley at the         !
  !  Department of Physical Oceanography, Woods Hole Oceanographic Institution   !
  !  (WHOI). This code was rewritten in Fortran 90/2K, modularized, rendered     ! 
  !  somewhat understandable, and parallelized by Geoff Cowles at SMAST/UMASSD.  !
  !                                                                              !
  !    The Development was initially supported by the Georgia Sea Grant College  !
  !  Program for the study of the complex dynamic system in Georgia estuaries.   !
  !  The code improvement has been supported by Dr. Chen's research grants       !
  !  received from NSF and NOAA Coastal Ocean research grants and SMAST-NASA     !
  !  fishery research grants. The University of Massachusetts-Dartmouth          !
  !                                                                              !
  !    FVCOM is a three dimensional,time dependent,primitive equations,          !
  !  coastal ocean circulation model. The model computes the momentum,           !
  !  continuity, temperature, salinity, and density equations and is closed      !
  !  physically and mathematically using the Mellor and Yamada level-2.5         !
  !  turbulent closure submodel. The irregular bottom slope is represented       !
  !  using a sigma-coordinate transformation, and the horizontal grids           !
  !  comprise unstructured triangular cells. The finite-volume method (FVM)      !
  !  used in this model combines the advantages of the finite-element            !
  !  method (FEM) for geometric flexibility and the finite-difference            !
  !  method (FDM) for simple discrete computation. Current, temperature,         !
  !  and salinity in the model are computed in the integral form of the          !
  !  equations, which provides a better representation of the conservation       !
  !  laws for mass, momentum, and heat in the coastal region with complex        !
  !  geometry.                                                                   !
  !                                                                              !
  !    All users should read this agreement carefully.  A user, who receives any !  
  !  version of the source code of FVCOM, must accept all the terms and          !
  !  conditions of this agreement and also agree that this agreement is like any !
  !  written negotiated agreement signed by you. You may be required to have     !
  !  another written agreement directly with Dr. Changsheng Chen at SMAST/UMASS-D!
  !  or Dr. Brian Rothschild, Director of the SMAST/UMASS-D that supplements     !
  !  all or portions of this agreement. Dr. Changsheng Chen, leader of the       !
  !  FVCOM development team, owns all intellectual property rights to the        !
  !  software. The University of Massachusetts-Dartmouth and the Georgia Sea     !
  !  Grant Program share the copyright of the software. All copyrights are       !
  !  reserved. Unauthorized reproduction and re-distribution of this program     !
  !  are expressly prohibited. This program is only permitted for use in         !
  !  non-commercial academic research and education.  Commercial use must be     !
  !  approved by Dr. Chen and is subject to a license fee. Registration is       !
  !  required for all new users.  Users should realize that this model software  !
  !  is a research product without any warranty. Users must take full            !
  !  responsibility for any mistakes that might be caused by any inappropriate   !
  !  modification of the source code.  Modification is not encouraged for users  !
  !  who do not have a deep understanding of the finite-volume numerical methods !
  !  used in FVCOM. Contributions made to correcting and modifying the programs  !
  !  will be credited, but will not affect copyrights. No duplicate              !
  !  configurations of FVCOM are allowed in the same geographical region,        !
  !  particularly in the regions where FVCOM has been already been applied.      !
  !  Users who want to use FVCOM in a region that the SMAST/UMASS Marine         !
  !  Ecosystem Dynamics Modeling (MEDM) group (led by Dr. Chen) is working on    !
  !  must request permission from Dr. Chen. No competition is allowed in the     !
  !  same region using FVCOM, especially with Dr. Chen's group. FVCOM has been   !
  !  validated for many standard model test cases.  Users are welcome to do any  !
  !  further model validation experiments. These experiments shall be carried    !
  !  out in collaboration with the SMAST/UMASSD model development team. To avoid !
  !  or reduce deriving any incorrect conclusions due to an inappropriate use of !
  !  FVCOM, users are required to contact the scientific leaders of the FVCOM    !
  !  development team (Dr. Chen at SMAST/UMASS-D and Dr. Beardsley at WHOI)      !
  !  before any formal publications are prepared for model validation.           !
  !                                                                              !
  !    For public use, all users should name this model as "FVCOM". In any       !
  !  publications with the use of FVCOM, acknowledgement must be included. The   !
  !  rationale behind this FVCOM distribution policy is straightforward.  New    !
  !  researchers and educators who want to use FVCOM and agree to the above      !
  !  requirements get free access to the latest version of FVCOM and the         !
  !  collective wisdom and experience of the FVCOM development team and existing !
  !  users. Problems arising in new FVCOM applications, both related to          !
  !  conceptual as well as numerical and coding issues, can be shared with the   !
  !  development team and other users who can work together on physics and code  !
  !  improvements that over time will lead to a better FVCOM.                    !
  !                                                                              !
  !    FVCOM has been developed to date with state and federal funding with the  !
  !  idea that FVCOM will become a community model that new users start to use   !
  !  the model and its scientific usefulness and numerical accuracy and          !
  !  efficiency continue to improve.  The FVCOM distribution policy is designed  !
  !  to encourage this transition while maintaining a central core group         !
  !  responsible for overall FVCOM development direction, implementing official!
  !  code improvements, and maintaining well tested and documented updated code  !
  !  versions.                                                                   !       
  !                                                                              !
  !                                                                              !
  !  External forces used to drive this model:                                   !
  !                                                                              !
  !  1) Tidal amplitudes and phases at open boundaries (initial designs          !
  !         include 6 tidal consituents, more can be added as needed);           !
  !  2) Wind Stress [3 ways: a) uniform wind speed and direction, b)             !
  !         spatially distributed wind velocity field, and c)the MM5 model-out   !
  !         wind fields]                                                         !
  !  3) Surface heat flux [3 ways: a) uniform heat flux, b) spatially            !
  !         distributed heat flux, and c) the MM5-output heat flux fields        !
  !         All the surface net heat flux and short-wave radiation are needed    !
  !         in the input file                                                    ! 
  !  4) River discharges: specify the location and discharge volume,             !
  !         temperature, and salinity                                            !
  !  5) Groundwater input: currently diffused bottom flux only                   !
  !                                                                              !
  !  Initial conditions:                                                         !
  !                                                                              !
  !  The model can be prognostically run for both barotropic and baroclinic      !
  !  cases.                                                                      !
  !                                                                              !
  !  Tidal forcing can be added into the system with zero velocity               !
  !  field at initial or specified the 3-D tidal initial velocity field          !
  !  using the model-predicted harmonic tidal currents.                          !
  !                                                                              !
  !  Initial fields of temperature and salinity needed to be specified           !
  !  by using either climatological field, real-time observed field or           !
  !  idealized functions. The model has included Gregorian time for the          !
  !  time simulation for tidal currents.                                         !
  !                                                                              !
  !  For the purpose of interdisciplinary studies, biological, chemical, and     !
  !  sediment suspension models are available for FVCOM.  These submodels are    !
  !  directly driven by the FVCOM physical model. A description of these         !
  !  submodels follows.                                                          !
  !                                                                              !
  !  Generalized biological modules-a software platform that allows users to     !
  !             build his own biological model in FVCOM                          !
  !                                                                              !
  !  NPZ model--a 3 component nutrient-phytoplankton-zooplankton model           !
  !                                                                              !
  !  NPZD model--an 8 component nutrient-phytolankton-zooplankton-detritus       !
  !              model;                                                          !
  !                                                                              !
  !  NPZDB-model-a 9 phosphorus-controlled component nutrient-                   !
  !               phytoplankton-zooplankton-detritus-bacteria model;             !
  !                                                                              !
  !  Water quality model with inclusion of the benthic process                   !
  !                                                                              !
  !  Sediment model--A new module that was developed by Dr. Geoff Cowles         !
  !                                                                              !
  !  Lagrangian particle tracking:                                               !
  !                                                                              !
  !  A bilinear interpolation scheme is used to determine the particle           !
  !  velocity for the Lagrangian particle tracking. A random walk process        !
  !  also could be included with a specified function related to horizontal      !
  !  and vertical diffusion coefficients                                         !
  !                                                                              !
  !  Key reference:                                                              !
  !                                                                              !
  !   Chen, C., H. Liu, and R. C. Beardsley, 2003. An unstructured grid,         !
  !       finite-volume, three-dimensional, primitive equations ocean            !
  !       model: application to coastal ocean and estuaries, Journal             !
  !       of Atmospheric and Oceanic Technology,  20, 159-186.                   !
  !                                                                              !
  !                                                                              !
  !                                                                              !
  !  Please direct criticisms and suggestions to                                 !
  !                                                                              !
  !               Changsheng Chen                                                !
  !               School for Marine Science and Technology                       !
  !               University of Massachusetts-Dartmouth                          !
  !               New Bedford, MA 02742                                          !
  !               Phone: 508-910-6388, Fax: 508-910-6371                         !
  !               E-mail: c1chen@umassd.edu                                      !
  !               Web: http://fvcom.smast.umassd.edu                             !
  !                                                                              !
  ! What are new for version 2.7 Beta?                                           !
  !    1) Multiple choice to set up the coordinate in the vertical is developed  !
  !       by Dr. Qi                                                              !
  !       a) sigma levels                                                        !
  !       b) general vertical levels                                             ! 
  !       c) constant layer transformation                                       !
  !                                                                              !
  ! What are new for version 2.6?                                                !
  !    1) A new Lagrangian particle tracking module is added with multiprocessor !
  !       formulation, restart capability and general tracking et al.            !
  !       (by Dr. Cowles)                                                        !
  !                                                                              !
  ! What are new for version 2.5?                                                !
  !    1) A new spherical coordinate version is added with an accurate treatment !
  !       of the north pole (Arctic Ocedan) by Dr. Chen et al. at UMASSD         !
  !    2) Spherical and x-y coordinate versions was merged into a single code    !
  !       with users' choice for either coordiante for their application         !
  !    3) Multiple choices of open radiation boundary conditions are added       !
  !    4) General turbulence modules were linked by Dr. Cowles                   !
  !    5) The selection for a 2-D barotrophic application is added               !
  !    6) bugs in paralleziation and wet/dry point treatments are corrected      !
  ! For more detailed information, please read the upgrade user manaul           !
  !                                                                              !
  ! What will be included in version 2.4 (under validation tests)                !
  !    1) Generalized Biological Modules (developed by Dr. Tian and Chen)        !
  !    2) 3-D sediment model(developed by Dr. Cowles)                            !
  !    3) Reduced Kalman Filter, Ensemble Kalman Filter, and Ensemble Tranistion !
  !       Kalman Filter (implented into FVCOM with UMASSD and MIT groups led     !
  !       by Chen and Roziili                                                    !
  !    4) Full nonlinear ice models (developed orignally to FVCOM by Dr.Dupont   !
  !       and modified to be implemented to a parallelized new spherical         !
  !       coordinate version by the MEDM group at UMASSD                         !
  !                                                                              !   
  ! Enjoy!                                                                       !
  !==============================================================================!

  !==============================================================================!
  !  INCLUDE MODULES                                                             !
  !==============================================================================!

  USE MOD_UTILS
  USE CONTROL

  USE MOD_PAR  

  USE MOD_STARTUP

  USE MOD_TIME
  USE MOD_CLOCK

  USE MOD_INPUT
  USE MOD_NCDIO
  USE MOD_NCLL

  USE MOD_SETUP
  USE MOD_SET_TIME

  USE MOD_FORCE
  USE MOD_OBCS

  USE MOD_ASSIM
  USE MOD_NESTING

  USE MOD_REPORT
  USE MOD_LAG
  USE PROBES
# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT
# endif
# if defined (BALANCE_2D)
  USE MOD_BALANCE_2D
# endif  
# if defined (ONE_D_MODEL)
  USE MOD_ONEDTIDE
# endif
# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif    
# if defined (BioGen)
  USE MOD_BIO_3D
# endif    
!--------------------------------------------------------------
# if defined (ICE)
      use mod_ice
# endif
!--------------------------------------------------------------
# if defined (NH)
  USE NON_HYDRO
# endif

# if defined (SEMI_IMPLICIT) || (NH)
  USE MOD_PETSc, ONLY : PETSc_SET, PETSc_CLEANUP
# endif
  !------------------------------------------------------------------------------|
  IMPLICIT NONE

  integer  :: i, i1, slice_node

  character(len=*),parameter::CVS_Id="$Id: fvcom.F,v 1.1.2.29 2008/06/09 22:37:09 gao Exp $" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date: 2008/06/09 22:37:09 $" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name: New_Input $" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision: 1.1.2.29 $" ! [sng] File revision string

# if defined (RRK_PRE)
  CHARACTER(LEN=120):: RRKFILES
  CHARACTER(LEN=8)  :: RRKINP1
  CHARACTER(LEN=4)  :: RRKINP2
# endif

  type(watch) Timer


  !==============================================================================!
  ! INITIALIZE ALL CONTROL VARIABLES
  !==============================================================================!
  CALL INITIALIZE_CONTROL("FVCOM")

#  if defined (MULTIPROCESSOR)
  ! INTIALIZE MPI CONTROL VARIABLES
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
#  endif

  !==============================================================================!
  !   INITIALIZE A STOP WATCH TIMER FOR TESTING SUBROUTINE EFFICENCY             !
  !==============================================================================!
  CALL WATCH_INIT(TIMER)

  !==============================================================================!
  !   IMPORT CASENAME AND COMMAND LINE ARGUMENTS AND START LOG FILE              !
  !==============================================================================!
  CALL COMMANDLINEIO(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)       
  if(DBG_SET(dbg_log)) Call WRITE_BANNER(PAR,NPROCS,MYID)

  !==============================================================================!
  ! SET DEFAULT VALUES AND READ NAME LISTS                                            
  !==============================================================================!

  CALL NAMELIST

  !==============================================================================!
  !   SET MODEL CONTROL PARAMTERS BASED ON NAME LIST HERE                        !
  !==============================================================================!
  CALL CNTRL_PRMTRS

  !==============================================================================!
  !   SET THE STARTUP TYPE TO BE USED!                                           !
  !==============================================================================!
  CALL SET_STARTUP_TYPE ! see: startup_type.F

  !==============================================================================!
  !   OPEN ALL FILES NEEDED BASED ON THE RUN PARAMETERS                          !
  !==============================================================================!
  CALL OPEN_ALL

  !==============================================================================!
  !   SET MODEL TIME BASED ON THE NAMELIST TIME STRINGS OR RESTART FILE          !
  !==============================================================================!
  CALL SETUP_TIME

  !==============================================================================!
  !   LOAD GRID CONNECTIVITY AND OBC LIST FOR METIS DECOMPOSITION                !
  !==============================================================================!
  CALL LOAD_GRID

  !==============================================================================!
  !   SETUP THE DOMAIN FOR PARALLEL OR SERIAL RUNNING                            !
  !==============================================================================!
  CALL SETUP_DOMAIN

  !==============================================================================!
  !   ALLOCATE ALL DOMAIN SIZE VARIABLES HERE                                    !
  !==============================================================================!
  CALL ALLOCATE_ALL

  !==============================================================================!
  !   LOAD/SETUP PHYSICAL QUANTITIES (CORIOLIS, GRAVITY, SPONGE LAYER, XY/LATLON)!
  !==============================================================================!
  CALL COORDS_N_CONST

  !==============================================================================!
  ! CALCULATE GRID METRICS - NEIGHBORS, GRADIENTS, CELL AREA, INTERP COEFF'S     !
  !==============================================================================!
  CALL GRID_METRICS

  !==============================================================================!
  !  SETUP THE MODEL FORCING                                                     !
  !==============================================================================!
  CALL SETUP_FORCING

  !==============================================================================!
  !  GET THE PARAMETERS OF BIOLOGICAL MODEL                                      !
  !==============================================================================!
#  if defined (BioGen)
   KBV=KB
   CALL GET_PARAMETER
!   RIVER_TS_SETTING = 'NONE' !No bio river yet
#  endif

  !==============================================================================!
  !  SETUP OTHER TOOLS, MODELS AND DATA ASSIMILATION                             !
  !==============================================================================!

# if defined (ICE)
  IF(ICE_MODEL) CALL ICE_INIT_0
# endif

  IF(NESTING_ON) CALL SETUP_NEST

  !==============================================================================!
  !  SETUP PETSc FOR SEMI_IMPLICIT AND NON-HYDROSTATIC MODULE                    !
  !==============================================================================!
# if defined (SEMI_IMPLICIT) || (NH)
  CALL PETSc_SET
# endif

  !==============================================================================!
  !  SET THE INITIAL CONDITIONS FOR THE MODEL                                    !
  !==============================================================================!
  CALL STARTUP

  !==============================================================================!
  !  CALL ARCHIVE TO SETUP THE OUTPUT AND DUMP CONSTANT VALUES                   !
  !==============================================================================!
  CALL ARCHIVE

  IF(NESTING_ON .and. NESTING_MODE == "large_domain") CALL ARCHIVE_NEST

  CALL SET_PROBES(PROBES_ON,PROBES_NUMBER,PROBES_FILE)

  CALL SET_LAG
  
  !==============================================================================!
  !  SELECT THE RUN MODE AND EXECUTE THE MAIN LOOP
  !==============================================================================!
  SELECT CASE(FVCOM_RUN_MODE)
     ! RUN MODE SET IN mod_assim.F(set_assim_param)

     ! =============================================================================!
     ! == PURE SIMULATION MODE - Instantanious data assimilation only ==============!
  CASE(FVCOM_PURE_SIM)
     ! =============================================================================!


     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "============== STARTING MAIN LOOP AT:==============="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     !////////////////////////// MAIN LOOP //////////////////////////////////////////
     DO IINT=ISTART,IEND

        IntTime=IntTime + IMDTI

        CALL INTERNAL_STEP

        !==============================================================================!
        !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
        !==============================================================================!
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        
        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

        !==============================================================================!
        !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
        !==============================================================================!
        CALL ARCHIVE

        CALL DUMP_PROBE_DATA
        !==============================================================================!
        !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
        !==============================================================================!
        CALL SHUTDOWN_CHECK(D1)

        !==============================================================================!
        !    LAGRANGIAN PARTICLE TRACKING                                              |
        !==============================================================================!
        CALL LAG_UPDATE

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
        IF(NESTING_ON .and. NESTING_MODE == "large_domain") CALL ARCHIVE_NEST

     END DO
     !////////////////////////// END MAIN LOOP //////////////////////////////////////

  ! ================================================================================!
  ! == THIS MAIN LOOP USES NUDGING METHODS TO ASSIMILATE ===========================!
  CASE(FVCOM_NUDGE_AVG_SST)
  ! ================================================================================!

     ! SET UP ASSIMILATION VARIABLES
     CALL LOAD_ASSIM_DATA

     CALL SET_SSTGRD_ASSIM_DATA

     CALL ALLOC_BUFFER

     ! SET THE ASSIMILATION/SIMULATION RESET TIME
     ASSIM_RESTART_TIME = StartTime

     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
     INT_START = ISTART
     INT_COUNT = ISTART
     INT_END   = ISTART + 2*(IEND-ISTART)

     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
        Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
        write(IPT,*) "===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')
     
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
        write(IPT,*) "===================================================="
     end if

     DO WHILE(IntTime < EndTime)

        ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + MuSpD  ! ADD ONE MORE DAY
        SST_SAVE_INDEX = 0
        SST_SAVE_TIME  = IntTime + sst_save_interval

        CALL SST_SAVE_STATE


        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Simulation  ===================="
        end if
        !==============================================================================!
        !    RUN PURE SIMULATION MODE:
        !==============================================================================!

        ASSIM_FLAG = 0
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)

        END DO


        !==============================================================================!
        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
        !==============================================================================!
        CALL SST_RESTORE_STATE
        CALL SST_OBSERVATION_UPDATE


        if(DBG_SET(dbg_log)) THEN
           Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
           write(IPT,*) "======= Start 1 Day Assimilation  ===================="
        end if

        !==============================================================================!
        !    RUN DATA ASSIMILATION MODE:
        !==============================================================================!
        ASSIM_FLAG = 1
        DO WHILE(IntTime < ASSIM_RESTART_TIME)

           IINT = IINT + 1
           IntTime=IntTime + IMDTI

           CALL INTERNAL_STEP

           !==============================================================================!
           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
           !==============================================================================!
           if(DBG_SET(dbg_log)) &
                & Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)

           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

           !==============================================================================!
           !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
           !==============================================================================!
           CALL ARCHIVE

           CALL DUMP_PROBE_DATA
           !==============================================================================!
           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
           !==============================================================================!
           CALL SHUTDOWN_CHECK(D1)


           !==============================================================================!
           !    LAGRANGIAN PARTICLE TRACKING                                              |
           !==============================================================================!
           CALL LAG_UPDATE

        END DO

        ! RUN THE NEXT DAY

     END DO


  CASE(FVCOM_KALMAN_1)
  CASE(FVCOM_KALMAN_2)
  CASE(FVCOM_KALMAN_3)

  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN FVCOM_RUN_MODE :'"//TRIM(FVCOM_RUN_MODE),&
          & "Options are the following: '"//TRIM(FVCOM_PURE_SIM)//"&
          &' or '"//TRIM(FVCOM_NUDGE_AVG_SST) )
  END SELECT

#  if defined (SEMI_IMPLICIT) || (NH)
   CALL PETSc_CLEANUP
#  endif
  
  if(DBG_SET(dbg_log)) write(IPT,*)"TADA!"
  CALL PSHUTDOWN

END PROGRAM FVCOM
