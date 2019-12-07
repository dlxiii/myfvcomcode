! $Id: mod_ncdio.F,v 1.20.2.36 2008/06/09 22:26:11 gao Exp $
! $Name: New_Input $
! $Revision: 1.20.2.36 $
MODULE MOD_NCDIO
  !==============================================================================!
  !  NetCDF Io for FVCOM using CF Metadata Convention                            !
  !                                                                              !
  !    see: http://www.cgd.ucar.edu/cms/eaton/cf-metadata/ for info              !
  !                                                                              !
  !    current time dependent variables set up                                   !
  !         el:    surface elevation                                             !
  !          u:    x-velocity. In spherical coordinate,lon-velocity              !                         
  !          v:    y-velocity. In spherical coordinate,lat-velocity              !                        
  !         ww:    z-velocity                                                    !
  !         kh:    turbulent diffusivity                                         !
  !         km:    turbulent viscosity                                           !
  !         t1:    temperature                                                   !
  !         s1:    salinity                                                      !
  !         ua:    vertically-averaged x-velocity                                !
  !                In spherical coordinate,vertically-averaged lon-velocity      !
  !         va:    vertically-averaged y-velocity                                !
  !                In spherical coordinate,vertically-averaged lat-velocity      !
  !          d:    depth at procs                                                !
  !        dye:    dye at procs                                                  !
  !       aice:    ice concentration on procs                                    !
  !       vice:    ice thichness on procs                                        !
  !      uuice:    ice x-velocity                                                !
  !      vvice:    ice y-velocity                                                !
  !     uuwind:    wind speed in x direction                                     !
  !     vvwind:    wind speed in y direction                                     !
  !                                                                              !
  !       wd:      wet/dry flag (0 or 1)                                         !
  !                                                                              !
  !       vort:    vorticity                                                     !
  !    to add additional variables:                                              !
  !      1.) add to list above                                                   !
  !      2.) add *_vid to variables vid in section "new variable vid"            !
  !      3.) go to definition section "new variable definition"                  !
  !      4.) add io section "new variable io"                                    !
  !==============================================================================!

  USE ALL_VARS
  USE MOD_PREC
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_TIME
  USE MOD_INPUT
  
#  if defined (SEDIMENT)
  USE mod_sed
#  endif

#  if defined (ICE)
  USE mod_ICE
  USE mod_ICE2D
#  endif
# if defined (BioGen)
  USE MOD_BIO_3D
# endif    
# if defined (WATER_QUALITY)
  USE MOD_WQM, ONLY : NB,WQM,WQM_NAME  
# endif  
# if defined (NH)
  USE NON_HYDRO, ONLY: W4ZT, RHS, QP, NHQDRX, NHQDRY, NHQDRZ, NHQ2DX, NHQ2DY
# endif
  implicit none

  Character(LEN=50) :: CoordVar 

  LOGICAL, private :: FOUND
  logical, private :: NEED_INIT = .TRUE.

  TYPE(NCDIM), POINTER :: DIM_nele
  TYPE(NCDIM), POINTER :: DIM_node
  TYPE(NCDIM), POINTER :: DIM_three
  TYPE(NCDIM), POINTER :: DIM_four

  TYPE(NCDIM), POINTER :: DIM_siglay
  TYPE(NCDIM), POINTER :: DIM_siglev

  TYPE(NCDIM), POINTER :: DIM_time
  TYPE(NCDIM), POINTER :: DIM_DateStrLen

  TYPE(NCDIM), POINTER :: DIM_nobc
  TYPE(NCDIM), POINTER :: DIM_nlsf

  TYPE(NCDIM), POINTER :: DIM_MaxNode
  TYPE(NCDIM), POINTER :: DIM_MaxElem
#  if defined (ICE)
  TYPE(NCDIM), POINTER :: DIM_GRID
  TYPE(NCDIM), POINTER :: DIM_NCAT
  TYPE(NCDIM), POINTER :: DIM_ntilay
#  endif


  save
  
CONTAINS
!=============================================================  
  SUBROUTINE ARCHIVE
    IMPLICIT NONE
    INTEGER :: STATUS
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ARCHIVE"

    if(NEED_INIT) then
       IF(USE_MPI_IO_MODE) THEN
          CALL MPI_IO_SYNCHRONIZE(INIT_CODE)
       ELSE
          CALL CALL_FUNC(INIT_CODE,status)
          IF (status/=0) call fatal_error("ARCHIVE:: Bad INIT_CODE",&
               & "Could not retrieve valid function pointer?")
       END IF

       NEED_INIT = .FALSE.

    end if

    IF(NC_ON)THEN
       IF(NC_DAT%FTIME%NEXT_IO <= IntTime) THEN
        
          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(NC_CODE)
          ELSE
             CALL CALL_FUNC(NC_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad NC_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF
       END IF
    END IF

    
    IF(NCAV_ON)THEN
       IF(NC_AVG%FTIME%NEXT_IO <= IntTime) THEN
          
          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(NCAV_CODE)
          ELSE
             CALL CALL_FUNC(NCAV_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad NCAV_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF
       END IF
    END IF


    IF(RST_ON)THEN
       IF(NC_RST%FTIME%NEXT_IO <= IntTime) then
          
          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(RESTART_CODE)
          ELSE
             CALL CALL_FUNC(RESTART_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad RESTART_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF
       END IF
    END IF
        
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ARCHIVE"
  END SUBROUTINE ARCHIVE
!=============================================================  
  SUBROUTINE INIT_NCDIO
    IMPLICIT NONE

    IF(DBG_SET(DBG_LOG)) &
         & write(IPT,*)"! SETTING UP NCDIO: CREATING AND DUMPING OUTPUT FILE META DATA"

    CALL DEFINE_DIMENSIONS

#  if defined (SPHERICAL)
    CoordVar="lat lon"
# else
    CoordVar="x y"
#endif


    IF(NC_ON) then
       
       IF(.not. ASSOCIATED(NC_DAT)) Call Fatal_Error &
            & ("INIT_NCDIO: THE DATA FILE OBJECT IS NOT ASSOCIATED ")
           
       CALL SETUP_DATFILE


    END IF
  
    IF(NCAV_ON) then
       
       IF(.not. ASSOCIATED(NC_AVG)) Call Fatal_Error &
            & ("INIT_NCDIO: THE AVERAGE FILE OBJECT IS NOT ASSOCIATED ")


       CALL SETUP_AVGFILE
    END IF

    IF(RST_ON) then


       IF(.not. ASSOCIATED(NC_RST)) Call Fatal_Error &
            & ("INIT_NCDIO: THE RESTART FILE OBJECT IS NOT ASSOCIATED ")

       CALL SETUP_RSTFILE
    END IF

  END SUBROUTINE INIT_NCDIO
!=============================================================  
! NEED AN INTERFACE WITH NO ARGS FOR FUNCTION POINTERS
!=============================================================  
  SUBROUTINE DUMP_NC_DAT
    IMPLICIT NONE

    IF (ICING_MODEL .AND. .NOT. IOPROC) CALL ICING(IntTime)

    CALL DUMP_DATA(NC_DAT)
  END SUBROUTINE DUMP_NC_DAT
  !=============================================================  
  SUBROUTINE DUMP_NC_AVG
    IMPLICIT NONE
    CALL DUMP_DATA(NC_AVG)
  END SUBROUTINE DUMP_NC_AVG
  !=============================================================  
  SUBROUTINE DUMP_NC_RST
    IMPLICIT NONE
    CALL DUMP_DATA(NC_RST)
  END SUBROUTINE DUMP_NC_RST
!=============================================================  
  SUBROUTINE SETUP_DATFILE
    IMPLICIT NONE
    character(len=80) :: tmp
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_DATAFILE"

    ! if USE_MPIO

    NC_DAT => ADD(NC_DAT,GRID_FILE_OBJECT() )
    
    IF(NC_FILE_DATE) THEN
       NC_DAT => ADD(NC_DAT,FILE_DATE_OBJECT() )
    END IF

    IF(NC_GRID_METRICS) THEN
       NC_DAT => ADD(NC_DAT,GRID_METRICS_FILE_OBJECT() )
    END IF

    IF(NC_VELOCITY) THEN
       NC_DAT => ADD(NC_DAT,VELOCITY_FILE_OBJECT() )
    END IF

    IF(NC_VERTICAL_VEL) THEN
       NC_DAT => ADD(NC_DAT,VERTICAL_VEL_FILE_OBJECT() )
    END IF

    IF(NC_AVERAGE_VEL) THEN
       NC_DAT => ADD(NC_DAT,AVERAGE_VEL_FILE_OBJECT() )
    END IF
    
    IF(NC_VORTICITY) THEN
       NC_DAT => ADD(NC_DAT,VORTICITY_FILE_OBJECT() )
    END IF
    
    IF(NC_SALT_TEMP) THEN
       NC_DAT => ADD(NC_DAT,SALT_TEMP_FILE_OBJECT() )
    END IF

    IF(NC_TURBULENCE) THEN
       NC_DAT => ADD(NC_DAT,TURBULENCE_FILE_OBJECT() )
    END IF

    IF (NC_SURFACE_HEAT .and. HEATING_ON) THEN
       NC_DAT => ADD(NC_DAT,SURFACE_HEATING_FILE_OBJECT() )
    END IF

    IF (NC_WIND_VEL) THEN
       NC_DAT => ADD(NC_DAT,WIND_VELOCITY_FILE_OBJECT() )
    END IF

    IF (NC_WIND_STRESS .and. WIND_ON) THEN
       NC_DAT => ADD(NC_DAT,WIND_STRESS_FILE_OBJECT() )
    END IF

    IF (NC_EVAP_PRECIP .and. PRECIPITATION_ON) THEN
       NC_DAT => ADD(NC_DAT,PRECIPITATION_FILE_OBJECT() )
    END IF

    IF(WETTING_DRYING_ON) THEN
       NC_DAT => ADD(NC_DAT, WET_DRY_FILE_OBJECT() )
    END IF

    IF(ICING_MODEL) THEN
       NC_DAT => ADD(NC_DAT, ICING_FILE_OBJECT() )
    END IF

    IF (GROUNDWATER_ON .and. NC_GROUNDWATER) THEN
       NC_DAT => ADD(NC_DAT,GROUNDWATER_FILE_OBJECT() )
    END IF

#   if defined (BioGen)
    IF(NC_BIO)THEN
      NC_DAT => ADD(NC_DAT,BIO_FILE_OBJECT() )
    END IF
#   endif      

#   if defined (WATER_QUALITY)
    IF(NC_WQM)THEN
      NC_DAT => ADD(NC_DAT,WQM_FILE_OBJECT() )
    END IF
#   endif      

#   if defined (NH)
    IF(NC_NH_QP) THEN
      NC_DAT => ADD(NC_DAT,QP_FILE_OBJECT() )
    END IF
    IF(NC_NH_RHS) THEN
      NC_DAT => ADD(NC_DAT,RHS_FILE_OBJECT() )
    END IF
#   endif      
    
#   if defined (ICE)
    !-----------------------------------------------------------------
    ! state variables
    NC_DAT => ADD(NC_DAT,ICE_DATA_STATE_FILE_OBJECT() )
    !-----------------------------------------------------------------
    ! velocity
    !-----------------------------------------------------------------
    NC_DAT => ADD(NC_DAT,ICE_VEL_FILE_OBJECT() )
    !-----------------------------------------------------------------
# endif

    IF (STARTUP_TYPE /= "crashrestart") THEN
       CALL NC_WRITE_FILE(NC_DAT)
       NC_DAT%FTIME%NEXT_STKCNT = 1

    ELSE
       NC_DAT%CONNECTED = .TRUE.
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_DATAFILE"
  END SUBROUTINE SETUP_DATFILE
!=============================================================  
  SUBROUTINE SETUP_AVGFILE
    IMPLICIT NONE
    character(len=80) :: tmp
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_AVGFILE"


    NC_AVG => ADD(NC_AVG,GRID_FILE_OBJECT() )

    IF (STARTUP_TYPE /= "crashrestart") THEN

       CALL NC_WRITE_FILE(NC_AVG)
       NC_AVG%FTIME%NEXT_STKCNT = 1

    ELSE
       NC_AVG%CONNECTED = .TRUE.
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_AVGFILE"
  END SUBROUTINE SETUP_AVGFILE
!=============================================================  
  SUBROUTINE SETUP_RSTFILE
    IMPLICIT NONE
    character(len=80) :: tmp
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_RSTFILE"


    NC_RST => ADD(NC_RST,GRID_FILE_OBJECT() )

    NC_RST => ADD(NC_RST,FILE_DATE_OBJECT() )

    NC_RST => ADD(NC_RST,VELOCITY_FILE_OBJECT() )

    NC_RST => ADD(NC_RST,AVERAGE_VEL_FILE_OBJECT() )

    NC_RST => ADD(NC_RST,VERTICAL_VEL_FILE_OBJECT() )

    NC_RST => ADD(NC_RST,TURBULENCE_FILE_OBJECT() )
    
    NC_RST => ADD(NC_RST,SALT_TEMP_FILE_OBJECT() )

    NC_RST => ADD(NC_RST,RESTART_EXTRAS_FILE_OBJECT() )
    
    IF(WETTING_DRYING_ON) THEN
       NC_RST => ADD(NC_RST, WET_DRY_FILE_OBJECT() )
    END IF

#   if defined (BioGen)
    NC_RST => ADD(NC_RST,BIO_FILE_OBJECT() )
#   endif      

#   if defined (WATER_QUALITY)
    NC_RST => ADD(NC_RST,WQM_FILE_OBJECT() )
#   endif      

#   if defined (NH)
    NC_RST => ADD(NC_RST,NH_RST_FILE_OBJECT() )
#   endif

# if defined (ICE)
      !-----------------------------------------------------------------
      ! state variables
    NC_RST => ADD(NC_RST,ICE_RESTART_STATE_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
    NC_RST => ADD(NC_RST,ICE_VEL_FILE_OBJECT() )
      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------
    NC_RST => ADD(NC_RST,ICE_EXTRA_FILE_OBJECT() )
# endif

    IF (STARTUP_TYPE /= "crashrestart") THEN
       CALL NC_WRITE_FILE(NC_RST)
       NC_RST%FTIME%NEXT_STKCNT = 1
    ELSE
       NC_RST%CONNECTED = .TRUE.
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_RSTFILE"
  END SUBROUTINE SETUP_RSTFILE
!=============================================================  
  SUBROUTINE DUMP_DATA(NCF)
    IMPLICIT NONE  
    TYPE(NCFILE), POINTER ::NCF
    TYPE(NCFTIME), POINTER :: FTM

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START DUMP_DATA"

    IF(DBG_SET(DBG_IO)) CALL PRINT_FILE(NCF)

    FTM => NCF%FTIME


    IF (FTM%MAX_STKCNT .NE. 0 .AND. &
         & FTM%NEXT_STKCNT > FTM%MAX_STKCNT) THEN

       FTM%NEXT_STKCNT=0
       CALL INCRIMENT_FNAME(NCF%FNAME)
       NCF%CONNECTED=.FALSE.
   
       ! WRITE NEW FILE'S CONSTANT DATA (GRID ETC)
       CALL NC_WRITE_FILE(NCF)

       ! INCRIMENT THE STACK COUNT
       FTM%NEXT_STKCNT = 1

    END IF

    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! DUMPING DATA FROM FILE OBJECT: &
         &"//trim(NCF%FNAME)

    ! IF UPDATE IODATA BECOMES SPECIFIC TO DIFFERENT DATA SETS IT
    ! WILL HAVE TO BE MOVED INSIDE OF THE PARTICULAR OUTPUT STATEMENTS
    CALL UPDATE_IODATA(NCF,IntTime)

    CALL NC_WRITE_FILE(NCF)


    ! ONCE THE FILE IS WRITEN INCRIMENT THE FILE OBJECT TIME
    FTM%PREV_IO = IntTime
    FTM%NEXT_IO = FTM%NEXT_IO + FTM%INTERVAL

    ! INCRIMENT THE STACK COUNT
    FTM%PREV_STKCNT = FTM%NEXT_STKCNT
    FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END DUMP_DATA"
    
  END SUBROUTINE DUMP_DATA
!=============================================================  
  SUBROUTINE DEFINE_DIMENSIONS
    USE BCS
    IMPLICIT NONE
    INTEGER :: SBUF, SOURCE
    INTEGER :: RBUF, DEST
    

    DIM_nele       => NC_MAKE_DIM(name='nele',len=ngl)
    DIM_node       => NC_MAKE_DIM(name='node',len=mgl)

    DIM_three      => NC_MAKE_DIM(name='three',len=3)
    DIM_four       => NC_MAKE_DIM(name='four',len=4)

    DIM_siglev     => NC_MAKE_DIM(name='siglev',len=kb)
    DIM_siglay     => NC_MAKE_DIM(name='siglay',len=kbm1)
    
    DIM_DateStrLen => NC_MAKE_DIM(name='DateStrLen',len=DateStrLen)
    DIM_time       => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)

# if defined (ICE)
    DIM_GRID       => NC_MAKE_DIM(name='dimension2',len=JMT_LOCAL)
    DIM_NCAT       => NC_MAKE_DIM(name='category',len=NCAT)
    DIM_Ntilay     => NC_MAKE_DIM(name='totallayer',len=Ntilay)
# endif   
   
    ! ONLY USED IF OBC IS ON
    IF(OBC_ON) THEN
       DIM_nobc       => NC_MAKE_DIM(name='nobc',len=IOBCN_GL)
    ELSE
       nullify(DIM_NOBC)
    END IF
    ! ONLY USED IF LOND SHORE FLOW BOUNDARY IS ON
    IF(OBC_LONGSHORE_FLOW_ON) THEN
       DIM_nlsf       => NC_MAKE_DIM(name='nlsf',len=nobclsf_GL)
    ELSE
       nullify(DIM_NLSF)
    END IF

    DIM_MaxNode    => NC_MAKE_RUNTIME_DIM(name='maxnode',len=MX_NBR_ELEM+3)

    DIM_MaxElem    => NC_MAKE_RUNTIME_DIM(name='maxelem',len=MX_NBR_ELEM+1)


  END SUBROUTINE DEFINE_DIMENSIONS
!=============================================================  
  FUNCTION GRID_FILE_OBJECT() RESULT(NCF)
    USE MOD_CLOCK
    USE MOD_FORCE
    IMPLICIT NONE
    
    INTEGER, ALLOCATABLE, SAVE :: partition(:)
    INTEGER :: status, I
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START GRID_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(partition(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:partition")
       partition = 0

       allocate(xm(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XM")
       xm = 0.0_SP

       allocate(ym(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YM")
       ym = 0.0_SP

       allocate(LON(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LON")
       lon = 0.0_SP

       allocate(LAT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LAT")
       lat = 0.0_SP

       allocate(xmc(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XMC")
       xmc = 0.0_SP

       allocate(ymc(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YMC")
       ymc = 0.0_SP

       allocate(LONC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LONC")
       lonc = 0.0_SP

       allocate(LATC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LATC")
       latc = 0.0_SP

       allocate(zz(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ZZ")
       zz = 0.0_SP

       allocate(z(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Z")
       z = 0.0_SP

       allocate(h(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:H")
       h = 0.0_SP
       
       allocate(nvgl(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NVGL")
       nvgl = 0

       allocate(EL(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL")
       EL = 0.0_SP


!!$       ! FILL WITH DUMMY VALUES...
!!$       DO I = 1 , CHAR_MAX_ATTLEN
!!$          TIDE_FORCING_COMMENTS(I:I) = 'X'
!!$          
!!$          RIVER_FORCING_COMMENTS(I:I) = 'X'
!!$          
!!$          GWATER_FORCING_COMMENTS(I:I) = 'X'
!!$          
!!$          HEAT_FORCING_COMMENTS(I:I) = 'X'
!!$          
!!$          WINDS_FORCING_COMMENTS(I:I) = 'X'
!!$          
!!$          PRECIP_FORCING_COMMENTS(I:I) = 'X'
!!$       END DO
       
       
    END IF


    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values=trim(case_title)) 
    NCF => ADD(NCF,ATT)


    ATT => NC_MAKE_ATT(name='institution',values=trim(institution)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='source',values=trim(fvcom_version)) 
    NCF => ADD(NCF,ATT)

    call get_timestamp(temp)
    timestamp = 'model started at: '//trim(temp)

    ATT => NC_MAKE_ATT(name='history',values=trim(timestamp)) 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='references',values=trim(fvcom_website)) 
    NCF => ADD(NCF,ATT)

    netcdf_convention = 'CF-1.0'
    ATT => NC_MAKE_ATT(name='Conventions',values=trim(netcdf_convention)) 
    NCF => ADD(NCF,ATT)

#  if defined (SPHERICAL)
    ATT => NC_MAKE_ATT(name='CoordinateSystem',values="GeoReferenced" ) 
    NCF => ADD(NCF,ATT)
#endif

    ATT => NC_MAKE_ATT(name='CoordinateProjection',values=PROJECTION_REFERENCE ) 
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    ATT=> NC_Make_Runtime_Att_CHR(name='Tidal_Forcing',values=TIDE_FORCING_COMMENTS)
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
    
    ATT=> NC_Make_Runtime_Att_CHR(name='River_Forcing',values=RIVER_FORCING_COMMENTS)
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    ATT=> NC_Make_Runtime_Att_CHR(name='GroundWater_Forcing',values=GWATER_FORCING_COMMENTS)
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
    
    ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Heat_Forcing',values=HEAT_FORCING_COMMENTS)
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    ATT=> NC_Make_Runtime_Att_CHR(name='Surface_Wind_Forcing',values=WINDS_FORCING_COMMENTS)
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    ATT=> NC_Make_Runtime_Att_CHR(name='Surface_PrecipEvap_Forcing',values=PRECIP_FORCING_COMMENTS)
    IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)

    IF (ICING_MODEL) THEN
       ATT=> NC_Make_Runtime_Att_CHR(name='Icing_Model_Forcing',values=ICING_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
    END IF

    IF (ICE_MODEL) THEN
       ATT=> NC_Make_Runtime_Att_CHR(name='Ice_Model_Forcing',values=ICE_FORCING_COMMENTS)
       IF(ASSOCIATED(ATT)) NCF => ADD(NCF,ATT)
    END IF


    IF(OBC_LONGSHORE_FLOW_ON) THEN
       aTT=> NC_MAKE_ATT(name='Special_Physical_processes',&
            & values='long shore flow adjustment for thermal wind and win&
            &d driven setup')
       NCF => ADD(NCF,ATT)
    END IF

    ! ADD THE VARIABLES

    ! NPROCS
    VAR  => NC_MAKE_AVAR(name='nprocs',values=nprocs)

    ATT  => NC_MAKE_ATT(name='long_name',values='number of processors') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

#  if defined (MULTIPROCESSOR)
    ! PARTITION
    IF(.not. ALLOCATED(partition)) ALLOCATE(partition(NT)); partition=myid
    VAR  => NC_MAKE_AVAR(name='partition',values=partition,DIM1=DIM_nele)


    ATT  => NC_MAKE_ATT(name='long_name',values='partition') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)
# endif

    ! X
    VAR  => NC_MAKE_AVAR(name='x',values=xm,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal x-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    ! Y
    VAR  => NC_MAKE_AVAR(name='y',values=ym,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal y-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! LON
    VAR  => NC_MAKE_AVAR(name='lon',values=lon,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    ! LAT
    VAR  => NC_MAKE_AVAR(name='lat',values=LAT,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF (ALLOCATED(XMC)) THEN
       ! XMC
       VAR  => NC_MAKE_AVAR(name='xc',values=xmc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal x-coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(YMC)) THEN
       ! YMC
       VAR  => NC_MAKE_AVAR(name='yc',values=ymc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal y-coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(LONC)) THEN
       ! LONC
       VAR  => NC_MAKE_AVAR(name='lonc',values=lonc,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal longitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
       VAR  => ADD(VAR,ATT)
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(LATC)) THEN
       ! LATC
       VAR  => NC_MAKE_AVAR(name='latc',values=LATC,DIM1=DIM_nele)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='zonal latitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF


    IF (ALLOCATED(zz)) THEN
       ! siglay
       VAR  => NC_MAKE_AVAR(name='siglay',&
            & values=zz,&
            & DIM1= DIM_node,&
            & DIM2= DIM_siglay )
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Sigma Layers') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='ocean_sigma/general_coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_min',values=-1.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_max',values=0.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='formula_terms',values='sigma: siglay eta: zeta depth: h') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(Z)) THEN
       ! siglev
       VAR  => NC_MAKE_AVAR(name='siglev',&
            & values=z, DIM1= DIM_node, DIM2= DIM_siglev )
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Sigma Levels') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='ocean_sigma/general_coordinate') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_min',values=-1.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_max',values=0.0_SPA) 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='formula_terms',values='sigma:siglay eta: zeta depth: h') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(h)) THEN
       ! h
       VAR  => NC_MAKE_AVAR(name='h',&
            & values=h, DIM1= DIM_node)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Bathymetry') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='depth') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='units',values='meters')
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='positive',values='down') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='Bathymetry_Mesh') 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    ! nv
    VAR  => NC_MAKE_AVAR(name='nv',&
         & values=nvgl, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding element') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! IINT
    VAR  => NC_MAKE_AVAR(name='iint',&
         & values=IO_IINT, DIM1= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='internal mode iteration number') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! time
    VAR => FLOAT_TIME_OBJECT &
         &(Data=IO_DAYS, &
         & USE_MJD=use_real_world_time, &
         & DIM=DIM_TIME)

    NCF  => ADD(NCF,VAR)


    ! Itime
    VAR  => ITIME_OBJECT &
         &(Data=IO_MJD, &
         & Use_MJD=use_real_world_time, &
         & DIM=DIM_TIME)

    NCF  => ADD(NCF,VAR)

    ! Itime2
    VAR => ITIME2_OBJECT &
         &(Data=IO_MSEC, &
         & Use_MJD=use_real_world_time, &
         & DIM=DIM_TIME)

    NCF => ADD(NCF,VAR)

    IF (use_real_world_time) THEN

       VAR => DATETIME_OBJECT &
            &(Data=IO_timestr,&
            & DIMSTR=DIM_DateStrLen,&
            & DIMTIME=DIM_TIME,&
            TIMEZONE=TIMEZONE)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF (ALLOCATED(EL)) THEN
       ! zeta
       VAR  => NC_MAKE_AVAR(name='zeta',&
            & values=EL, DIM1= DIM_node, DIM2= DIM_time)

       ATT  => NC_MAKE_ATT(name='long_name',values='Water Surface Elevation') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='meters') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='positive',values='up') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='standard_name',values='sea_surface_elevation') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='SSH_Mesh') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
       VAR  => ADD(VAR,ATT)

       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GRID_FILE_OBJECT"

  END FUNCTION GRID_FILE_OBJECT
!=============================================================  
  FUNCTION GRID_METRICS_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START GRID_METRICS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(NBEGL(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NBEGL")
       NBEGL = 0

       allocate(NTSN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NTSN")
       NTSN = 0

       allocate(NBSNGL(MGL,DIM_MaxNode%DIM),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NBSNGL")
       NBSNGL = 0

       allocate(NTVE(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NTVE")
       NTVE = 0

       allocate(NBVEGL(MGL,DIM_MaxElem%DIM),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NBVEGL")
       NBVEGL = 0

       allocate(A1U(NGL,4),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:A1U")
       A1U = 0.0_SP

       allocate(A2U(NGL,4),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:A2U")
       A2U = 0.0_SP

       allocate(AWX(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:AWX")
       AWX = 0.0_SP

       allocate(AWY(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:AWY")
       AWY = 0.0_SP

       allocate(AW0(NGL,3),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:AW0")
       AW0 = 0.0_SP

       allocate(ART2(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ART2")
       ART2 = 0.0_SP

       allocate(ART1(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ART1")
       ART1 = 0.0_SP


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! NBE
    VAR  => NC_MAKE_PVAR(name='nbe',&
         & values=NBEGL, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='elements surrounding each element') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NTSN
    VAR  => NC_MAKE_AVAR(name='ntsn',&
         & values=NTSN, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='#nodes surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NBSN
    VAR  => NC_MAKE_PVAR(name='nbsn',&
         & values=NBSNGL, DIM1= DIM_node, DIM2= DIM_MaxNode)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NTVE
    VAR  => NC_MAKE_AVAR(name='ntve',&
         & values=NTVE, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='#elems surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NBVE
    VAR  => NC_MAKE_PVAR(name='nbve',&
         & values=NBVEGL, DIM1= DIM_node, DIM2= DIM_MaxElem)

    ATT  => NC_MAKE_ATT(name='long_name',values='elems surrounding each node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! A1U
    VAR  => NC_MAKE_AVAR(name='a1u',&
         & values=a1u, DIM1= DIM_nele, DIM2= DIM_four)

    ATT  => NC_MAKE_ATT(name='long_name',values='a1u') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! A2U
    VAR  => NC_MAKE_AVAR(name='a2u',&
         & values=a2u, DIM1= DIM_nele, DIM2= DIM_four)

    ATT  => NC_MAKE_ATT(name='long_name',values='a2u') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! AW0
    VAR  => NC_MAKE_AVAR(name='aw0',&
         & values=aw0, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='aw0') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! AWX
    VAR  => NC_MAKE_AVAR(name='awx',&
         & values=awx, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='awx') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! AWY
    VAR  => NC_MAKE_AVAR(name='awy',&
         & values=awy, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='awy') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ART2
    VAR  => NC_MAKE_AVAR(name='art2',&
         & values=art2, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='Area of elements around a node') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ART1
    VAR  => NC_MAKE_AVAR(name='art1',&
         & values=art1, DIM1= DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='Area of Node-Base Con&
         &trol volume') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GRID_METRICS_FILE_OBJECT"

  END FUNCTION GRID_METRICS_FILE_OBJECT

!=============================================================  
  FUNCTION FILE_DATE_OBJECT() RESULT(NCF)
    USE MOD_CLOCK
    IMPLICIT NONE

    INTEGER :: status
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START FILE_DATE_OBJECT"
    
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! FILE_DATE
    VAR  => NC_MAKE_AVAR(name='file_date',&
         & values=IO_FILE_DATE, DIM1= DIM_DateStrLen, DIM2 = DIM_time)

    IF(USE_REAL_WORLD_TIME) THEN
       ATT  => NC_MAKE_ATT(name='time_zone',values=TRIM(TIMEZONE))
    ELSE
       ATT  => NC_MAKE_ATT(name='time_zone',values="UTC")
    END IF

    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

  END FUNCTION FILE_DATE_OBJECT


!=============================================================  
  FUNCTION VELOCITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VELOCITY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(U(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:U")
       U = 0.0_sp

       allocate(V(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:V")
       V = 0.0_sp

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! U
    VAR  => NC_MAKE_AVAR(name='u',&
         & values=U, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! V
    VAR  => NC_MAKE_AVAR(name='v',&
         & values=V, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VELOCITY_FILE_OBJECT"

  END FUNCTION VELOCITY_FILE_OBJECT
!=============================================================  
  FUNCTION VORTICITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VORTIICTY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(VORT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VORT")
       VORT = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! VORTICITY
    VAR  => NC_MAKE_AVAR(name='vorticity',&
         & values=VORT, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Ertels 2d potential vorticity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VORTICITY_FILE_OBJECT"

  END FUNCTION VORTICITY_FILE_OBJECT
!=============================================================  
  FUNCTION AVERAGE_VEL_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START AVERAGE_VEL_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(UA(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UA")
       UA = 0.0_SP

       allocate(VA(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VA")
       VA = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UA
    VAR  => NC_MAKE_AVAR(name='ua',&
         & values=UA, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertically Averaged x-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! VA
    VAR  => NC_MAKE_AVAR(name='va',&
         & values=VA, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertically Averaged y-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END AVERAGE_VEL_FILE_OBJECT"

  END FUNCTION AVERAGE_VEL_FILE_OBJECT
  
!=============================================================  
  FUNCTION VERTICAL_VEL_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START VERTICAL_VEL_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(WW(NGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WW")
       WW = 0.0_SP

       allocate(WTS(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WTS")
       WTS = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! WTS
    VAR  => NC_MAKE_AVAR(name='omega',&
         & values=WTS, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertical Sigma Coordinate Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)



    ! WW
    VAR  => NC_MAKE_AVAR(name='ww',&
         & values=WW, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Upward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END VERTICAL_VEL_FILE_OBJECT"

  END FUNCTION VERTICAL_VEL_FILE_OBJECT
!=============================================================  
  FUNCTION SALT_TEMP_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SALT_TEMP_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(T1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:T1")
       T1 = 0.0_SP

       allocate(S1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:S1")
       S1 = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! T
    VAR  => NC_MAKE_AVAR(name='temp',&
         & values=T1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)


    ! T
    VAR  => NC_MAKE_AVAR(name='salinity',&
         & values=S1, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)    

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SALT_TEMP_FILE_OBJECT"

  END FUNCTION SALT_TEMP_FILE_OBJECT
!=============================================================  
  FUNCTION TURBULENCE_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    TYPE(NCDIM),  POINTER :: DIM3

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START TURBULENCE_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(KM(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:KM")
       KM = 0.0_SP

       allocate(KH(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:KH")
       KH = 0.0_SP

       allocate(KQ(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:KQ")
       KQ = 0.0_SP
       
# if defined (GOTM)
       allocate(TKE(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:TKE")
       TKE = 0.0_SP

       allocate(TEPS(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:TEPS")
       TEPS = 0.0_SP
# else
       allocate(Q2(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Q2")
       Q2 = 0.0_SP

       allocate(Q2L(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Q2L")
       Q2L = 0.0_SP
       
       allocate(L(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:L")
       L = 0.0_SP
# endif


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! KM
    VAR  => NC_MAKE_AVAR(name='km',&
         & values=km, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Eddy Viscosity For Momentum') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! KH
    VAR  => NC_MAKE_AVAR(name='kh',&
         & values=kh, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Eddy Viscosity For Scalars') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! KQ
    VAR  => NC_MAKE_AVAR(name='kq',&
         & values=kq, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Eddy Viscosity For Q2/Q2L') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m 2 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

# if defined (GOTM)
    ! TKE
    VAR  => NC_MAKE_AVAR(name='tke',&
         & values=tke, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Energy') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m2 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! TEPS
    VAR  => NC_MAKE_AVAR(name='teps',&
         & values=teps, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Energy Dissipation Rate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='w kg-1 ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)
# else
    ! Q2
    VAR  => NC_MAKE_AVAR(name='q2',&
         & values=q2, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Energy') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m2 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! Q2L
    VAR  => NC_MAKE_AVAR(name='q2l',&
         & values=q2l, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Kinetic Ene&
         &rgy X Turbulent Macroscale') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m3 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! L
    VAR  => NC_MAKE_AVAR(name='l',&
         & values=l, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Turbulent Macroscale') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m3 s-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

# endif

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END TURBULENCE_FILE_OBJECT"

  END FUNCTION TURBULENCE_FILE_OBJECT
!=============================================================  
  FUNCTION SURFACE_HEATING_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SURFACE_HEATING_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(SWRAD_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:SWRAD")
       swrad_watts = 0.0_SP

       allocate(WTSURF_WATTS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WTSURF")
       WTSURF_watts = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! SWRAD
    VAR  => NC_MAKE_AVAR(name='short_wave',&
         & values=swrad_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Short Wave Radiation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WTSURF - NET HEAT FLUX
    VAR  => NC_MAKE_AVAR(name='net_heat_flux',&
         & values=WTSURF_WATTS, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Surface Net Heat Flux') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W m-2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SURFACE_HEATING_FILE_OBJECT"

  END FUNCTION SURFACE_HEATING_FILE_OBJECT
!=============================================================  
  FUNCTION WIND_VELOCITY_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WIND_VELOCITY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(UUWIND(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UUWIND")
       UUWIND = 0.0_SP

       allocate(VVWIND(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VVWIND")
       VVWIND = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UUWIND
    VAR  => NC_MAKE_AVAR(name='uwind_speed',&
         & values=uuwind, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VVWIND
    VAR  => NC_MAKE_AVAR(name='vwind_speed',&
         & values=vvwind, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WIND_VELOCITY_FILE_OBJECT"
    
  END FUNCTION WIND_VELOCITY_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION WIND_STRESS_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WIND_STRESS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(WUSURF_save(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WUSURF")
       WUSURF_save = 0.0_SP

       allocate(WVSURF_save(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WVSURF")
       WVSURF_save = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UUWIND
    VAR  => NC_MAKE_AVAR(name='uwind_stress',&
         & values=WUSURF_save, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! VVWIND
    VAR  => NC_MAKE_AVAR(name='vwind_stress',&
         & values=WVSURF_save, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='Wind Stress') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Pa') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WIND_STRESS_FILE_OBJECT"
    
  END FUNCTION WIND_STRESS_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION PRECIPITATION_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    TYPE(NCDIM),  POINTER :: DIM1
    TYPE(NCDIM),  POINTER :: DIM2
    TYPE(NCDIM),  POINTER :: DIM3

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START PRECIPITATION_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(QPREC(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QPREC2")
       QPREC = 0.0_SP

       allocate(QEVAP(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QEVAP2")
       QEVAP = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! PRECIPITATION
    VAR  => NC_MAKE_AVAR(name='precip',&
         & values=qprec, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Precipitation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='description',values='Precipitation, ocean &
         &lose water is negative') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! EVAPORATION
    VAR  => NC_MAKE_AVAR(name='evap',&
         & values=QEVAP, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Evaporation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='description',values='Evaporation, ocean &
         &lose water is negative') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END PRECIPITATION_FILE_OBJECT"
    
  END FUNCTION PRECIPITATION_FILE_OBJECT

!=============================================================  
!=============================================================  
!=============================================================  
  FUNCTION RESTART_EXTRAS_FILE_OBJECT() RESULT(NCF)
    USE BCS
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START RESTART_EXTRAS_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(COR(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:COR")
       COR = 0.0_SP

       allocate(CC_SPONGE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:CC_SPONGE")
       CC_SPONGE = 0.0_SP

       allocate(ET(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL")
       ET = 0.0_SP


! OPEN BOUNDARY SETTINGS
       IF (OBC_ON) THEN
          ! NEED SPECIAL VARIABLE FOR GLOBAL NODE NUMBER OF LOCAL NODES...
          allocate(I_OBC_N_OUTPUT(IOBCN_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:I_OBC_N_OUTPUT")
          I_OBC_N_OUTPUT = 0
          
          
          allocate(type_obc(IOBCN_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:type_obc")
          type_OBC = 0
       END IF

! LONG SHORE FLOW BOUNDARY SETTINGS
       IF(OBC_LONGSHORE_FLOW_ON) THEN
          allocate(IBCLSF_OUTPUT(NOBCLSF_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:IBCLSF_OUTPUT")
          IBCLSF_OUTPUT = 0
          
          allocate(RBC_GEO(NOBCLSF_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RBC_GEO")
          RBC_GEO = 0.0_SP
          
          allocate(RBC_WDF(NOBCLSF_GL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RBC_WDF")
          RBC_WDF = 0.0_SP
       END IF
       
          
       allocate(TMEAN1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Tmean1")
       TMEAN1 = 0.0_SP
       
       allocate(SMEAN1(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Smean1")
       SMEAN1 = 0.0_SP

# if defined (EQUI_TIDE)
       ALLOCATE(EL_EQI(MGL), stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL_EQI")
       EL_EQI = 0.0_SP
# endif

# if defined (ATMO_TIDE)
       ALLOCATE(EL_ATMO(MGL), stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL_ATMO")
       EL_ATMO = 0.0_SP
# endif

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! COR
    VAR  => NC_MAKE_AVAR(name='cor',&
         & values=COR, DIM1= DIM_nele)

    ATT  => NC_MAKE_ATT(name='long_name',values='Coriolis Parameter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! CC_SPONGE
    VAR  => NC_MAKE_AVAR(name='cc_sponge',&
         & values=cc_sponge, DIM1= DIM_nele)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sponge Layer Parameter') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='nd') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! et
    VAR  => NC_MAKE_AVAR(name='et',&
         & values=et, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Water Surface Elevation At Last Timestep') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='positive',values='up') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_surface_elevation') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SSH_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! TMEAN1
    VAR  => NC_MAKE_AVAR(name='tmean1',&
         & values=tmean1, DIM1= DIM_node, DIM2= DIM_siglay )

    ATT  => NC_MAKE_ATT(name='long_name',values='mean initial temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! SMEAN1
    VAR  => NC_MAKE_AVAR(name='smean1',&
         & values=smean1, DIM1= DIM_node, DIM2= DIM_siglay )

    ATT  => NC_MAKE_ATT(name='long_name',values='mean initial salinity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(OBC_ON)THEN
       ! OBC_GRID
       VAR  => NC_MAKE_AVAR(name='obc_nodes',&
            & values=I_OBC_N_OUTPUT, DIM1= DIM_nobc)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Node Number') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
       
       ! OBC_TYPE
       VAR  => NC_MAKE_AVAR(name='obc_type',&
            & values=type_obc, DIM1= DIM_nobc)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Open Boundary Type') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='obc_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    IF(OBC_LONGSHORE_FLOW_ON)THEN
       ! long shore flow grid
       VAR  => NC_MAKE_AVAR(name='lsf_nodes',&
            & values=ibclsf_output, DIM1= DIM_nlsf)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Longshore Flow Node Number') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='lsf_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
       
       ! 
       VAR  => NC_MAKE_AVAR(name='wdf',&
            & values=RBC_WDF, DIM1= DIM_nlsf)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Wind Driven Flow Adjustment Scaling') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_range',values='[0 1]') 
       VAR  => ADD(VAR,ATT)       

       ATT  => NC_MAKE_ATT(name='grid',values='lsf_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)

       VAR  => NC_MAKE_AVAR(name='geo',&
            & values=RBC_GEO, DIM1= DIM_nlsf)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='Thermal Wind Flow Adjustment Scaling') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='valid_range',values='[0 1]') 
       VAR  => ADD(VAR,ATT)       

       ATT  => NC_MAKE_ATT(name='grid',values='lsf_grid') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF


# if defined (EQUI_TIDE)
    VAR  => NC_MAKE_AVAR(name='el_eqi',&
         & values=el_eqi, DIM1= DIM_node)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Equilibrium tide adjustment height') 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    
    NCF  => ADD(NCF,VAR)
# endif
    
# if defined (ATMO_TIDE)
    VAR  => NC_MAKE_AVAR(name='el_atmo',&
         & values=el_atmo, DIM1= DIM_node)
    
    ATT  => NC_MAKE_ATT(name='long_name',values='Atmospheric tide adjustment height') 
    VAR  => ADD(VAR,ATT)
     
    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
       
    NCF  => ADD(NCF,VAR)
# endif


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END RESTART_EXTRAS_FILE_OBJECT"
    
  END FUNCTION RESTART_EXTRAS_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION WET_DRY_FILE_OBJECT() RESULT(NCF)
    USE MOD_WD
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WET_DRY_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(ISWETN(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETN")
       ISWETN = 0

       allocate(ISWETC(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETC")
       ISWETC = 0 

       allocate(ISWETNT(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETNT")
       ISWETN = 0

       allocate(ISWETCT(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETCT")
       ISWETCT = 0 

       allocate(ISWETCE(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ISWETCE")
       ISWETCE = 0


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! WET NODES
    VAR  => NC_MAKE_AVAR(name='wet_nodes',&
         & values=ISWETN, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Nodes') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET CELLS
    VAR  => NC_MAKE_AVAR(name='wet_cells',&
         & values=ISWETC, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Cells') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET NODES AT LAST INT STEP
    VAR  => NC_MAKE_AVAR(name='wet_nodes_prev_int',&
         & values=ISWETNT, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Nodes_At_Previous_Internal_Step') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET CELLS AT LAST EXT STEP
    VAR  => NC_MAKE_AVAR(name='wet_cells_prev_int',&
         & values=ISWETCT, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Cells_At_Previous_Internal_Step') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! WET CELLS AT LAST EXT STEP
    VAR  => NC_MAKE_AVAR(name='wet_cells_prev_ext',&
         & values=ISWETCE, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Wet_Cells_At_Previous_External_Step') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WET_DRY_FILE_OBJECT"
    
  END FUNCTION WET_DRY_FILE_OBJECT

!=============================================================  
!=============================================================  
  FUNCTION ICING_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICING_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(ICING_0kts(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_0kts")
       ICING_0kts = 0.0_SP

       allocate(ICING_10kts(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_10kts")
       ICING_10kts = 0.0_SP

       allocate(icing_wndX(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_wndY")
       ICING_wndX = 0.0_SP

       allocate(icing_wndY(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_wndX")
       ICING_wndY = 0.0_SP

       allocate(icing_satmp(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ICING_satmp")
       ICING_satmp = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! ICING_0KTS
    VAR  => NC_MAKE_AVAR(name='icing_0kts',&
         & values=icing_0kts, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Hazard@0knots') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m C s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_10KTS
    VAR  => NC_MAKE_AVAR(name='icing_10kts',&
         & values=icing_10kts, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Hazard@10knots') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m C s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_WNDX
    VAR  => NC_MAKE_AVAR(name='icing_wndx',&
         & values=icing_wndx, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Wind x-direction') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_WNDY
    VAR  => NC_MAKE_AVAR(name='icing_wndy',&
         & values=icing_wndy, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Wind y-direction') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m s^-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! ICING_SATMP
    VAR  => NC_MAKE_AVAR(name='icing_satmp',&
         & values=icing_satmp, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Icing Surface Air Temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICING_FILE_OBJECT"
    
  END FUNCTION ICING_FILE_OBJECT
!=============================================================  
!=============================================================  
  FUNCTION GROUNDWATER_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: GROUNDWATER_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(BFWDIS(MGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BFWDIS")
       BFWDIS = 0.0_SP

       IF(GROUNDWATER_TEMP_ON) THEN
          allocate(BFWTMP(MGL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BFWTMP")
          BFWTMP = 0.0_SP
       END IF

       IF(GROUNDWATER_SALT_ON) THEN
          allocate(BFWSLT(MGL),stat=status)
          IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BFWSLT")
          BFWSLT = 0.0_SP
       END IF


    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! GROUNDWATER VOLUME FLUX
    VAR  => NC_MAKE_AVAR(name='groundwater_flux',&
         & values=bfwdis, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='groundwater volume flux') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m3 s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! GROUNDWATER INFLOW TEMPERATURE
    IF(GROUNDWATER_TEMP_ON) THEN
       VAR  => NC_MAKE_AVAR(name='groundwater_temp',&
            & values=bfwdis, DIM1= DIM_node, DIM2= DIM_time)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='groundwater inflow temperature') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='degrees_C') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF

    ! GROUNDWATER INFLOW SALINITY
    IF(GROUNDWATER_SALT_ON) THEN
       VAR  => NC_MAKE_AVAR(name='groundwater_salt',&
            & values=bfwdis, DIM1= DIM_node, DIM2= DIM_time)
       
       ATT  => NC_MAKE_ATT(name='long_name',values='groundwater inflow salinity') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='units',values='1e-3') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
       VAR  => ADD(VAR,ATT)
       
       ATT  => NC_MAKE_ATT(name='type',values='data') 
       VAR  => ADD(VAR,ATT)
       
       NCF  => ADD(NCF,VAR)
    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END GROUND_FILE_OBJECT"
    
  END FUNCTION GROUNDWATER_FILE_OBJECT
!=============================================================  

# if defined (NH)
!=============================================================
  FUNCTION QP_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START QP_FILE_OBJECT"

    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(QP(MGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:QP")
       QP = 0.0_sp

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! QP
    VAR  => NC_MAKE_AVAR(name='qp',&
         & values=QP, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pertubation Pressure')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END QP_FILE_OBJECT"

  END FUNCTION QP_FILE_OBJECT
!=============================================================
  FUNCTION RHS_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START RHS_FILE_OBJECT"

    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(RHS(MGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:RHS")
       RHS = 0.0_sp

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! RHS
    VAR  => NC_MAKE_AVAR(name='rhs',&
         & values=RHS, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Right Hand Side Term')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='Kg M-3 S-2')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END RHS_FILE_OBJECT"

  END FUNCTION RHS_FILE_OBJECT
!=============================================================
  FUNCTION NH_RST_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START NH_RST_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(W4ZT(MGL,KB),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:W4ZT")
       W4ZT = 0.0_SP

       allocate(NHQDRX(NGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQDRX")
       NHQDRX = 0.0_SP

       allocate(NHQDRY(NGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQDRY")
       NHQDRY = 0.0_SP

       allocate(NHQDRZ(MGL,KBM1),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQDRZ")
       NHQDRZ = 0.0_SP

       allocate(NHQ2DX(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQ2DX")
       NHQ2DX = 0.0_SP

       allocate(NHQ2DY(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:NHQ2DY")
       NHQ2DY = 0.0_SP

    END IF

  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ! W4ZT
    VAR  => NC_MAKE_AVAR(name='w4zt',&
         & values=W4ZT, DIM1= DIM_node, DIM2= DIM_siglev, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Temp Vertical Velocity in Z')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQDRX
    VAR  => NC_MAKE_AVAR(name='nhqdrx',&
         & values=NHQDRX, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient XCor')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQDRY
    VAR  => NC_MAKE_AVAR(name='nhqdry',&
         & values=NHQDRY, DIM1= DIM_nele, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient YCor')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQDRZ
    VAR  => NC_MAKE_AVAR(name='nhqdrz',&
         & values=NHQDRZ, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient ZCor')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQ2DX
    VAR  => NC_MAKE_AVAR(name='nhq2dx',&
         & values=NHQ2DX, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient XCor 2D')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NHQ2DY
    VAR  => NC_MAKE_AVAR(name='nhq2dy',&
         & values=NHQ2DY, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nonhydrostatic Pressure Gradient YCor 2D')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N M-3')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END NH_RST_FILE_OBJECT"
    
  END FUNCTION NH_RST_FILE_OBJECT
# endif

!!! ggao/0104/2008  !! restart file for ice model
!!---------------------------------------------------------------------------
# if defined (ICE)
!=============================================================  
  FUNCTION ICE_DATA_STATE_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    REAL(SP),SAVE, POINTER :: AICE_P(:)
    REAL(SP),SAVE, POINTER :: VICE_P(:)

    character(len=100)    :: timestamp, temp, netcdf_convention

    integer :: NI,K,I,J

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_STATE_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
       
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


   ! ICE AREA 

    AICE_P => aice(1:m,1)

!    VAR  => NC_MAKE_AVAR(name='aice',&
!         & values=aice, DIM1= DIM_node, DIM2= Dim_grid, DIM3 = DIM_time)

    VAR  => NC_MAKE_PVAR(name='aice',&
         & values=aice_P, DIM1= DIM_node, DIM2 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_area_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
   ! ICE VOLUME
    VICE_P => vice(1:m,1)

!    VAR  => NC_MAKE_AVAR(name='vice',&
!         & values=vice, DIM1= DIM_node, DIM2= Dim_grid, DIM3 = DIM_time)

    VAR  => NC_MAKE_PVAR(name='vice',&
         & values=vice_P, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice volume per unit grid area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_volume_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_STATE0_FILE_OBJECT"

  END FUNCTION ICE_DATA_STATE_FILE_OBJECT
!=============================================================


!=============================================================  
  FUNCTION ICE_RESTART_STATE_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    REAL(SP),SAVE, POINTER :: AICEN_P(:,:)
    REAL(SP),SAVE, POINTER :: VICEN_P(:,:)
    REAL(SP),SAVE, POINTER :: VSNON_P(:,:)
    REAL(SP),SAVE, POINTER :: Esnon_P(:,:)
    REAL(SP),SAVE, POINTER :: Tsfcn_P(:,:)
    REAL(SP),SAVE, POINTER :: Eicen_P(:,:)

    integer :: NI,K,I,J

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_STATE_FILE_OBJECT"
    
    
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

   ! ICE AREA 
   
    AICEN_P => aicen(1:imt_local,1,1:ncat)
 
    VAR  => NC_MAKE_PVAR(name='AICEN',&
         & values=AICEN_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_area_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
   ! ICE VOLUME
      VICEN_P => vicen(1:imt_local,1,1:ncat)
    
    VAR  => NC_MAKE_PVAR(name='VICEN',&
         & values=Vicen_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice volume per unit grid area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_volume_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
   ! SNOW VOLUME
       Vsnon_P =>  vsnon(1:imt_local,1,1:ncat)
    
    VAR  => NC_MAKE_PVAR(name='VSNON',&
         & values=Vsnon_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='snow volume per unit grid area') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='snow_volume_category') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
     

   ! ice/snow surface temperature
   
    Tsfcn_P =>  Tsfcn(1:imt_local,1,1:ncat)
    
    
    VAR  => NC_MAKE_PVAR(name='TSFCN',&
         & values=Tsfcn_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice/snow surface temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_surface_temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='deg C') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
     
   !  snow energy 
       Esnon_P => esnon(1:imt_local,1,1:ncat)
    
    VAR  => NC_MAKE_PVAR(name='ESNON',&
         & values=Esnon_P, DIM1= DIM_node, DIM2= DIM_NCAT, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='snow energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='snow_energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='J/m$^2$') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
!!-------------------------------------------------------------------    
   !  ice energy category and layer

    Eicen_P => eicen(1:imt_local,1,1:ntilay)
   
    VAR  => NC_MAKE_PVAR(name='EICEN',&
         & values=Eicen_P, DIM1= DIM_node, DIM2= DIM_ntilay, DIM3 = DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_energy ') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='J/m$^2$') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_STATE_FILE_OBJECT"

  END FUNCTION ICE_RESTART_STATE_FILE_OBJECT
!=============================================================

!=============================================================  
  FUNCTION ICE_VEL_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_VEL_FILE_OBJECT"
    
    ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
    ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
    ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

       IOPROC_ALLOCATED = .true.

       allocate(UICE2(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UICE")
       UICE2 = 0.0_SP

       allocate(VICE2(NGL),stat=status)
       IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VICE")
       VICE2 = 0.0_SP

    END IF


  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()


    ! UICE2
    VAR  => NC_MAKE_AVAR(name='uuice',&
         & values=Uice2, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Ice x-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! VICE2
    VAR  => NC_MAKE_AVAR(name='vvice',&
         & values=VICE2, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Ice y-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_VEL_FILE_OBJECT"

  END FUNCTION ICE_VEL_FILE_OBJECT

!============================================================= 
  FUNCTION ICE_EXTRA_FILE_OBJECT() RESULT(NCF)
    IMPLICIT NONE

    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    character(len=100)    :: timestamp, temp, netcdf_convention

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ICE_EXTRA_FILE_OBJECT"
    
  ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

      !-----------------------------------------------------------------
      ! fresh water, salt, and heat flux
      !-----------------------------------------------------------------

   ! fresh water
    
    VAR  => NC_MAKE_AVAR(name='fresh',&
         & values=fresh, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='flux of water ice to ocean') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='fresh_water') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m^2/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

   ! salt
    
    VAR  => NC_MAKE_AVAR(name='fsalt',&
         & values=fsalt, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='flux of salt ice to ocean') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='flux_salt') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='kg/m^2/s') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)


   ! net heat
    
    VAR  => NC_MAKE_AVAR(name='fhnet',&
         & values=fhnet, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='heat flux ice to ocean') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='flux_heat') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='W/m^2') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)
    
      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------
   ! ice strength
    
    VAR  => NC_MAKE_AVAR(name='strength',&
         & values=strength, DIM1= DIM_node,Dim2 =DIM_GRID, DIM3= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice strength') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='ice_strength') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='N/m') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)
    
    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)    

    NCF  => ADD(NCF,VAR)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
    ! ice cell mask
    VAR  => NC_MAKE_AVAR(name='isicec',&
         & values=isicec, DIM1= DIM_nele, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='ice cell mask') 
    VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='units',values='') 
!    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ICE_EXTRA_FILE_OBJECT"

  END FUNCTION ICE_EXTRA_FILE_OBJECT
!=============================================================
# endif
!!---------------------------------------------------------------------------
!!---------------------------------------------------------------------------
!! ggao/0104/2008  !! restart file for ice model
!=============================================================  
# if defined (BioGen)
  FUNCTION BIO_FILE_OBJECT() RESULT(NCF)
   IMPLICIT NONE

   INTEGER :: status
   LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II
   REAL(SP), ALLOCATABLE :: BIO_ALL_TMP(:,:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START BIO_FILE_OBJECT"
   
   DO II=1,NTT
     ALLOCATE(BIO_ALL_TMP(0:MT,KB))
     BIO_ALL_TMP(:,:) = BIO_ALL(:,:,II)
    
   ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
   ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
   ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
   IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

     IOPROC_ALLOCATED = .true.

     ALLOCATE(BIO_ALL_TMP(MGL,KB),STAT=STATUS)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:BIO_ALL")
     BIO_ALL_TMP = 0.0_SP

   END IF


  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()

  ! BIO
     VAR  => NC_MAKE_AVAR(name=TRIM(BIO_NAME(II,1)),       &
           & values=BIO_ALL_TMP, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

     ATT  => NC_MAKE_ATT(name='long_name',values=TRIM(BIO_NAME(II,3))) 
     VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
!    VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='units',values=TRIM(BIO_NAME(II,2)) )
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
     VAR  => ADD(VAR,ATT)
    
     ATT  => NC_MAKE_ATT(name='type',values='data') 
     VAR  => ADD(VAR,ATT)    

     NCF  => ADD(NCF,VAR)
  
     DEALLOCATE(BIO_ALL_TMP)

   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END BIO_FILE_OBJECT"

  END FUNCTION BIO_FILE_OBJECT
# endif
!=============================================================
# if defined (WATER_QUALITY)
  FUNCTION WQM_FILE_OBJECT() RESULT(NCF)
   IMPLICIT NONE

   INTEGER :: status
   LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
   TYPE(NCFILE), POINTER :: NCF
   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCATT),  POINTER :: ATT
   INTEGER :: II
   REAL(SP),ALLOCATABLE :: WQM_TMP(:,:)

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START WQM_FILE_OBJECT"
   
   DO II = 1,NB
     ALLOCATE(WQM_TMP(0:MT,KB))
     WQM_TMP(:,:) = WQM(:,:,II)
    
   ! IO PROC MUST ALLOCATE SPACE FOR THE ARRAYS 
   ! THESE ARRAYS MUST HAVE THE ATTRIBUTE SAVE AND FOR CLARITY
   ! SHOULD HAVE THE SAME NAME AS THOSE USED ON THE OTHER PROCESSORS!
   IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED) THEN

     IOPROC_ALLOCATED = .true.

     ALLOCATE(WQM_TMP(MGL,KB),STAT=STATUS)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:WQM")
     WQM_TMP = 0.0_SP

   END IF


  ! ALLOCATE THE NEW FILE OBJECT
   NCF => NEW_FILE()

  ! WQM
    
     VAR  => NC_MAKE_AVAR(name=TRIM(WQM_NAME(II,1)),&
           & values=WQM_TMP, DIM1= DIM_node, DIM2= DIM_siglay, DIM3 = DIM_time)

     ATT  => NC_MAKE_ATT(name='long_name',values=TRIM(WQM_NAME(II,3))) 
     VAR  => ADD(VAR,ATT)

!    ATT  => NC_MAKE_ATT(name='standard_name',values='sea_water_temperature') 
!    VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='units',values=TRIM(WQM_NAME(II,2))) 
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
     VAR  => ADD(VAR,ATT)

     ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
     VAR  => ADD(VAR,ATT)
    
     ATT  => NC_MAKE_ATT(name='type',values='data') 
     VAR  => ADD(VAR,ATT)    

     NCF  => ADD(NCF,VAR)
  
     DEALLOCATE(WQM_TMP)
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END WQM_FILE_OBJECT"

  END FUNCTION WQM_FILE_OBJECT
# endif
!=============================================================  

SUBROUTINE UPDATE_IODATA(NCF,NOW)
  IMPLICIT NONE
  TYPE(TIME), INTENT(IN) :: NOW
  TYPE(NCVAR), POINTER :: VAR1,VAR2
  TYPE(NCFILE), POINTER :: NCF
  LOGICAL :: FOUND
!====================================================================
! THIS SUBROUTINE IS IN CHARGE OF UPDATING ANY VARIABLES FOR IO WHICH
! ARE NOT ALREADY UPDATED BY FVCOM DURING THE MAIN LOOP:
!====================================================================

  if(.not. Associated(NCF)) CALL FATAL_ERROR&
       &("UPDATE_IODATA: THE FILE OBJECT IS NOT ASSOCIATED!")
  
  VAR1 => FIND_VAR(NCF,"time",FOUND)
  IF(FOUND) CALL UPDATE_FLOAT_TIME(VAR1,NOW)

  VAR1 => FIND_VAR(NCF,"Itime",FOUND)
  IF(FOUND) THEN
     VAR2 => FIND_VAR(NCF,"Itime2",FOUND)
     IF (.NOT.FOUND) THEN
        CALL WARNING&
             & ("FOUND ONLY PART OF INTEGER TIME VARIABLE IN OUT PUT FILE!")
     ELSE
        CALL UPDATE_ITIME(VAR1,VAR2,NOW)
     END IF
  END IF

  VAR1 => FIND_VAR(NCF,"Times",FOUND)
  IF(FOUND) CALL UPDATE_DATETIME(VAR1,NOW)


  ! IINT IS A LONG LONG INTEGER BUT WE CAN'T WRITE LONG INTEGERS...
  VAR1 => FIND_VAR(NCF,"iint",FOUND)
  IF(FOUND) VAR1%scl_int = IINT

  
  VAR1 => FIND_VAR(NCF,"file_date",FOUND)
  IF(FOUND) CALL UPDATE_DATETIME(VAR1,GET_NOW())

END SUBROUTINE UPDATE_IODATA
!=============================================================  
!=============================================================  
!=============================================================  
!=============================================================      
  SUBROUTINE SETUP_MPI_IO_MODE(TF,COMMGRP) 
    !===================================================================================|
    !  INITIALIZE MPI ENVIRONMENT                                                       |
    !===================================================================================|
    LOGICAL, INTENT(INOUT) :: TF
    INTEGER, INTENT(OUT)   :: COMMGRP
    INTEGER :: total_group  
    INTEGER :: fvcom_group
    INTEGER :: SBUF,RBUF, trueval, i
    INTEGER, allocatable :: fvcom_subset(:)
    INTEGER IERR


    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "STARTING SETUP_MPI_IO_MODE"


    if (NPROCS .LE. 3) THEN
       if (TF) &
            & CALL WARNING("FVCOM CAN NOT USE MPI IO MODE WHEN RUN &
            &ON LESS THAN 4 PROCESSORS", "CONTINUEING WITHOUT THIS OPTION!")
       TF = .false.
    end if

    IF (TF) THEN

# if defined (MULTIPROCESSOR)
       ! MODIFY CONTROL VARIABLE EFFECTED BY THE USE OF MPI IO MODE
       NPROCS_TOTAL = NPROCS
       IOPROCID = NPROCS

       IOPROC=.FALSE.
       IF(MYID==NPROCS) IOPROC=.TRUE.

       ! NOW REDUCE THE NUMBER OF PROCESSORS BY ONE
       NPROCS= NPROCS - 1

       ! RETURN HANDLE TO LIST OF GROUPS
       call mpi_comm_group(mpi_comm_world,total_group,ierr)

       !first: define subset
       !comp procs will have process_id/myid:  0-->(nprocs-1)/1-->nprocs
       !mr_printy  will have process_id/myid:  nprocs/nprocs_tot
       allocate(fvcom_subset(nprocs)) ; fvcom_subset = 0
       do i=1,nprocs
          fvcom_subset(i) = i-1
       end do
       call mpi_group_incl(total_group,nprocs,fvcom_subset,fvcom_group,ierr)
       call mpi_comm_create(mpi_comm_world,fvcom_group,COMMGRP,ierr)
       deallocate(fvcom_subset)


       call mpi_barrier(mpi_comm_world,ierr)

       IF(.not. IOPROC) then
          sbuf = myid
          rbuf = 0
          call mpi_allreduce(sbuf,rbuf,1,mpi_integer,mpi_sum,COMMGRP,ierr)

          trueval = 0
          do i=1,nprocs
             trueval = trueval + i
          end do
          
          if( trueval .NE. rbuf) &
               & CALL FATAL_ERROR("TESTING GROUP COMMUNICATION FOR MPI &
               &IO MODE FAILED")
       END IF
       
       
       sbuf = myid 
       rbuf = 0
       call mpi_allreduce(sbuf,rbuf,1,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)

       trueval = 0
       do i=1,nprocs_total
          trueval = trueval + i
       end do

       if( trueval .NE. rbuf ) &
            & CALL FATAL_ERROR("TESTING WORLD COMMUNICATION FOR MPI &
            &IO MODE FAILED")

       IF (DBG_SET(dbg_log)) &
            & write(IPT,*) "!  MPI IO MODE IS ACTIVE"
# endif

    ELSE 

       NPROCS_TOTAL = NPROCS
       IOPROCID = -1
       IOPROC=.FALSE.

       ! MUST SET INTENT(OUT) COMMGRP
# if defined (MULTIPROCESSOR)
       ! IF MULTIPROCESSOR BUT NOT MPI IO SET IT TO COMM WORLD!
       COMMGRP = MPI_COMM_WORLD
# else
       ! IF NOT MULTIPROCESSOR SET A DUMMY VALUE
       COMMGRP = -999

# endif


       IF (DBG_SET(dbg_log)) &
            & write(IPT,*) "!  MPI IO MODE IS NOT ACTIVE"

    END IF

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "END SETUP_MPI_IO_MODE"
  END SUBROUTINE SETUP_MPI_IO_MODE


# if defined (MULTIPROCESSOR)

  SUBROUTINE MPI_IO_LOOP
    implicit none
    integer RBUF, J, IERR
    integer SOURCE
    integer status
    INTEGER STAT(MPI_STATUS_SIZE)
    character(len=5) :: lpc

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START MPI_IO_LOOP"

    IN_MPI_IO_LOOP = .TRUE.

    if (MYID .NE. IOPROCID) then
       IF(DBG_SET(DBG_LOG))THEN
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++"
          WRITE(IPT,*) "! I AM NOT THE IO PROC:"
          WRITE(IPT,*) "! SETTING IN_MPI_IO_LOOP = T"
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++"
          if(DBG_SET(dbg_sbr)) &
               & write(IPT,*) "END MPI_IO_LOOP"
       END IF
       
       RETURN

    else   
       IF(DBG_SET(DBG_LOG))THEN
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++++++++++++"
          WRITE(IPT,*) "! I AM THE IO PROC: STARTING IO LOOP"
          WRITE(IPT,*) "!++++++++++++++++++++++++++++++++++++++++++++"
       END IF
    end if


    J = 1
    DO 
       write(lpc,'(I5.5)') J
       IF(DBG_SET(DBG_IO)) write(IPT,*)"========IO LOOP COUNT "//lpc//"==============="

       !      if (J==5) then
       !         call system("sleep 5")
       !         call report_error("testing!")
       !      end if

       IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "WAITING FOR IO CODE:"

       STAT=0
!       SOURCE = MPI_ANY_SOURCE ! ALLOW ANY PROCESSOR TO COMMUNICATE WITH THE IO PROC
       SOURCE = 0 ! ONLY THE MASTER CAN COMMMUNCATE WITH IO PROC
       CALL MPI_RECV(RBUF,1,MPI_INTEGER,SOURCE&
            &,SYNC_TAG,MPI_COMM_WORLD,STAT,IERR)

       select case(RBUF)
       case(EXT_CODE)
          ! JUST CALL STOP...
          IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "IO PROC RECIEVED EXIT CODE:"
          IF(DBG_SET(DBG_IO)) WRITE(IPT,*) "THATS ALL FOLKS!"

          Call mpi_finalize(IERR)
          stop

       case(WAIT_CODE)

          IF(DBG_SET(DBG_IO)) write(IPT,*) "OTHER PROCS WAITING FOR ME:"

          Call mpi_barrier(mpi_comm_world,ierr)

          IF(DBG_SET(DBG_IO)) write(IPT,*) "FINISHED WRITING LAST FILE: GO"

       case default

          IF(DBG_SET(DBG_IO)) write(IPT,*) "IO PROC GOT GO FOR CALL:",RBUF

          Call mpi_barrier(mpi_comm_world,ierr)

          IF(DBG_SET(DBG_IO)) write(IPT,*) "IO PROC SINKED FOR GO"

          call CALL_FUNC(RBUF,status)
          IF (status/=0) call fatal_error("IO PROC Recieved bad go code",&
               & "Could not retrieve valid function pointer?")

          IF(DBG_SET(DBG_IO)) write(IPT,*) "FINISHED MPI_IO CALL!"

       end select

       J = J + 1

    END DO

  END SUBROUTINE MPI_IO_LOOP



  SUBROUTINE MPI_IO_SYNCHRONIZE(CODE)
    implicit none
    INTEGER, INTENT(IN) ::  CODE
    integer :: ierr, STATUS
    INTEGER STAT(MPI_STATUS_SIZE)

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "START MPI_IO_SYNCHRONIZE"



    IF (MYID .EQ. IOPROCID) CALL FATAL_ERROR("IO PROC SHOULD NEVER&
         & BE IN MPI_IO_SYNCHRONIZE")


    IF(DBG_SET(DBG_IO)) write(IPT,*) "SYNCHRONIZE FVCOM_GROUP WITH IOPROC:"

    select CASE(CODE)
    case(WAIT_CODE)

       if (MYID==1) then
          IF(DBG_SET(DBG_IO)) write(IPT,*) "MASTER SENDING: WAIT_CODE"
          CALL MPI_SEND(CODE,1,MPI_INTEGER,IOPROCID-1,SYNC_TAG&
               &,MPI_COMM_WORLD,IERR)
       else
          IF(DBG_SET(DBG_IO)) write(IPT,*) "WAITING FOR IO PROC"
       end if

       call mpi_barrier(mpi_comm_world,ierr)

       IF(DBG_SET(DBG_IO)) write(IPT,*) "IO PROC HAS FINISHED LAS&
            &T WRITE COMMAND"

    case default

       if (MYID==1) then
          IF(DBG_SET(DBG_IO)) write(IPT,*) "MASTER SENDING GO CODE:",CODE
          CALL MPI_SEND(CODE,1,MPI_INTEGER,IOPROCID-1,SYNC_TAG,MPI_COMM_WORLD,IERR)
       else
          IF(DBG_SET(DBG_IO)) write(IPT,*) "WAITING FOR IO PROC"
       end if

       call mpi_barrier(mpi_comm_world,ierr)

          call CALL_FUNC(CODE,status)
          IF (status/=0) call fatal_error("MPI_IO_SYNCHRONIZE: Bad go code",&
               & "Could not retrieve valid function pointer?")

          IF(DBG_SET(DBG_IO)) write(IPT,*) "FINISHED MPI_IO CALL!"

    end select


  END SUBROUTINE MPI_IO_SYNCHRONIZE

# else
    SUBROUTINE MPI_IO_SYNCHRONIZE(CODE)
    implicit none
    integer code

    CALL FATAL_ERROR("THIS IS A DUMMY ROUTINE:","MPI_IO_SYNCHRONIZ&
         &E only exists when FVCOM is compiled with multiprocessor",&
         & "USE_MPI_IO_MODE should not be .true.")
  END SUBROUTINE MPI_IO_SYNCHRONIZE

# endif


END MODULE MOD_NCDIO
