MODULE MOD_SST2GRD
  USE MOD_NCTOOLS
  USE MOD_UTILS
  USE MOD_INPUT
  USE MOD_TIME
  USE CONTROL
  IMPLICIT NONE



  REAL(SP) :: CVAL
  LOGICAL :: CONSTANT

  INTEGER, PARAMETER       :: NX = 600
  INTEGER, PARAMETER       :: NY = 320
  INTEGER, PARAMETER       :: FUNIT = 46

  ! MICRO_SECONDS PER DAY
  INTEGER(ITIME), PARAMETER :: SST_INTERVAL = 86400*1000000

  INTEGER nsite   !number of valid measurements in regular grid

  INTEGER i,j,k,num,im,iday,iyear,CNT

  INTEGER :: IOS

  REAL DEP,zval,x0,y0, mean

  INTEGER, ALLOCATABLE :: NV(:,:)
  REAL(SP), DIMENSION(NX,NY) ::OSST,OX,OY, lon, lat
  INTEGER, DIMENSION(NX,NY) :: MASK
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: XX,YY,SST

  REAL(DP), ALLOCATABLE,DIMENSION(:)   ::SSTD
  REAL(SP), ALLOCATABLE, DIMENSION(:)    :: X,Y,Z,TMP
  CHARACTER CH*3,month(12)*3,yn*1,cyear*4


  CHARACTER(len=80) :: GRID_DIR
  CHARACTER(len=80) :: Lat_data_name
  CHARACTER(len=80) :: Lon_data_name
  CHARACTER(len=80) :: Data_prefix

  NAMELIST /NML_INPUT/         &
       & INPUT_DIR,            &
       & LAT_DATA_NAME,        &
       & LON_DATA_NAME,        &
       & DATA_PREFIX,          &
       & OUTPUT_DIR,           & 
       & GRID_DIR,             &
       & GRID_FILE,            &
       & PROJECTION_REFERENCE, &
       & START_DATE,           &
       & TIMEZONE           

  CHARACTER(len=120) :: FNAME
  LOGICAL           :: FEXIST


  TYPE(INTERP_WEIGHTS) :: WEIGHTS

  TYPE(NCFILE), POINTER :: NC_OUT
  TYPE(NCFTIME), POINTER :: FTM

  TYPE(TIME) :: NOW


CONTAINS



  SUBROUTINE GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)
    use mod_sng

    character(len=*), INTENT(IN)::CVS_Id  ! [sng] CVS Identification
    character(len=*), INTENT(IN)::CVS_Date ! [sng] Date string
    character(len=*), INTENT(IN)::CVS_Name ! [sng] File name string
    character(len=*), INTENT(IN)::CVS_Revision ! [sng] File revision string

    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
    ! Command-line parsing
    character(80)::arg_val ! [sng] command-line argument value
    character(200)::cmd_ln ! [sng] command-line
    character(80)::opt_sng ! [sng] Option string
    character(2)::dsh_key ! [sng] command-line dash and switch
    character(200)::prg_ID ! [sng] Program ID

    integer::arg_idx ! [idx] Counting index
    integer::arg_nbr ! [nbr] Number of command-line arguments
    integer::opt_lng ! [nbr] Length of option

    CONSTANT = .false.
    ! Main code
    call ftn_strini(cmd_ln) ! [sng] sng(1:len)=NUL

    call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command-line into single string
    call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID

    arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments

    if (arg_nbr .LE. 0 ) then
       if(MSR) WRITE(IPT,*) "You must specify an arugument:"
       if(MSR) Call MYHelpTxt
       call PSHUTDOWN
    end if

    arg_idx=1 ! [idx] Counting index
    do while (arg_idx <= arg_nbr)
       call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
       dsh_key=arg_val(1:2) ! [sng] First two characters of option
       if (dsh_key == "--") then
          opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
          if (opt_lng <= 0) then
             if(MSR) write(IPT,*) "Long option has no name"
             call PSHUTDOWN
          end if

          opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
          if (dbg_lvl >= dbg_io) then
             if(MSR) write (6,"(5a,i3)") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": DEBUG Double hyphen indicates multi-character option: ", &
                  "opt_sng = ",opt_sng(1:ftn_strlen(opt_sng)),", opt_lng = ",opt_lng
          end if
          if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
             call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level


          else if (opt_sng == "constant" .or. opt_sng == "Constant" ) then
             call ftn_arg_get(arg_idx,arg_val,cval) ! [enm] Debugging level

             write(ipt,*) "Seeting a constant value:", cval
             CONSTANT = .TRUE.

          else if (opt_sng == "Fileame" .or.opt_sng == "filename"&
               & .or.opt_sng == "FILENAME") then

             call ftn_arg_get(arg_idx,arg_val,FName) ! [sng] Input file
             FName=FName(1:ftn_strlen(FName))
             ! Convert back to a fortran string!

          else if (opt_sng == "help" .or.opt_sng == "HELP" .or. opt_sng&
               & == "Help") then

             if(MSR) call MYHelpTxt
             call PSHUTDOWN

          else ! Option not recognized
             arg_idx=arg_idx-1 ! [idx] Counting index
             if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
          endif ! endif option is recognized
          ! Jump to top of while loop
          cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
       endif ! endif long option

       if (dsh_key == "-V" .or.dsh_key == "-v" ) then

          if(MSR) write(IPT,*) prg_id
          call PSHUTDOWN

       else if (dsh_key == "-H" .or.dsh_key == "-h" ) then

          if(MSR) Call MYHelpTxt
          Call PSHUTDOWN

       else ! Option not recognized
          arg_idx=arg_idx-1 ! [idx] Counting index
          if(MSR) call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
       endif ! endif arg_val


    end do ! end while (arg_idx <= arg_nbr)

    CALL dbg_init(IPT_BASE,.false.)

  END SUBROUTINE GET_COMMANDLINE

  SUBROUTINE MYHELPTXT
    IMPLICIT NONE

    WRITE(IPT,*) "Add better help here!"
    WRITE(IPT,*) "! OPTIONS:"
    WRITE(IPT,*) "! --filename=XXX"
    WRITE(IPT,*) "! --constant=X.X (Set constant value)"
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   The namelist runfile for the program! "
    WRITE(IPT,*) "!   "
    WRITE(IPT,*) "!   Namelist OPTIONS: "
    WRITE(IPT,*) "!    INPUT_DIR"
    WRITE(IPT,*) "!    LAT_DATA_NAME"
    WRITE(IPT,*) "!    LON_DATA_NAME"
    WRITE(IPT,*) "!    DATA_PREFIX"
    WRITE(IPT,*) "!    OUTPUT_DIR"
    WRITE(IPT,*) "!    GRID_DIR"
    WRITE(IPT,*) "!    GRID_FILE"
    WRITE(IPT,*) "!    PROJECTION_REFERENCE"
    WRITE(IPT,*) "!    START_DATE"
    WRITE(IPT,*) "!    TIMEZONE"
    WRITE(IPT,*) "!    "
    WRITE(IPT,*) "!    "


    WRITE(IPT,*) "! NOTES: Do not run this program in parallel!"


  END SUBROUTINE MYHELPTXT

  SUBROUTINE READ_NAMELIST
    IMPLICIT NONE
    integer :: ios, i
    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Begins: Read_Name_List;"


    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Read_Name_List: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read IO Information
    READ(UNIT=NMLUNIT, NML=NML_INPUT,IOSTAT=ios)
   if(ios .NE. 0 ) THEN
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_INPUT)
       CALL FATAL_ERROR("Can Not Read NameList NML_INPUT from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_INPUT)

    CLOSE(NMLUNIT)


  END SUBROUTINE READ_NAMELIST


  SUBROUTINE CREATE_INTERP
    IMPLICIT NONE


    !READ FIRST SST MAP TO GET MASK OF INPUT DATA  
    write(CH,'(i3.3)') 1  ! 

    FNAME = TRIM(INPUT_DIR)//TRIM(data_prefix)//CH
    CALL FOPEN(FUNIT,TRIM(FNAME),'cfr')

    DO J=1,NY
       READ(FUNIT,*) (OSST(I,J),I=1,NX)   
    ENDDO
    CLOSE(FUNIT)


    ! BUILD LIST OF VALID DATA LOCATIONS
    mask = 0
    WHERE(OSST .GT. -90) mask = 1

    NSITE = SUM(MASK)


    ALLOCATE(X(NSITE))
    ALLOCATE(Y(NSITE))     
    ALLOCATE(Z(NSITE))
    ALLOCATE(TMP(NSITE))


    CNT =0 
    DO J=1,NY 
       DO I=1,NX
          IF(mask(I,J) == 1) then
             CNT = CNT +1
             X(CNT) = OX(I,J)
             Y(CNT) = OY(I,J)
          ENDIF
       ENDDO
    ENDDO


    WEIGHTS%Nin  = NSITE
    WEIGHTS%Nout = M

    ALLOCATE(WEIGHTS%INDEX(NSITE))
    ALLOCATE(WEIGHTS%PTW(M))

    CALL SORTRX(NSITE,X,WEIGHTS%INDEX)

    DO I = 1,WEIGHTS%Nin
       TMP(I) = X(WEIGHTS%INDEX(I))
    END DO
    X = TMP

    DO I = 1,WEIGHTS%Nin
       TMP(I) = Y(WEIGHTS%INDEX(I))
    END DO
    Y = TMP

    WRITE(*,*) "min/max(X)",minval(X),maxval(X)
    WRITE(*,*) "min/max(Y)",minval(Y),maxval(Y)

    DO I = 1, WEIGHTS%Nout
       CALL GEN_WTS(X,Y,XX(I),YY(I),WEIGHTS%INDEX,WEIGHTS%PTW(I))
    END DO

  END SUBROUTINE CREATE_INTERP



  SUBROUTINE MY_OUTFILE
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    TYPE(NCDIM),  POINTER :: DIM_node
    TYPE(NCDIM),  POINTER :: DIM_nele
    TYPE(NCDIM),  POINTER :: DIM_three
    TYPE(NCDIM),  POINTER :: DIM_DateStrLen
    TYPE(NCDIM),  POINTER :: DIM_time

    ALLOCATE(SST(M))

    DIM_node => NC_MAKE_DIM(name='node',len=M)
    DIM_nele => NC_MAKE_DIM(name='nele',len=N)
    DIM_three => NC_MAKE_DIM(name='three',len=3)

    DIM_DateStrLen => NC_MAKE_DIM(name='DateStrLen',len=20)
    DIM_time => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    NC_OUT => NCF

    ALLOCATE(NCF%FTIME)


    NCF%FNAME = TRIM(OUTPUT_DIR)//'sst.nc'

    ! ADD THE FILE ATTRIBUTES
    ATT => NC_MAKE_ATT(name='title',values="some sst data") 
    NCF => ADD(NCF,ATT)

    ATT => NC_MAKE_ATT(name='source',values='FVCOM') 
    NCF => ADD(NCF,ATT)

    ! X
    VAR  => NC_MAKE_AVAR(name='x',values=xx,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal x-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    ! Y
    VAR  => NC_MAKE_AVAR(name='y',values=yy,DIM1=DIM_node)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal x-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! nv
    VAR  => NC_MAKE_AVAR(name='nv',&
         & values=nv, DIM1= DIM_nele, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding element') 
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


    ! SST
    VAR  => NC_MAKE_AVAR(name='sst',&
         & values=SST, DIM1= DIM_node, DIM2= DIM_time)

    ATT  => NC_MAKE_ATT(name='long_name',values='Sea Surface Temperature') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='celcius') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


  END SUBROUTINE MY_OUTFILE



END MODULE MOD_SST2GRD
