Module Mod_Nesting
  use all_vars
  use mod_utils
  use mod_ncdio
  use mod_par
  implicit none


  integer :: MT_NEST,MGL_NEST
  integer :: NT_NEST,NGL_NEST


  !  NODES          -- Map Global
  INTEGER, POINTER :: MID_NEST_GL(:)

  !  ELEMENTS:      --  Map Global 
  INTEGER, POINTER :: NID_NEST_GL(:)

  !  NODES          -- MAP2Local -- Map local2global
  INTEGER, POINTER :: MID_NEST(:),MID_NEST_OUTPUT(:)
  !  ELEMENTS:      -- MAP2Local -- Map local2global
  INTEGER, POINTER :: NID_NEST(:),NID_NEST_OUTPUT(:)

  INTEGER,  ALLOCATABLE :: NV_NEST(:,:)

  ! For output only
  REAL(SP), ALLOCATABLE :: X_NEST(:),Y_NEST(:)
  REAL(SP), ALLOCATABLE :: LON_NEST(:),LAT_NEST(:)
  REAL(SP), ALLOCATABLE :: XMC_NEST(:),YMC_NEST(:)
  REAL(SP), ALLOCATABLE :: LONC_NEST(:),LATC_NEST(:)
  REAL(SP), ALLOCATABLE :: ZZ_NEST(:,:),Z_NEST(:,:)
  REAL(SP), ALLOCATABLE :: H_NEST(:)


  ! Used in code for subdomain mode only
  REAL(SP), ALLOCATABLE :: UA_NEST(:),VA_NEST(:),EL_NEST(:)
  REAL(SP), ALLOCATABLE :: U_NEST(:,:),V_NEST(:,:)
  REAL(SP), ALLOCATABLE :: S1_NEST(:,:),T1_NEST(:,:)

  ! For data IO
  REAL(SP), ALLOCATABLE :: UA_BLK(:,:),VA_BLK(:,:),EL_BLK(:,:)
  REAL(SP), ALLOCATABLE :: U_BLK(:,:,:),V_BLK(:,:,:)
  REAL(SP), ALLOCATABLE :: S1_BLK(:,:,:),T1_BLK(:,:,:)

  TYPE(NCVAR), POINTER :: VAR_UA
  TYPE(NCVAR), POINTER :: VAR_VA
  TYPE(NCVAR), POINTER :: VAR_EL
  TYPE(NCVAR), POINTER :: VAR_U
  TYPE(NCVAR), POINTER :: VAR_V
  TYPE(NCVAR), POINTER :: VAR_S1
  TYPE(NCVAR), POINTER :: VAR_T1

  INTEGER :: BLK_INDEX
  
  LOGICAL, PRIVATE :: NEED_INIT_NEST = .TRUE.


  ! VARAIBLE FOR NESTING OUTPUT
  TYPE(NCFILE), POINTER :: NESTING_FILE
  TYPE(NCFILE), POINTER :: NC_NEST

  TYPE(NCDIM), POINTER  :: DIM_nele_nest
  TYPE(NCDIM), POINTER  :: DIM_node_nest
  TYPE(NCDIM), POINTER  :: DIM_time_nest
  TYPE(NCDIM), POINTER  :: DIM_BLK_nest

  ! VARIABLES FOR NESTING INPUT
  REAL(DP) :: LD_SEC
  TYPE(TIME) :: LD_TSTEP
  TYPE(TIME) :: BLOCK_ENDS
  TYPE(TIME) :: BLOCK_STARTS
  INTEGER :: NBLOCKS


  !--Parameters in NameList NML_NESTING
  LOGICAL NESTING_ON                   !!TRUE IF OUTPUT RESART FILES
  CHARACTER(LEN=80) NESTING_NODES_FILE !!FILE WITH NODE LIST FOR NESTING
  CHARACTER(LEN=80) NESTING_ELEMS_FILE !!FILE WITH NODE LIST FOR NESTING
  INTEGER :: NESTING_BLOCKSIZE         !!SIZE OF DATA BLOCKS IN FILE
  CHARACTER(LEN=80) NESTING_FILE_NAME  !!NAME OF RESTART FILE
  CHARACTER(LEN=80) NESTING_MODE       !!LARGE DOMAIN or SUBDOMAIN

  NAMELIST /NML_NESTING/     &
       & NESTING_ON,         &
       & NESTING_NODES_FILE, &
       & NESTING_ELEMS_FILE, &
       & NESTING_BLOCKSIZE,  &
       & NESTING_FILE_NAME,  &
       & NESTING_MODE
  
  INTERFACE SET_VAR
     MODULE PROCEDURE SET_VAR_1D
     MODULE PROCEDURE SET_VAR_2D
  END INTERFACE


CONTAINS
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE DEFINE_EXTRA_DIMENSIONS
    IMPLICIT NONE

    DIM_nele_nest => NC_MAKE_DIM(name='nele',len=ngl_nest)
    DIM_node_nest => NC_MAKE_DIM(name='node',len=mgl_nest)
    DIM_time_nest => NC_MAKE_DIM(name='time',len=NF90_UNLIMITED)
    DIM_BLK_nest  => NC_MAKE_DIM(NAME='nesting_block',len=NESTING_BLOCKSIZE)

    RETURN
  END SUBROUTINE DEFINE_EXTRA_DIMENSIONS
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE NAME_LIST_INITIALIZE_NEST
    USE CONTROL

    IMPLICIT NONE

    !--Parameters in NameList NML_NESTING
    NESTING_ON = .FALSE.
    NESTING_NODES_FILE = trim(casename)//"nesting_nodes.dat"
    NESTING_ELEMS_FILE = trim(casename)//"nesting_cells.dat"
    NESTING_BLOCKSIZE = -1
    NESTING_FILE_NAME = trim(casename)//"nesting.dat"
    NESTING_MODE ="'large_domain','subdomain'"

    RETURN
  END SUBROUTINE NAME_LIST_INITIALIZE_NEST
  !==============================================================================!
  !  
  !==============================================================================!  
  SUBROUTINE NAME_LIST_PRINT_NEST
    USE CONTROL

    IMPLICIT NONE

    WRITE(UNIT=IPT,NML=NML_NESTING)

    RETURN
  END SUBROUTINE NAME_LIST_PRINT_NEST
  !==============================================================================!
  !
  !==============================================================================!
  SUBROUTINE NAME_LIST_READ_NEST
    USE MOD_UTILS
    USE CONTROL

    IMPLICIT NONE

    INTEGER :: ios,I
    CHARACTER(LEN=120) :: FNAME
    IF(DBG_SET(dbg_sbr)) write(IPT,*) "Subroutine Begins: name_list_read_nest;"

    ios = 0
    FNAME = "./"//trim(casename)//"_run.nml"
    IF(DBG_SET(dbg_io)) write(IPT,*) "Get_nestpar: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    !READ NESTING FLAG
    READ(UNIT=NMLUNIT, NML=NML_NESTING,IOSTAT=ios)  
    IF(ios /= 0)THEN
       IF(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_NESTING)
       CALL FATAL_ERROR("Can Not Read NameList NML_NESTING from file: "//trim(FNAME))
    END IF

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_NESTING)

    CLOSE(NMLUNIT)

    ! DO SOME BASIC CHECKING ON THE NESTING NAMELIST
    IF(NESTING_ON) THEN
       IF(LEN_TRIM(NESTING_NODES_FILE) == 0 ) CALL FATAL_ERROR &
            & ("THE NESTING NODES FILE NAME IS EMPTY IN THE NAME LIST")

       IF(LEN_TRIM(NESTING_ELEMS_FILE) == 0 ) CALL FATAL_ERROR &
            & ("THE NESTING ELEMS FILE NAME IS EMPTY IN THE NAME LIST")

       IF( NESTING_BLOCKSIZE < 1) CALL FATAL_ERROR &
            & ("THE NESTING_BLOCKSIZE IS LESS THAN ONE IN THE NESTING NAME LIST")

       IF(LEN_TRIM(NESTING_FILE_NAME) == 0 ) CALL FATAL_ERROR &
            & ("THE NESTING FILE NAME IS EMPTY IN THE NAME LIST")

    END IF

    if(DBG_SET(dbg_sbr)) &
         & write(IPT,*) "Subroutine Ends: name_list_read_nest;"    

  END SUBROUTINE NAME_LIST_READ_NEST
  !============================================================================!
  !   
  !============================================================================!   
  SUBROUTINE LOAD_NESTING_OPB(mode,nnode,nodes,nelem,elems)
    USE CONTROL
    IMPLICIT NONE
    CHARACTER(len=80), INTENT(IN) :: MODE
    INTEGER, INTENT(OUT) :: nnode, nelem
    INTEGER, INTENT(OUT), POINTER :: nodes(:),elems(:)

    INTEGER, POINTER :: e1_dummy(:), e2_dummy(:)
    INTEGER, POINTER :: n1_dummy(:), n2_dummy(:)

    INTEGER :: I,J,n1,n2,n3
    CHARACTER(LEN=80) :: PATHNFILE, TEMP,TEMP2
    INTEGER :: ISCAN,IOS,sender,IERR

    IF(MSR)THEN

       !============================================
       !
       !  READ THE NESTING NODES
       !
       !============================================

       PATHNFILE = TRIM(INPUT_DIR)//TRIM(NESTING_NODES_FILE)
       CALL FOPEN(MNESTUNIT,trim(pathnfile),'cfr')

       ISCAN = SCAN_FILE(MNESTUNIT,"Node_Nest Number", ISCAL = nnode)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of NESTING NODES FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'The header must contain: "Node_Nest Number = "', &
               & 'Followed by an integer number of Nesting Nodes')
       END IF


       REWIND MNESTUNIT

       DO WHILE(.TRUE.)
          READ(MNESTUNIT,*,IOSTAT=IOS,END=199)N1,N2,N3
          if (IOS == 0) then
             BackSpace MNESTUNIT
             exit
          end if

          CYCLE

199       Call FATAL_ERROR('Improper formatting of NODE NESTING FILE:',&
               &'Reached end of file with out finding NESTING data?',&
               &'FORMAT: NEST# SUBDOMAIN# LARGEDOMAIN# (ALL INTEGERS)')

       END DO
       
       ALLOCATE(n1_dummy(nnode))
       ALLOCATE(n2_dummy(nnode))

       I = 0
       DO WHILE(.TRUE.)

          READ(MNESTUNIT,*,IOSTAT=IOS) N1,N2,N3
          IF(IOS < 0) exit

          I = I + 1
          IF(I > NNODE) CALL FATAL_ERROR('Number of rows of data in the NODE NESTING file &
               &exceeds the stated number of nodes in the header ?')


          n1_dummy(i) = N2
          n2_dummy(i) = N3

       END DO



       IF( 1 > N2 .or. N2 > MGL) then
          write(temp,'(I8)') I
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('OPEN BOUNDARY NODE NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE AND ENSURE OPEN BOUNDARY NODES <= '//trim(temp2))
       END IF

       IF( 1 > N3 .or. N3 > NGL) then
          write(temp,'(I8)') I
          CALL FATAL_ERROR('OPEN BOUNDARY NODE NUMBER'//trim(temp)//&
               & ' IS NOT IN THE VALID RANGE',&
               & 'THE OPEN BOUNDARY NODE TYPE MUST BE GREATER THAN 0',&
               & 'AND LESS THAN 11. SEE MOD_OBC.F FOR DESCRIPTION')
       END IF

       !============================================
       !
       !  READ THE NESTING ELEMENTS
       !
       !============================================

       PATHNFILE = TRIM(INPUT_DIR)//TRIM(NESTING_ELEMS_FILE)
       CALL FOPEN(NNESTUNIT,trim(pathnfile),'cfr')


       ISCAN = SCAN_FILE(NNESTUNIT,"Elem_Nest Number", ISCAL = nelem)
       IF(ISCAN /= 0) then
          write(temp,'(I2)') ISCAN
          call fatal_error('Improper formatting of NESTING ELEMENTS FILE: ISCAN ERROR&
               &# '//trim(temp),&
               & 'The header must contain: "Elem_Nest Number = "', &
               & 'Followed by an integer number of Nesting Elements')
       END IF

       REWIND NNESTUNIT

       DO WHILE(.TRUE.)
          READ(NNESTUNIT,*,IOSTAT=IOS,END=299)N1,N2,N3
          if (IOS == 0) then
             BackSpace NNESTUNIT
             exit
          end if

          CYCLE

299       Call FATAL_ERROR('Improper formatting of ELEMENT NESTING FILE:',&
               &'Reached end of file with out finding NESTING data?',&
               &'FORMAT: NEST# SUBDOMAIN# LARGEDOMAIN# (ALL INTEGERS)')

       END DO

       ALLOCATE(e1_dummy(NELEM))
       ALLOCATE(e2_dummy(NELEM))

       I = 0
       DO WHILE(.TRUE.)

          READ(NNESTUNIT,*,IOSTAT=IOS) N1,N2,N3
          IF(IOS < 0) exit

          I = I + 1
          IF(I > NELEM) CALL FATAL_ERROR('Number of rows of data in the ELEMENT NESTING file &
               &exceeds the stated number of elements in the header ?')


          e1_dummy(i) = N2
          e2_dummy(i) = N3

       END DO



       IF(MODE == 'large_domain') THEN

          NODES => n2_dummy
          ELEMS => e2_dummy

          DEALLOCATE(n1_dummy)
          DEALLOCATE(e1_dummy)

       ELSEIF(MODE == 'subdomain') THEN

          NODES => n1_dummy
          ELEMS => e1_dummy

          DEALLOCATE(n2_dummy)
          DEALLOCATE(e2_dummy)

       ELSE

          CALL FATAL_ERROR("INVALID NESTING MODE - PLEASE EXAMINE YOUR RUNFILE!")
       END IF


       ! CHECK TO MAKE SURE VALUES ARE REASONABLE!
       IF( 1 > MINVAL(NODES)) THEN
          write(temp,'(I8)') MINLOC(NODES)
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('NESTING NODE NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING NODES ARE 1 <= '//trim(temp2))
       END IF

       IF(  MAXVAL(NODES) > MGL) THEN
          write(temp,'(I8)') MAXLOC(NODES)
          write(temp2,'(I8)') MGL
          CALL FATAL_ERROR('NESTING NODE NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING NODES ARE 1 <= '//trim(temp2))
       END IF

       IF( 1 > MINVAL(ELEMS)) THEN
          write(temp,'(I8)') MINLOC(ELEMS)
          write(temp2,'(I8)') NGL
          CALL FATAL_ERROR('NESTING ELEMENT NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING ELEMENTS ARE 1 <= '//trim(temp2))
       END IF

       IF(  MAXVAL(ELEMS) > NGL) THEN
          write(temp,'(I8)') MAXLOC(ELEMS)
          write(temp2,'(I8)') NGL
          CALL FATAL_ERROR('NESTING ELEMENTS NUMBER'//trim(temp)//&
               & 'IS NOT IN THE GLOBAL DOMAIN',&
               & 'CHECK INPUT FILE NESTING ELEMENTS ARE 1 <= '//trim(temp2))
       END IF


    END IF

    ! BROADCAST TO OTHER PROCS
    IF(PAR)THEN
# if defined (MULTIPROCESSOR)
       SENDER = MSRID -1 ! SEND FROM MASTER
       CALL MPI_BCAST(NNODE,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(NELEM,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)

       IF( .NOT. MSR) THEN

          ALLOCATE(NODES(NNODE))
          ALLOCATE(ELEMS(NELEM))

       END IF


       CALL MPI_BCAST(NODES,NNODE,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
       CALL MPI_BCAST(ELEMS,NELEM,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
# endif
    END IF

    RETURN
  END SUBROUTINE LOAD_NESTING_OPB
  !============================================================================!
  !
  !============================================================================!   
  SUBROUTINE ALLOC_NEST ! CALL IN ALLOCATE ALL - USE THE MEM_COUNT!
    IMPLICIT NONE
    INTEGER :: STATUS
    INTEGER NDB
    !==============================================================================!
# if !defined (DOUBLE_PRECISION)
    NDB = 1       !!GWC BASE THIS ON KIND
# else
    NDB = 2
# endif

    IF(NESTING_MODE == "large_domain")THEN
       ALLOCATE(X_NEST(MT_NEST)      , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:X_NEST")

       ALLOCATE(Y_NEST(MT_NEST)      , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:Y_NEST")

       ALLOCATE(LON_NEST(MT_NEST)    , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:LON_NEST")

       ALLOCATE(LAT_NEST(MT_NEST)    , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:LAT_NEST")

       ALLOCATE(XMC_NEST(NT_NEST)    , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:XMC_NEST")

       ALLOCATE(YMC_NEST(NT_NEST)    , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:YMC_NEST")

       ALLOCATE(LONC_NEST(NT_NEST)   , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:LONC_NEST")

       ALLOCATE(LATC_NEST(NT_NEST)   , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:LATC_NEST")

       ALLOCATE(ZZ_NEST(MT_NEST,KB), stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:ZZ_NEST")

       ALLOCATE(Z_NEST(MT_NEST,KB) , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:Z_NEST")

       ALLOCATE(H_NEST(MT_NEST)      , stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:H_NEST")

       MEMCNT = MEMCNT + NT_NEST * 4 * NDB +  MT_NEST * (5 + 2*KB)* NDB
    END IF

    ! ALWAYS ALLOCATE NV_NEST
    ALLOCATE(NV_NEST(NT_NEST,3),stat=status)
    if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:VA_NEST")
    ! SET IN SETUP_NEST
    MEMCNT = MEMCNT + NT_NEST * 3

    IF(NESTING_MODE == "subdomain") THEN
       ! ALLOCATE SPACE FOR INTERPOLATED
       ALLOCATE(UA_NEST(NT_NEST),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:UA_NEST")

       ALLOCATE(VA_NEST(NT_NEST),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:VA_NEST")

       ALLOCATE(EL_NEST(MT_NEST),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:EL_NEST")

       ALLOCATE(U_NEST(NT_NEST,KB),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:U_NEST")

       ALLOCATE(V_NEST(NT_NEST,KB),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:V_NEST")

       ALLOCATE(S1_NEST(MT_NEST,KB),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:S1_NEST")

       ALLOCATE(T1_NEST(MT_NEST,KB),stat=status)
       if(status /= 0) call fatal_error("ALLOC_NEST: can not allocate:T1_NEST")

       MEMCNT = MEMCNT + NT_NEST * (2 + 2*KB)* NDB +  MT_NEST * (1 + 2*KB)* NDB

    END IF


    ! ALLOCATE SPACE FOR DATA BLOCK IO FILE
    ALLOCATE(UA_BLK(NT_NEST,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:UA_BLK")

    ALLOCATE(VA_BLK(NT_NEST,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:VA_BLK")

    ALLOCATE(EL_BLK(MT_NEST,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:EL_BLK")

    ALLOCATE(U_BLK(NT_NEST,KB,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:U_BLK")

    ALLOCATE(V_BLK(NT_NEST,KB,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:V_BLK")

    ALLOCATE(S1_BLK(MT_NEST,KB,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:S1_BLK")

    ALLOCATE(T1_BLK(MT_NEST,KB,NESTING_BLOCKSIZE),stat=status)
    if(status /= 0) call fatal_error("ALLOC_BLK: can not allocate:T1_BLK")

    MEMCNT = MEMCNT + NT_NEST * (2 + 2*KB)* NDB * NESTING_BLOCKSIZE &
         & +  MT_NEST * (1 + 2*KB)* NDB * NESTING_BLOCKSIZE

    CALL DEFINE_EXTRA_DIMENSIONS

    RETURN
  END SUBROUTINE ALLOC_NEST

  !==================================================================================!
  !
  !==================================================================================!
  Subroutine SETUP_NEST

    IMPLICIT NONE

    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCATT), POINTER :: ATT
    LOGICAL :: FOUND

    TYPE(TIME):: TIMETEST

    INTEGER, ALLOCATABLE :: EFOUND(:),NFOUND(:),TMPFOUND(:)
    INTEGER :: nTRI(3), I, J
    CHARACTER(LEN=10) :: STRNG


    IF(NESTING_MODE == "large_domain")THEN
    ! PUT DATA IN NEST GRID VARAIBLES FOR OUTPUT
       X_NEST=XM(MID_NEST(:))
       
       Y_NEST=YM(MID_NEST(:))

       LON_NEST= LON(MID_NEST(:))

       LAT_NEST= LAT(MID_NEST(:))

       XMC_NEST=XMC(NID_NEST(:))

       YMC_NEST=YMC(NID_NEST(:))

       LONC_NEST=LONC(NID_NEST(:))

       LATC_NEST=LATC(NID_NEST(:))

       ZZ_NEST=ZZ(MID_NEST(:),:)

       Z_NEST=Z(MID_NEST(:),:)

       H_NEST=H(MID_NEST(:))

!JQI       IF(MOD(NSTEPS,NESTING_BLOCKSIZE)/=0) THEN
       IF(MOD(NSTEPS-NESTING_BLOCKSIZE,NESTING_BLOCKSIZE-1)/=0) THEN
          
          IF(DBG_SET(DBG_LOG)) THEN
             WRITE(IPT,*) "MODEL STEPS=",nsteps
             WRITE(IPT,*) "NESTING BLOCKSIZE=",NESTING_BLOCKSIZE
             WRITE(IPT,*) "REMAINDER = ",MOD(NSTEPS,NESTING_BLOCKSIZE)
          END IF
          CALL FATAL_ERROR&
            &("THE NUMBER OF STEPS IN A NESTING MODEL RUN MUST BE",&
            & "MOD(NSTEPS-NESTING_BLOCKSIZE,NESTING_BLOCKSIZE-1)=0")
       END IF
!JQI       NESTING_BLOCKSIZE= NESTING_BLOCKSIZE +1
       


    ELSE
       ! SETUP SUBDOMAIN MODEL FOR INPUT
       CALL OPEN_NESTOUTPUT

       !=======================
       ! SANITY CHECK
       !=======================
       DIM => FIND_UNLIMITED(NC_NEST,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT",&
            & "FILE NAME: "//TRIM(NC_NEST%FNAME),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")


       NBLOCKS = DIM%DIM ! THE NUMBER OF BLOCKS IN THE FILE (TIME STEPS)


       ! GET THE FILE BLOCK SIZE
       DIM => FIND_DIM(NC_NEST,'nesting_block',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT",&
            & "FILE NAME: "//TRIM(NC_NEST%FNAME),&
            &"COULD NOT FIND THE DIMENSION 'nesting_block'")

       NESTING_BLOCKSIZE = DIM%DIM



       ! CHECK THE NUMBER OF NODES AND ELEMENTS
       DIM => FIND_DIM(NC_NEST,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT",&
            & "FILE NAME: "//TRIM(NC_NEST%FNAME),&
            &"COULD NOT FIND THE DIMENSION 'nele'")

       if(ngl_nest /= DIM%DIM) call fatal_error&
            ("The Nesting Element file does not match the data from the l&
            &arge domain netcdf file?")
       

       DIM => FIND_DIM(NC_NEST,'node',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT",&
            & "FILE NAME: "//TRIM(NC_NEST%FNAME),&
            &"COULD NOT FIND THE DIMENSION 'node'")

       if(mgl_nest /= DIM%DIM) call fatal_error&
            ("The Nesting Node file does not match the data from the l&
            &arge domain netcdf file?")

       ! GET THE TIME STEP FROM THE LARGE DOMAIN RUN
       VAR => FIND_VAR(NC_NEST,"delta_T",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT",&
            & "FILE NAME: "//TRIM(NC_NEST%FNAME),&
            &"COULD NOT FIND THE VARIABLE 'delta_T'")
       CALL NC_CONNECT_AVAR(VAR, LD_SEC)
       CALL NC_READ_VAR(VAR)
       write(ipt,*) "LD_SEC:",LD_SEC
       LD_TSTEP = SECONDS2TIME(LD_SEC)

       ! TEST THE START TIME


       TIMETEST = get_file_time(NC_NEST,1)

       CALL PRINT_TIME(TIMETEST,IPT,"BLOCKENDS")
       BLOCK_STARTS = TIMETEST - (NESTING_BLOCKSIZE -1)*LD_TSTEP

       CALL PRINT_TIME(BLOCK_STARTS,IPT,"BLOCKSTARTS")
       CALL PRINT_TIME(STARTTIME,IPT,"STARTTIME")


       CALL PRINT_TIME(LD_TSTEP,IPT,"STEP")
       WRITE(IPT,*) "NB:",NESTING_BLOCKSIZE

       IF(StartTime < (TIMETEST - (NESTING_BLOCKSIZE)*LD_TSTEP)) CALL FATAL_ERROR&
            &("THE SUBDOMAIN START TIME DOES NOT MATCH THE NESTING FILE START TIME")

       ! SET THE CURRENT TIME FOR THE END OF THE FIRST BLOCK!
       BLOCK_ENDS = TIMETEST
!JQI       BLOCK_STARTS = TIMETEST - (NESTING_BLOCKSIZE)*LD_TSTEP
       BLOCK_STARTS = TIMETEST - (NESTING_BLOCKSIZE-1)*LD_TSTEP

       ! TEST THE END TIME
       TIMETEST = get_file_time(NC_NEST,NBLOCKS)
!JQI       IF(TIMETEST < EndTime) CALL FATAL_ERROR &
!JQI            &("THE SUBDOMAIN END TIME DOES NOT MATCH THE NESTING FILE END TIME")
       



       !!========================================
       !! CONNECT MEMORY AND READ FIRST TIME STEP
       !!=======================================
       BLK_INDEX = 1

       VAR => FIND_VAR(NC_NEST,"zeta",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'zeta' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, EL_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_EL => VAR

       VAR => FIND_VAR(NC_NEST,"ua",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'ua' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, UA_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_UA => VAR
       
       VAR => FIND_VAR(NC_NEST,"va",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'va' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, VA_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_VA=> VAR

       VAR => FIND_VAR(NC_NEST,"u",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'u' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, U_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_U=> VAR

       VAR => FIND_VAR(NC_NEST,"v",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'v' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, V_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_V=> VAR

       VAR => FIND_VAR(NC_NEST,"temp",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'temp' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, T1_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_T1=> VAR

       VAR => FIND_VAR(NC_NEST,"salinity",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN NESTING FILE OBJECT NO 'salinity' FOUND?")
       CALL NC_CONNECT_AVAR(VAR, S1_BLK)
       CALL NC_READ_VAR(VAR,BLK_INDEX)
       VAR_S1=> VAR

    END IF


    ! CREATE A TEMPRORARY LIST OF NODES TO CHECK AGAINST
    ALLOCATE(EFOUND(NT_NEST),NFOUND(MT_NEST),TMPFOUND(MT_NEST))
    EFOUND=0
    NFOUND=0
    TMPFOUND=0

    ! MAKE SURE NO NODE OR ELEMENT IS NAMED TWICE!
    DO I =1,NT_NEST ! ELEMS
       WHERE(NID_NEST==NID_NEST(I))
          EFOUND = EFOUND +1
       END WHERE
    END DO

    DO I =1,NT_NEST
       IF(EFOUND(I)/=1) THEN
          WRITE(STRNG,'(I8)') EGID_X(NID_NEST(I))
          CALL FATAL_ERROR("THE GLOBAL ELEMENT # "//TRIM(STRNG)//&
               & "APPEARS MORE THAN ONCE IN THE NESTING LIST?")
       END IF
    END DO


    DO I =1,MT_NEST ! NODES
       WHERE(MID_NEST==MID_NEST(I))
          NFOUND = NFOUND +1
       END WHERE
    END DO

    DO I =1,MT_NEST
       IF(NFOUND(I)/=1) THEN
          WRITE(STRNG,'(I8)') NGID_X(MID_NEST(I))
          CALL FATAL_ERROR("THE GLOBAL NODE # "//TRIM(STRNG)//&
               & "APPEARS MORE THAN ONCE IN THE NESTING LIST?")
       END IF
    END DO

    NFOUND = 0
    EFOUND = 0

    DO I =1,NT_NEST
       ! PUT THE NODES FOR THAT CELL IN THE NV_NEST LIST
       NTRI(1:3) = NV(NID_NEST(I),1:3)

       ! MAKE SURE THOSE NODES ARE IN THE NESTING LIST
       DO J =1,3
          TMPFOUND = 0
          WHERE(MID_NEST==NTRI(J))
             TMPFOUND = 1

             NFOUND = 1
          END WHERE

          IF(SUM(TMPFOUND)/=1)THEN
             WRITE(STRNG,'(I8)') NGID_X(nTRI(J))
             CALL FATAL_ERROR("The global node # "//TRIM(STRNG)//" is in a nesting Element,",&
                  &"but not in the nesting node list?")
          END IF

       END DO

       NV_NEST(I,:) = NGID_X(NTRI)

    END DO


    DO I =1,MT_NEST
       IF(NFOUND(I)/=1) THEN

          WRITE(STRNG,'(I8)') NGID_X(MID_NEST(I))
          CALL FATAL_ERROR("The global node # "//TRIM(STRNG)//" is the nesting node list,",&
               &"but not in any nesting element?")
       END IF
    END DO



  END Subroutine SETUP_NEST
  !==================================================================================!
  !
  !==================================================================================!
  FUNCTION NESTING_GRID_OBJECT() RESULT(NCF)
    USE MOD_CLOCK
    IMPLICIT NONE
    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT
    character(len=100)    :: timestamp, temp

    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED)THEN
       IOPROC_ALLOCATED = .TRUE.

       ALLOCATE(X_NEST(MGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:X_NEST")
       X_NEST = 0.0_SP

       ALLOCATE(Y_NEST(MGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Y_NEST")
       Y_NEST = 0.0_SP

       ALLOCATE(LON_NEST(MGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LON_NEST")
       LON_NEST = 0.0_SP

       ALLOCATE(LAT_NEST(MGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LAT_NEST")
       LAT_NEST = 0.0_SP

       ALLOCATE(XMC_NEST(NGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:XMC_NEST")
       XMC_NEST = 0.0_SP

       ALLOCATE(YMC_NEST(NGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:YMC_NEST")
       YMC_NEST = 0.0_SP

       ALLOCATE(LONC_NEST(NGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LONC_NEST")
       LONC_NEST = 0.0_SP

       ALLOCATE(LATC_NEST(NGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:LATC_NEST")
       LATC_NEST = 0.0_SP

       ALLOCATE(ZZ_NEST(MGL_NEST,KB),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:ZZ_NEST")
       ZZ_NEST = 0.0_SP

       ALLOCATE(Z_NEST(MGL_NEST,KB),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:Z_NEST")
       Z_NEST = 0.0_SP

       ALLOCATE(H_NEST(MGL_NEST),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:H_NEST")
       H_NEST = 0.0_SP

    END IF

    ! ALLOCATE THE NEW FILE OBJECT
    NCF => NEW_FILE()

    ATT => NC_MAKE_ATT(name='source',values=trim(fvcom_version)) 
    NCF => ADD(NCF,ATT)

    call get_timestamp(temp)
    timestamp = 'model started at: '//trim(temp)

    ATT => NC_MAKE_ATT(name='history',values=trim(timestamp)) 
    NCF => ADD(NCF,ATT)

    ! NEST NODES
    VAR  => NC_MAKE_PVAR(name='nesting_nodes', values=MID_NEST_OUTPUT&
         &, DIM1= DIM_node_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nesting Node Number of Large Domain') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='nest_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NEST ELEMS
    VAR  => NC_MAKE_PVAR(name='nesting_elems', values=NID_NEST_OUTPUT&
         &, DIM1= DIM_nele_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='Nesting Element Number of Large Domain') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='nest_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NEST X
    VAR  => NC_MAKE_AVAR(name='x', values=X_nest&
         &, DIM1= DIM_node_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal x-coordinate')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='nest_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NEST Y
    VAR  => NC_MAKE_AVAR(name='y', values=Y_nest&
         &, DIM1= DIM_node_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal y-coordinate')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='nest_grid') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NEST LON
    VAR  => NC_MAKE_AVAR(name='lon',values=lon_nest,DIM1=DIM_node_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)

    ! NEST LAT
    VAR  => NC_MAKE_AVAR(name='lat',values=lat_nest,DIM1=DIM_node_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodal latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! NEST XMC
    VAR  => NC_MAKE_AVAR(name='xc',values=xmc_nest,DIM1=DIM_nele_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='zonal x-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)


    ! NEST YMC
    VAR  => NC_MAKE_AVAR(name='yc',values=ymc_nest,DIM1=DIM_nele_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='zonal y-coordinate') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! NEST LONC
    VAR  => NC_MAKE_AVAR(name='lonc',values=lonc_nest,DIM1=DIM_nele_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='zonal longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='longitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_east') 
    VAR  => ADD(VAR,ATT)
    NCF  => ADD(NCF,VAR)


    ! NEST LATC
    VAR  => NC_MAKE_AVAR(name='latc',values=latc_nest,DIM1=DIM_nele_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='zonal latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='standard_name',values='latitude') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='degrees_north') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! nest siglay
    VAR  => NC_MAKE_AVAR(name='siglay', values=zz_nest, DIM1= DIM_node_nest,&
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


    ! nest siglev
    VAR  => NC_MAKE_AVAR(name='siglev',&
         & values=z_nest, DIM1= DIM_node_nest, DIM2= DIM_siglev )

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


    ! nest h
    VAR  => NC_MAKE_AVAR(name='h',&
         & values=h_nest, DIM1= DIM_node_nest)

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

    ! nest nv
    VAR  => NC_MAKE_AVAR(name='nv',&
         & values=nv_nest, DIM1= DIM_nele_nest, DIM2= DIM_three)

    ATT  => NC_MAKE_ATT(name='long_name',values='nodes surrounding element') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! IINT
    VAR  => NC_MAKE_AVAR(name='iint',&
         & values=IO_IINT, DIM1= DIM_time_nest)

    ATT  => NC_MAKE_ATT(name='long_name',values='internal mode iteration number') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! time
    VAR => FLOAT_TIME_OBJECT &
         &(Data=IO_DAYS, &
         & USE_MJD=use_real_world_time, &
         & DIM=DIM_TIME_nest)

    NCF  => ADD(NCF,VAR)


    ! Itime
    VAR  => ITIME_OBJECT &
         &(Data=IO_MJD, &
         & Use_MJD=use_real_world_time, &
         & DIM=DIM_TIME_nest)

    NCF  => ADD(NCF,VAR)

    ! Itime2
    VAR => ITIME2_OBJECT &
         &(Data=IO_MSEC, &
         & Use_MJD=use_real_world_time, &
         & DIM=DIM_TIME_nest)

    NCF => ADD(NCF,VAR)

    IF (use_real_world_time) THEN

       VAR => DATETIME_OBJECT &
            &(Data=IO_timestr,&
            & DIMSTR=DIM_DateStrLen,&
            & DIMTIME=DIM_TIME_nest,&
            TIMEZONE=TIMEZONE)

       NCF  => ADD(NCF,VAR)
    END IF


    LD_SEC = seconds(imdti)
    VAR  => NC_MAKE_AVAR(name='delta_T', values=LD_SEC)   
    ATT  => NC_MAKE_ATT(name='long_name',values='Large Domain Time Step')
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='seconds')
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END NESTING_GRID_OBJECT"

  END FUNCTION NESTING_GRID_OBJECT
  !==================================================================================!
  !
  !==================================================================================!
  FUNCTION NESTING_FILE_OBJECT() RESULT(NCF)
    USE MOD_CLOCK
    IMPLICIT NONE
    INTEGER :: status
    LOGICAL, SAVE :: IOPROC_ALLOCATED = .FALSE.
    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR),  POINTER :: VAR
    TYPE(NCATT),  POINTER :: ATT

    IF(IOPROC .AND. .NOT. IOPROC_ALLOCATED)THEN
       IOPROC_ALLOCATED = .TRUE.

       ALLOCATE(EL_BLK(MGL_NEST,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL_BLK")
       EL_BLK = 0.0_SP

       ALLOCATE(UA_BLK(NGL_NEST,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:UA_BLK")
       UA_BLK = 0.0_SP

       ALLOCATE(VA_BLK(NGL_NEST,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:VA_BLK")
       VA_BLK = 0.0_SP

       ALLOCATE(U_BLK(NGL_NEST,KB,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:U_BLK")
       U_BLK = 0.0_SP

       ALLOCATE(V_BLK(NGL_NEST,KB,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:V_BLK")
       V_BLK = 0.0_SP

       ALLOCATE(S1_BLK(MGL_NEST,KB,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:S1_BLK")
       S1_BLK = 0.0_SP

       ALLOCATE(T1_BLK(MGL_NEST,KB,NESTING_BLOCKSIZE),STAT=STATUS) 
       IF (STATUS /=0 )  &
            & CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:T1_BLK")
       T1_BLK = 0.0_SP

    END IF

    NCF => NEW_FILE()

    ! nest zeta
    VAR  => NC_MAKE_AVAR(name='zeta',&
         & values=EL_BLK, DIM1= DIM_node_nest, DIM2 = DIM_BLK_NEST, &
	 & DIM3= DIM_TIME_NEST)

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

    ! UA_BLK
    VAR  => NC_MAKE_AVAR(name='ua',&
         & values=UA_BLK, DIM1= DIM_nele_nest, DIM2= DIM_BLK_NEST, &
	 & DIM3 = DIM_TIME_NEST)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertically Averaged x-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)


    ! VA_BLK
    VAR  => NC_MAKE_AVAR(name='va',&
         & values=VA_BLK, DIM1= DIM_nele_NEST, DIM2= DIM_BLK_NEST, &
	 & DIM3 = DIM_TIME_NEST)

    ATT  => NC_MAKE_ATT(name='long_name',values='Vertically Averaged y-velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! U_BLK
    VAR  => NC_MAKE_AVAR(name='u', values=U_BLK, DIM1= DIM_nele_nest,&
         & DIM2= DIM_siglay, DIM3= DIM_BLK_NEST, DIM4= DIM_TIME_NEST)
    ATT  => NC_MAKE_ATT(name='long_name',values='Eastward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! V_BLK
    VAR  => NC_MAKE_AVAR(name='v', values=V_BLK, DIM1= DIM_nele_nest, &
         & DIM2= DIM_siglay, DIM3 = DIM_BLK_NEST,DIM4 = DIM_TIME_NEST)

    ATT  => NC_MAKE_ATT(name='long_name',values='Northward Water Velocity') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='units',values='meters s-1') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid',values='fvcom_grid') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)

    NCF  => ADD(NCF,VAR)

    ! T_BLK
    VAR  => NC_MAKE_AVAR(name='temp',&
         & values=T1_BLK, DIM1= DIM_node_NEST, DIM2= DIM_siglay, DIM3 = DIM_BLK_NEST, &
	 & DIM4 = DIM_time_NEST)

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


    ! S_BLK
    VAR  => NC_MAKE_AVAR(name='salinity',&
         & values=S1_BLK, DIM1= DIM_node_NEST, DIM2= DIM_siglay, DIM3 = DIM_BLK_NEST, &
	 & DIM4 = DIM_time_NEST)

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

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END NESTING_FILE_OBJECT"

  END FUNCTION NESTING_FILE_OBJECT
  !===========================================================================!
  !
  !===========================================================================!  
  SUBROUTINE ARCHIVE_NEST
    IMPLICIT NONE
!JQI    logical, save :: NEED_INIT = .TRUE.
    INTEGER :: STATUS
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START ARCHIVE_NEST"

    IF(NESTING_ON)THEN

       ! SETUP THE OUTPUT IF THIS HAS NOT BEEN DONE YET
!JQI       if(NEED_INIT) then
       if(NEED_INIT_NEST) then
          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(INITNEST_CODE)
          ELSE
             CALL CALL_FUNC(INITNEST_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad INITNEST_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF

          BLK_INDEX = 0

!JQI          NEED_INIT = .FALSE.
          NEED_INIT_NEST = .FALSE.
       end if

       ! STORE DATA EVERYTIME
       CALL ASSIGN2BLOCK

       ! WRITE THE FILE IF IT IS TIME
       IF(BLK_INDEX == NESTING_BLOCKSIZE) then

          IF(USE_MPI_IO_MODE) THEN
             CALL MPI_IO_SYNCHRONIZE(NESTING_CODE)
          ELSE
             CALL CALL_FUNC(NESTING_CODE,status)
             IF (status/=0) call fatal_error("ARCHIVE:: Bad NESTING_CODE",&
                  & "Could not retrieve valid function pointer?")
          END IF

          CALL RESET_BLOCK

       END IF

    END IF


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END ARCHIVE_NEST"
  END SUBROUTINE ARCHIVE_NEST
  !==========================================================================!
  !
  !==========================================================================!
  SUBROUTINE OPEN_NESTOUTPUT
    USE CONTROL
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF

    character(len=160) :: pathnfile

    Nullify(NCF)

    IF(NESTING_ON .AND. NESTING_MODE == "subdomain") then
       !NESTING OUTPUT FILE
       ! LOOK IN OUTPUT DIR FOR PREVIOUS OUTPUT
       pathnfile= trim(OUTPUT_DIR)//TRIM(NESTING_FILE_NAME)

       ! OPEN THE FILE AND LOAD 
       NCF => NEW_FILE()
       NCF%FNAME=trim(pathnfile)
       Call NC_OPEN(NCF)
       CALL NC_LOAD(NCF)
       NC_NEST => NCF

    END IF

  END SUBROUTINE OPEN_NESTOUTPUT
  !========================================================================!
  !
  !========================================================================!
  SUBROUTINE DUMP_NC_NEST
    IMPLICIT NONE
    CALL DUMP_DATA(NC_NEST)
  END SUBROUTINE DUMP_NC_NEST
  !========================================================================!
  !
  !========================================================================!  
  SUBROUTINE SETUP_NESTFILE
    USE CONTROL
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    TYPE(TIME) :: DUMMY
    character(len=80) :: tmp
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START SETUP_NESTFILE"

    ! ALLOCATE THE NEW FILE OBJECT
    NC_NEST => NEW_FILE()
    NC_NEST%FNAME = NESTING_FILE_NAME


    NC_NEST => ADD(NC_NEST,NESTING_GRID_OBJECT() )

    NC_NEST => ADD(NC_NEST,NESTING_FILE_OBJECT() )


    ALLOCATE(NC_NEST%FTIME)
    CALL NC_WRITE_FILE(NC_NEST)
    NC_NEST%FTIME%NEXT_STKCNT=1

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END SETUP_NESTFILE"
  END SUBROUTINE SETUP_NESTFILE
  !========================================================================!
  !
  !========================================================================!
  SUBROUTINE ASSIGN2BLOCK
    IMPLICIT NONE

    BLK_INDEX = BLK_INDEX + 1

    UA_BLK(:,BLK_INDEX) = UA(NID_NEST(:))
    VA_BLK(:,BLK_INDEX) = VA(NID_NEST(:))
    EL_BLK(:,BLK_INDEX) = EL(MID_NEST(:))
    U_BLK(:,:,BLK_INDEX) = U(NID_NEST(:),:)
    V_BLK(:,:,BLK_INDEX) = V(NID_NEST(:),:)
    S1_BLK(:,:,BLK_INDEX) = S1(MID_NEST(:),:)
    T1_BLK(:,:,BLK_INDEX) = T1(MID_NEST(:),:)

  END SUBROUTINE ASSIGN2BLOCK
  !========================================================================!
  SUBROUTINE RESET_BLOCK
    IMPLICIT NONE

    IF(BLK_INDEX /= NESTING_BLOCKSIZE) CALL FATAL_ERROR("NESTING : WHAT THE HECK IS GOING ON?")

    UA_BLK(:,1) =  UA_BLK(:,BLK_INDEX)
    VA_BLK(:,1) =  VA_BLK(:,BLK_INDEX)
    EL_BLK(:,1) =  EL_BLK(:,BLK_INDEX)
    U_BLK(:,:,1) =  U_BLK(:,:,BLK_INDEX)
    V_BLK(:,:,1) =  V_BLK(:,:,BLK_INDEX)
    S1_BLK(:,:,1) =  S1_BLK(:,:,BLK_INDEX)
    T1_BLK(:,:,1) =  T1_BLK(:,:,BLK_INDEX)


    BLK_INDEX =1

  END SUBROUTINE RESET_BLOCK
  !========================================================================!  
  !
  !========================================================================!    
  SUBROUTINE SET_VAR_1D(NOW,BLK_DAT,DAT)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW    
    REAL(SP), ALLOCATABLE :: BLK_DAT(:,:),DAT(:)
    REAL(DP)     :: denom, numer
    REAL(DP)     :: prev_w, next_w
    INTEGER      :: PREV_I, NEXT_I, STEP

    TYPE(TIME) :: TDIFF

    ! SEE IF WE NEED NEW DATA
    IF(NOW > BLOCK_ENDS) THEN
       CALL  LOAD_NEXT_BLOCK
    ELSEIF(NOW < BLOCK_STARTS) THEN
       CALL FATAL_ERROR("THIS IS UNEXPECTED: MODIFY MOD_NESTING TO HANDLE THIS IF YOU NEED IT!")
    END IF

    ! FIND WHICH INDEX IN THE BLOCK BRACKETS NOW
    TDIFF = NOW - BLOCK_STARTS

    STEP = INT(seconds(TDIFF)/seconds(LD_TSTEP))
!JQI    PREV_I = STEP + 1 ! INDEX STARTS FROM ONE
    PREV_I = STEP         ! INDEX STARTS FROM ONE
    NEXT_I = PREV_I + 1

    ! GET THE INTERPOLATION WEIGHTS

    ! THE DIFFERENCE BETWEEN THE CURRENT TIME IN THE MODEL AND THE PREVIOUS DATA
    NUMER = SECONDS(NOW - (BLOCK_STARTS + STEP * LD_TSTEP))
    DENOM = SECONDS(LD_TSTEP)

    ! TAKE THE RATIO IN DOUBLE PRECISION AND CONVERT IF MODEL IS NOT DOUBLE
    NEXT_W = NUMER/DENOM
    PREV_W = 1.0_DP - NEXT_W

    ! GET THE DATA

    DAT = NEXT_W * BLK_DAT(:,NEXT_I) + PREV_W * BLK_DAT(:,PREV_I)
    
  END SUBROUTINE SET_VAR_1D
  !========================================================================!  
  !
  !========================================================================!    
  SUBROUTINE SET_VAR_2D(NOW,BLK_DAT,DAT)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW    
    REAL(SP), ALLOCATABLE :: BLK_DAT(:,:,:),DAT(:,:)
    REAL(DP)     :: denom, numer
    REAL(DP)     :: prev_w, next_w
    INTEGER      :: PREV_I, NEXT_I, STEP

    TYPE(TIME) :: TDIFF

    ! SEE IF WE NEED NEW DATA
    IF(NOW > BLOCK_ENDS) THEN
       CALL  LOAD_NEXT_BLOCK
    ELSEIF(NOW < BLOCK_STARTS) THEN
       CALL FATAL_ERROR("THIS IS UNEXPECTED: MODIFY MOD_NESTING TO HANDLE THIS IF YOU NEED IT!")
    END IF

    ! FIND WHICH INDEX IN THE BLOCK BRACKETS NOW
    TDIFF = NOW - BLOCK_STARTS

    STEP = INT(seconds(TDIFF)/seconds(LD_TSTEP))
!JQI    PREV_I = STEP + 1 ! INDEX STARTS FROM ONE
    PREV_I = STEP         ! INDEX STARTS FROM ONE
    NEXT_I = PREV_I + 1

    ! GET THE INTERPOLATION WEIGHTS

    ! THE DIFFERENCE BETWEEN THE CURRENT TIME IN THE MODEL AND THE PREVIOUS DATA
    NUMER = SECONDS(NOW - (BLOCK_STARTS + STEP * LD_TSTEP))
    DENOM = SECONDS(LD_TSTEP)

    ! TAKE THE RATIO IN DOUBLE PRECISION AND CONVERT IF MODEL IS NOT DOUBLE
    NEXT_W = NUMER/DENOM
    PREV_W = 1.0_DP - NEXT_W

    ! GET THE DATA

    DAT = NEXT_W * BLK_DAT(:,:,NEXT_I) + PREV_W * BLK_DAT(:,:,PREV_I)
    
  END SUBROUTINE SET_VAR_2D

  SUBROUTINE LOAD_NEXT_BLOCK
    IMPLICIT NONE
    BLK_INDEX = BLK_INDEX + 1

    BLOCK_ENDS   = get_file_time(NC_NEST,BLK_INDEX)
!JQI    BLOCK_STARTS = BLOCK_ENDS - (NESTING_BLOCKSIZE)*LD_TSTEP
    BLOCK_STARTS = BLOCK_ENDS - (NESTING_BLOCKSIZE-1)*LD_TSTEP
    CALL NC_READ_VAR(VAR_EL,BLK_INDEX)
    CALL NC_READ_VAR(VAR_UA,BLK_INDEX)
    CALL NC_READ_VAR(VAR_VA,BLK_INDEX)
    CALL NC_READ_VAR(VAR_U,BLK_INDEX)
    CALL NC_READ_VAR(VAR_V,BLK_INDEX)
    CALL NC_READ_VAR(VAR_T1,BLK_INDEX)
    CALL NC_READ_VAR(VAR_S1,BLK_INDEX)
   
  END SUBROUTINE LOAD_NEXT_BLOCK


END Module Mod_Nesting
