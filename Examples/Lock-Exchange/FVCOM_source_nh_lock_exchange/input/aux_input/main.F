PROGRAM AUXILIARY_INPUT
  use mod_input
  use mod_auxiliary
  use mod_setup
  USE MOD_SET_TIME
  USE MOD_REPORT
  IMPLICIT NONE
  INTEGER :: STATUS, I
  CHARACTER(LEN=4) :: BFLAG
  character(len=*),parameter::CVS_Id="$Id: main.F,v 1.1.2.1 2008/04/11 15:24:41 dstuebe Exp $" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date: 2008/04/11 15:24:41 $" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name: New_Input $" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision: 1.1.2.1 $" ! [sng] File revision string



  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("Other Input")

# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")


  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  IF (DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*) "! ========================="
     WRITE(IPT,*) "!   BEGIN RESTART CREATE"
     WRITE(IPT,*) "! ========================="
  END IF

  CALL NAME_LIST_INITIALIZE

  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS

  IF(USE_REAL_WORLD_TIME)THEN

     NOW = READ_DATETIME(START_DATE,DATE_FORMAT,TIMEZONE,status)
     IF(.not. STATUS) CALL FATAL_ERROR&
          & ('Bad Start Date format!', &
          & TRIM(START_DATE))
  ELSE
     CALL IDEAL_TIME_STRING2TIME(START_DATE,BFLAG,NOW,IINT)
     IF(BFLAG == 'step') CALL FATAL_ERROR&
          &("You must specify a time, not a step, for this restart file", &
          & "The Step will be set by the old restart file...")

  END IF


  CALL PRINT_REAL_TIME(NOW,IPT,"Restart Time Is:")

  CALL CHECK_IO_DIRS

  CALL OPEN_FILES

  CALL READ_COLDSTART_SIGMA
  CLOSE(SIGMAUNIT)

  CALL READ_COLDSTART_GRID(GRIDUNIT,MGL,NGL,NV)
  m = MGL
  mt = MGL
  n = ngl
  nt = ngl

  ALLOCATE(NVGL(0:NGL,3)); NVGL=0
  NVGL(1:NGL,1:3) = NV(1:NGL,1:3)

  CALL ALLOCATE_SPACE

  ALLOCATE(X_LCL(0:MGL),Y_LCL(0:MGL))
  CALL READ_COLDSTART_COORDS(GRIDUNIT,MGL,X_LCL,Y_LCL)
  CLOSE(GRIDUNIT)

  CALL COORDINATE_UNITS(X_LCL,Y_LCL)
  CALL SETUP_CENTER_COORDS

  CALL READ_COLDSTART_DEPTH(DEPTHUNIT,MGL,X_LCL,Y_LCL,H)
  CALL N2E2D(H,H1)
  CLOSE(DEPTHUNIT)

  KBM1 = KB - 1
  KBM2 = KB - 2
  CALL Setup_Sigma
  CALL SETUP_SIGMA_DERIVATIVES

  DEALLOCATE(X_LCL,Y_LCL)
!====================================================================
! CREATE YOUR SPECIALIZED INPUT DATA HERE AND SAVE TO A FILE!

  CALL CREATE_HVC(NN_HVC,CC_HVC)

  CALL CREATE_BOTTOM_ROUGHNESS(CC_Z0B)

  CALL DUMP_AUX

!====================================================================
  IF (DBG_SET(DBG_LOG)) WRITE(IPT,*) "! TADA!!!"
  CALL PSHUTDOWN


END PROGRAM AUXILIARY_INPUT
