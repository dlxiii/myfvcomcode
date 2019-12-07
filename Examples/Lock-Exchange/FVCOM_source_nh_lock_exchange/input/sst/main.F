PROGRAM sst2grd
  USE MOD_NCDIO, only : update_iodata
  USE MOD_SST2GRD
  USE LIMS
  USE CONTROL, only : IPT
  USE MOD_UTILS
  USE MOD_INTERP
  USE MOD_INPUT
  USE mod_set_time

  ! SEE PARAMETER SETTINGS IN MOD_SST2GRD.F90
  IMPLICIT NONE
  INTEGER :: STATUS,value
  character(len=4) :: bflag
  

  character(len=*),parameter::CVS_Id="$Id: main.F,v 1.1.2.3 2008/01/11 12:42:46 necofs Exp $" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date: 2008/01/11 12:42:46 $" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name: New_Input $" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision: 1.1.2.3 $" ! [sng] File revision string



  WRITE(IPT,*) "BEGIN SST INTERP:"

  ! INTIALIZE MEMORY FROM libfvcom.a
  CALL INITIALIZE_CONTROL("SST2GRID")


# if defined (MULTIPROCESSOR)
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
# endif

  IF (PAR) CALL FATAL_ERROR &
       & ("THIS PROGRAM IS WAY TO SHORT TO BOTHER WITH WRITING AS A PARALLE CODE.",&
       & "STOP AND RERUN WITH JUST ONE PROCESSOR!")


  CALL GET_COMMANDLINE(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)

  WRITE(IPT,*) "BEGIN SST INTERP:"

  CALL READ_NAMELIST

  CALL CNTRL_PRMTRS
     
  if(USE_REAL_WORLD_TIME) then
       
       NOW = READ_DATETIME(Start_Date,"YMD",TIMEZONE,status)
       IF(.not.status) call fatal_error &
            &("could not parse time_origin or time_zone passed for spectral tidal forcing file?")
       
  CALL PRINT_REAL_TIME(NOW,IPT,'START DATE')

    else

       CALL IDEAL_TIME_STRING2TIME(Start_date,BFLAG,NOW,IINT)
       IF(BFLAG == 'step') CALL FATAL_ERROR&
            &("You must specify a time, not a step, for this restart file", &
            & "The Step will be set by the old restart file...")

       CALL PRINT_TIME(NOW,IPT,'START DATE')

    END if




  
! OPEN AND READ THE FVCOM GRID FILE
  FNAME = TRIM(GRID_DIR)//TRIM(GRID_FILE)
  WRITE(IPT,*) "OPENING GRIDFILE: "//TRIM(FNAME)
  Call FOPEN(GRIDUNIT,TRIM(FNAME),'cfr')
  CALL READ_COLDSTART_GRID(GRIDUNIT,M,N,NV)
  ALLOCATE(XX(M),YY(M))
  CALL READ_COLDSTART_COORDS(GRIDUNIT,M,XX,YY)

  write(*,*) "READ FVCOM GRID xx, yy: "
  WRITE(*,*) "min/max(xx)",minval(xx),maxval(xx)
  WRITE(*,*) "min/max(lat)",minval(yy),maxval(yy)

  ! CREATE THE OUTPUT FILE TREE
  CALL MY_OUTFILE

  CALL PRINT_FILE(NC_OUT)
  CALL NC_WRITE_FILE(NC_OUT)

  IF (.not. CONSTANT) THEN
     
     ! open the grid file for the sst data longitude
     CALL FOPEN(FUNIT,TRIM(INPUT_DIR)//TRIM(LON_DATA_NAME),'cfr')
     DO J=1,NY
        READ(FUNIT,*) (LON(I,J), I=1,NX)
     ENDDO
     CLOSE(FUNIT)
     
     ! open the grid file for the sst data latitude
     CALL FOPEN(FUNIT,TRIM(INPUT_DIR)//TRIM(LAT_DATA_NAME),'cfr')
     
     DO J=1,NY
        READ(FUNIT,*) (LAT(I,J), I=1,NX)
     ENDDO
     CLOSE(FUNIT)
     
     
     write(*,*) "READ SST DATA lon and lat: "
     WRITE(*,*) "min/max(lon)",minval(lon),maxval(lon)
     WRITE(*,*) "min/max(lat)",minval(lat),maxval(lat)
     
     CALL DEGREES2METERS(LON,LAT,PROJECTION_REFERENCE,OX,OY,nx,ny)
     
     write(*,*) "CONVERTED TO METERS: "  
     WRITE(*,*) "min/max(OX)",minval(OX),maxval(OX)
     WRITE(*,*) "min/max(OY)",minval(OY),maxval(OY)
     
     
     
     ! READ IN FIRST SST FILE AND CREATE INTERP WEIGHTS
     CALL CREATE_INTERP
     
     
     
     ! NOW LOAD EACH DATA FILE AND WRITE OUT RESULTS
     IDAY = 0
     DO   !365 or 366
        IDAY = IDAY + 1
        WRITE(*,*) "PROCESSING IDAY: ",IDAY
        
        ! input data
        write(CH,'(i3.3)') iday  !
        FNAME = TRIM(INPUT_DIR)//TRIM(DATA_PREFIX)//CH
        INQUIRE(FILE=FNAME,EXIST=FEXIST)
        IF(.NOT. FEXIST) THEN
           CALL WARNING("COULD NOT FIND FILE: "//TRIM(FNAME),&
                & "Incriment Year!:")
           
           iday=1
           write(CH,'(i3.3)') iday
           read(DATA_PREFIX(4:7),*) value
           value = value + 1
           write(ipt,*) "Year=",value
           write(DATA_PREFIX(4:7),'(I4.4)') value


           FNAME = TRIM(INPUT_DIR)//TRIM(DATA_PREFIX)//CH
           INQUIRE(FILE=FNAME,EXIST=FEXIST)
           IF(.NOT. FEXIST)  THEN
              WRITE(IPT,*) "COULD NOT FIND ANY MORE DATA!"
              Write(IPT,*) "Last File Name: ",TRIM(FNAME)
              CALL PSHUTDOWN
           END IF
           
        END IF
        
        CALL FOPEN(FUNIT,TRIM(FNAME),'cfr')
        
        
        DO J=1,NY
           READ(FUNIT,*) (OSST(I,J),I=1,NX)   
        ENDDO
        CLOSE(FUNIT)
        
        ! USE MASK TO CREATE LINEAR ARRAY
        CNT =0 
        DO J=1,NY 
           DO I=1,NX
              IF(mask(I,J) == 1) then
                 CNT = CNT +1
                 Z(CNT) = OSST(I,J)
              ENDIF
           ENDDO
        ENDDO
        
        DO I = 1,WEIGHTS%Nin
           TMP(I) = Z(WEIGHTS%INDEX(I))
        END DO
        Z = TMP
        
        DO I = 1,WEIGHTS%Nout
           CALL INTERP_WEIGH(Z,WEIGHTS%PTW(I),SST(I))       
        END DO
        
        mean = sum(sst,1)/size(sst,1)

        CALL UPDATE_IODATA(NC_OUT,NOW)
        write(*,*) "TIME IS: "//TRIM(IO_TIMESTR)
        PRINT*,"sst: min/max/mean", minval(sst),maxval(sst),mean
        
        
        FTM => NC_OUT%FTIME
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL NC_WRITE_FILE(NC_OUT)
        
        NOW = NOW + SST_INTERVAL
        
        !---------------------------------------------------
        ! output sst in grid 
        !---------------------------------------------------
        
        
        
        
     ENDDO   !
     
  ELSE ! IF IT SHOULD JUST BE A CONSTANT VALUE...

     SST = cval

     DO I = 1,365
     
        call update_ioData(NC_OUT,NOW)
        
        FTM => NC_OUT%FTIME
        FTM%NEXT_STKCNT = FTM%NEXT_STKCNT + 1
        CALL NC_WRITE_FILE(NC_OUT)
        
        NOW = NOW + SST_INTERVAL
        
     END DO
  END IF

     
END PROGRAM sst2grd

