MODULE MOD_binary
  USE MOD_BULK
  USE MOD_SURFACEFORCE
  USE MOD_SET_TIME
  implicit none
  
  TYPE bin_data
     type(TIME) :: dtm
     REAL(SP), POINTER :: data(:)
  END TYPE bin_data
  
  TYPE(bin_data), POINTER :: SWR_NEXT, SWR_PREV
  TYPE(bin_data), POINTER :: NHF_NEXT, NHF_PREV
  TYPE(bin_data), POINTER :: WNDX_NEXT, WNDX_PREV
  TYPE(bin_data), POINTER :: WNDY_NEXT, WNDY_PREV
  TYPE(bin_data), POINTER :: EVP_NEXT, EVP_PREV
  TYPE(bin_data), POINTER :: PRC_NEXT, PRC_PREV

  !  STORAGE FOR OUTPUT OF BINARY FILES
  REAL(SP), POINTER :: WNDYGL(:), WNDXGL(:)
  REAL(SP), POINTER :: SWRGL(:),  NHFGL(:)
  REAL(SP), POINTER :: PRCGL(:),  EVPGL(:)

  
CONTAINS

  SUBROUTINE UPDATE_BINARY(NOW)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: UPDATE_BINARY"

    IF(PRECIPITATION_ON) THEN
       CALL UPDATE_EVP(NOW,EVP,PRC)
    END IF

    IF(WIND_ON) THEN

       IF(WIND_TYPE == 'stress') THEN
          CALL UPDATE_WND(NOW,STRX,STRY)
       ELSE IF(WIND_TYPE == 'speed') THEN
          CALL UPDATE_WND(NOW,SPDX,SPDY)
          CALL SIMPLE_DRAG(SPDX,SPDY,STRX,STRY)
       END IF
       
!    WRITE(IPT,*) "MIN/MAX(SPDX)::",minval(spdx),maxval(spdx)
!    WRITE(IPT,*) "MIN/MAX(SPDy)::",minval(spdy),maxval(spdy)

!    WRITE(IPT,*) "MIN/MAX(STRX)::",minval(strx),maxval(strx)
!    WRITE(IPT,*) "MIN/MAX(STRy)::",minval(stry),maxval(stry)

    END IF

    IF(HEATING_ON) THEN
       CALL UPDATE_HFX(NOW,SWR,NHF)
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: UPDATE_BINARY"
    
  END SUBROUTINE UPDATE_BINARY

  SUBROUTINE DUMP_BINARY(NOW)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_BINARY"

    IF(PRECIPITATION_ON) THEN
       CALL DUMP_EVP(NOW)
    END IF

    IF(WIND_ON) THEN
       CALL DUMP_WND(NOW)
    END IF

    IF(HEATING_ON) THEN
       CALL DUMP_HFX(NOW)
    END IF

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_BINARY"

  END SUBROUTINE DUMP_BINARY
  
  FUNCTION NEW_DATA(dims)
    IMPLICIT NONE
    TYPE(bin_data), POINTER :: NEW_DATA
    integer, intent(IN) :: DIMS
    INTEGER :: STATUS

    ALLOCATE(NEW_DATA,stat=status)
    IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE TYPE POINTER?")
    
    ALLOCATE(NEW_DATA%DATA(0:DIMS), STAT=STATUS)
    IF(status /=0) CALL FATAL_ERROR("NEW_DATA: COULD NOT ALLOCATE DATA POINTER?")
    
  END FUNCTION NEW_DATA


  SUBROUTINE LOAD_BINARY(WND,HFX,EVP)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WND,HFX,EVP
    INTEGER :: STATUS
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:LOAD_BINARY"
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Wind Stress file:"&
         &//TRIM(WND)
    inquire(file=trim(WND),exist=WIND_ON)
    IF(WIND_ON) THEN
       
       IF(WIND_TYPE /= 'stress' .and. WIND_TYPE/='speed')THEN
          CALL FATAL_ERROR("To convert a binary wind file,",&
               &"you must specify 'binary speed' or 'binary stress'")
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND WIND FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(WNDUNIT,WND,'cur')
       
       WNDX_NEXT => NEW_DATA(N)
       WNDX_PREV => NEW_DATA(N)
       WNDY_NEXT => NEW_DATA(N)
       WNDY_PREV => NEW_DATA(N)
       
       CALL READ_WND(WNDX=WNDX_PREV,WNDY=WNDY_PREV)
       CALL READ_WND(WNDX=WNDX_NEXT,WNDY=WNDY_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(WNDX_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(WNDX_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND WIND FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO WIND FILE FOUND"
       
    END IF
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for HeatFlux file:"&
         &//TRIM(HFX)
    inquire(file=trim(HFX),exist=HEATING_ON)
    IF(HEATING_ON) THEN
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HEATING FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(HFXUNIT,HFX,'cur')
       
       SWR_NEXT => NEW_DATA(M)
       SWR_PREV => NEW_DATA(M)
       NHF_NEXT => NEW_DATA(M)
       NHF_PREV => NEW_DATA(M)

       CALL READ_HFX(SWR=SWR_PREV,NHF=NHF_PREV)
       CALL READ_HFX(SWR=SWR_NEXT,NHF=NHF_NEXT)

       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SWR_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND HEATING FILE: READ FIRST DATA POINTS"
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO HEATING FILE FOUND"
       
       
    END IF
    
    IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "Looking for Evaporation file:"&
         &//TRIM(EVP)
    inquire(file=trim(EVP),exist=PRECIPITATION_ON)
    IF(PRECIPITATION_ON) THEN
 
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND PRECIPITATION FILE: OPEN AND READ"
       
       IF(MSR) CALL FOPEN(EVPUNIT,EVP,'cur')

       EVP_NEXT => NEW_DATA(M)
       EVP_PREV => NEW_DATA(M)
       PRC_NEXT => NEW_DATA(M)
       PRC_PREV => NEW_DATA(M)

       CALL READ_EVP(EVP=EVP_PREV,PRC=PRC_PREV)
       CALL READ_EVP(EVP=EVP_NEXT,PRC=PRC_NEXT)
       
       IF(DBG_SET(DBG_LOG)) THEN
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"FIRST TIME POINT",timezone)
          CALL PRINT_REAL_TIME(SWR_NEXT%dtm,IPT,"SECOND TIME POINT",timezone)
       END IF
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "FOUND PRECIPITATION FILE: READ FIRST DATA POINTS"
       
    ELSE
       
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "! NO PRECIPITATION FILE FOUND"
       
    END IF

    IF (.not. PRECIPITATION_ON .and. &
         & .not. HEATING_ON .and. &
         & .not. WIND_ON) CALL FATAL_ERROR("FOUND NO BINARY FORCING INPUT FILES?")
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:LOAD_BINARY"
    
  END SUBROUTINE LOAD_BINARY
!------------------------------------------------------------------
  SUBROUTINE IOERROR(IOS,MSG)
    IMPLICIT NONE
    INTEGER IOS
    CHARACTER(LEN=*) MSG
    CHARACTER(LEN=4) IOSC
    
    IF(IOS ==0) RETURN

    WRITE(IOSC,'(I4)') IOS
    
    CALL FATAL_ERROR("ERROR DURING FILE IO:"//TRIM(IOSC),&
         TRIM(MSG))
  
  END SUBROUTINE IOERROR

!-------------------------------------------------------------------
  SUBROUTINE READ_WND(WNDX,WNDY)
    IMPLICIT NONE
    TYPE(bin_data) :: WNDX,WNDY
    REAL(SP) :: hour
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: WNDYGL(:),WNDXGL(:)
   
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_WND"


    IF(MSR) THEN
 
       IF(PAR) THEN
          ALLOCATE(WNDYGL(NGL))
          ALLOCATE(WNDXGL(NGL))
       ELSE
          WNDXGL => WNDX%DATA(1:NGL)
          WNDYGL => WNDY%DATA(1:NGL)
       END IF


       READ(UNIT=wndunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from wind file")
       READ(UNIT= wndunit,IOSTAT=ios)(WNDXGL(i),WNDYGL(i),i=1,ngl)
       CALL IOERROR(IOS,"Can't read data from wind file")


       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(WNDX)::",MINVAL(WNDXgl),maxval(WNDXgl),sum(WNDXgl)/real(Ngl)
          WRITE(IPT,*) "MIN/MAX/MEAN(WNDY)::",MINVAL(WNDYgl),maxval(WNDYgl),sum(WNDYgl)/real(Ngl)
       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL PDEAL(MYID,MSRID,NPROCS,EMAP,WNDXGL,WNDX%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,EMAP,WNDYGL,WNDY%DATA)

       IF(MSR) THEN
          DEALLOCATE(WNDXGL)
          DEALLOCATE(WNDYGL)
       END IF
# endif
    END IF

    NULLIFY(WNDXGL,WNDYGL)

    WNDX%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    WNDY%dtm = WNDX%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_WND"
  END SUBROUTINE READ_WND

!-------------------------------------------------------------------
  SUBROUTINE READ_HFX(SWR,NHF)
    IMPLICIT NONE
    TYPE(bin_data) :: SWR,NHF
    REAL(SP) :: hour
    integer :: i, SOURCE, ios
    Real(SP), POINTER :: NHFGL(:),SWRGL(:)
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_HFX"
    

    IF(MSR) THEN
       IF(PAR) THEN
          ALLOCATE(NHFGL(MGL))
          ALLOCATE(SWRGL(MGL))
       ELSE
          SWRGL => SWR%DATA(1:MGL)
          NHFGL => NHF%DATA(1:MGL)
       END IF



       READ(UNIT=hfxunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from heatflux file")
       READ(UNIT=hfxunit,IOSTAT=ios)(NHFGL(i),SWRGL(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from heatflux file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(SWR)::",MINVAL(evpgl),maxval(evpgl),sum(evpgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(NHF)::",MINVAL(NHFgl),maxval(NHFgl),sum(NHFgl)/real(Mgl)
       END IF


    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1

       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,SWRGL,SWR%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,NHFGL,NHF%DATA)

       IF(MSR) THEN
          DEALLOCATE(SWRGL)
          DEALLOCATE(NHFGL)
       END IF
# endif
    END IF

    NULLIFY(SWRGL,NHFGL)

    SWR%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    NHF%dtm = SWR%dtm
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_HFX"

  END SUBROUTINE READ_HFX

!-------------------------------------------------------------------
  SUBROUTINE READ_EVP(EVP,PRC)
    IMPLICIT NONE
    TYPE(bin_data) :: EVP,PRC
    REAL(SP) :: hour
    integer :: i, SOURCE,ios
    Real(SP), POINTER :: EVPGL(:),PRCGL(:)
  
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START: READ_EVP"

    IF(MSR) THEN

         IF(PAR) THEN
          ALLOCATE(EVPGL(MGL))
          ALLOCATE(PRCGL(MGL))
       ELSE
          EVPGL => EVP%DATA(1:MGL)
          PRCGL => PRC%DATA(1:MGL)
       END IF

       READ(UNIT=evpunit,IOSTAT=ios) hour
       CALL IOERROR(IOS,"Can't read hour from evap/prec file")
       READ(UNIT=evpunit,IOSTAT=ios)(evpgl(i),prcgl(i),i=1,mgl)
       CALL IOERROR(IOS,"Can't read data from evap/prec file")

       IF(DBG_SET(DBG_IO)) THEN
          WRITE(IPT,*) "MIN/MAX/MEAN(EVAP)::",MINVAL(evpgl),maxval(evpgl),sum(evpgl)/real(Mgl)
          WRITE(IPT,*) "MIN/MAX/MEAN(PREC)::",MINVAL(prcgl),maxval(prcgl),sum(prcgl)/real(Mgl)
       END IF

    END IF


    ! IF NOT PAR, VARIABLES ARE ALREADY POINTED CORRECTLY
    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       
       SOURCE = MSRID -1
       CALL MPI_BCAST(hour,1,MPI_F,SOURCE,MPI_FVCOM_GROUP)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,prcgl,PRC%DATA)
       CALL PDEAL(MYID,MSRID,NPROCS,NMAP,evpgl,EVP%DATA)
       
       IF(MSR) THEN
          DEALLOCATE(PRCGL)
          DEALLOCATE(EVPGL)
       END IF

# endif
    END IF

    NULLIFY(EVPGL,PRCGL)
    
    PRC%dtm = seconds2time(hour*3600.0_SP) + ZEROTIME

    EVP%dtm = PRC%dtm

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END: READ_EVP"

  END SUBROUTINE READ_EVP

  SUBROUTINE UPDATE_WND(NOW,WNDX,WNDY)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: WNDX(:), WNDY(:)
    TYPE(BIN_DATA), POINTER :: A, B
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_WND "

    DO       
       IF(NOW .LT. WNDX_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(WNDX_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE WIND FILE")
          
       ELSE IF(NOW .gt. WNDX_NEXT%dtm)THEN
          
          A=> WNDX_PREV
          WNDX_PREV => WNDX_NEXT
          
          B => WNDY_PREV
          WNDY_PREV => WNDY_NEXT
          
          CALL READ_WND(WNDX=A,WNDY=B)

          WNDX_NEXT => A
          WNDY_NEXT => B
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO

    NUMER = SECONDS(NOW - WNDX_PREV%dtm)

    DENOM = SECONDS(WNDY_NEXT%dtm - WNDX_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    WNDX= WNDX_NEXT%data *fw + WNDX_PREV%data *bw 

    WNDY= WNDY_NEXT%data *fw + WNDY_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_WND "

  END SUBROUTINE UPDATE_WND

  SUBROUTINE UPDATE_HFX(NOW,SWR,NHF)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: SWR(:), NHF(:)
    TYPE(BIN_DATA), POINTER :: A, B
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_HFX "

    DO       
       IF(NOW .LT. SWR_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(SWR_PREV%dtm,IPT,"DATA TIME",timezone)

          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE HEATING FILE")
          
       ELSE IF(NOW .gt. SWR_NEXT%dtm)THEN
          
          A=> SWR_PREV
          SWR_PREV => SWR_NEXT
          
          B => NHF_PREV
          NHF_PREV => NHF_NEXT
          
          CALL READ_HFX(SWR=A,NHF=B)
          SWR_NEXT => A
          NHF_NEXT => B
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - SWR_PREV%dtm)

    DENOM = SECONDS(SWR_NEXT%dtm - SWR_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    SWR= SWR_NEXT%data *fw + SWR_PREV%data *bw 

    NHF= NHF_NEXT%data *fw + NHF_PREV%data *bw 


    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_HFX "

  END SUBROUTINE UPDATE_HFX

  SUBROUTINE UPDATE_EVP(NOW,EVP,PRC)
    IMPLICIT NONE
    TYPE(TIME) :: NOW
    REAL(SP), POINTER :: EVP(:), PRC(:)
    TYPE(BIN_DATA), POINTER :: A, B
    REAL(DP)     :: denom, numer
    REAL(SP)     :: fw, bw
    
    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:UPDATE_EVP "

    DO       
       IF(NOW .LT. EVP_PREV%dtm) THEN

          CALL PRINT_REAL_TIME(NOW,IPT,"OUTPUT TIME",timezone)
          CALL PRINT_REAL_TIME(EVP_PREV%dtm,IPT,"DATA TIME",timezone)


          CALL FATAL_ERROR("CAN NOT REWIND BINARY FILES",&
               & "SOMETHING IS WRONG WITH TIME IN THE EVAP/PREC FILE")
          
       ELSE IF(NOW .gt. EVP_NEXT%dtm)THEN
          
          A=> EVP_PREV
          EVP_PREV => EVP_NEXT
          
          B => PRC_PREV
          PRC_PREV => PRC_NEXT
          
          CALL READ_EVP(EVP=A,PRC=B)
          EVP_NEXT => A
          PRC_NEXT => B
          
       ELSE
          
          EXIT
          
       END IF
       
    END DO
       
    NUMER = SECONDS(NOW - EVP_PREV%dtm)

    DENOM = SECONDS(PRC_NEXT%dtm - EVP_PREV%dtm)
    
    fw = NUMER/DENOM
    bw = 1.0_DP - NUMER/DENOM

    EVP= EVP_NEXT%data *fw + EVP_PREV%data *bw 

    PRC= PRC_NEXT%data *fw + PRC_PREV%data *bw 

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:UPDATE_EVP "

  END SUBROUTINE UPDATE_EVP



!/////////////////////////////////////////////////////////////
  Subroutine MAKE_BINARY_OUTFILE(WNDNM,HFXNM,EVPNM)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WNDNM,HFXNM,EVPNM
    integer status
    INTEGER(ITIME) :: dummy
    CHARACTER(LEN=4) :: FLAG

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:MAKE_BINARY_OUTFILE "

    IF (.NOT. MSR) RETURN

    IF(WIND_ON) THEN
       CALL FOPEN(WNDUNIT,WNDNM,'our')


       IF(PAR) THEN
          ALLOCATE(WNDXGL(0:NGL))
          ALLOCATE(WNDYGL(0:NGL))
       ELSE
          WNDXGL => STRX
          WNDYGL => STRY
       END IF
    END IF

    IF (HEATING_ON) THEN
       CALL FOPEN(HFXUNIT,HFXNM,'our')

       IF(PAR) THEN
          ALLOCATE(SWRGL(0:MGL))
          ALLOCATE(NHFGL(0:MGL))
       ELSE
          swrgl  => SWR
          nhfgl  => NHF
       END IF
    END IF

    IF(PRECIPITATION_ON) THEN
       CALL FOPEN(EVPUNIT,EVPNM,'our')

       IF(PAR) THEN
          ALLOCATE(evpgl(0:MGL))
          ALLOCATE(prcgl(0:MGL))
       ELSE
          evpgl => EVP
          prcgl => PRC
       END IF
    END IF


    ! CHECK TO MAKE SURE ATLEAST ONE IS ACTIVE
    IF (.not. PRECIPITATION_ON .and. &
         & .not. HEATING_ON .and. &
         & .not. WIND_ON) CALL FATAL_ERROR("FOUND NO BINARY FORCING INPUT FILES?")

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:MAKE_BINARY_OUTFILE "

  END SUBROUTINE MAKE_BINARY_OUTFILE
! -------------------------------------------------------------------  
  SUBROUTINE DUMP_HFX(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_HFX "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,SWR,SWRGL)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,NHF,NHFGL)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(swrgl(1:mGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("swrgl is nan")
    write(ipt,*) "SWRGL: min/max/mean",minval(swrgl(1:MGL)),maxval(swrgl(1:MGL)),mean

    mean = sum(nhfgl(1:mGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("nhfgl is nan")
    write(ipt,*) "NHFGL: min/max/mean",minval(nhfgl(1:MGL)),maxval(nhfgl(1:MGL)),mean

    WRITE(hfxunit) hour
    write(hfxunit)(NHFGL(i),SWRGL(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_HFX "

  END SUBROUTINE DUMP_HFX


  SUBROUTINE DUMP_EVP(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_EVP "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,PRC,prcgl)
       CALL PCOLLECT(MYID,MSRID,NPROCS,NMAP,EVP,evpgl)
# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(evpgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Evpgl is nan")
    write(ipt,*) "EVPGL: min/max/mean",minval(Evpgl(1:MGL)),maxval(Evpgl(1:MGL)),mean

    mean = sum(Prcgl(1:mGL))/real(MGL,SP)
    if(isnan(mean)) call fatal_error("Prcgl is nan")
    write(ipt,*) "PRCGL: min/max/mean",minval(Prcgl(1:MGL)),maxval(Prcgl(1:MGL)),mean

    WRITE(evpunit) hour
    write(evpunit)(evpgl(i),prcgl(i),i=1,mgl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_EVP "

  END SUBROUTINE DUMP_EVP


  SUBROUTINE DUMP_WND(NOW)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP) :: hour, mean
    INTEGER :: STATUS, I

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "START:DUMP_WND "

    IF (PAR)THEN
# if defined(MULTIPROCESSOR)

       CALL PCOLLECT(MYID,MSRID,NPROCS,EMAP,STRX,WNDXGL)
       CALL PCOLLECT(MYID,MSRID,NPROCS,EMAP,STRY,WNDYGL)

# endif
    END IF

    IF (.NOT. MSR) RETURN

    hour = DAYS(NOW - ZEROTIME) * 24.0_SP

    write(ipt,*) "time", hour

    mean = sum(wndxgl(1:NGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("wndxgl is nan")
    write(ipt,*) "WNDXGL: min/max/mean",minval(wndxgl(1:NGL)),maxval(wndxgl(1:NGL)),mean


    mean = sum(wndygl(1:NGL))/real(NGL,SP)
    if(isnan(mean)) call fatal_error("wndygl is nan")
    write(ipt,*) "WNDYGL: min/max/mean",minval(wndygl(1:NGL)),maxval(wndygl(1:NGL)),mean

    WRITE(wndunit) hour
    write(wndunit)(WNDXGL(i),WNDYGL(i),i=1,ngl)

    IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "END:DUMP_WND "

  END SUBROUTINE DUMP_WND



END MODULE MOD_BINARY
