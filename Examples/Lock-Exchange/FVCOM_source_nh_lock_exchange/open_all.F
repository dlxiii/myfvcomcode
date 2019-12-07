! $Id: open_all.F,v 1.1.2.6 2008/04/03 17:05:10 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.6 $
SUBROUTINE OPEN_ALL
!===============================================================================!
! OPEN FILES
! UNLESS OTHERWISE SPECIFED THE ROUTINES CALLED HERE ARE IN mod_input.F
!===============================================================================!
  USE CONTROL
  USE MOD_INPUT
  IMPLICIT NONE
  CHARACTER(LEN=160) :: FNAME
  LOGICAL :: FEXIST
  

  CALL NULLIFY_FILE_POINTERS

  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE(STARTUP_TYPE_HOTSTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               OPEN INIT FILES FOR HOTSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     CALL CHECK_IO_DIRS
     
     CALL OPEN_STARTUP_FILE
     
     CALL OPEN_FORCING

     CALL OPEN_NEW_OUTPUT
          
          
     !=================================================
     ! CRASHSTART 
  CASE(STARTUP_TYPE_CRASHRESTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              OPENING FILES FOR CRASHRESTART                    !'
        WRITE(IPT,*)'!                                                                !'
     end if
 
     CALL CHECK_IO_DIRS

     CALL OPEN_CRASHSTART
     
     CALL OPEN_FORCING
     
     !=================================================
     ! COLDSTART 
  CASE(STARTUP_TYPE_COLDSTART)
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              OPENING FILES FOR COLDSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
          

     CALL CHECK_IO_DIRS

        ! MAKE SURE THE RUN FILE DOES NOT REQUEST A START FILE
        

     IF (STARTUP_TS_TYPE .eq. STARTUP_TYPE_OBSERVED) THEN
        CALL OPEN_STARTUP_FILE

     ELSE IF (STARTUP_TS_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
        CALL OPEN_STARTUP_FILE
        
     ELSE IF (STARTUP_UV_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
        CALL OPEN_STARTUP_FILE
        
     ELSE IF (STARTUP_TURB_TYPE .eq. STARTUP_TYPE_SETVALUES) THEN
        CALL OPEN_STARTUP_FILE
     ELSE 
        if(dbg_set(dbg_log)) write(ipt,*) "! No Startup file needed fo&
             &r this cold start"

     END IF
    
     ! OPEN THE OTHER COLD START FILES (GRID,DEPTH SPONGE, ETC)
     IF (MSR) CALL OPEN_COLDSTART ! ONLY MASTER READS THESE FILES
     
     CALL OPEN_FORCING

     CALL OPEN_NEW_OUTPUT 

  END SELECT

END SUBROUTINE OPEN_ALL
