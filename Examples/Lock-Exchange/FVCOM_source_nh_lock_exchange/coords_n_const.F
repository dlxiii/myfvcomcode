! $Id: coords_n_const.F,v 1.1.2.7 2008/04/11 14:17:01 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.7 $
SUBROUTINE COORDS_N_CONST

  USE MOD_INPUT
  USE MOD_SETUP
  use non_hydro
  IMPLICIT NONE
  INTEGER STATUS, I

!===============================================================================!
! OPEN FILES
! UNLESS OTHERWISE SPECIFED THE ROUTINES CALLED HERE ARE IN mod_input.F
!===============================================================================!
    
  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE("hotstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING FILES FOR HOTSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
               
     ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
     ALLOCATE(Y_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
     Y_LCL=0.0_SP
     ALLOCATE(X_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")
     X_LCL=0.0_SP
     ALLOCATE(H_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE H_LCL")
     H_LCL=0.0_SP

     CALL LOAD_RESTART_COORDS(X_LCL,Y_LCL)
     CALL COORDINATE_UNITS(X_LCL,Y_LCL)
     CALL SETUP_CENTER_COORDS
     
     DEALLOCATE(X_LCL)
     DEALLOCATE(Y_LCL)


     CALL LOAD_RESTART_DEPTH(H_LCL)
     CALL SETUP_DEPTH
     DEALLOCATE(H_LCL) ! COULD BE LOADED DIRECTLY - MUST SET MAX/MIN

     CALL LOAD_RESTART_CORIOLIS(COR)! LOAD DIRECTLY

     CALL SETUP_GRAVITY
     
     CALL LOAD_RESTART_SPONGE(CC_SPONGE) ! LOAD DIRECTLY
     
     STYPE = STYPE_RESTART
     CALL LOAD_RESTART_SIGMA(Z,Z1) ! LOAD DIRECTLY TO ALL_VARS:Z,Z1
     CALL SETUP_SIGMA_DERIVATIVES

     CALL SETUP_OBCTYPES

     CALL SETUP_HORIZONTAL_MIXING_COEFFICIENT

     CALL SETUP_BOTTOM_ROUGHNESS

     !=================================================
     ! CRASHSTART 
  CASE("crashrestart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              READING FILES FOR CRASHRESTART                    !'
        WRITE(IPT,*)'!                                                                !'
     end if
     ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
     ALLOCATE(Y_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
     Y_LCL=0.0_SP
     ALLOCATE(X_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")
     X_LCL=0.0_SP
     ALLOCATE(H_LCL(0:MT),stat=status)
     IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE H_LCL")
     H_LCL=0.0_SP

          
     CALL LOAD_RESTART_COORDS(X_LCL,Y_LCL)
     CALL COORDINATE_UNITS(X_LCL,Y_LCL)
     CALL SETUP_CENTER_COORDS

     DEALLOCATE(X_LCL)
     DEALLOCATE(Y_LCL)

     
     CALL LOAD_RESTART_DEPTH(H_LCL)
     CALL SETUP_DEPTH
     DEALLOCATE(H_LCL) ! COULD BE LOADED DIRECTLY - MUST SET MAX/MIN

     CALL LOAD_RESTART_CORIOLIS(COR)

     CALL SETUP_GRAVITY
     
     CALL LOAD_RESTART_SPONGE(CC_SPONGE) ! LOAD DIRECTLY TO ALL_VARS:CC_SPONGE


     STYPE = STYPE_RESTART
     CALL LOAD_RESTART_SIGMA(Z,Z1) ! LOAD DIRECTLY TO ALL_VARS:Z,Z1
     CALL SETUP_SIGMA_DERIVATIVES

     CALL SETUP_OBCTYPES

     CALL SETUP_HORIZONTAL_MIXING_COEFFICIENT

     CALL SETUP_BOTTOM_ROUGHNESS

     !=================================================
     ! COLDSTART 
  CASE("coldstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!              READING FILES FOR COLDSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     IF(MSR) THEN
        ! ALLOCATE SPACE FOR THE GLOBAL GRID DATA
        ALLOCATE(Y_GBL(0:MGL),stat=status)
        IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_GBL")
        Y_GBL=0.0_SP
        ALLOCATE(X_GBL(0:MGL),stat=status)
        IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_GBL")
        X_GBL=0.0_SP
     END IF

    ! ALLOCATE SPACE FOR THE LOCAL GRID DATA
    ALLOCATE(Y_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE Y_LCL")
    Y_LCL=0.0_SP
    ALLOCATE(X_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE X_LCL")
    X_LCL=0.0_SP

    ALLOCATE(H_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE H_LCL")
    H_LCL=0.0_SP

    ALLOCATE(C_LCL(0:MT),stat=status)
    IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE C_LCL")
    C_LCL=0.0_SP

     
     CALL LOAD_COLDSTART_COORDS(X_GBL,Y_GBL,X_LCL,Y_LCL)
     CALL COORDINATE_UNITS(X_LCL,Y_LCL)
     CALL SETUP_CENTER_COORDS


     DEALLOCATE(X_LCL)
     DEALLOCATE(Y_LCL)

     CALL LOAD_COLDSTART_DEPTH(X_GBL,Y_GBL,H_LCL)
     CALL SETUP_DEPTH

     DEALLOCATE(H_LCL)

     CALL LOAD_COLDSTART_CORIOLIS(X_GBL,Y_GBL,C_LCL)
     CALL SETUP_CORIOLIS

     DEALLOCATE(C_LCL)

     CALL SETUP_GRAVITY

     ! THESE VARIABLES ARE NOT ALLOCATED YET
     CALL LOAD_COLDSTART_SPONGE(X_GBL,Y_GBL,NSPONGE,N_SPG,R_SPG,C_SPG,X_SPG,Y_SPG)
     CALL SETUP_SPONGE

 
!     CALL LOAD_COLDSTART_SIGMA
     CALL SETUP_SIGMA
     CALL SETUP_SIGMA_DERIVATIVES

     IF(MSR) THEN
        XG = X_GBL
        DEALLOCATE(X_GBL)
        DEALLOCATE(Y_GBL)
     END IF

     CALL SETUP_OBCTYPES

     CALL SETUP_HORIZONTAL_MIXING_COEFFICIENT

     CALL SETUP_BOTTOM_ROUGHNESS

  END SELECT
  
  
END SUBROUTINE COORDS_N_CONST
