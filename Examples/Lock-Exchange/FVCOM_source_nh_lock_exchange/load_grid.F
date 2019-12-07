! $Id: load_grid.F,v 1.1.2.6 2008/06/05 14:53:09 jqi Exp $
! $Name: New_Input $
! $Revision: 1.1.2.6 $
SUBROUTINE LOAD_GRID
!===============================================================================!
! GET THE GRID CONNECTIVITY AND OBC NODE LIST FOR DOMAIN
!    DECOMPOSITION AND GENMAP
!===============================================================================!
  USE CONTROL
  USE MOD_INPUT
  USE MOD_OBCS
  USE ALL_VARS
  USE MOD_NESTING
  IMPLICIT NONE
  integer i


  SELECT CASE(STARTUP_TYPE)
     !=================================================
     ! HOTSTART 
  CASE("hotstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING GRID FOR HOTSTART                        !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     CALL LOAD_RESTART_GRID(NVG) ! SET DIMENSIONS HERE
          
     CALL LOAD_RESTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
        
     CALL LOAD_RESTART_LSF_GRID(NOBCLSF_GL,IBCLSF_GL, RBC_GEO_GL,RBC_WDF_GL)
     

     !=================================================
     ! CRASHSTART 
  CASE("crashrestart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING GRID FOR CRASHSTART                      !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     CALL LOAD_RESTART_GRID(NVG)  ! SET DIMENSIONS HERE
 
     CALL LOAD_RESTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)
     
     CALL LOAD_RESTART_LSF_GRID(NOBCLSF_GL,IBCLSF_GL, RBC_GEO_GL,RBC_WDF_GL)  
     
     !=================================================
     ! COLDSTART 
  CASE("coldstart")
     !=================================================
     if(DBG_SET(dbg_log)) then 
        WRITE(IPT,*)'!               READING GRID FOR COLDSTART                       !'
        WRITE(IPT,*)'!                                                                !'
     end if
     
     ! SET DIMENSIONS HERE
     CALL LOAD_COLDSTART_SIGMA 

     CALL LOAD_COLDSTART_GRID(NVG)
     
     CALL LOAD_COLDSTART_OBC_GRID(IOBCN_GL,I_OBC_GL, TYPE_OBC_GL)

     CALL LOAD_COLDSTART_LSF(NOBCLSF_GL,IBCLSF_GL, RBC_GEO_GL,RBC_WDF_GL)

     IF(NESTING_ON) &
          & CALL LOAD_NESTING_OPB (NESTING_MODE, MGL_NEST,MID_NEST_GL,NGL_NEST,NID_NEST_GL)     
  END SELECT
  

END SUBROUTINE LOAD_GRID
