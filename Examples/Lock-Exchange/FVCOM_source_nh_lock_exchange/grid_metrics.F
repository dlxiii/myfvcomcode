! $Id: grid_metrics.F,v 1.1.2.9 2008/04/03 16:30:46 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.9 $
SUBROUTINE GRID_METRICS
  USE MOD_PAR
  USE MOD_OBCS
  USE ALL_VARS
  USE MOD_NORTHPOLE
  IMPLICIT NONE
  

  ! THESE SUBROUTINES SHOULD BE FURTHER BROKEN DOWN INTO THEIR
  ! COMPONENT OPERATIONS AND THE WHOLE THING SHOULD BE MODULAR
  
  
!============================================
!Set up fluxes and control Volumes
!============================================
  CALL TRIANGLE_GRID_EDGE
    
!============================================
!Calculate Element and Control Volume Areas
!============================================
  CALL CELL_AREA            


!====================================================
! Calculate Shape Coefficients for Flux Construction
!====================================================
# if defined (GCN)
  CALL SHAPE_COEF_GCN       
# else
  CALL SHAPE_COEF_GCY  
# endif
  
#  if defined (SPHERICAL)
  CALL FIND_NORTHPOLE
  CALL SHAPE_COEF_XY 
  CALL FIND_CELLSIDE
#  endif

!============================================
!Calculate Node Control Volume Edge Lengths
!============================================
  CALL EDGE_LEN


!====================================================
!  EXCHANGE SHAPE FACTOR INFORMATION
!====================================================
# if defined (MULTIPROCESSOR)
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,A1U,A2U) 
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,AWX,AWY,AW0) 
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ALPHA) 
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ART)
# endif

!====================================================
! Calculate the horizontal gradient of the bathymetry
!====================================================
 CALL DEPTH_GRADIENT


!====================================================
! SETUP BOUNDARY NEIGHBORS AND SUCH
!====================================================
  CALL SETUP_OBC

!====================================================
! SETUP THE LONGSHORE FLOW INDEX AND ANGLES
!====================================================
  IF(OBC_LONGSHORE_FLOW_ON) CALL longshore_flow




END SUBROUTINE GRID_METRICS

