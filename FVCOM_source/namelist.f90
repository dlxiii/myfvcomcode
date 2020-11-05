










!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

SUBROUTINE NAMELIST
  USE MOD_UTILS
  USE CONTROL
  USE MOD_INPUT
  USE MOD_NESTING
  USE MOD_STATION_TIMESERIES  
  USE MOD_SPARSE_TIMESERIES
  USE MOD_DYE
  USE MOD_HEATFLUX

  IMPLICIT NONE


  !==============================================================================!
  ! SET DEFAULT VALUES IN NAME LIST                                                   
  !==============================================================================!
  CALL NAME_LIST_INITIALIZE

  CALL NAME_LIST_INITIALIZE_NEST

  CALL NAME_LIST_INITIALIZE_DYE
  CALL HEATING_CALCULATE_NAMELIST_INITIALIZE

  CALL STATION_NAME_LIST_INITIALIZE 

  ! IF FVCOM IS ONLY PRINTING A BLANK NAME LIST FOR A NEW CASE:
  if (BLANK_NAMELIST) then
     CALL NAME_LIST_PRINT

     ! NESTING ONLY WORKS IN PARALLEL
     CALL NAME_LIST_PRINT_NEST

     CALL NAME_LIST_PRINT_DYE
     CALL HEATING_CALCULATE_NAMELIST_PRINT

     CALL STATION_NAME_LIST_PRINT 

     CALL PSHUTDOWN
  end if

  !==============================================================================!
  !   SETUP MODEL RUN PARAMETERS                                                 !
  !==============================================================================!

  !READ DATA IN THE NAME LIST FILE
  CALL NAME_LIST_READ ! ALL PROCS READ THIS

  CALL NAME_LIST_READ_NEST
  IF(NESTING_ON .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN NESTING. STOP RUNNING...'
    CALL PSTOP
  END IF
  IF(NCNEST_ON .AND. SERIAL)THEN
    IF(MSR) WRITE(*,*) 'PLEASE USE MORE THAN ONE PROCESSOR TO RUN NCNEST. STOP RUNNING...'
    CALL PSTOP
  END IF


  CALL NAME_LIST_READ_DYE
  CALL HEATING_CALCULATE_NAMELIST_READ

  CALL STATION_NAME_LIST_READ

  !PRINT THE NAME LIST DATA TO THE SCREEN FOR THE LOG
  IF(DBG_SET(DBG_LOG)) CALL NAME_LIST_PRINT 


END SUBROUTINE NAMELIST
