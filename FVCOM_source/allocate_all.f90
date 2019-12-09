










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

SUBROUTINE ALLOCATE_ALL
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_OBCS
  USE MOD_WD
  USE MOD_SPHERICAL
  USE MOD_NESTING
  IMPLICIT NONE
  INTEGER IERR
  
  IF(DBG_SET(DBG_LOG)) THEN
     WRITE(IPT,*  )'!'
     WRITE(IPT,*)'!                ALLOCATING MEMORY    '
     WRITE(IPT,*  )'!'
  END IF
  
  memcnt=0.0_SP

  CALL ALLOC_VARS(DBG_SET(DBG_LOG))
  
  CALL ALLOC_OBC_DATA

  ! ALWAYS INITIALIZE WET DRY VARIABLE IF DEFINED
  CALL ALLOC_WD_DATA





   MEMTOT = MEMCNT*4
   IF(PAR)CALL MPI_REDUCE(MEMCNT,MEMTOT,1,MPI_F,MPI_SUM,0,MPI_FVCOM_GROUP,IERR)

   IF(DBG_SET(DBG_LOG)) THEN
      WRITE(IPT,*)'!  # MBYTES TOTAL MEM    :',MEMTOT/1E+6
      IF(PAR)WRITE(IPT,*)'!  # AVERAGE MBYTES/PROC :',MEMTOT/(1E+6*NPROCS)
   END IF


END SUBROUTINE ALLOCATE_ALL
  
