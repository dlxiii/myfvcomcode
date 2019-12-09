










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

MODULE RRKVAL
END MODULE RRKVAL

MODULE MOD_RRK
   USE MOD_INPUT, only : NC_START
   USE RRKVAL
   USE CONTROL
   USE MOD_UTILS
   USE MOD_NCTOOLS
   IMPLICIT NONE
   SAVE
   
   INTEGER  ::  RRK_EOFCONTR
   CHARACTER(LEN=80) :: REF_START_DATE
   CHARACTER(LEN=80) :: REF_END_DATE
   CHARACTER(LEN=80) :: RRK_START_DATE
   CHARACTER(LEN=80) :: RRK_END_DATE

   TYPE(TIME) :: REF_START_TIME
   TYPE(TIME) :: REF_END_TIME
   TYPE(TIME) :: RRK_START_TIME
   TYPE(TIME) :: RRK_END_TIME

   TYPE(TIME) :: RRK_CYC

   CHARACTER(LEN=80) :: RRK_ASSIM_INTERVAL
   TYPE(TIME) :: RRK_INTERVAL


   INTEGER      REF_INT           !!GLOBAL NUMBER OF THE READING FILE INTERVALS 
   INTEGER      RRK_NOBSMAX
   INTEGER      RRK_OPTION        !!OPTION 1 FOR BAROTROPIC CASE; OPTION 2 FOR BAROCLINIC CASE
   INTEGER      RRK_NEOF          !!NUMBER OF THE EOF  
   REAL(SP) ::  RRK_PSIZE         !!PERTURBATION SIZE  
   REAL(SP) ::  RRK_PSCALE        !!PSEUDO MODEL ERROR    
   REAL(SP) ::  RRK_RSCALE        !!SCALE FACTOR APPLIED TO ONE STANDARD DEVIATION FOR R

   LOGICAL  :: RRK_ON

END MODULE MOD_RRK

