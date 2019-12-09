










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

!==============================================================================|
!   CALCULATE GHOST VELOCITY FOR EXTERNAL MODE                                 !
!==============================================================================|
   SUBROUTINE GHOSTUV2(I,JJ,UAKK,VAKK)
!==============================================================================|

   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: I,JJ
   INTEGER              :: J1,J2
   REAL(SP)             :: DELTX,DELTY,ALPHA1
   REAL(SP)             :: UTMP,VTMP
   REAL(SP), INTENT(OUT):: UAKK,VAKK

   UAKK = 0.0_SP; VAKK = 0.0_SP

   IF(ISBCE(I) /= 2)THEN
     J1 = JJ+1-INT((JJ+1)/4)*3
     J2 = JJ+2-INT((JJ+2)/4)*3
     DELTX = VX(NV(I,J1))-VX(NV(I,J2))
     DELTY = VY(NV(I,J1))-VY(NV(I,J2))

     ALPHA1 = ATAN2(DELTY,DELTX)

     UTMP = UA(I)*COS(ALPHA1)+VA(I)*SIN(ALPHA1)
     VTMP = -UA(I)*SIN(ALPHA1)+VA(I)*COS(ALPHA1)

!     VTMP = -VTMP
     VTMP = 0.0_SP

     UAKK = UTMP*COS(ALPHA1)-VTMP*SIN(ALPHA1)
     VAKK = UTMP*SIN(ALPHA1)+VTMP*COS(ALPHA1)
   ELSE IF(ISBCE(I) == 2)THEN
     UAKK = UA(I)
     VAKK = VA(I)
   END IF

   RETURN
   END SUBROUTINE GHOSTUV2


!==============================================================================|
!   CALCULATE GHOST VELOCITY FOR INTERNAL MODE                                 !
!==============================================================================|
   SUBROUTINE GHOSTUV3(I,JJ,UAKK,VAKK)
!==============================================================================|

   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: I,JJ
   INTEGER              :: J1,J2,K
   REAL(SP)             :: DELTX,DELTY,ALPHA1
   REAL(SP)             :: UTMP,VTMP
   REAL(SP), INTENT(OUT):: UAKK(KB),VAKK(KB)

   UAKK = 0.0_SP; VAKK = 0.0_SP

   IF(ISBCE(I) /= 2)THEN
     J1 = JJ+1-INT((JJ+1)/4)*3
     J2 = JJ+2-INT((JJ+2)/4)*3
     DELTX = VX(NV(I,J1))-VX(NV(I,J2))
     DELTY = VY(NV(I,J1))-VY(NV(I,J2))

     ALPHA1 = ATAN2(DELTY,DELTX)

     DO K = 1,KBM1
       UTMP = U(I,K)*COS(ALPHA1)+V(I,K)*SIN(ALPHA1)
       VTMP = -U(I,K)*SIN(ALPHA1)+V(I,K)*COS(ALPHA1)

!       VTMP = -VTMP
       VTMP = 0.0_SP

       UAKK(K) = UTMP*COS(ALPHA1)-VTMP*SIN(ALPHA1)
       VAKK(K) = UTMP*SIN(ALPHA1)+VTMP*COS(ALPHA1)
     END DO
     
   ELSE IF(ISBCE(I) == 2)THEN
     UAKK = U(I,:)
     VAKK = V(I,:)
   END IF

   RETURN
   END SUBROUTINE GHOSTUV3
