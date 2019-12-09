










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
! PURPOSE ARE  DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/
SUBROUTINE ADCOR

   USE ALL_VARS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE MOD_WD


   IMPLICIT NONE
   REAL(SP) :: UFC(0:NT,KB),VFC(0:NT,KB)
   REAL(SP),PARAMETER :: BETA0=0.5_SP
   REAL(SP) ::CURCOR,PRECOR
   INTEGER :: I,K  
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP  


   UFC=0.0_SP
   VFC=0.0_SP

   DO I = 1, N

       DO K = 1, KBM1
         CURCOR=BETA0*COR(I)*VF(I,K)
         PRECOR=(1._SP-BETA0)*COR(I)*V(I,K)
         UFC(I,K)=UBETA(I,K)-(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       END DO
   
   END DO

   DO I = 1, N

       DO K = 1, KBM1
         CURCOR=BETA0*COR(I)*UF(I,K)
         PRECOR=(1._SP-BETA0)*COR(I)*U(I,K)
         VFC(I,K)=VBETA(I,K)+(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       END DO

   END DO

   DO I=1,N
       DO K=1,KBM1
         UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*UFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*VFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       END DO
   END DO

   DO I =1,N
      IF(ISWETCT(I)*ISWETC(I) .NE. 1)THEN
         DO K=1,KBM1
            UF(I,K)=0.0_SP
            VF(I,K)=0.0_SP
         END DO
      END IF
   END DO

   RETURN
   END SUBROUTINE ADCOR
