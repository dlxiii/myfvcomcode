










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
!   CALCULATE CONVECTION AND DIFFUSION FLUXES FOR EXTERNAL MODE                !
!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_GCN(XFLUX,YFLUX)
!==============================================================================|

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE BCS
   USE MOD_OBCS
   USE MOD_WD





   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2,II
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN_TMP
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: FACT,FM1,ISWETTMP
   INTEGER  :: STAT_MF    
   REAL(SP) :: TPA,TPB


   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY

!#  if defined (THIN_DAM)
   REAL(SP) :: A1UIA1,A1UIA2,A1UIA3,A1UIA4,A2UIA1,A2UIA2,A2UIA3,A2UIA4
   REAL(SP) :: A1UIB1,A1UIB2,A1UIB3,A1UIB4,A2UIB1,A2UIB2,A2UIB3,A2UIB4
   INTEGER  :: J11,J12,J21,J22,E1,E2,ISBCE1,ISBC_TMP,IB_TMP
   LOGICAL  :: ISMATCH
!#  endif

   REAL(SP) :: BTPS
   REAL(SP) :: U_TMP,V_TMP,UAC_TMP,VAC_TMP,WUSURF_TMP,WVSURF_TMP,WUBOT_TMP,WVBOT_TMP,UAF_TMP,VAF_TMP 

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advave_edge_gcn.F"

!------------------------------------------------------------------------------!

   SELECT CASE(HORIZONTAL_MIXING_TYPE)
   CASE ('closure')
      FACT = 1.0_SP
      FM1  = 0.0_SP
   CASE('constant')
      FACT = 0.0_SP
      FM1  = 1.0_SP
   CASE DEFAULT
      CALL FATAL_ERROR("UNKNOW HORIZONTAL MIXING TYPE:",&
           & TRIM(HORIZONTAL_MIXING_TYPE) )
   END SELECT

!
!-------------------------INITIALIZE FLUXES------------------------------------!
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP
   PSTX  = 0.0_SP
   PSTY  = 0.0_SP


!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
   
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)

     J1=IENODE(I,1)
     J2=IENODE(I,2)

     DIJ=0.5_SP*(D(J1)+D(J2))
     ELIJ=0.5_SP*(EL(J1)+EL(J2))




     IF(ISWETCE(IA)*ISWETC(IA) == 1 .OR. ISWETCE(IB)*ISWETC(IB) == 1)THEN
!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)

     A1UIA1 = A1U(IA,1)
     A1UIA2 = A1U(IA,2)
     A1UIA3 = A1U(IA,3)
     A1UIA4 = A1U(IA,4)
     A2UIA1 = A2U(IA,1)
     A2UIA2 = A2U(IA,2)
     A2UIA3 = A2U(IA,3)
     A2UIA4 = A2U(IA,4)
!---------------------------------------------------------------
     COFA1=A1UIA1*UA(IA)+A1UIA2*UA(K1)+A1UIA3*UA(K2)+A1UIA4*UA(K3)
     COFA2=A2UIA1*UA(IA)+A2UIA2*UA(K1)+A2UIA3*UA(K2)+A2UIA4*UA(K3)
     COFA5=A1UIA1*VA(IA)+A1UIA2*VA(K1)+A1UIA3*VA(K2)+A1UIA4*VA(K3)
     COFA6=A2UIA1*VA(IA)+A2UIA2*VA(K1)+A2UIA3*VA(K2)+A2UIA4*VA(K3)

!     COFA1=A1U(IA,1)*UA(IA)+A1U(IA,2)*UA(K1)+A1U(IA,3)*UA(K2)+A1U(IA,4)*UA(K3)
!     COFA2=A2U(IA,1)*UA(IA)+A2U(IA,2)*UA(K1)+A2U(IA,3)*UA(K2)+A2U(IA,4)*UA(K3)
!     COFA5=A1U(IA,1)*VA(IA)+A1U(IA,2)*VA(K1)+A1U(IA,3)*VA(K2)+A1U(IA,4)*VA(K3)
!     COFA6=A2U(IA,1)*VA(IA)+A2U(IA,2)*VA(K1)+A2U(IA,3)*VA(K2)+A2U(IA,4)*VA(K3)
     
     XIJ=XIJC(I)-XC(IA)
     YIJ=YIJC(I)-YC(IA)
     UIJ1=UA(IA)+COFA1*XIJ+COFA2*YIJ
     VIJ1=VA(IA)+COFA5*XIJ+COFA6*YIJ

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
     IB_TMP = IB

     A1UIB1 = A1U(IB_TMP,1)
     A1UIB2 = A1U(IB_TMP,2)
     A1UIB3 = A1U(IB_TMP,3)
     A1UIB4 = A1U(IB_TMP,4)
     A2UIB1 = A2U(IB_TMP,1)
     A2UIB2 = A2U(IB_TMP,2)
     A2UIB3 = A2U(IB_TMP,3)
     A2UIB4 = A2U(IB_TMP,4)

     COFA3=A1UIB1*UA(IB_TMP)+A1UIB2*UA(K1)+A1UIB3*UA(K2)+A1UIB4*UA(K3)
     COFA4=A2UIB1*UA(IB_TMP)+A2UIB2*UA(K1)+A2UIB3*UA(K2)+A2UIB4*UA(K3)
     COFA7=A1UIB1*VA(IB_TMP)+A1UIB2*VA(K1)+A1UIB3*VA(K2)+A1UIB4*VA(K3)
     COFA8=A2UIB1*VA(IB_TMP)+A2UIB2*VA(K1)+A2UIB3*VA(K2)+A2UIB4*VA(K3)

!     COFA3=A1U(IB,1)*UA(IB)+A1U(IB,2)*UA(K1)+A1U(IB,3)*UA(K2)+A1U(IB,4)*UA(K3)
!     COFA4=A2U(IB,1)*UA(IB)+A2U(IB,2)*UA(K1)+A2U(IB,3)*UA(K2)+A2U(IB,4)*UA(K3)
!     COFA7=A1U(IB,1)*VA(IB)+A1U(IB,2)*VA(K1)+A1U(IB,3)*VA(K2)+A1U(IB,4)*VA(K3)
!     COFA8=A2U(IB,1)*VA(IB)+A2U(IB,2)*VA(K1)+A2U(IB,3)*VA(K2)+A2U(IB,4)*VA(K3)
     
     XIJ=XIJC(I)-XC(IB_TMP)
     YIJ=YIJC(I)-YC(IB_TMP)
     UIJ2=UA(IB_TMP)+COFA3*XIJ+COFA4*YIJ
     VIJ2=VA(IB_TMP)+COFA7*XIJ+COFA8*YIJ

!    NORMAL VELOCITY
     UIJ=0.5_SP*(UIJ1+UIJ2)
     VIJ=0.5_SP*(VIJ1+VIJ2)
     UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB_TMP)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!     VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)
     ! David moved HPRNU and added HVC
     VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB_TMP)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB_TMP)))/HPRNU

!    SHEAR STRESSES
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
     FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))

!    ADD CONVECTIVE AND VISCOUS FLUXES
     XADV=DIJ*UN_TMP*&
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
     YADV=DIJ*UN_TMP* &
          ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP

!    ACCUMULATE FLUX
     ISBC_TMP = ISBC(I)
     XFLUX(IA)=XFLUX(IA)+(XADV+FXX*EPOR(IA))*(1.0_SP-ISBC_TMP)*IUCP(IA)
     YFLUX(IA)=YFLUX(IA)+(YADV+FYY*EPOR(IA))*(1.0_SP-ISBC_TMP)*IUCP(IA)
     XFLUX(IB)=XFLUX(IB)-(XADV+FXX*EPOR(IB))*(1.0_SP-ISBC_TMP)*IUCP(IB)
     YFLUX(IB)=YFLUX(IB)-(YADV+FYY*EPOR(IB))*(1.0_SP-ISBC_TMP)*IUCP(IB)

     END IF



!    ACCUMULATE BAROTROPIC FLUX
!for spherical coordinator and domain across 360^o latitude         
     PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC(I)
     PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC(I)
     PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC(I)
     PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC(I)

   END DO

!#  if !defined (SEMI_IMPLICIT)

   DO I = 1,N
     ISWETTMP = ISWETCE(I)*ISWETC(I)
     XFLUX(I) = XFLUX(I)*ISWETTMP
     YFLUX(I) = YFLUX(I)*ISWETTMP
   END DO


!
!-------------------------SET BOUNDARY VALUES----------------------------------!
!

!  MODIFY BOUNDARY FLUX
      DO I=1,N
        IF(ISBCE(I) == 2) THEN
          XFLUX(I)=(XFLUX(I)+Fluxobn(I)*UA(I))*IUCP(I)
          YFLUX(I)=(YFLUX(I)+Fluxobn(I)*VA(I))*IUCP(I)
        ENDIF
      END DO


!  ADJUST FLUX FOR RIVER INFLOW
   IF(NUMQBC > 0) THEN
     IF(RIVER_INFLOW_LOCATION == 'node')THEN
       DO K=1,NUMQBC
         J=INODEQ(K)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         XFLUX(I1)=XFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
         XFLUX(I2)=XFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*COS(ANGLEQ(K))
         YFLUX(I2)=YFLUX(I2)-0.5_SP*QDIS(K)*VLCTYQ(K)*SIN(ANGLEQ(K))
       END DO
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       DO K=1,NUMQBC
         I1=ICELLQ(K)
         VLCTYQ(K)=QDIS(K)/QAREA(K)
         TEMP=QDIS(K)*VLCTYQ(K)
         XFLUX(I1)=XFLUX(I1)-TEMP*COS(ANGLEQ(K))
         YFLUX(I1)=YFLUX(I1)-TEMP*SIN(ANGLEQ(K))
       END DO
     END IF
   END IF

!  ADJUST FLUX FOR OPEN BOUNDARY MEAN FLOW



   if(dbg_set(dbg_sbr)) write(ipt,*) "End: advave_edge_gcn.F"

   END SUBROUTINE ADVAVE_EDGE_GCN
!==============================================================================|
