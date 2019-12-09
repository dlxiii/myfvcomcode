










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
   SUBROUTINE ADVECTION_EDGE_GCY(XFLUX,YFLUX)
!==============================================================================|
!   Calculate the Advection and Diffusion Terms of 3D Velocity Field           |
!   These Terms will be vertically integrated to form the Mean Terms in        |
!   the Gx and Gy Terms of the External Mode Equation                          |
!   Ghost cell boundary conditions are used here                               |
!==============================================================================|

   USE MOD_UTILS
   USE ALL_VARS
   USE BCS
   USE MOD_SPHERICAL
   USE MOD_NORTHPOLE
   USE MOD_WD

   IMPLICIT NONE
   REAL(SP), INTENT(OUT), DIMENSION(0:NT,KB) :: XFLUX,YFLUX
   REAL(SP) :: DIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2
   REAL(SP) :: ISWETTMP

   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP) :: UK1(KB),UK2(KB),UK3(KB),UK4(KB),UK5(KB),UK6(KB), &
               VK1(KB),VK2(KB),VK3(KB),VK4(KB),VK5(KB),VK6(KB)
!------------------------------------------------------------------------------|

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advection_edge_gcy.F"


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
!--Initialize Variables--------------------------------------------------------|
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP

!
!--Loop Over Edges and Accumulate Fluxes-For Each Element----------------------|
!

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
    IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
     J1=IENODE(I,1)
     J2=IENODE(I,2)

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)
     XIJA=XIJC(I)-XC(IA)
     YIJA=YIJC(I)-YC(IA)
     XIJB=XIJC(I)-XC(IB)
     YIJB=YIJC(I)-YC(IB)

     UK1 = U(K1,:)
     UK2 = U(K2,:)
     UK3 = U(K3,:)
     UK4 = U(K4,:)
     UK5 = U(K5,:)
     UK6 = U(K6,:)
     VK1 = V(K1,:)
     VK2 = V(K2,:)
     VK3 = V(K3,:)
     VK4 = V(K4,:)
     VK5 = V(K5,:)
     VK6 = V(K6,:)

     IF(K1 == 0) CALL GHOSTUV3(IA,1,UK1,VK1)
     IF(K2 == 0) CALL GHOSTUV3(IA,2,UK2,VK2)
     IF(K3 == 0) CALL GHOSTUV3(IA,3,UK3,VK3)
     IF(K4 == 0) CALL GHOSTUV3(IB,1,UK4,VK4)
     IF(K5 == 0) CALL GHOSTUV3(IB,2,UK5,VK5)
     IF(K6 == 0) CALL GHOSTUV3(IB,3,UK6,VK6)

     DO K=1,KBM1

       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))

       !!FORM THE LEFT FLUX
       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*UK1(K)+A1U(IA,3)*UK2(K)+A1U(IA,4)*UK3(K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*UK1(K)+A2U(IA,3)*UK2(K)+A2U(IA,4)*UK3(K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*VK1(K)+A1U(IA,3)*VK2(K)+A1U(IA,4)*VK3(K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*VK1(K)+A2U(IA,3)*VK2(K)+A2U(IA,4)*VK3(K)
       UIJ1=U(IA,K)+COFA1*XIJA+COFA2*YIJA
       VIJ1=V(IA,K)+COFA5*XIJA+COFA6*YIJA

       !!FORM THE RIGHT FLUX
       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*UK4(K)+A1U(IB,3)*UK5(K)+A1U(IB,4)*UK6(K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*UK4(K)+A2U(IB,3)*UK5(K)+A2U(IB,4)*UK6(K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*VK4(K)+A1U(IB,3)*VK5(K)+A1U(IB,4)*VK6(K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*VK4(K)+A2U(IB,3)*VK5(K)+A2U(IB,4)*VK6(K)
       UIJ2=U(IB,K)+COFA3*XIJB+COFA4*YIJB
       VIJ2=V(IB,K)+COFA7*XIJB+COFA8*YIJB

       !!COMPUTE THE NORMAL VELOCITY ACROSS THE EDGE
       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       UN=VIJ*DLTXC(I) - UIJ*DLTYC(I)

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2) + FM1)
       ! David moved HPRNU and added HVC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))


       !!UPWIND THE ADVECTIVE FLUX
       XADV=DIJ*UN*((1.0_SP-SIGN(1.0_SP,UN))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN))*UIJ1)*0.5_SP
       YADV=DIJ*UN*((1.0_SP-SIGN(1.0_SP,UN))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN))*VIJ1)*0.5_SP


       !!COMPUTE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)


       !!ACCUMULATE THE FLUX
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB


     END DO
    END IF
   END DO


   DO I=1,N
     ISWETTMP = ISWETCT(I)*ISWETC(I)
     DO K=1,KBM1
       XFLUX(I,K) = XFLUX(I,K)*ISWETTMP
       YFLUX(I,K) = YFLUX(I,K)*ISWETTMP
     END DO
   END DO


!
!--Boundary Conditions on Flux-------------------------------------------------|
!
   DO I=1,N
     IF(ISBCE(I) == 2)THEN
       DO K=1,KBM1
         XFLUX(I,K)=0.0_SP
         YFLUX(I,K)=0.0_SP
       END DO
     END IF
   END DO

!
!--Adjust Flux for Fresh Water Inflow------------------------------------------|
!

   IF(NUMQBC > 0) THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
       DO II=1,NUMQBC
         J=INODEQ(II)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ(J,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           XFLUX(I2,K)=XFLUX(I2,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
!           YFLUX(I2,K)=YFLUX(I2,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           XFLUX(I2,K)=XFLUX(I2,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
           YFLUX(I2,K)=YFLUX(I2,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       DO II=1,NUMQBC
         I1=ICELLQ(II)
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ1(I1,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ1(I1,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ1(I1,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     END IF
   END IF

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: advection_edge_gcy.F"

   RETURN
   END SUBROUTINE ADVECTION_EDGE_GCY
!==============================================================================|



