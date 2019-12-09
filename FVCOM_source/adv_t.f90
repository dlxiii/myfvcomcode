










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
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

SUBROUTINE ADV_T               

  !------------------------------------------------------------------------------|

  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_PAR
  USE BCS
  USE MOD_OBCS
  USE MOD_WD 
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE

  IMPLICIT NONE
  REAL(SP), DIMENSION(0:MT,KB)      :: XFLUX,XFLUX_ADV,RF
  REAL(SP), DIMENSION(0:MT)         :: PUPX,PUPY,PVPX,PVPY  
  REAL(SP), DIMENSION(0:MT)         :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
  REAL(SP), DIMENSION(3*(NT),KBM1)  :: DTIJ 
  REAL(SP), DIMENSION(3*(NT),KBM1)  :: UVN
  REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 !,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
  REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
  REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
  REAL(SP) :: FACT,FM1
  INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
  REAL(SP) :: T1MIN, T1MAX, T2MIN, T2MAX

!!$#  if defined (SPHERICAL)
!!$  REAL(DP) TY,TXPI,TYPI
!!$  REAL(DP) :: XTMP1,XTMP
!!$  REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
!!$  REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
!!$  REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
!!$  REAL(DP) :: TXPI_TMP,TYPI_TMP
!!$#  endif




  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: adv_t"
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
  !--Initialize Fluxes-----------------------------------------------------------!
  !
  XFLUX     = 0.0_SP
  XFLUX_ADV = 0.0_SP

  !
  !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
  !
!!#  if !defined (1)
  DO I=1,NCV
     I1=NTRG(I)
     !     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
        DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
        ! USE U,V
        UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I) 
     END DO
  END DO
!!#  else
!!  DO I=1,NCV
!!     I1=NTRG(I)
!!     !     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!        DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!        ! USE US,VS
!!        UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!  END DO
!!#  endif

  !
  !--Add the Shortwave Radiation Body Force--------------------------------------!
  !
  RF = 0.0_sp
  IF(HEATING_ON) THEN
     SELECT CASE(HEATING_TYPE)
     CASE('body')

        ! REMOVED START TIME FOR BODY HEATING
        !     TTIME=THOUR
        !     if(ttime < thour_hs) then
        !       RF=0.0_SP
        !     else

        DO  K=1,KBM1
           DO  I=1,M
              !         IF(TTIME < THOUR_HS)THEN
              !           RF(I,K)=0.0_SP
              !         ELSE
              ! REMOVED DEPTH CHECK FROM ECOM-SI, FVCOM DOES NOT NEED IT!
              !          ZEDP = 0.0_SP
              !           IF(DT(I) > 0.0_SP) THEN
              ZDEP=0.5_SP*(Z(I,K)+Z(I,K+1))*DT(I)
              !           END IF
              RF(I,K)=-SWRAD(I)*((RHEAT/ZETA1)*EXP(ZDEP/ZETA1) &
                   +((1-RHEAT)/ZETA2)*EXP(ZDEP/ZETA2))*DT(I)
              !         END IF
           END DO
        END DO
        !     endif

     CASE('flux')
        RF = 0.0_SP
     CASE DEFAULT
        CALL FATAL_ERROR('The surface heating type is set incorrectly:',&
             & TRIM(HEATING_TYPE))
     END SELECT
  END IF

!--ADJUST VOLUME'S HEAT CONTENT FOR EVAPORATION AND PRECIPITATION ------------------!
  IF (PRECIPITATION_ON) THEN
     RF(:,1)=RF(:,1)+ROFVROS*(QEVAP+QPREC)*T1(:,1)
  END IF

  !
  !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
  !

  DO K=1,KBM1
     PTPX  = 0.0_SP
     PTPY  = 0.0_SP
     PTPXD = 0.0_SP
     PTPYD = 0.0_SP
     DO I=1,M
        DO J=1,NTSN(I)-1
           I1=NBSN(I,J)
           I2=NBSN(I,J+1)
	 
!J. Ge for tracer advection
         IF(BACKWARD_ADVECTION==.FALSE.)THEN
           IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
            FFD=0.5_SP*(T1(I,K)+T1(I2,K)-TMEAN1(I,K)-TMEAN1(I2,K))
            FF1=0.5_SP*(T1(I,K)+T1(I2,K))
       	   ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
            FFD=0.5_SP*(T1(I1,K)+T1(I,K)-TMEAN1(I1,K)-TMEAN1(I,K))
            FF1=0.5_SP*(T1(I1,K)+T1(I,K))
	   ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
            FFD=0.5_SP*(T1(I,K)+T1(I,K)-TMEAN1(I,K)-TMEAN1(I,K))
            FF1=0.5_SP*(T1(I,K)+T1(I,K))
	   ELSE
            FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
            FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
	   END IF 
         ELSE
           IF(BACKWARD_STEP==1)THEN
              IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
               FFD=0.5_SP*((T0(I,K)+T1(I,K))*0.5+(T0(I2,K)+T1(I2,K))*0.5-TMEAN1(I,K)-TMEAN1(I2,K))
               FF1=0.5_SP*((T0(I,K)+T1(I,K))*0.5+(T0(I2,K)+T1(I2,K))*0.5)
              ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((T0(I1,K)+T1(I1,K))*0.5+(T0(I,K)+T1(I,K))*0.5-TMEAN1(I1,K)-TMEAN1(I,K))
               FF1=0.5_SP*((T0(I1,K)+T1(I1,K))*0.5+(T0(I,K)+T1(I,K))*0.5)
	      ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((T0(I,K)+T1(I,K))*0.5+(T0(I,K)+T1(I,K))*0.5-TMEAN1(I,K)-TMEAN1(I,K))
               FF1=0.5_SP*((T0(I,K)+T1(I,K))*0.5+(T0(I,K)+T1(I,K))*0.5)
	      ELSE
               FFD=0.5_SP*((T0(I1,K)+T1(I1,K))*0.5+(T0(I2,K)+T1(I2,K))*0.5-TMEAN1(I1,K)-TMEAN1(I2,K))
               FF1=0.5_SP*((T0(I1,K)+T1(I1,K))*0.5+(T0(I2,K)+T1(I2,K))*0.5)
	      END IF 
           ELSEIF(BACKWARD_STEP==2)THEN
              IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
               FFD=0.5_SP*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP+(T2(I2,K)+T0(I2,K)+T1(I2,K))/3.0_SP-TMEAN1(I,K)-TMEAN1(I2,K))
               FF1=0.5_SP*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP+(T2(I2,K)+T0(I2,K)+T1(I2,K))/3.0_SP)
              ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((T2(I1,K)+T0(I1,K)+T1(I1,K))/3.0_SP+(T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP-TMEAN1(I1,K)-TMEAN1(I,K))
               FF1=0.5_SP*((T2(I1,K)+T0(I1,K)+T1(I1,K))/3.0_SP+(T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP)
	      ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
               FFD=0.5_SP*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP+(T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP-TMEAN1(I,K)-TMEAN1(I,K))
               FF1=0.5_SP*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP+(T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP)
	      ELSE
               FFD=0.5_SP*((T2(I1,K)+T0(I1,K)+T1(I1,K))/3.0_SP+(T2(I2,K)+T0(I2,K)+T1(I2,K))/3.0_SP-TMEAN1(I1,K)-TMEAN1(I2,K))
               FF1=0.5_SP*((T2(I1,K)+T0(I1,K)+T1(I1,K))/3.0_SP+(T2(I2,K)+T0(I2,K)+T1(I2,K))/3.0_SP)
	      END IF 
            ENDIF
        ENDIF
!J. Ge for tracer advection
	 
!!$#        if defined (SPHERICAL)
!!$           XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$           XTMP1 = VX(I2)-VX(I1)
!!$           IF(XTMP1 >  180.0_SP)THEN
!!$              XTMP = -360.0_SP*TPI+XTMP
!!$           ELSE IF(XTMP1 < -180.0_SP)THEN
!!$              XTMP =  360.0_SP*TPI+XTMP
!!$           END IF
!!$           TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$           TYPI=(VY(I1)-VY(I2))*tpi
!!$           ! ERROR HERE
!!$           !#    if defined (NORTHPOLE)
!!$           IF(NODE_NORTHAREA(I) == 1)THEN
!!$              VX1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * COS(VX(I1)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I1)*DEG2RAD))
!!$              VY1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * SIN(VX(I1)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I1)*DEG2RAD))
!!$
!!$              VX2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * COS(VX(I2)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I2)*DEG2RAD))
!!$              VY2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * SIN(VX(I2)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I2)*DEG2RAD))
!!$
!!$              TXPI = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+SIN(VY(I)*DEG2RAD)))
!!$              TYPI = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+SIN(VY(I)*DEG2RAD)))
!!$              IF(I /= NODE_NORTHPOLE)THEN
!!$                 TXPI_TMP = TYPI*COS(VX(I)*DEG2RAD)-TXPI*SIN(VX(I)*DEG2RAD)
!!$                 TYPI_TMP = TXPI*COS(VX(I)*DEG2RAD)+TYPI*SIN(VX(I)*DEG2RAD)
!!$                 TYPI_TMP = -TYPI_TMP
!!$
!!$                 TXPI = TXPI_TMP
!!$                 TYPI = TYPI_TMP
!!$              END IF
!!$           END IF
!!$           ! END ERROR
!!$           
!!$           PTPX(I)=PTPX(I)+FF1*TYPI
!!$           PTPY(I)=PTPY(I)+FF1*TXPI
!!$           PTPXD(I)=PTPXD(I)+FFD*TYPI
!!$           PTPYD(I)=PTPYD(I)+FFD*TXPI
!!$#        else
!!$           PTPX(I)=PTPX(I)+FF1*(VY(I1)-VY(I2))
!!$           PTPY(I)=PTPY(I)+FF1*(VX(I2)-VX(I1))
!!$           PTPXD(I)=PTPXD(I)+FFD*(VY(I1)-VY(I2))
!!$           PTPYD(I)=PTPYD(I)+FFD*(VX(I2)-VX(I1))
!!$#        endif

           
           PTPX(I)=PTPX(I)+FF1*DLTYTRIE(i,j)
           PTPY(I)=PTPY(I)+FF1*DLTXTRIE(i,j)
           PTPXD(I)=PTPXD(I)+FFD*DLTYTRIE(i,j)
           PTPYD(I)=PTPYD(I)+FFD*DLTXTRIE(i,j)

        END DO
! gather all neighboring control volumes connecting at dam node 
        PTPX(I)=PTPX(I)/ART2(I)
        PTPY(I)=PTPY(I)/ART2(I)
        PTPXD(I)=PTPXD(I)/ART2(I)
        PTPYD(I)=PTPYD(I)/ART2(I)

     END DO

     IF(K == KBM1)THEN
        DO I=1,M
           PFPXB(I) = PTPX(I)
           PFPYB(I) = PTPY(I)
        END DO
     END IF

     DO I=1,M

        VISCOFF(I) = VISCOFH(I,K)

     END DO
     IF(K == KBM1) THEN
        AH_BOTTOM(1:M) = (FACT*VISCOFF(1:M) + FM1) * NN_HVC(1:M)
     END IF


     DO I=1,NCV_I
        IA=NIEC(I,1)
        IB=NIEC(I,2)


!!$        XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!!$        YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!$#      if defined (SPHERICAL)
!!$        X1_DP=XIJE(I,1)
!!$        Y1_DP=YIJE(I,1)
!!$        X2_DP=XIJE(I,2)
!!$        Y2_DP=YIJE(I,2)
!!$        CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!!$        XI=XII		
!!$        XTMP  = XI*TPI-VX(IA)*TPI
!!$        XTMP1 = XI-VX(IA)
!!$        IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$        END IF
!!$        DXA=XTMP*COS(DEG2RAD*VY(IA))  
!!$        DYA=(YI-VY(IA))*TPI
!!$
!!$        XTMP  = XI*TPI-VX(IB)*TPI
!!$        XTMP1 = XI-VX(IB)
!!$        IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$        END IF
!!$
!!$        DXB=XTMP*COS(DEG2RAD*VY(IB))
!!$        DYB=(YI-VY(IB))*TPI
!!$#      else
!!$        DXA=XI-VX(IA)
!!$        DYA=YI-VY(IA)
!!$        DXB=XI-VX(IB)
!!$        DYB=YI-VY(IB)
!!$#      endif
!!$        FIJ1=T1(IA,K)+DXA*PTPX(IA)+DYA*PTPY(IA)
!!$        FIJ2=T1(IB,K)+DXB*PTPX(IB)+DYB*PTPY(IB)

!J. Ge for tracer advection
        IF(BACKWARD_ADVECTION==.FALSE.)THEN
          FIJ1=T1(IA,K)+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
          FIJ2=T1(IB,K)+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)
        ELSE
          IF(BACKWARD_STEP==1)THEN
            FIJ1=(T0(IA,K)+T1(IA,K))*0.5+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
            FIJ2=(T0(IB,K)+T1(IB,K))*0.5+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)
          ELSEIF(BACKWARD_STEP==2)THEN
            FIJ1=(T2(IA,K)+T0(IA,K)+T1(IA,K))/3.0_SP+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
            FIJ2=(T2(IA,K)+T0(IB,K)+T1(IB,K))/3.0_SP+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)
          ENDIF
        ENDIF
!J. Ge for tracer advection

        T1MIN=MINVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
        T1MIN=MIN(T1MIN, T1(IA,K))
        T1MAX=MAXVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
        T1MAX=MAX(T1MAX, T1(IA,K))
        T2MIN=MINVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
        T2MIN=MIN(T2MIN, T1(IB,K))
        T2MAX=MAXVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
        T2MAX=MAX(T2MAX, T1(IB,K))
        IF(FIJ1 < T1MIN) FIJ1=T1MIN
        IF(FIJ1 > T1MAX) FIJ1=T1MAX
        IF(FIJ2 < T2MIN) FIJ2=T2MIN
        IF(FIJ2 > T2MAX) FIJ2=T2MAX

        UN=UVN(I,K)

!        VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
        ! David moved HPRNU and added HVC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))


        TXX=0.5_SP*(PTPXD(IA)+PTPXD(IB))*VISCOF
        TYY=0.5_SP*(PTPYD(IA)+PTPYD(IB))*VISCOF

        FXX=-DTIJ(I,K)*TXX*DLTYE(I)
        FYY= DTIJ(I,K)*TYY*DLTXE(I)

        EXFLUX=-UN*DTIJ(I,K)* &
             ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

        XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
        XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

        XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
        XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)


     END DO


  END DO !! K LOOP

  IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,XFLUX,XFLUX_ADV)

  DO K=1,KBM1
     IF(IOBCN > 0) THEN
        DO I=1,IOBCN
           I1=I_OBC_N(I)
           XFLUX_OBC(I,K)=XFLUX_ADV(I1,K)
        END DO
     END IF
  END DO



  !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
  !




  !--------------------------------------------------------------------
  !   The central difference scheme in vertical advection
  !--------------------------------------------------------------------
  DO I=1, M
     IF(ISWETN(I)*ISWETNT(I) == 1) THEN


        DO K=1, KBM1



!J. Ge for tracer advection
           IF(BACKWARD_ADVECTION==.FALSE.)THEN
             IF(K == 1) THEN
              TEMP=-WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
             ELSE IF(K == KBM1) THEN
              TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
             ELSE
              TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
             END IF
           ELSE
             IF(BACKWARD_STEP==1)THEN
               IF(K == 1) THEN
                TEMP=-WTS(I,K+1)*((T0(I,K)+T1(I,K))*0.5*DZ(I,K+1)+(T0(I,K+1)+T1(I,K+1))*0.5*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
               ELSE IF(K == KBM1) THEN
                TEMP= WTS(I,K)*((T0(I,K)+T1(I,K))*0.5*DZ(I,K-1)+(T0(I,K-1)+T1(I,K-1))*0.5*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
               ELSE
                 TEMP= WTS(I,K)*((T0(I,K)+T1(I,K))*0.5*DZ(I,K-1)+(T0(I,K-1)+T1(I,K-1))*0.5*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*((T0(I,K)+T1(I,K))*0.5*DZ(I,K+1)+(T0(I,K+1)+T1(I,K+1))*0.5*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
               END IF
             ELSEIF(BACKWARD_STEP==2)THEN
               IF(K == 1) THEN
                TEMP=-WTS(I,K+1)*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP*DZ(I,K+1)+(T2(I,K+1)+T0(I,K+1)+T1(I,K+1))/3.0_SP*DZ(I,K))/   &
                   (DZ(I,K)+DZ(I,K+1))
               ELSE IF(K == KBM1) THEN
                TEMP= WTS(I,K)*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP*DZ(I,K-1)+(T2(I,K-1)+T0(I,K-1)+T1(I,K-1))/3.0_SP*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
               ELSE
                 TEMP= WTS(I,K)*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP*DZ(I,K-1)+(T2(I,K-1)+T0(I,K-1)+T1(I,K-1))/3.0_SP*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
                   WTS(I,K+1)*((T2(I,K)+T0(I,K)+T1(I,K))/3.0_SP*DZ(I,K+1)+(T2(I,K+1)+T0(I,K+1)+T1(I,K+1))/3.0_SP*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
               END IF
             ENDIF
           ENDIF
!J. Ge for tracer advection


           IF(ISONB(I) == 2) THEN
              !         XFLUX(I,K)=TEMP*ART1(I)/DZ(I,K)
              XFLUX(I,K)=TEMP*ART1(I)
           ELSE
              !         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(I,K)
              XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)

           END IF
        ENDDO
     END IF
  END DO  !!K LOOP
  !
  !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
  !
  IF(RIVER_TS_SETTING == 'calculated') THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              JJ=INODEQ(J)
              STPOINT=TDIS(J)  
              DO K=1,KBM1
              !   STPOINT    = T1(JJ,K)
                 !             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
                 XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT
              END DO
           END DO
        END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              J1=N_ICELLQ(J,1)
              J2=N_ICELLQ(J,2)
              STPOINT=TDIS(J) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
              DO K=1,KBM1
                 !STPOINT1 = T1(J1,K)
                 !STPOINT2 = T1(J2,K)
                 !             XFLUX(J1,K)=XFLUX(J1,K)-  &
                 !                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT1/DZ1(J1,K)
                 !             XFLUX(J2,K)=XFLUX(J2,K)-  &
                 !                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT2/DZ1(J2,K)
                 XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT   !1
                 XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT   !2
              END DO
           END DO
        END IF
     END IF
  END IF
  !---------------------------------------------------------------------


  ! APPLY GROUND WATER TEMPERATURE FORCING
  IF(GROUNDWATER_ON .and. GROUNDWATER_TEMP_ON)THEN
     DO I=1,M
        XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*BFWTMP(I)
     END DO
  ELSEIF(GROUNDWATER_ON) THEN
     DO I=1,M
 !J. Ge for tracer advection
        IF(BACKWARD_ADVECTION==.FALSE.)THEN
          XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*T1(I,KBM1)
        ELSE
          IF(BACKWARD_STEP==1)THEN
            XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*(T0(I,KBM1)+T1(I,KBM1))*0.5
          ELSEIF(BACKWARD_STEP==2)THEN
            XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*(T2(I,KBM1)+T0(I,KBM1)+T1(I,KBM1))/3.0_SP
          ENDIF
        ENDIF
!J. Ge for tracer advection
     END DO
  END IF


  !
  !--Update Temperature----------------------------------------------------------!
  !


  DO I=1,M
     IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
        DO K=1,KBM1
           XFLUX(I,K) = XFLUX(I,K) - RF(I,K)*ART1(I)    !/DZ(I,K)
           !       TF1(I,K)=(T1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/DTFA(I))

           TF1(I,K)=(T1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))

        END DO
     ELSE
        DO K=1,KBM1
           TF1(I,K)=T1(I,K)
        END DO
     END IF
  END DO



  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: adv_t"

END SUBROUTINE ADV_T
!==============================================================================|
