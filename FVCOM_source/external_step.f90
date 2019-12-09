










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

SUBROUTINE EXTERNAL_STEP

  USE MOD_NESTING
  USE MOD_UTILS
  USE ALL_VARS
  USE MOD_TIME
  USE MOD_OBCS
  USE MOD_WD

  USE MOD_PAR

  USE MOD_ICE
  USE MOD_ICE2D




  IMPLICIT NONE
  REAL(SP) :: TMP
  INTEGER :: K, I, J, JN, J1,i1,i2

!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: external_step"
 
  !! David for VISIT
  !!=======================================================
  !!=======================================================


!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  IF(IRAMP /= 0) THEN
     TMP = real(IINT-1,sp)+real(IEXT,sp)/real(ISPLIT,sp)
     RAMP=TANH(TMP/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF

!
!------SURFACE BOUNDARY CONDITIONS FOR EXTERNAL MODEL--------------------------!
!
  CALL BCOND_GCN(9,0)



!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
  ELRK1 = EL1
  ELRK  = EL
  UARK  = UA
  VARK  = VA


! New Open Boundary Condition ----5
  
  
!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!

  DO K=1,4
     
     RKTIME = ExtTime + IMDTE * (ALPHA_RK(K) - 1.0_DP)

!     CALL PRINT_REAL_TIME(RKTIME,IPT,"RUNGE-KUTTA")


! New Open Boundary Condition ----6
     
     
!FREE SURFACE AMPLITUDE UPDATE  --> ELF
     CALL EXTEL_EDGE(K)
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)
     
     
       
       

     ! New Open Boundary Condition ----7

     ! VALUES FOR THE OPEN BOUNDARY ARE ONLY UPDATED IN THE LOCAL DOMAIN
     ! THE HALO IS NOT SET HERE
     CALL BCOND_GCN(1,K)
 
     IF(NESTING_ON )THEN
        CALL SET_VAR(ExtTime,EL=ELF)
     END IF
     
     DO I=1,IBCN(1)
        JN = OBC_LST(1,I)
        J=I_OBC_N(JN)
        ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
     END DO

     
     ! DAVID ADDED THIS EXCHANGE CALL:
     ! IT SEEMS LIKELY THAT THE HALO VALUES OF ELF WILL BE USED
     ! BEFORE THEY ARE SET CORRECTLY OTHERWISE
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)

!---------------For Dam Model-----------------------
! Jadon


     CALL N2E2D(ELF,ELF1)
          
     IF(WETTING_DRYING_ON)CALL WET_JUDGE

     CALL FLUX_OBN(K)

     !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff


     CALL EXTUV_EDGE(K)


     CALL BCOND_GCN(2,K)


     IF(NESTING_ON )THEN
       CALL SET_VAR(ExtTime,UA=UAF)
       CALL SET_VAR(ExtTime,VA=VAF)
     END IF

     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ELF)


     
     !UPDATE WATER SURFACE ELEVATION
     CALL ASSIGN_ELM1_TO_ELM2

     EL  = ELF
     EL1 = ELF1


     
     !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
     CALL N2E2D(EL,EL1)
     
     !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
     D   = H + EL
     D1  = H1 + EL1
     UA  = UAF
     VA  = VAF
     DTFA = D

     
     ! New Open Boundary Condition ----8

     
     !!ENSURE ALL CELLS ARE WET IN NO FLOOD/DRY CASE  
     
     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA,D1)
!!#   if defined (1)
!!     IF(PAR .AND. K==3)CALL AEXCHANGE(EC,MYID,NPROCS,UAS,VAS)
!!#   endif

     !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
     IF(K == 3)THEN
        UARD = UARD + UA*D1
        VARD = VARD + VA*D1
        EGF  = EGF  + EL/ISPLIT
        


!!#   if defined (1)
!!        UARDS = UARDS + UAS*D1
!!        VARDS = VARDS + VAS*D1
!!#   endif
     END IF
     
     !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
     IF(K == 4.AND.IOBCN > 0) THEN
        DO I=1,IOBCN
           J=I_OBC_N(I)
           TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
           UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
        END DO
     END IF
!end !defined (TWO_D_MODEL)

     
     !UPDATE WET/DRY FACTORS
     IF(WETTING_DRYING_ON)CALL WD_UPDATE(1)

  END DO     !! END RUNGE-KUTTA LOOP
  

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: external_step"

END SUBROUTINE EXTERNAL_STEP
