! $Id: external_step.F,v 1.1.2.14 2008/05/27 14:53:48 jqi Exp $
! $Name: New_Input $
! $Revision: 1.1.2.14 $
# if !defined(ONE_D_MODEL) && !defined (SEMI_IMPLICIT)
SUBROUTINE EXTERNAL_STEP

  USE MOD_UTILS
  USE ALL_VARS
  USE MOD_TIME
  USE MOD_OBCS
  USE MOD_WD
# if defined (VISIT)
  USE MOD_VISIT
# endif

# if defined(MULTIPROCESSOR)
  USE MOD_PAR
# endif

  USE MOD_ICE
  USE MOD_ICE2D

  IMPLICIT NONE
  REAL(SP) :: TMP
  INTEGER :: K, I, J, JN, J1
  integer  slice_node, i1

# if defined (ICE)
    real(SP) :: DUVI,ice_ocnDX,ice_ocnDY
# endif
!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: external_step"
 
  !! David for VISIT
  !!=======================================================
# if defined (VISIT)
  visit_cycle=IINT-1 
  visit_time_int = DAYS(IntTime)
  visit_time_ext = DAYS(ExtTime)
  Call VisitCheck
# endif
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
# if defined (GCN)
  CALL BCOND_GCN(9)
# else
  CALL BCOND_GCY(9)
# endif

# if defined (ICE)
       !  ggao 01042008 for ice ocean couple

       DO I=1,NT
          IF(ISICEC(I)==1) THEN
       !! (magnitude of relative ocean current)*rhow*drag*aice
         DUVI = DRAGW*SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2)  ! m/s
       !! ice/ocean stress
       ice_ocnDX = DUVI*((U(I,1)*COSW-V(I,1)*SINW)-(uice2(i)*cosw-vice2(i)*sinw))
       ice_ocnDY = DUVI*((V(I,1)*COSW+U(I,1)*SINW)-(vice2(i)*cosw+uice2(i)*sinw))

        wusurf2(I)=wusurf2(I)*(1.0_SP-AIU(I))-ice_ocnDX*aiu(I)*1.0E-3 !*RAMP
        wvsurf2(I)=wvsurf2(I)*(1.0_SP-AIU(I))-ice_ocnDY*aiu(I)*1.0E-3 !*RAMP
        ENDIF
        END DO
#   if defined (MULTIPROCESSOR)
        IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUSURF2,WVSURF2)
# endif

#endif


!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
  ELRK1 = EL1
  ELRK  = EL
  UARK  = UA
  VARK  = VA

# if defined (EQUI_TIDE)
  ELRK_EQI = EL_EQI
# endif
# if defined (ATMO_TIDE)
  ELRK_ATMO = EL_ATMO
# endif 
# if defined (AIR_PRESSURE)
  ELRK_AIR = EL_AIR
# endif

! New Open Boundary Condition ----5
# if defined (MEAN_FLOW)
  ELRKT  = ELT
  ELRKP  = ELP
  UARKNT = UANT   !change BKI
  VARKNT = VANT   !change BKI
  UARKN  = UAN    !change BKI
  VARKN  = VAN    !change BKI
# endif 
  
  
!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!

  DO K=1,4
     
     RKTIME = ExtTime + IMDTE * (ALPHA_RK(K) - 1.0_DP)

!     CALL PRINT_REAL_TIME(RKTIME,IPT,"RUNGE-KUTTA")


! New Open Boundary Condition ----6
# if defined (MEAN_FLOW)
     CALL BCOND_TIDE_2D
     !       CALL BCOND_NG_2D        ! change BKI
     CALL FLUX_OBN2D(K)
     CALL FLUX_OBC2D
# endif 
     
     
!FREE SURFACE AMPLITUDE UPDATE  --> ELF
     CALL EXTEL_EDGE(K)
# if defined (MULTIPROCESSOR)
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)
# endif
     
     
# if defined (EQUI_TIDE)
     CALL ELEVATION_EQUI
     ELF_EQI = ELRK_EQI +ALPHA_RK(K)*(ELF_EQI-ELRK_EQI) 
# endif
# if defined (ATMO_TIDE)
     CALL ELEVATION_ATMO
     ELF_ATMO = ELRK_ATMO +ALPHA_RK(K)*(ELF_ATMO-ELRK_ATMO) 
# endif
       
#       if defined (AIR_PRESSURE)
        CALL BCOND_PA_AIR 
        ELF_AIR = ELRK_AIR +ALPHA_RK(K)*(ELF_AIR-ELRK_AIR) 
#       endif
       

     ! New Open Boundary Condition ----7
# if defined (MEAN_FLOW)
     IF (ntidenode > 0) THEN
        DO I = 1, ntidenode     ! need to calculate NEXT_OBC column of ELPF
           J = I_TIDENODE_N(I)
           ELPF(I) = ELF(J) - ELTF(I)
        END DO
     END IF
     
     CALL EXTELPF_EDGE(K)
     
     IF (IOBCN > 0) THEN
        DO I = 1, IOBCN
           J = I_OBC_N(I)
           J1= I_OBC_NODE(J)
           ELF(J) = ELTF(J1) + ELPF(J1)
        END DO
     END IF

# else

     ! VALUES FOR THE OPEN BOUNDARY ARE ONLY UPDATED IN THE LOCAL DOMAIN
     ! THE HALO IS NOT SET HERE
#   if defined (GCN)
     CALL BCOND_GCN(1)
#   else
     CALL BCOND_GCY(1)
#   endif
     
     DO I=1,IBCN(1)
        JN = OBC_LST(1,I)
        J=I_OBC_N(JN)
        ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
     END DO
# endif

     
     ! DAVID ADDED THIS EXCHANGE CALL:
     ! IT SEEMS LIKELY THAT THE HALO VALUES OF ELF WILL BE USED
     ! BEFORE THEY ARE SET CORRECTLY OTHERWISE
# if defined (MULTIPROCESSOR)
     IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ELF)
# endif

     CALL N2E2D(ELF,ELF1)
          
     IF(WETTING_DRYING_ON)CALL WET_JUDGE

     CALL FLUX_OBN(K)

     !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
# if defined (GCN)
     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
# else
     CALL ADVAVE_EDGE_GCY(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
# endif


     CALL EXTUV_EDGE(K)

# if defined (GCN)
     CALL BCOND_GCN(2)
# else
     CALL BCOND_GCY(2)
# endif

       i1= nt/8
       do i=1, i1
         slice_node = (i-1)*8 + 5
         uaf(slice_node+2) = uaf(slice_node)

         slice_node = (i-1)*8 + 6
         uaf(slice_node+2) = uaf(slice_node)

         slice_node = (i-1)*8 + 4
         uaf(slice_node-2) = uaf(slice_node)

         slice_node = (i-1)*8 + 3
         uaf(slice_node-2) = uaf(slice_node)
       enddo
       vaf = 0.0_sp
           
# if defined (MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ELF)
# endif

!======== symmetric treatment ===================
   i1 = mt/5
   do i=1, i1
     slice_node = (i-1)*5 + 3
     elf(slice_node-2) = elf(slice_node)
     elf(slice_node-1) = elf(slice_node)
     elf(slice_node+1) = elf(slice_node)
     elf(slice_node+2) = elf(slice_node)
   enddo
!================================================
     
     !UPDATE WATER SURFACE ELEVATION
     CALL ASSIGN_ELM1_TO_ELM2

     EL  = ELF
     EL1 = ELF1
     
# if defined (EQUI_TIDE)
     EL_EQI = ELF_EQI
# endif

# if defined (ATMO_TIDE)
     EL_ATMO = ELF_ATMO
# endif       
# if defined (AIR_PRESSURE)
     EL_AIR = ELF_AIR
# endif       
     
     !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
     CALL N2E2D(EL,EL1)
     
     !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
     D   = H + EL
     D1  = H1 + EL1
     UA  = UAF
     VA  = VAF
     DTFA = D

     
     ! New Open Boundary Condition ----8
# if defined (MEAN_FLOW)
     ELT = ELTF
     UAT = UATF
     VAT = VATF
     ! The part below is equivalent to do both NODE MATCH and EXCHANGE for ELPF
#   if defined (MULTIPROCESSOR)
     IF (ntidenode > 0) THEN
        DO I = 1, ntidenode
           J = I_TIDENODE_N(I)
           ELPF(I) = ELF(J) - ELTF(I)
        END DO
     END IF
#   endif
     ELP = ELPF
     IF(K == 4 .and. nmfcell > 0)THEN
        DO I = 1, nmfcell
           J = I_MFCELL_N(I)
           OBC2D_X_TIDE(I) = OBC2D_X_TIDE(I) + UANT(I) * D1(J)
           OBC2D_Y_TIDE(I) = OBC2D_Y_TIDE(I) + VANT(I) * D1(J)
        ENDDO
     END IF
     CALL BCOND_BKI_2D(K)            ! change BKI
# endif

     
     !!ENSURE ALL CELLS ARE WET IN NO FLOOD/DRY CASE  
# if !defined (WET_DRY)
     CALL DEPTH_CHECK
# endif
     
     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
# if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UA,VA,D1)
#   if defined (WET_DRY)
     IF(PAR .AND. K==3)CALL AEXCHANGE(EC,MYID,NPROCS,UAS,VAS)
#   endif
# endif



# if !defined (TWO_D_MODEL)
     !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
     IF(K == 3)THEN
        UARD = UARD + UA*D1
        VARD = VARD + VA*D1
        EGF  = EGF  + EL/ISPLIT
        
#   if defined (EQUI_TIDE)
        EGF_EQI = EGF_EQI + EL_EQI/ISPLIT
#   endif

#   if defined (ATMO_TIDE)
        EGF_ATMO = EGF_ATMO + EL_ATMO/ISPLIT
#   endif       
#   if defined (AIR_PRESSURE)
        EGF_AIR = EGF_AIR + EL_AIR/ISPLIT
#   endif

#   if defined (WET_DRY)
        UARDS = UARDS + UAS*D1
        VARDS = VARDS + VAS*D1
#   endif
     END IF
     
     !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
     IF(K == 4.AND.IOBCN > 0) THEN
        DO I=1,IOBCN
           J=I_OBC_N(I)
           TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
           UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
        END DO
     END IF
# endif    
!end !defined (TWO_D_MODEL)

     
     !UPDATE WET/DRY FACTORS
     IF(WETTING_DRYING_ON)CALL WD_UPDATE(1)


  END DO     !! END RUNGE-KUTTA LOOP
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: external_step"

END SUBROUTINE EXTERNAL_STEP
# endif
