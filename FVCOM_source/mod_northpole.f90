










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

MODULE MOD_NORTHPOLE
   USE ALL_VARS 
   USE MOD_SPHERICAL
   USE MOD_PAR
   USE MOD_UTILS

   IMPLICIT NONE
   SAVE
   INTEGER :: NODE_NORTHPOLE            !Node index at the north pole point
   INTEGER :: MP,NP,NPE,NPCV
   INTEGER, ALLOCATABLE :: NODE_NORTHAREA(:)
   INTEGER, ALLOCATABLE :: CELL_NORTHAREA(:)   
   INTEGER, ALLOCATABLE :: NPEDGE_LST(:)   
   INTEGER, ALLOCATABLE :: NCEDGE_LST(:)   
   INTEGER, ALLOCATABLE :: MP_LST(:),NP_LST(:)   
   REAL(DP), ALLOCATABLE :: A1U_XY(:,:),A2U_XY(:,:)
   REAL(DP), ALLOCATABLE :: AW0_XY(:,:),AWX_XY(:,:),AWY_XY(:,:)


   CONTAINS
!==============================================================================|
     SUBROUTINE FIND_NORTHPOLE
     IMPLICIT NONE
     INTEGER :: I,ITMP,NODE_NORTHPOLE_GL,IERR,NDOM
     INTEGER,ALLOCATABLE :: TMP(:)
     
     NODE_NORTHPOLE = 0
!     NODE_NORTHPOLE_GL = 0
     NP = 0
     MP = 0
     NPE = 0
     NPCV = 0

     NDOM=0
     DO I = 1,MT
!        IF(ABS(VY(I)-90.0_SP) .LE. (NEAREST(90.0_SP,1.0_sp) - 90.0_SP))THEN
        IF(ABS(VY(I)-90.0_SP) .LE. 10E-4)THEN
           NODE_NORTHPOLE = I
           ndom=ndom+1
        END IF
     END DO
     
     IF (nDOM .GT. 1) CALL Fatal_Error('Found more than one north pole node in the grid',&
          &'This should never happen, TGE should crash first?')


     ! SINCE SPHERICAL NOW ALWAYS LOOKS FOR A NORTH POLE, DO NOT CALL
     ! FATAL_ERROR IF WE DON'T FIND A NORTH POLE!

     NODE_NORTHPOLE_GL = NGID_X(NODE_NORTHPOLE)


     IF (PAR) THEN
        ! GATHER THE NORTH POLE NODE GLOBAL ID TO ALL PROCS AND COMPARE RESULTS
        ALLOCATE(TMP(NPROCS)); tmp = 0
        CALL MPI_ALLGATHER(NODE_NORTHPOLE_GL,1,MPI_INTEGER,TMP,1,MPI_INTEGER,MPI_FVCOM_GROUP,IERR)

        NDOM=0
        DO I=1,NPROCS
           IF( (0 .LT. TMP(I)) .and. (TMP(I) .LE. MGL)) THEN
              NDOM=NDOM+1
              
              IF(NODE_NORTHPOLE_GL.EQ.0) THEN
                 ! Set the north pole global node
                 NODE_NORTHPOLE_GL = TMP(I)
              ELSE
                 ! Make sure it is the same as you global node number
                 IF (NODE_NORTHPOLE_GL.NE.TMP(I)) CALL FATAL_ERROR &
                      &("TWO DOMAINS REPORT DIFFERENT GLOBAL NORTH POLE NODES?")
              END IF
           END IF
        END DO
        DEALLOCATE(TMP)
     END IF

     if(dbg_set(dbg_log)) THEN
        WRITE(IPT,*) "! //////////////////////////////////////////////"
        IF(ndom .eq.0) THEN
           WRITE(IPT,*) "! NO NORTH POLE FOUND; PROCEED WITH SPHERICAL MODEL!"
        ELSE
           WRITE(IPT,*) "! FOUND THE GLOBAL NORTH POLE NODE::", NODE_NORTHPOLE_GL
           IF(NDOM .gt. 1)THEN
              WRITE(IPT,*) "THE NORTH POLE IS ON A PROCESSOR BOUNDARY"
              WRITE(IPT,*) "Be afraid, be very afraid...."
              CALL FATAL_ERROR("THE NORTH POLE IS ON A PROCESSOR BOUNDARY")
           END IF
        END IF
        WRITE(IPT,*) "! //////////////////////////////////////////////"
     end if

!     ! ALL OTHER PROCESSORS ESCAPE HERE
!     IF (NODE_NORTHPOLE .EQ. 0) RETURN
     
     ALLOCATE(NODE_NORTHAREA(0:MT)); NODE_NORTHAREA = 0
     ALLOCATE(CELL_NORTHAREA(0:NT)); CELL_NORTHAREA = 0

     ! ALL OTHER PROCESSORS ESCAPE HERE
     IF (NODE_NORTHPOLE .EQ. 0) RETURN
     

     ITMP = NODE_NORTHPOLE
     ALLOCATE(TMP(MT)); TMP = 0     
     MP = 0
     DO I=1,NTSN(ITMP)-1
        MP = MP + 1
        TMP(MP) = NBSN(ITMP,I)
        NODE_NORTHAREA(NBSN(ITMP,I)) = 1
     END DO
     MP = MP + 1
     TMP(MP) = ITMP
     NODE_NORTHAREA(ITMP) = 1
     
     ALLOCATE(MP_LST(MP))
     MP_LST(1:MP) = TMP(1:MP)
     DEALLOCATE(TMP)
     
     ALLOCATE(TMP(NT)); TMP = 0          
     NP = 0 
     DO I=1,NTVE(ITMP)
        NP = NP + 1
        TMP(NP) = NBVE(ITMP,I)
        CELL_NORTHAREA(NBVE(ITMP,I)) = 1
     END DO
     
     ALLOCATE(NP_LST(NP))
     NP_LST(1:NP) = TMP(1:NP)
     DEALLOCATE(TMP)
     
     
     RETURN
     END SUBROUTINE FIND_NORTHPOLE
!==============================================================================|
     
!==============================================================================|

     SUBROUTINE FIND_CELLSIDE
     
     IMPLICIT NONE
     INTEGER  ::  I,IA,IB
     INTEGER, ALLOCATABLE :: TEMP(:)
     
     ! ALL OTHER PROCESSORS ESCAPE HERE
     IF (NODE_NORTHPOLE .EQ. 0) RETURN


    ALLOCATE(TEMP(NE));  TEMP = ZERO
     NPE = 0
     
     DO I=1,NE
       IA = IEC(I,1)
       IB = IEC(I,2)
       IF(CELL_NORTHAREA(IA) == 1 .OR. CELL_NORTHAREA(IB) == 1)THEN
         NPE = NPE + 1
	 TEMP(NPE) = I
       END IF
     END DO
     
     ALLOCATE(NPEDGE_LST(NPE))
     NPEDGE_LST(1:NPE) = TEMP(1:NPE)
     DEALLOCATE(TEMP)
     
     ALLOCATE(TEMP(NCV));  TEMP = ZERO
     NPCV = 0
     
     DO I=1,NCV
       IA = NIEC(I,1)
       IB = NIEC(I,2)
       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         NPCV = NPCV + 1
	 TEMP(NPCV) = I
       END IF
     END DO
     
     ALLOCATE(NCEDGE_LST(NPCV))
     NCEDGE_LST(1:NPCV) = TEMP(1:NPCV)
     DEALLOCATE(TEMP)
     
     RETURN
     END SUBROUTINE FIND_CELLSIDE
       	 
!==============================================================================|

!==============================================================================|
   SUBROUTINE ADVAVE_EDGE_XY(XFLUX,YFLUX,IFCETA)

   USE MOD_OBCS
   IMPLICIT NONE
   INTEGER  :: I,J,K,IA,IB,J1,J2,K1,K2,K3,I1,I2,II
   REAL(SP) :: DIJ,ELIJ,XIJ,YIJ,UIJ,VIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ 
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   REAL(SP) :: FACT,FM1,ISWETTMP

   REAL(SP) :: TPA,TPB

   REAL(SP) :: UIJ1_TMP,VIJ1_TMP,UIJ2_TMP,VIJ2_TMP,TXXIJ_TMP,TYYIJ_TMP
   REAL(SP) :: XADV_TMP,YADV_TMP,PSTX_TMP,PSTY_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,UN_TMP
   REAL(SP) :: DLTXC_TMP,DLTYC_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP
   REAL(SP) :: UAIA,VAIA,UAIB,VAIB,UAK1,VAK1,UAK2,VAK2,UAK3,VAK3
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP
   REAL(SP) :: XIJ_TMP,YIJ_TMP
   
   REAL(SP) :: IFCETA
  
   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
!------------------------------------------------------------------------------!

   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
 
   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advave_edge_XY"

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
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     IF(CELL_NORTHAREA(IA) == 1)THEN
       XFLUX(IA) = 0.0_SP
       YFLUX(IA) = 0.0_SP
       PSTX(IA)  = 0.0_SP
       PSTY(IA)  = 0.0_SP
     END IF  
     IF(CELL_NORTHAREA(IB) == 1)THEN  
       XFLUX(IB) = 0.0_SP
       YFLUX(IB) = 0.0_SP
       PSTX(IB)  = 0.0_SP
       PSTY(IB)  = 0.0_SP
     END IF  
   END DO  
!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!

   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

     DIJ=0.5_SP*(D(J1)+D(J2))
     ELIJ=0.5_SP*(EL(J1)+EL(J2))

!   ggao 0103-2007   consider the sea surface pressure change
!     change end



!    FLUX FROM LEFT
     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
         
     UAIA = -VA(IA)*COS(XC(IA)*DEG2RAD)-UA(IA)*SIN(XC(IA)*DEG2RAD)
     VAIA = -VA(IA)*SIN(XC(IA)*DEG2RAD)+UA(IA)*COS(XC(IA)*DEG2RAD)
     UAK1 = -VA(K1)*COS(XC(K1)*DEG2RAD)-UA(K1)*SIN(XC(K1)*DEG2RAD)
     VAK1 = -VA(K1)*SIN(XC(K1)*DEG2RAD)+UA(K1)*COS(XC(K1)*DEG2RAD)
     UAK2 = -VA(K2)*COS(XC(K2)*DEG2RAD)-UA(K2)*SIN(XC(K2)*DEG2RAD)
     VAK2 = -VA(K2)*SIN(XC(K2)*DEG2RAD)+UA(K2)*COS(XC(K2)*DEG2RAD)
     UAK3 = -VA(K3)*COS(XC(K3)*DEG2RAD)-UA(K3)*SIN(XC(K3)*DEG2RAD)
     VAK3 = -VA(K3)*SIN(XC(K3)*DEG2RAD)+UA(K3)*COS(XC(K3)*DEG2RAD)
     
     COFA1=A1U_XY(IA,1)*UAIA+A1U_XY(IA,2)*UAK1   &
          +A1U_XY(IA,3)*UAK2+A1U_XY(IA,4)*UAK3
     COFA2=A2U_XY(IA,1)*UAIA+A2U_XY(IA,2)*UAK1   &
          +A2U_XY(IA,3)*UAK2+A2U_XY(IA,4)*UAK3
     COFA5=A1U_XY(IA,1)*VAIA+A1U_XY(IA,2)*VAK1   &
          +A1U_XY(IA,3)*VAK2+A1U_XY(IA,4)*VAK3
     COFA6=A2U_XY(IA,1)*VAIA+A2U_XY(IA,2)*VAK1   &
          +A2U_XY(IA,3)*VAK2+A2U_XY(IA,4)*VAK3
     
     XIJC_TMP = REARTH * COS(YIJC(I)*DEG2RAD) * COS(XIJC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*DEG2RAD))
     YIJC_TMP = REARTH * COS(YIJC(I)*DEG2RAD) * SIN(XIJC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*DEG2RAD))
     XCIA_TMP = REARTH * COS(YC(IA)*DEG2RAD) * COS(XC(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IA)*DEG2RAD))
     YCIA_TMP = REARTH * COS(YC(IA)*DEG2RAD) * SIN(XC(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IA)*DEG2RAD))
		  
     XIJ_TMP = XIJC_TMP-XCIA_TMP
     YIJ_TMP = YIJC_TMP-YCIA_TMP

     UIJ1=UAIA+COFA1*XIJ_TMP+COFA2*YIJ_TMP
     VIJ1=VAIA+COFA5*XIJ_TMP+COFA6*XIJ_TMP

!    FLUX FROM RIGHT
     K1=NBE(IB,1)
     K2=NBE(IB,2)
     K3=NBE(IB,3)
          
     UAIB = -VA(IB)*COS(XC(IB)*DEG2RAD)-UA(IB)*SIN(XC(IB)*DEG2RAD)
     VAIB = -VA(IB)*SIN(XC(IB)*DEG2RAD)+UA(IB)*COS(XC(IB)*DEG2RAD)
     UAK1 = -VA(K1)*COS(XC(K1)*DEG2RAD)-UA(K1)*SIN(XC(K1)*DEG2RAD)
     VAK1 = -VA(K1)*SIN(XC(K1)*DEG2RAD)+UA(K1)*COS(XC(K1)*DEG2RAD)
     UAK2 = -VA(K2)*COS(XC(K2)*DEG2RAD)-UA(K2)*SIN(XC(K2)*DEG2RAD)
     VAK2 = -VA(K2)*SIN(XC(K2)*DEG2RAD)+UA(K2)*COS(XC(K2)*DEG2RAD)
     UAK3 = -VA(K3)*COS(XC(K3)*DEG2RAD)-UA(K3)*SIN(XC(K3)*DEG2RAD)
     VAK3 = -VA(K3)*SIN(XC(K3)*DEG2RAD)+UA(K3)*COS(XC(K3)*DEG2RAD)
     
     COFA3=A1U_XY(IB,1)*UAIB+A1U_XY(IB,2)*UAK1   &
          +A1U_XY(IB,3)*UAK2+A1U_XY(IB,4)*UAK3
     COFA4=A2U_XY(IB,1)*UAIB+A2U_XY(IB,2)*UAK1   &
          +A2U_XY(IB,3)*UAK2+A2U_XY(IB,4)*UAK3
     COFA7=A1U_XY(IB,1)*VAIB+A1U_XY(IB,2)*VAK1   &
          +A1U_XY(IB,3)*VAK2+A1U_XY(IB,4)*VAK3
     COFA8=A2U_XY(IB,1)*VAIB+A2U_XY(IB,2)*VAK1   &
          +A2U_XY(IB,3)*VAK2+A2U_XY(IB,4)*VAK3
     
     XCIB_TMP = REARTH * COS(YC(IB)*DEG2RAD) * COS(XC(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IB)*DEG2RAD))
     YCIB_TMP = REARTH * COS(YC(IB)*DEG2RAD) * SIN(XC(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IB)*DEG2RAD))
		  
     XIJ_TMP = XIJC_TMP-XCIB_TMP
     YIJ_TMP = YIJC_TMP-YCIB_TMP
     UIJ2=UAIB+COFA3*XIJ_TMP+COFA4*YIJ_TMP
     VIJ2=VAIB+COFA7*XIJ_TMP+COFA8*YIJ_TMP

!    VISCOSITY COEFFICIENT
     VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
     VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

     !VISCOF=HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
     ! David moved HPRNU and added VHC
     VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU
     

     VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))
     VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))

     VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))
     VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD) &
               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))

     DLTXC_TMP = VX2_TMP-VX1_TMP
     DLTYC_TMP = VY2_TMP-VY1_TMP

!    SHEAR STRESSES         
     TXXIJ=(COFA1+COFA3)*VISCOF
     TYYIJ=(COFA6+COFA8)*VISCOF
     TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
     FXX=DIJ*(TXXIJ*DLTYC_TMP-TXYIJ*DLTXC_TMP)
     FYY=DIJ*(TXYIJ*DLTYC_TMP-TYYIJ*DLTXC_TMP)

     UIJ1_TMP = UIJ1
     VIJ1_TMP = VIJ1
     UIJ2_TMP = UIJ2
     VIJ2_TMP = VIJ2

     UIJ_TMP=0.5_SP*(UIJ1+UIJ2)
     VIJ_TMP=0.5_SP*(VIJ1+VIJ2)
     UN_TMP=-UIJ_TMP*DLTYC_TMP + VIJ_TMP*DLTXC_TMP

!    ADD CONVECTIVE AND VISCOUS FLUXES

!    ACCUMULATE FLUX
     XADV_TMP=DIJ*UN_TMP*&
              ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2_TMP                    &
              +(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1_TMP)*0.5_SP  
     YADV_TMP=DIJ*UN_TMP* &
              ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2_TMP                    &
	      +(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1_TMP)*0.5_SP  

     IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
       XFLUX(IA)=XFLUX(IA)+(XADV_TMP+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
       YFLUX(IA)=YFLUX(IA)+(YADV_TMP+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
       XFLUX(IB)=XFLUX(IB)-(XADV_TMP+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
       YFLUX(IB)=YFLUX(IB)-(YADV_TMP+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
       XFLUX(IA)=XFLUX(IA)+(XADV_TMP+FXX*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
       YFLUX(IA)=YFLUX(IA)+(YADV_TMP+FYY*EPOR(IA))*(1.0_SP-ISBC(I))*IUCP(IA)
     ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN
       XFLUX(IB)=XFLUX(IB)-(XADV_TMP+FXX*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
       YFLUX(IB)=YFLUX(IB)-(YADV_TMP+FYY*EPOR(IB))*(1.0_SP-ISBC(I))*IUCP(IB)
     END IF 


!    ACCUMULATE BAROTROPIC FLUX

!     VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))
!     VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))

!     VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))
!     VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD) &
!               * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))

!     DLTXC_TMP = VX2_TMP-VX1_TMP
!     DLTYC_TMP = VY2_TMP-VY1_TMP


     IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
      PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
      PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
      PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
      PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
     ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
      PSTX(IA)=PSTX(IA)-GRAV_E(IA)*D1(IA)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
      PSTY(IA)=PSTY(IA)+GRAV_E(IA)*D1(IA)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
     ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
      PSTX(IB)=PSTX(IB)+GRAV_E(IB)*D1(IB)*ELIJ*DLTYC_TMP/(2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
      PSTY(IB)=PSTY(IB)-GRAV_E(IB)*D1(IB)*ELIJ*DLTXC_TMP/(2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
     END IF

   END DO

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: advave_edge_XY"

   RETURN
   END SUBROUTINE ADVAVE_EDGE_XY
!==============================================================================|

              	 

!==============================================================================|
   SUBROUTINE ADVECTION_EDGE_XY(XFLUX,YFLUX)

   IMPLICIT NONE
!   REAL(SP), INTENT(OUT), DIMENSION(0:NT,KB) :: XFLUX,YFLUX
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: DIJ
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

   REAL(SP) :: UIJ1_TMP,VIJ1_TMP,UIJ2_TMP,VIJ2_TMP,TXXIJ_TMP,TYYIJ_TMP
   REAL(SP) :: XADV_TMP,YADV_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,UN_TMP
   REAL(SP) :: DLTXC_TMP,DLTYC_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP
   REAL(SP) :: UIA,VIA,UIB,VIB,UK1,VK1,UK2,VK2,UK3,VK3,UK4,VK4,UK5,VK5,UK6,VK6
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP

   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
!------------------------------------------------------------------------------|
   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: advection_edge_XY"

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
   DO K = 1,KBM1
     DO II=1,NPE
       I=NPEDGE_LST(II)
       IA=IEC(I,1)
       IB=IEC(I,2)
       IF(CELL_NORTHAREA(IA) == 1)THEN
         XFLUX(IA,K) = 0.0_SP
         YFLUX(IA,K) = 0.0_SP
       END IF
       IF(CELL_NORTHAREA(IB) == 1)THEN	 
         XFLUX(IB,K) = 0.0_SP
         YFLUX(IB,K) = 0.0_SP
       END IF	 
     END DO
   END DO    

!
!--Loop Over Edges and Accumulate Fluxes-For Each Element----------------------|
!


   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
!     DIJ= 0.5_SP*(DT(J1)+DT(J2))

     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)

     XIJC_TMP = REARTH * COS(YIJC(I)*DEG2RAD) * COS(XIJC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*DEG2RAD))
     YIJC_TMP = REARTH * COS(YIJC(I)*DEG2RAD) * SIN(XIJC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*DEG2RAD))
     XCIA_TMP = REARTH * COS(YC(IA)*DEG2RAD) * COS(XC(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IA)*DEG2RAD))
     YCIA_TMP = REARTH * COS(YC(IA)*DEG2RAD) * SIN(XC(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IA)*DEG2RAD))
     XCIB_TMP = REARTH * COS(YC(IB)*DEG2RAD) * COS(XC(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IB)*DEG2RAD))
     YCIB_TMP = REARTH * COS(YC(IB)*DEG2RAD) * SIN(XC(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IB)*DEG2RAD))
		  
     XIJA = XIJC_TMP-XCIA_TMP
     YIJA = YIJC_TMP-YCIA_TMP
     XIJB = XIJC_TMP-XCIB_TMP
     YIJB = YIJC_TMP-YCIB_TMP

     DO K=1,KBM1

       DIJ= 0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))

       UIA = -V(IA,K)*COS(XC(IA)*DEG2RAD)-U(IA,K)*SIN(XC(IA)*DEG2RAD)
       VIA = -V(IA,K)*SIN(XC(IA)*DEG2RAD)+U(IA,K)*COS(XC(IA)*DEG2RAD)
       UIB = -V(IB,K)*COS(XC(IB)*DEG2RAD)-U(IB,K)*SIN(XC(IB)*DEG2RAD)
       VIB = -V(IB,K)*SIN(XC(IB)*DEG2RAD)+U(IB,K)*COS(XC(IB)*DEG2RAD)
       UK1 = -V(K1,K)*COS(XC(K1)*DEG2RAD)-U(K1,K)*SIN(XC(K1)*DEG2RAD)
       VK1 = -V(K1,K)*SIN(XC(K1)*DEG2RAD)+U(K1,K)*COS(XC(K1)*DEG2RAD)
       UK2 = -V(K2,K)*COS(XC(K2)*DEG2RAD)-U(K2,K)*SIN(XC(K2)*DEG2RAD)
       VK2 = -V(K2,K)*SIN(XC(K2)*DEG2RAD)+U(K2,K)*COS(XC(K2)*DEG2RAD)
       UK3 = -V(K3,K)*COS(XC(K3)*DEG2RAD)-U(K3,K)*SIN(XC(K3)*DEG2RAD)
       VK3 = -V(K3,K)*SIN(XC(K3)*DEG2RAD)+U(K3,K)*COS(XC(K3)*DEG2RAD)
       UK4 = -V(K4,K)*COS(XC(K4)*DEG2RAD)-U(K4,K)*SIN(XC(K4)*DEG2RAD)
       VK4 = -V(K4,K)*SIN(XC(K4)*DEG2RAD)+U(K4,K)*COS(XC(K4)*DEG2RAD)
       UK5 = -V(K5,K)*COS(XC(K5)*DEG2RAD)-U(K5,K)*SIN(XC(K5)*DEG2RAD)
       VK5 = -V(K5,K)*SIN(XC(K5)*DEG2RAD)+U(K5,K)*COS(XC(K5)*DEG2RAD)
       UK6 = -V(K6,K)*COS(XC(K6)*DEG2RAD)-U(K6,K)*SIN(XC(K6)*DEG2RAD)
       VK6 = -V(K6,K)*SIN(XC(K6)*DEG2RAD)+U(K6,K)*COS(XC(K6)*DEG2RAD)
       !!FORM THE LEFT FLUX
       COFA1=A1U_XY(IA,1)*UIA+A1U_XY(IA,2)*UK1   &
            +A1U_XY(IA,3)*UK2+A1U_XY(IA,4)*UK3
       COFA2=A2U_XY(IA,1)*UIA+A2U_XY(IA,2)*UK1   &
            +A2U_XY(IA,3)*UK2+A2U_XY(IA,4)*UK3
       COFA5=A1U_XY(IA,1)*VIA+A1U_XY(IA,2)*VK1   &
            +A1U_XY(IA,3)*VK2+A1U_XY(IA,4)*VK3
       COFA6=A2U_XY(IA,1)*VIA+A2U_XY(IA,2)*VK1   &
            +A2U_XY(IA,3)*VK2+A2U_XY(IA,4)*VK3
       UIJ1=UIA+COFA1*XIJA+COFA2*YIJA
       VIJ1=VIA+COFA5*XIJA+COFA6*YIJA

       !!FORM THE RIGHT FLUX
       COFA3=A1U_XY(IB,1)*UIB+A1U_XY(IB,2)*UK4   &
            +A1U_XY(IB,3)*UK5+A1U_XY(IB,4)*UK6
       COFA4=A2U_XY(IB,1)*UIB+A2U_XY(IB,2)*UK4   &
            +A2U_XY(IB,3)*UK5+A2U_XY(IB,4)*UK6
       COFA7=A1U_XY(IB,1)*VIB+A1U_XY(IB,2)*VK4   &
            +A1U_XY(IB,3)*VK5+A1U_XY(IB,4)*VK6
       COFA8=A2U_XY(IB,1)*VIB+A2U_XY(IB,2)*VK4   &
            +A2U_XY(IB,3)*VK5+A2U_XY(IB,4)*VK6
       UIJ2=UIB+COFA3*XIJB+COFA4*YIJB
       VIJ2=VIB+COFA7*XIJB+COFA8*YIJB

       !!    VISCOSITY COEFFICIENT EDGE
       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)
        
!       VISCOF = HORCON*(FACT*0.5_SP*(VISCOF1+VISCOF2)/HPRNU + FM1)
       ! David moved HPRNU and added VHC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD)       &
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*DEG2RAD))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*DEG2RAD))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*DEG2RAD))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*DEG2RAD))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP
       
       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC_TMP-TXYIJ*DLTXC_TMP)
       FYY=DIJ*(TXYIJ*DLTYC_TMP-TYYIJ*DLTXC_TMP)

       UIJ1_TMP = UIJ1
       VIJ1_TMP = VIJ1
       UIJ2_TMP = UIJ2
       VIJ2_TMP = VIJ2
     
       UIJ_TMP=0.5_SP*(UIJ1+UIJ2)
       VIJ_TMP=0.5_SP*(VIJ1+VIJ2)
       UN_TMP=-UIJ_TMP*DLTYC_TMP + VIJ_TMP*DLTXC_TMP

       !!COMPUTE BOUNDARY FLUX AUGMENTERS   
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)


       !!ACCUMULATE THE FLUX
       XADV_TMP=DIJ*UN_TMP*&
                ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2_TMP                 &
	        +(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1_TMP)*0.5_SP
       YADV_TMP=DIJ*UN_TMP* &
                ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2_TMP                 &
	          +(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1_TMP)*0.5_SP 
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN 
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       END IF
     END DO
   END DO

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: advection_edge_XY"

   RETURN
   END SUBROUTINE ADVECTION_EDGE_XY
!==============================================================================|

!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_XY(XFLUX,YFLUX,CETA,STG,K_STG)

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ,UN
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: SITA,DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2

!   REAL(SP) :: TY
   REAL(SP) :: UIJ1_TMP,VIJ1_TMP,UIJ2_TMP,VIJ2_TMP,UIJ3_TMP,VIJ3_TMP
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP
   REAL(SP) :: TXXIJ_TMP,TYYIJ_TMP
   REAL(SP) :: XADV_TMP,YADV_TMP,PSTX_TMP,PSTY_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP
   REAL(SP) :: DLTXC_TMP,DLTYC_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP
   REAL(SP) :: UIA,VIA,UIB,VIB,UK1,VK1,UK2,VK2,UK3,VK3,UK4,VK4,UK5,VK5,UK6,VK6 
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP

   REAL(SP) :: CETA

   REAL(SP) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY

   INTEGER :: STG, K_STG
!------------------------------------------------------------------------------!
  ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
 
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: adv_uv_edge_xy"

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
!-----Initialize Flux Variables------------------------------------------------!
!

   DO K = 1,KBM1
     DO II=1,NPE
       I=NPEDGE_LST(II)
       IA=IEC(I,1)
       IB=IEC(I,2)
       IF(CELL_NORTHAREA(IA) == 1)THEN
         XFLUX(IA,K) = 0.0_SP
         YFLUX(IA,K) = 0.0_SP
         PSTX_TM(IA,K) = 0.0_SP
         PSTY_TM(IA,K) = 0.0_SP
       END IF
       IF(CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IB,K) = 0.0_SP
         YFLUX(IB,K) = 0.0_SP
         PSTX_TM(IB,K) = 0.0_SP
         PSTY_TM(IB,K) = 0.0_SP
       END IF 
     END DO
   END DO         

!
!-----Loop Over Edges and Accumulate Flux--------------------------------------!
!

   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
!     DIJ=0.5_SP*(DT(J1)+DT(J2))


     ELIJ=0.5_SP*(EGF(J1)+EGF(J2))
!   ggao 0103-2007   consider the sea surface pressure change
!     change end



     K1=NBE(IA,1)
     K2=NBE(IA,2)
     K3=NBE(IA,3)
     K4=NBE(IB,1)
     K5=NBE(IB,2)
     K6=NBE(IB,3)

     XIJC_TMP = REARTH * COS(YIJC(I)*DEG2RAD) * COS(XIJC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*DEG2RAD))
     YIJC_TMP = REARTH * COS(YIJC(I)*DEG2RAD) * SIN(XIJC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*DEG2RAD))
		  
     XCIA_TMP = REARTH * COS(YC(IA)*DEG2RAD) * COS(XC(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IA)*DEG2RAD))
     YCIA_TMP = REARTH * COS(YC(IA)*DEG2RAD) * SIN(XC(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IA)*DEG2RAD))
     XCIB_TMP = REARTH * COS(YC(IB)*DEG2RAD) * COS(XC(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IB)*DEG2RAD))
     YCIB_TMP = REARTH * COS(YC(IB)*DEG2RAD) * SIN(XC(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+SIN(YC(IB)*DEG2RAD))
		  
     XIJA = XIJC_TMP-XCIA_TMP
     YIJA = YIJC_TMP-YCIA_TMP
     XIJB = XIJC_TMP-XCIB_TMP
     YIJB = YIJC_TMP-YCIB_TMP

     DO K=1,KBM1

       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))

       UIA = -V(IA,K)*COS(XC(IA)*DEG2RAD)-U(IA,K)*SIN(XC(IA)*DEG2RAD)
       VIA = -V(IA,K)*SIN(XC(IA)*DEG2RAD)+U(IA,K)*COS(XC(IA)*DEG2RAD)
       UIB = -V(IB,K)*COS(XC(IB)*DEG2RAD)-U(IB,K)*SIN(XC(IB)*DEG2RAD)
       VIB = -V(IB,K)*SIN(XC(IB)*DEG2RAD)+U(IB,K)*COS(XC(IB)*DEG2RAD)
       UK1 = -V(K1,K)*COS(XC(K1)*DEG2RAD)-U(K1,K)*SIN(XC(K1)*DEG2RAD)
       VK1 = -V(K1,K)*SIN(XC(K1)*DEG2RAD)+U(K1,K)*COS(XC(K1)*DEG2RAD)
       UK2 = -V(K2,K)*COS(XC(K2)*DEG2RAD)-U(K2,K)*SIN(XC(K2)*DEG2RAD)
       VK2 = -V(K2,K)*SIN(XC(K2)*DEG2RAD)+U(K2,K)*COS(XC(K2)*DEG2RAD)
       UK3 = -V(K3,K)*COS(XC(K3)*DEG2RAD)-U(K3,K)*SIN(XC(K3)*DEG2RAD)
       VK3 = -V(K3,K)*SIN(XC(K3)*DEG2RAD)+U(K3,K)*COS(XC(K3)*DEG2RAD)
       UK4 = -V(K4,K)*COS(XC(K4)*DEG2RAD)-U(K4,K)*SIN(XC(K4)*DEG2RAD)
       VK4 = -V(K4,K)*SIN(XC(K4)*DEG2RAD)+U(K4,K)*COS(XC(K4)*DEG2RAD)
       UK5 = -V(K5,K)*COS(XC(K5)*DEG2RAD)-U(K5,K)*SIN(XC(K5)*DEG2RAD)
       VK5 = -V(K5,K)*SIN(XC(K5)*DEG2RAD)+U(K5,K)*COS(XC(K5)*DEG2RAD)
       UK6 = -V(K6,K)*COS(XC(K6)*DEG2RAD)-U(K6,K)*SIN(XC(K6)*DEG2RAD)
       VK6 = -V(K6,K)*SIN(XC(K6)*DEG2RAD)+U(K6,K)*COS(XC(K6)*DEG2RAD)

       COFA1=A1U_XY(IA,1)*UIA+A1U_XY(IA,2)*UK1   &
            +A1U_XY(IA,3)*UK2+A1U_XY(IA,4)*UK3
       COFA2=A2U_XY(IA,1)*UIA+A2U_XY(IA,2)*UK1   &
            +A2U_XY(IA,3)*UK2+A2U_XY(IA,4)*UK3
       COFA5=A1U_XY(IA,1)*VIA+A1U_XY(IA,2)*VK1   &
            +A1U_XY(IA,3)*VK2+A1U_XY(IA,4)*VK3
       COFA6=A2U_XY(IA,1)*VIA+A2U_XY(IA,2)*VK1   &
            +A2U_XY(IA,3)*VK2+A2U_XY(IA,4)*VK3

       UIJ1=UIA+COFA1*XIJA+COFA2*YIJA
       VIJ1=VIA+COFA5*XIJA+COFA6*YIJA

       COFA3=A1U_XY(IB,1)*UIB+A1U_XY(IB,2)*UK4   &
            +A1U_XY(IB,3)*UK5+A1U_XY(IB,4)*UK6
       COFA4=A2U_XY(IB,1)*UIB+A2U_XY(IB,2)*UK4   &
            +A2U_XY(IB,3)*UK5+A2U_XY(IB,4)*UK6
       COFA7=A1U_XY(IB,1)*VIB+A1U_XY(IB,2)*VK4   &
            +A1U_XY(IB,3)*VK5+A1U_XY(IB,4)*VK6
       COFA8=A2U_XY(IB,1)*VIB+A2U_XY(IB,2)*VK4   &
            +A2U_XY(IB,3)*VK5+A2U_XY(IB,4)*VK6

       UIJ2=UIB+COFA3*XIJB+COFA4*YIJB
       VIJ2=VIB+COFA7*XIJB+COFA8*YIJB

!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!
       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

       !VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
     
       ! David moved HPRNU and added VHC
       VISCOF=(FACT*0.5_SP*(VISCOF1*CC_HVC(IA)+VISCOF2*CC_HVC(IB)) + FM1*0.5_SP*(CC_HVC(IA)+CC_HVC(IB)))/HPRNU

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD)       &
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*DEG2RAD))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,1))*DEG2RAD))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*DEG2RAD))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD)&
                  * 2._SP /(1._SP+SIN(VY(IENODE(I,2))*DEG2RAD))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP
       
       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX=DIJ*(TXXIJ*DLTYC_TMP-TXYIJ*DLTXC_TMP)
       FYY=DIJ*(TXYIJ*DLTYC_TMP-TYYIJ*DLTXC_TMP)


       UIJ1_TMP = UIJ1
       VIJ1_TMP = VIJ1
       UIJ2_TMP = UIJ2
       VIJ2_TMP = VIJ2

       UIJ_TMP=0.5_SP*(UIJ1+UIJ2)
       VIJ_TMP=0.5_SP*(VIJ1+VIJ2)
       EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYC_TMP + VIJ_TMP*DLTXC_TMP)

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
       XADV_TMP=EXFLUX_TMP*&
                ((1.0_SP-SIGN(1.0_SP,EXFLUX_TMP))*UIJ2_TMP                 &
                +(1.0_SP+SIGN(1.0_SP,EXFLUX_TMP))*UIJ1_TMP)*0.5_SP
       YADV_TMP=EXFLUX_TMP* &
                ((1.0_SP-SIGN(1.0_SP,EXFLUX_TMP))*VIJ2_TMP                 &
 	        +(1.0_SP+SIGN(1.0_SP,EXFLUX_TMP))*VIJ1_TMP)*0.5_SP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         XFLUX(IA,K)=XFLUX(IA,K)+XADV_TMP*TPA+(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IA)
         YFLUX(IA,K)=YFLUX(IA,K)+YADV_TMP*TPA+(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN 
         XFLUX(IB,K)=XFLUX(IB,K)-XADV_TMP*TPB-(FXX+3.0_SP*FXX*FLOAT(ISBC(I)))*EPOR(IB)
         YFLUX(IB,K)=YFLUX(IB,K)-YADV_TMP*TPB-(FYY+3.0_SP*FYY*FLOAT(ISBC(I)))*EPOR(IB)
       END IF

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * COS(VX(IENODE(I,1))*DEG2RAD) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*DEG2RAD) * SIN(VX(IENODE(I,1))*DEG2RAD) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,1))*DEG2RAD))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * COS(VX(IENODE(I,2))*DEG2RAD) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*DEG2RAD) * SIN(VX(IENODE(I,2))*DEG2RAD) &
                 * 2._SP /(1._SP+sin(VY(IENODE(I,2))*DEG2RAD))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP

       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/  &
	              (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/  &
	              (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
	              (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
	              (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC_TMP/   &
	              (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
         PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC_TMP/   &
	              (2._SP /(1._SP+sin(YC(IA)*DEG2RAD)))
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN       
         PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC_TMP/  &
	              (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
         PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC_TMP/  &
	              (2._SP /(1._SP+sin(YC(IB)*DEG2RAD)))
       END IF
     END DO
   END DO


   DO K = 1,KBM1
     DO II=1,NP
       I=NP_LST(II)
       XFLUX(I,K)=XFLUX(I,K)+PSTX_TM(I,K)
       YFLUX(I,K)=YFLUX(I,K)+PSTY_TM(I,K)
     END DO
   END DO

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO II=1,NP
     I=NP_LST(II)
     IF(CELL_NORTHAREA(I) == 1)THEN
       DO K=1,KBM1
         IF(K == 1) THEN
           UIJ1_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
           VIJ1_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)
           UIJ2_TMP = -V(I,K+1)*COS(XC(I)*DEG2RAD)-U(I,K+1)*SIN(XC(I)*DEG2RAD)
           VIJ2_TMP = -V(I,K+1)*SIN(XC(I)*DEG2RAD)+U(I,K+1)*COS(XC(I)*DEG2RAD)
           XFLUXV=-W(I,K+1)*(UIJ1_TMP*DZ1(I,K+1)+UIJ2_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
           YFLUXV=-W(I,K+1)*(VIJ1_TMP*DZ1(I,K+1)+VIJ2_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
         ELSE IF(K == KBM1) THEN
           UIJ1_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
           VIJ1_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)
           UIJ2_TMP = -V(I,K-1)*COS(XC(I)*DEG2RAD)-U(I,K-1)*SIN(XC(I)*DEG2RAD)
           VIJ2_TMP = -V(I,K-1)*SIN(XC(I)*DEG2RAD)+U(I,K-1)*COS(XC(I)*DEG2RAD)
           XFLUXV= W(I,K)*(UIJ1_TMP*DZ1(I,K-1)+UIJ2_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))
           YFLUXV= W(I,K)*(VIJ1_TMP*DZ1(I,K-1)+VIJ2_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))
         ELSE
           UIJ1_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
           VIJ1_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)
           UIJ2_TMP = -V(I,K-1)*COS(XC(I)*DEG2RAD)-U(I,K-1)*SIN(XC(I)*DEG2RAD)
           VIJ2_TMP = -V(I,K-1)*SIN(XC(I)*DEG2RAD)+U(I,K-1)*COS(XC(I)*DEG2RAD)
           UIJ3_TMP = -V(I,K+1)*COS(XC(I)*DEG2RAD)-U(I,K+1)*SIN(XC(I)*DEG2RAD)
           VIJ3_TMP = -V(I,K+1)*SIN(XC(I)*DEG2RAD)+U(I,K+1)*COS(XC(I)*DEG2RAD)
           XFLUXV= W(I,K)*(UIJ1_TMP*DZ1(I,K-1)+UIJ2_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))-  &
                   W(I,K+1)*(UIJ1_TMP*DZ1(I,K+1)+UIJ3_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
           YFLUXV= W(I,K)*(VIJ1_TMP*DZ1(I,K-1)+VIJ2_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))-  &
                   W(I,K+1)*(VIJ1_TMP*DZ1(I,K+1)+VIJ3_TMP*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K+1))
         END IF
         U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
         V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)
!         XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
!                    +DRHOX(I,K)-COR(I)*V_TMP*DT1(I)*ART(I)
!         YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
!                    +DRHOY(I,K)+COR(I)*U_TMP*DT1(I)*ART(I)
         XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                    +DRHOX(I,K)-COR(I)*V_TMP*DT1(I)*DZ1(I,K)*ART(I)
         YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                    +DRHOY(I,K)+COR(I)*U_TMP*DT1(I)*DZ1(I,K)*ART(I)

       END DO
     END IF
   END DO


   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: adv_uv_edge_xy"

   RETURN
   END SUBROUTINE ADV_UV_EDGE_XY
!==============================================================================!

!==============================================================================|
   SUBROUTINE EXTEL_EDGE_XY(K,XFLUX)       

   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:MT)
   REAL(SP) :: DIJ,UIJ,VIJ
   INTEGER  :: I,J,K,I1,IA,IB,JJ,J1,J2,II

   REAL(SP) :: UIJ_TMP,VIJ_TMP,EXFLUX_TMP
   REAL(SP) :: DLTXE_TMP,DLTYE_TMP
   REAL(SP) :: VX1_TMP,VX2_TMP,VY1_TMP,VY2_TMP

   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: extel_edge_XY"
!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA) = 0.0_SP
     END IF
     IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB) = 0.0_SP
     END IF  
   END DO  

!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     DIJ = D1(I1)

     UIJ = UA(I1)
     VIJ = VA(I1)

     IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
       UIJ_TMP = -VIJ*COS(XC(I1)*DEG2RAD)-UIJ*SIN(XC(I1)*DEG2RAD)
       VIJ_TMP = -VIJ*SIN(XC(I1)*DEG2RAD)+UIJ*COS(XC(I1)*DEG2RAD)
       
       VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

       DLTXE_TMP = VX2_TMP-VX1_TMP
       DLTYE_TMP = VY2_TMP-VY1_TMP
       
       EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
     END IF  
     
     IF(IA == NODE_NORTHPOLE) THEN    
       XFLUX(IA) = XFLUX(IA)-EXFLUX_TMP
     ELSE IF(IB == NODE_NORTHPOLE)THEN
       XFLUX(IB) = XFLUX(IB)+EXFLUX_TMP
     END IF    

   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: extel_edge_XY"

   RETURN
   END SUBROUTINE EXTEL_EDGE_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE EXTUV_EDGE_XY(K)       

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESX,RESY
   INTEGER  :: I,II

   REAL(SP) :: UARK_TMP,VARK_TMP,UAF_TMP,VAF_TMP,UA_TMP,VA_TMP

   REAL(SP) :: WUSURF2_TMP,WVSURF2_TMP,WUBOT_TMP,WVBOT_TMP
   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: extuv_edge_XY"

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   DO II=1,NP
     I=NP_LST(II)
     UA_TMP = -VA(I)*COS(XC(I)*DEG2RAD)-UA(I)*SIN(XC(I)*DEG2RAD)
     VA_TMP = -VA(I)*SIN(XC(I)*DEG2RAD)+UA(I)*COS(XC(I)*DEG2RAD)

     WUSURF2_TMP = -WVSURF2(I)*COS(XC(I)*DEG2RAD)-WUSURF2(I)*SIN(XC(I)*DEG2RAD)
     WVSURF2_TMP = -WVSURF2(I)*SIN(XC(I)*DEG2RAD)+WUSURF2(I)*COS(XC(I)*DEG2RAD)

     WUBOT_TMP = -WVBOT(I)*COS(XC(I)*DEG2RAD)-WUBOT(I)*SIN(XC(I)*DEG2RAD)
     WVBOT_TMP = -WVBOT(I)*SIN(XC(I)*DEG2RAD)+WUBOT(I)*COS(XC(I)*DEG2RAD)

     RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)                &
               -COR(I)*VA_TMP*D1(I)*ART(I)                       &
               -(WUSURF2_TMP+WUBOT_TMP)*ART(I)
     RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)                &
               +COR(I)*UA_TMP*D1(I)*ART(I)                       &
               -(WVSURF2_TMP+WVBOT_TMP)*ART(I)

!
!--UPDATE----------------------------------------------------------------------|
!
     UARK_TMP = -VARK(I)*COS(XC(I)*DEG2RAD)-UARK(I)*SIN(XC(I)*DEG2RAD)
     VARK_TMP = -VARK(I)*SIN(XC(I)*DEG2RAD)+UARK(I)*COS(XC(I)*DEG2RAD)

     UAF_TMP = (UARK_TMP*(H1(I)+ELRK1(I))              &
               -ALPHA_RK(K)*DTE*RESX(I)/ART(I))/(H1(I)+ELF1(I))
     VAF_TMP = (VARK_TMP*(H1(I)+ELRK1(I))              &
               -ALPHA_RK(K)*DTE*RESY(I)/ART(I))/(H1(I)+ELF1(I))
		       
     UAF(I)  = VAF_TMP*COS(XC(I)*DEG2RAD)-UAF_TMP*SIN(XC(I)*DEG2RAD)
     VAF(I)  = UAF_TMP*COS(XC(I)*DEG2RAD)+VAF_TMP*SIN(XC(I)*DEG2RAD)
     VAF(I)  = -VAF(I)			    

   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: extuv_edge_XY"

   RETURN
   END SUBROUTINE EXTUV_EDGE_XY
!==============================================================================|

!==============================================================================|

   SUBROUTINE VERTVL_EDGE_XY(XFLUX,CETA)         

!------------------------------------------------------------------------------|
   IMPLICIT NONE 
   REAL(SP) :: XFLUX(MT,KBM1)
   REAL(SP) :: DIJ,UIJ,VIJ,UN,EXFLUX,TMP1,DIJ1,UIJ1,VIJ1
   INTEGER  :: I,K,IA,IB,I1 ,J,JJ,J1,J2,II

   REAL(SP) :: UIJ_TMP,VIJ_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,UIJ1_TMP,VIJ1_TMP
   REAL(SP) :: DLTXE_TMP,DLTYE_TMP,EXFLUX_TMP

   REAL(SP) :: CETA
!------------------------------------------------------------------------------|
   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: Vertvl_edge_XY"

!----------------------INITIALIZE FLUX-----------------------------------------!

   DO K=1,KBM1
     DO II=1,NPCV
       I = NCEDGE_LST(II)
       IA  = NIEC(I,1)
       IB  = NIEC(I,2)
       IF(IA == NODE_NORTHPOLE)THEN
         XFLUX(IA,K) = 0.0_SP
       END IF
       IF(IB == NODE_NORTHPOLE)THEN  
         XFLUX(IB,K) = 0.0_SP
       END IF  
     END DO  
   END DO
!----------------------ACCUMULATE FLUX-----------------------------------------!

   DO II=1,NPCV
     I=NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
!     DIJ=DT1(I1)
     DO K=1,KBM1
       DIJ=DT1(I1)*DZ1(I1,K)
       UIJ=U(I1,K)
       VIJ=V(I1,K)

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -VIJ*COS(XC(I1)*DEG2RAD)-UIJ*SIN(XC(I1)*DEG2RAD)
         VIJ_TMP = -VIJ*SIN(XC(I1)*DEG2RAD)+UIJ*COS(XC(I1)*DEG2RAD)

       VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)&
                 * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         EXFLUX_TMP = DIJ*(-UIJ_TMP*DLTYE_TMP+VIJ_TMP*DLTXE_TMP)
       END IF  
     
       IF(IA == NODE_NORTHPOLE)THEN
         XFLUX(IA,K) = XFLUX(IA,K)-EXFLUX_TMP
       ELSE IF(IB == NODE_NORTHPOLE)THEN
         XFLUX(IB,K) = XFLUX(IB,K)+EXFLUX_TMP
       END IF
     END DO
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: Vertvl_edge_xy"

   RETURN
   END SUBROUTINE VERTVL_EDGE_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE ADV_S_XY(XFLUX,XFLUX_ADV,PSPX,PSPY,PSPXD,PSPYD,VISCOFF,K,CETA)               

!------------------------------------------------------------------------------|

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K

   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(M)           :: PSPX,PSPY,PSPXD,PSPYD,VISCOFF
!   REAL(SP), DIMENSION(3*(NT))      :: DTIJ 
   REAL(SP), DIMENSION(3*(NT),KBM1)      :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2 
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF   
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PSPX_TMP,PSPY_TMP,PSPXD_TMP,PSPYD_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2

!!  ggao edge calculation
   REAL(SP) :: XIJE1_TMP,YIJE1_TMP,XIJE2_TMP,YIJE2_TMP
   REAL(SP) :: S1MIN, S1MAX, S2MIN, S2MAX

   REAL(SP) :: CETA
!------------------------------------------------------------------------------!
   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: ADV_S_XY(K):",K

   
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
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA,K) = 0.0_SP
       XFLUX_ADV(IA,K) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB,K) = 0.0_SP
       XFLUX_ADV(IB,K) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!
   I = NODE_NORTHPOLE

   IF(I==0)  RETURN

   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
       
     VX_TMP = REARTH * COS(VY(I)*DEG2RAD) * COS(VX(I)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I)*DEG2RAD))
     VY_TMP = REARTH * COS(VY(I)*DEG2RAD) * SIN(VX(I)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I)*DEG2RAD))
		     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*DEG2RAD) * COS(VX(NV(I1,J1))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*DEG2RAD))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*DEG2RAD) * SIN(VX(NV(I1,J1))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*DEG2RAD))
		     
     VX2_TMP= REARTH * COS(YC(I1)*DEG2RAD) * COS(XC(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(YC(I1)*DEG2RAD))
     VY2_TMP= REARTH * COS(YC(I1)*DEG2RAD) * SIN(XC(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(YC(I1)*DEG2RAD))
		     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*DEG2RAD) * COS(VX(NV(I1,J2))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*DEG2RAD))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*DEG2RAD) * SIN(VX(NV(I1,J2))*DEG2RAD) &
                    * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*DEG2RAD))
		     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)
     
     U_TMP = -V(I1,K)*COS(XC(I1)*DEG2RAD)-U(I1,K)*SIN(XC(I1)*DEG2RAD)
     V_TMP = -V(I1,K)*SIN(XC(I1)*DEG2RAD)+U(I1,K)*COS(XC(I1)*DEG2RAD)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2
   VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

   IF(K == KBM1) THEN
     AH_BOTTOM(I) = (FACT*VISCOFF(I) + FM1)*NN_HVC(I)
   END IF

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     
     IF((IA <= M .AND. IB <= M) .AND. I1 <= N)THEN
!       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!  ggao edge calculation
       XIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       YIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       XIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       YIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       XI_TMP =0.5_SP*(XIJE1_TMP+XIJE2_TMP)
       YI_TMP =0.5_SP*(YIJE1_TMP+YIJE2_TMP)


       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
!         XI_TMP = REARTH * COS(YI*DEG2RAD) * COS(XI*DEG2RAD) &
!                  * 2._SP /(1._SP+sin(YI*DEG2RAD))
!         YI_TMP = REARTH * COS(YI*DEG2RAD) * SIN(XI*DEG2RAD) &
!                  * 2._SP /(1._SP+sin(YI*DEG2RAD))

         VXA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
         VYA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))

         VXB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
         VYB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))

!         IF(IA == NODE_NORTHPOLE)THEN
         DXA=XI_TMP-VXA_TMP
         DYA=YI_TMP-VYA_TMP
!         ELSE IF(IB == NODE_NORTHPOLE)THEN
         DXB=XI_TMP-VXB_TMP
         DYB=YI_TMP-VYB_TMP
!	 END IF
!       END IF

        IF(IA == NODE_NORTHPOLE)THEN
	  PSPX_TMP=-PSPY(IB)*COS(VX(IB)*DEG2RAD)-PSPX(IB)*SIN(VX(IB)*DEG2RAD)
          PSPY_TMP=-PSPY(IB)*SIN(VX(IB)*DEG2RAD)+PSPX(IB)*COS(VX(IB)*DEG2RAD)
   
	  PSPXD_TMP=-PSPYD(IB)*COS(VX(IB)*DEG2RAD)-PSPXD(IB)*SIN(VX(IB)*DEG2RAD)
          PSPYD_TMP=-PSPYD(IB)*SIN(VX(IB)*DEG2RAD)+PSPXD(IB)*COS(VX(IB)*DEG2RAD)
   
          FIJ1=S1(IA,K)+DXA*PSPX(IA)+DYA*PSPY(IA)
          FIJ2=S1(IB,K)+DXB*PSPX_TMP+DYB*PSPY_TMP

          !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB))) !/HPRNU

          TXX=0.5_SP*(PSPXD(IA)+PSPXD_TMP)*VISCOF
          TYY=0.5_SP*(PSPYD(IA)+PSPYD_TMP)*VISCOF
        ELSE IF(IB == NODE_NORTHPOLE)THEN
	  PSPX_TMP=-PSPY(IA)*COS(VX(IA)*DEG2RAD)-PSPX(IA)*SIN(VX(IA)*DEG2RAD)
          PSPY_TMP=-PSPY(IA)*SIN(VX(IA)*DEG2RAD)+PSPX(IA)*COS(VX(IA)*DEG2RAD)
   
	  PSPXD_TMP=-PSPYD(IA)*COS(VX(IA)*DEG2RAD)-PSPXD(IA)*SIN(VX(IA)*DEG2RAD)
          PSPYD_TMP=-PSPYD(IA)*SIN(VX(IA)*DEG2RAD)+PSPXD(IA)*COS(VX(IA)*DEG2RAD)
   
          FIJ1=S1(IA,K)+DXA*PSPX_TMP+DYA*PSPY_TMP
          FIJ2=S1(IB,K)+DXB*PSPX(IB)+DYB*PSPY(IB)

          !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

          TXX=0.5_SP*(PSPXD_TMP+PSPXD(IB))*VISCOF
          TYY=0.5_SP*(PSPYD_TMP+PSPYD(IB))*VISCOF
        END IF

       S1MIN=MINVAL(S1(NBSN(IA,1:NTSN(IA)-1),K))
       S1MIN=MIN(S1MIN, S1(IA,K))
       S1MAX=MAXVAL(S1(NBSN(IA,1:NTSN(IA)-1),K))
       S1MAX=MAX(S1MAX, S1(IA,K))
       S2MIN=MINVAL(S1(NBSN(IB,1:NTSN(IB)-1),K))
       S2MIN=MIN(S2MIN, S1(IB,K))
       S2MAX=MAXVAL(S1(NBSN(IB,1:NTSN(IB)-1),K))
       S2MAX=MAX(S2MAX, S1(IB,K))
       IF(FIJ1 < S1MIN) FIJ1=S1MIN
       IF(FIJ1 > S1MAX) FIJ1=S1MAX
       IF(FIJ2 < S2MIN) FIJ2=S2MIN
       IF(FIJ2 > S2MAX) FIJ2=S2MAX

!       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -V(I1,K)*COS(XC(I1)*DEG2RAD)-U(I1,K)*SIN(XC(I1)*DEG2RAD)
         VIJ_TMP = -V(I1,K)*SIN(XC(I1)*DEG2RAD)+U(I1,K)*COS(XC(I1)*DEG2RAD)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)
         VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)

         VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)
         VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         FXX=-DTIJ(I,K)*TXX*DLTYE_TMP
         FYY= DTIJ(I,K)*TYY*DLTXE_TMP

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP  
         EXFLUX_TMP = -UVN_TMP*DTIJ(I,K)*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP

         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX_TMP+FXX+FYY
           XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+EXFLUX_TMP
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX_TMP-FXX-FYY
           XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-EXFLUX_TMP
         END IF
       END IF
     END IF  
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: ADV_S_XY(K):",K

   RETURN
   END SUBROUTINE ADV_S_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE ADV_T_XY(XFLUX,XFLUX_ADV,PTPX,PTPY,PTPXD,PTPYD,VISCOFF,K,CETA)               

!------------------------------------------------------------------------------|

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(M)           :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT),KBM1)      :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PTPX_TMP,PTPY_TMP,PTPXD_TMP,PTPYD_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2

!!  ggao edge calculation
   REAL(SP) :: XIJE1_TMP,YIJE1_TMP,XIJE2_TMP,YIJE2_TMP
   REAL(SP) :: T1MIN, T1MAX, T2MIN, T2MAX

   REAL(SP) :: CETA
!------------------------------------------------------------------------------!

   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: ADV_T_XY(K):",K

   
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
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA,K) = 0.0_SP
       XFLUX_ADV(IA,K) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB,K) = 0.0_SP
       XFLUX_ADV(IB,K) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!
   I = NODE_NORTHPOLE

   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
       
     VX_TMP = REARTH * COS(VY(I)*DEG2RAD) * COS(VX(I)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I)*DEG2RAD))
     VY_TMP = REARTH * COS(VY(I)*DEG2RAD) * SIN(VX(I)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I)*DEG2RAD))
		     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*DEG2RAD) * COS(VX(NV(I1,J1))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*DEG2RAD))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*DEG2RAD) * SIN(VX(NV(I1,J1))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*DEG2RAD))
		     
     VX2_TMP= REARTH * COS(YC(I1)*DEG2RAD) * COS(XC(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(YC(I1)*DEG2RAD))
     VY2_TMP= REARTH * COS(YC(I1)*DEG2RAD) * SIN(XC(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(YC(I1)*DEG2RAD))
		     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*DEG2RAD) * COS(VX(NV(I1,J2))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*DEG2RAD))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*DEG2RAD) * SIN(VX(NV(I1,J2))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*DEG2RAD))
		     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)
     
     U_TMP = -V(I1,K)*COS(XC(I1)*DEG2RAD)-U(I1,K)*SIN(XC(I1)*DEG2RAD)
     V_TMP = -V(I1,K)*SIN(XC(I1)*DEG2RAD)+U(I1,K)*COS(XC(I1)*DEG2RAD)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2
   VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

   IF(K == KBM1) THEN
     AH_BOTTOM(I) = (FACT*VISCOFF(I) + FM1)*NN_HVC(I)
   END IF

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     IF(IA <= M .AND. IB <= M .AND. I1 <= N)THEN
!       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!  ggao edge calculation
       XIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       YIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       XIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       YIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       XI_TMP =0.5_SP*(XIJE1_TMP+XIJE2_TMP)
       YI_TMP =0.5_SP*(YIJE1_TMP+YIJE2_TMP)

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
!         XI_TMP = REARTH * COS(YI*DEG2RAD) * COS(XI*DEG2RAD) &
!                  * 2._SP /(1._SP+sin(YI*DEG2RAD))
!         YI_TMP = REARTH * COS(YI*DEG2RAD) * SIN(XI*DEG2RAD) &
!                  * 2._SP /(1._SP+sin(YI*DEG2RAD))

         VXA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
         VYA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))

         VXB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
         VYB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))

!        IF(IA == NODE_NORTHPOLE)THEN
           DXA=XI_TMP-VXA_TMP
           DYA=YI_TMP-VYA_TMP
!	 ELSE IF(IB == NODE_NORTHPOLE)THEN
           DXB=XI_TMP-VXB_TMP
           DYB=YI_TMP-VYB_TMP
!	 END IF
!       END IF

        IF(IA == NODE_NORTHPOLE)THEN
	  PTPX_TMP=-PTPY(IB)*COS(VX(IB)*DEG2RAD)-PTPX(IB)*SIN(VX(IB)*DEG2RAD)
          PTPY_TMP=-PTPY(IB)*SIN(VX(IB)*DEG2RAD)+PTPX(IB)*COS(VX(IB)*DEG2RAD)
   
	  PTPXD_TMP=-PTPYD(IB)*COS(VX(IB)*DEG2RAD)-PTPXD(IB)*SIN(VX(IB)*DEG2RAD)
          PTPYD_TMP=-PTPYD(IB)*SIN(VX(IB)*DEG2RAD)+PTPXD(IB)*COS(VX(IB)*DEG2RAD)
   
          FIJ1=T1(IA,K)+DXA*PTPX(IA)+DYA*PTPY(IA)
          FIJ2=T1(IB,K)+DXB*PTPX_TMP+DYB*PTPY_TMP

          !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB))) !/HPRNU

          TXX=0.5_SP*(PTPXD(IA)+PTPXD_TMP)*VISCOF
          TYY=0.5_SP*(PTPYD(IA)+PTPYD_TMP)*VISCOF
        ELSE IF(IB == NODE_NORTHPOLE)THEN
	  PTPX_TMP=-PTPY(IA)*COS(VX(IA)*DEG2RAD)-PTPX(IA)*SIN(VX(IA)*DEG2RAD)
          PTPY_TMP=-PTPY(IA)*SIN(VX(IA)*DEG2RAD)+PTPX(IA)*COS(VX(IA)*DEG2RAD)
   
	  PTPXD_TMP=-PTPYD(IA)*COS(VX(IA)*DEG2RAD)-PTPXD(IA)*SIN(VX(IA)*DEG2RAD)
          PTPYD_TMP=-PTPYD(IA)*SIN(VX(IA)*DEG2RAD)+PTPXD(IA)*COS(VX(IA)*DEG2RAD)
   
          FIJ1=T1(IA,K)+DXA*PTPX_TMP+DYA*PTPY_TMP
          FIJ2=T1(IB,K)+DXB*PTPX(IB)+DYB*PTPY(IB)

          !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

          TXX=0.5_SP*(PTPXD_TMP+PTPXD(IB))*VISCOF
          TYY=0.5_SP*(PTPYD_TMP+PTPYD(IB))*VISCOF
        END IF

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

!     IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -V(I1,K)*COS(XC(I1)*DEG2RAD)-U(I1,K)*SIN(XC(I1)*DEG2RAD)
         VIJ_TMP = -V(I1,K)*SIN(XC(I1)*DEG2RAD)+U(I1,K)*COS(XC(I1)*DEG2RAD)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)
         VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)

         VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)
         VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         FXX=-DTIJ(I,K)*TXX*DLTYE_TMP
         FYY= DTIJ(I,K)*TYY*DLTXE_TMP

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*DTIJ(I,K)*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
       
         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX_TMP+FXX+FYY
           XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+EXFLUX_TMP
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX_TMP-FXX-FYY
           XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-EXFLUX_TMP
         END IF
       END IF
     END IF  
   END DO 

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: ADV_T_XY(K):",K

   RETURN
   END SUBROUTINE ADV_T_XY
!==============================================================================|

   SUBROUTINE ADV_N_XY(XFLUX,PWPX,PWPY,ISS,ID,DEP2,SPCDIR,N32)
!  This subroutine is for wave advection at the north pole    yzhang
!------------------------------------------------------------------------------|
   USE SWCOMM3, ONLY :MDC,MSC

   IMPLICIT NONE

   SAVE

   REAL :: N32(MDC,MSC,0:MT),N32_TMP(MDC,MSC,0:MT)
   INTEGER  :: ID,ISS,IG,IG2,IDT,IDD ! LWU IG2 IS FOR PLBC
   REAL :: CANX,CANY,CANX_TMP,CANY_TMP,ADDEXFLUX2455
   REAL :: SPCDIR(MDC,6)
   REAL(SP) :: DEP2(MT)
   REAL    :: DEPLOC,KWAVELOC,CGLOC,NN,ND,SPCSIGL
   REAL(SP), DIMENSION(0:MT)     :: XFLUX
   REAL(SP), DIMENSION(M)           :: PWPX,PWPY
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,JTMP,JJ,II,L
   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UL_DEGREE,DL_DEGREE,CENTER_DEGREE,FF11
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PWPX_TMP,PWPY_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: ADDYIJE1,ADDYIJE2
   REAL(SP) ::  ADD_DLTXE ,ADD_DLTYE,ADD_DLTXTRIE,ADD_DLTYTRIE
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
   REAL(SP) :: XIJE1_TMP,YIJE1_TMP,XIJE2_TMP,YIJE2_TMP
   REAL(SP) :: AC1MIN, AC1MAX, AC2MIN, AC2MAX

!------------------------------------------------------------------------------!
   IF (NODE_NORTHPOLE .EQ. 0) RETURN

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: ADV_N_XY"

!   CAX = 0.0
!   CAY = 0.0

!
!--Initialize Fluxes-----------------------------------------------------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN
       XFLUX(IB) = 0.0_SP
     END IF
   END DO
!==========YZHANG_NORTHPOLE_TESTING================================
   I=NODE_NORTHPOLE
   PWPX(I)=0.0_SP
   PWPY(I)=0.0_SP
   ADD_DLTXTRIE=0.0_SP
   ADD_DLTYTRIE=0.0_SP
   DO J=1,NTSN(I)-1
     I1=NBSN(I,J)
     I2=NBSN(I,J+1)
     UL_DEGREE=22.5
     CENTER_DEGREE=0
     DL_DEGREE=337.5

     IF (VX(I1)<UL_DEGREE .OR. VX(I1)>DL_DEGREE)THEN
       IDT=MDC-(MDC*2/8)
       IDD=ID+IDT
!       IDD=ID
       IF(IDD>MDC)THEN
         IDD=IDD-MDC
       END IF
       N32_TMP(ID,ISS,I1)=N32(IDD,ISS,I1)
     END IF

     IF (VX(I2)<UL_DEGREE .OR. VX(I2)>DL_DEGREE)THEN
       IDT=MDC-(MDC*2/8)
       IDD=ID+IDT
!       IDD=ID
       IF(IDD>MDC)THEN
         IDD=IDD-MDC
       END IF
       N32_TMP(ID,ISS,I2)=N32(IDD,ISS,I2)
     END IF
     DO L=1,7
       CENTER_DEGREE=L*45.0
       UL_DEGREE=CENTER_DEGREE+22.5
       DL_DEGREE=CENTER_DEGREE-22.5

       IF(VX(I1)<UL_DEGREE .AND. VX(I1)>DL_DEGREE) THEN
         IDT=MDC-(MDC*(L+2)/8)
!         IDT=MDC-(MDC*L/8)
         IF(IDT<0)THEN
           IDT=IDT+MDC
         END IF
         IDD=ID+IDT
         IF(IDD>MDC)THEN
           IDD=IDD-MDC
         END IF
         N32_TMP(ID,ISS,I1)=N32(IDD,ISS,I1)
       END IF
       IF(VX(I2)<UL_DEGREE .AND. VX(I2)>DL_DEGREE) THEN
         IDT=MDC-(MDC*(L+2)/8)
!         IDT=MDC-(MDC*L/8)
         IF(IDT<0)THEN
           IDT=IDT+MDC
         END IF
         IDD=ID+IDT
         IF(IDD>MDC)THEN
           IDD=IDD-MDC
         END IF
         N32_TMP(ID,ISS,I2)=N32(IDD,ISS,I2)
       END IF
     END DO
     FF11=0.5*(N32_TMP(ID,ISS,I1)+N32_TMP(ID,ISS,I2))
     PWPX(I)=PWPX(I)+FF11*DLTYTRIE(I,J)
     PWPY(I)=PWPY(I)+FF11*DLTXTRIE(I,J)
     ADD_DLTXTRIE=ADD_DLTXTRIE+DLTXTRIE(I,J)
     ADD_DLTYTRIE=ADD_DLTYTRIE+DLTYTRIE(I,J)
   END DO

   PWPX(I)=PWPX(I)/ART2(I)
   PWPY(I)=PWPY(I)/ART2(I)
!============================================================================
   ADDYIJE1=0.0_SP
   ADDYIJE2=0.0_SP
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     ADDYIJE1=ADDYIJE1+YIJE(I,1)
     ADDYIJE2=ADDYIJE2+YIJE(I,2)
   END DO
   ADDYIJE1=ADDYIJE1/NPCV
   ADDYIJE2=ADDYIJE2/NPCV

   ADD_DLTXE=0.0_SP
   ADD_DLTYE=0.0_SP

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     IF(IA <= M .AND. IB <= M .AND. I1 <= N)THEN
!        XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!        YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!   ggao edge calculation
       XIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       YIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       XIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       YIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       XI_TMP =0.5_SP*(XIJE1_TMP+XIJE2_TMP)
       YI_TMP =0.5_SP*(YIJE1_TMP+YIJE2_TMP)

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         VXA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
         VYA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))

         VXB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
         VYB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))

!         IF(IA == NODE_NORTHPOLE)THEN
         DXA=XI_TMP-VXA_TMP
         DYA=YI_TMP-VYA_TMP
!         ELSE IF(IB == NODE_NORTHPOLE)THEN
         DXB=XI_TMP-VXB_TMP
         DYB=YI_TMP-VYB_TMP

         IF(IA == NODE_NORTHPOLE)THEN
           PWPX_TMP=-PWPY(IB)*COS(VX(IB)*DEG2RAD)-PWPX(IB)*SIN(VX(IB)*DEG2RAD)
           PWPY_TMP=-PWPY(IB)*SIN(VX(IB)*DEG2RAD)+PWPX(IB)*COS(VX(IB)*DEG2RAD)

!           PSPXD_TMP=-PSPYD(IB)*COS(VX(IB)*DEG2RAD)-PSPXD(IB)*SIN(VX(IB)*DEG2RAD)
!           PSPYD_TMP=-PSPYD(IB)*SIN(VX(IB)*DEG2RAD)+PSPXD(IB)*COS(VX(IB)*DEG2RAD)
!======================zhangyang======start=======================
           IF (VX(IB)<UL_DEGREE .OR. VX(IB)>DL_DEGREE)THEN
             IDT=MDC-(MDC*2/8)
             IDD=ID+IDT
!             IDD=ID
             IF(IDD>MDC)THEN
               IDD=IDD-MDC
             END IF
             N32_TMP(ID,ISS,IB)=N32(IDD,ISS,IB)
           END IF

           DO L=1,7
             CENTER_DEGREE=L*45.0
             UL_DEGREE=CENTER_DEGREE+22.5
             DL_DEGREE=CENTER_DEGREE-22.5
             IF(VX(IB)<UL_DEGREE .AND. VX(IB)>DL_DEGREE) THEN
               IDT=MDC-(MDC*(L+2)/8)
!               IDT=MDC-(MDC*L/8)
               IF(IDT<0)THEN
                 IDT=IDT+MDC
               END IF
               IDD=ID+IDT
               IF(IDD>MDC)THEN
                 IDD=IDD-MDC
               END IF
               N32_TMP(ID,ISS,IB)=N32(IDD,ISS,IB)
             END IF
           END DO
           FIJ1=N32(ID,ISS,IA)!+DXA*PWPX(IA)+DYA*PWPY(IA)
           FIJ2=N32_TMP(ID,ISS,IB)!+DXB*PWPX_TMP+DYB*PWPY_TMP
!================zhangyang======end==================================

          !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
!          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) +
!          FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB))) !/HPRNU

!          TXX=0.5_SP*(PSPXD(IA)+PSPXD_TMP)*VISCOF
!          TYY=0.5_SP*(PSPYD(IA)+PSPYD_TMP)*VISCOF

         ELSE IF(IB == NODE_NORTHPOLE)THEN
           PWPX_TMP=-PWPY(IA)*COS(VX(IA)*DEG2RAD)-PWPX(IA)*SIN(VX(IA)*DEG2RAD)
           PWPY_TMP=-PWPY(IA)*SIN(VX(IA)*DEG2RAD)+PWPX(IA)*COS(VX(IA)*DEG2RAD)

!           PSPXD_TMP=-PSPYD(IA)*COS(VX(IA)*DEG2RAD)-PSPXD(IA)*SIN(VX(IA)*DEG2RAD)
!           PSPYD_TMP=-PSPYD(IA)*SIN(VX(IA)*DEG2RAD)+PSPXD(IA)*COS(VX(IA)*DEG2RAD)

!======================zhangyang======start=======================
           IF (VX(IA)<UL_DEGREE .OR. VX(IA)>DL_DEGREE)THEN
             IDT=MDC-(MDC*2/8)
             IDD=ID+IDT
!             IDD=ID
             IF(IDD>MDC)THEN
               IDD=IDD-MDC
             END IF
             N32_TMP(ID,ISS,IA)=N32(IDD,ISS,IA)
           END IF

           DO L=1,7
             CENTER_DEGREE=L*45.0
             UL_DEGREE=CENTER_DEGREE+22.5
             DL_DEGREE=CENTER_DEGREE-22.5
             IF(VX(IA)<UL_DEGREE .AND. VX(IA)>DL_DEGREE) THEN
               IDT=MDC-(MDC*(L+2)/8)
!               IDT=MDC-(MDC*L/8)
               IF(IDT<0)THEN
                 IDT=IDT+MDC
               END IF
               IDD=ID+IDT
               IF(IDD>MDC)THEN
                 IDD=IDD-MDC
               END IF
               N32_TMP(ID,ISS,IA)=N32(IDD,ISS,IA)
             END IF
           END DO
           FIJ1=N32_TMP(ID,ISS,IA)!+DXA*PWPX_TMP+DYA*PWPY_TMP
           FIJ2=N32(ID,ISS,IB)!+DXB*PWPX(IB)+DYB*PWPY(IB)

!================zhangyang======end==================================
         END IF
         CALL SWAPAR1(I1,ISS,ID,DEP2(1),KWAVELOC,CGLOC)
         UIJ_TMP = UA(I1)
         VIJ_TMP = VA(I1)

         DO L=1,8
           CENTER_DEGREE=L*45.0-22.5
           UL_DEGREE=CENTER_DEGREE+22.5
           DL_DEGREE=CENTER_DEGREE-22.5

           IF(XC(I1)<UL_DEGREE .AND. XC(I1)>DL_DEGREE) THEN
             IDT=MDC-(MDC*(L+2)/8-5)
!             IDT=MDC-(MDC*L/8)
             IF(IDT<0)THEN
               IDT=IDT+MDC
             END IF
             IDD=ID+IDT
             IF(IDD>MDC)THEN
               IDD=IDD-MDC
             END IF
           END IF
         END DO

         CALL SPROXY(I1,ISS,IDD,CANX,CANY,CGLOC,SPCDIR(IDD,2),SPCDIR(IDD,3),UIJ_TMP,VIJ_TMP)

         VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
         VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
 
         VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
         VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP

         CANX_TMP = -CANY*COS(XC(I1)*DEG2RAD)-CANX*SIN(XC(I1)*DEG2RAD)
         CANY_TMP = -CANY*SIN(XC(I1)*DEG2RAD)+CANX*COS(XC(I1)*DEG2RAD)

         UVN_TMP = CANY_TMP*DLTXE_TMP - CANX_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP

         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA)=XFLUX(IA)+EXFLUX_TMP
!           XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+EXFLUX_TMP
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB)=XFLUX(IB)-EXFLUX_TMP
!           XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-EXFLUX_TMP
         END IF
       END IF
     END IF
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: ADV_N_XY(ID,ISS):",ID,ISS

   RETURN
   END SUBROUTINE ADV_N_XY

!==============================================================================|
!     CALCULATE THE BAROCLINIC PRESSURE GRADIENT IN SIGMA COORDINATES          |
!==============================================================================|

   SUBROUTINE BAROPG_XY(DRIJK1,DRIJK2) 

!==============================================================================|
   IMPLICIT NONE
   REAL(SP) :: DRIJK1(0:N,3,KBM1), DRIJK2(0:N,KBM1)
   REAL(SP) :: DIJ,DRHO1,DRHO2
   INTEGER  :: I,II,K,J,J1,J2,IJK
   REAL(SP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
!==============================================================================|

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: baropg_XY"


   DO II = 1, NP
     I=NP_LST(II)
     DO K=1,KBM1
       DRHOX(I,K)=0.0_SP
       DRHOY(I,K)=0.0_SP
     END DO
   END DO

   DO II = 1, NP
     I=NP_LST(II)
     DO K=1,KBM1
        DO J = 1, 3
          J1=J+1-INT((J+1)/4)*3
          J2=J+2-INT((J+2)/4)*3
          IJK=NBE(I,J)
          DIJ=0.5_SP*(DT(NV(I,J1))+DT(NV(I,J2)))

          VY1_TMP=REARTH*COS(VY(NV(I,J1))*DEG2RAD)*SIN(VX(NV(I,J1))*DEG2RAD)
          VY2_TMP=REARTH*COS(VY(NV(I,J2))*DEG2RAD)*SIN(VX(NV(I,J2))*DEG2RAD)

          DRHO1=(VY1_TMP-VY2_TMP)*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VY1_TMP-VY2_TMP)*DIJ*DRIJK2(I,K)

          DRHOX(I,K)=DRHOX(I,K)+DRHO1+DRHO2

          VX1_TMP=REARTH*COS(VY(NV(I,J1))*DEG2RAD)*COS(VX(NV(I,J1))*DEG2RAD)
          VX2_TMP=REARTH*COS(VY(NV(I,J2))*DEG2RAD)*COS(VX(NV(I,J2))*DEG2RAD)

	  DRHO1=(VX2_TMP-VX1_TMP)*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VX2_TMP-VX1_TMP)*DIJ*DRIJK2(I,K)

          DRHOY(I,K)=DRHOY(I,K)+DRHO1+DRHO2

       END DO
     END DO
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: baropg_XY"
   
   RETURN
   END SUBROUTINE BAROPG_XY
!==============================================================================|

!==============================================================================|
   SUBROUTINE SHAPE_COEF_XY

!----------------------------------------------------------------------!
!  This subrountine is used to calculate the coefficient for a linear  !
!  function on the x-y plane, i.e.:                                    !
!                     r(x,y;phai)=phai_c+cofa1*x+cofa2*y               !
!     innc(i)=0    cells on the boundary                               !
!     innc(i)=1    cells in the interior                               !
!----------------------------------------------------------------------!
     
   USE ALL_VARS
   IMPLICIT NONE
   REAL(DP) X1,X2,X3,Y1,Y2,Y3,DELT,AI1,AI2,AI3,BI1,BI2,BI3,CI1,CI2,CI3
   REAL(DP) DELTX,DELTY,TEMP1,ANG1,ANG2,B1,B2,ANGLE
   REAL(DP), ALLOCATABLE :: XC_TMP(:),YC_TMP(:),VX_TMP(:),VY_TMP(:)
   INTEGER  I,II,J,JJ,J1,J2
!
!---------------interior cells-----------------------------------------!
!
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: shape_coef_xy"

   ! ALL OTHER PROCESSORS ESCAPE HERE
   IF (NODE_NORTHPOLE .EQ. 0) RETURN


   ALLOCATE(A1U_XY(N,4)); A1U_XY = 0.0_SP
   ALLOCATE(A2U_XY(N,4)); A2U_XY = 0.0_SP
   ALLOCATE(AW0_XY(N,3)); AW0_XY = 0.0_SP
   ALLOCATE(AWX_XY(N,3)); AWX_XY = 0.0_SP
   ALLOCATE(AWY_XY(N,3)); AWY_XY = 0.0_SP
   
   ALLOCATE(XC_TMP(0:NT)); XC_TMP = 0.0_SP
   ALLOCATE(YC_TMP(0:NT)); YC_TMP = 0.0_SP
   ALLOCATE(VX_TMP(0:MT)); VX_TMP = 0.0_SP
   ALLOCATE(VY_TMP(0:MT)); VY_TMP = 0.0_SP
   
   DO I=1,NT
     XC_TMP(I) = REARTH * COS(YC(I)*DEG2RAD) * COS(XC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YC(I)*DEG2RAD))
     YC_TMP(I) = REARTH * COS(YC(I)*DEG2RAD) * SIN(XC(I)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YC(I)*DEG2RAD))
   END DO		  

   DO I=1,MT
     VX_TMP(I) = REARTH * COS(VY(I)*DEG2RAD) * COS(VX(I)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(I)*DEG2RAD))
     VY_TMP(I) = REARTH * COS(VY(I)*DEG2RAD) * SIN(VX(I)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(VY(I)*DEG2RAD))
   END DO		  

   DO I=1,N
     IF(ISBCE(I) == 0)THEN
       Y1 = YC_TMP(NBE(I,1))-YC_TMP(I)
       Y2 = YC_TMP(NBE(I,2))-YC_TMP(I)
       Y3 = YC_TMP(NBE(I,3))-YC_TMP(I)
       X1=XC_TMP(NBE(I,1))-XC_TMP(I)
       X2=XC_TMP(NBE(I,2))-XC_TMP(I)
       X3=XC_TMP(NBE(I,3))-XC_TMP(I)

       X1=X1/1000.0_SP
       X2=X2/1000.0_SP
       X3=X3/1000.0_SP
       Y1=Y1/1000.0_SP
       Y2=Y2/1000.0_SP
       Y3=Y3/1000.0_SP

       delt=(x1*y2-x2*y1)**2+(x1*y3-x3*y1)**2+(x2*y3-x3*y2)**2
       delt=delt*1000.0_SP

       a1u_XY(i,1)=(y1+y2+y3)*(x1*y1+x2*y2+x3*y3)- &
                (x1+x2+x3)*(y1**2+y2**2+y3**2)
       a1u_XY(i,1)=a1u_XY(i,1)/delt
       a1u_XY(i,2)=(y1**2+y2**2+y3**2)*x1-(x1*y1+x2*y2+x3*y3)*y1
       a1u_XY(i,2)=a1u_XY(i,2)/delt
       a1u_XY(i,3)=(y1**2+y2**2+y3**2)*x2-(x1*y1+x2*y2+x3*y3)*y2
       a1u_XY(i,3)=a1u_XY(i,3)/delt
       a1u_XY(i,4)=(y1**2+y2**2+y3**2)*x3-(x1*y1+x2*y2+x3*y3)*y3
       a1u_XY(i,4)=a1u_XY(i,4)/delt

       a2u_XY(i,1)=(x1+x2+x3)*(x1*y1+x2*y2+x3*y3)- &
                (y1+y2+y3)*(x1**2+x2**2+x3**2)
       a2u_XY(i,1)=a2u_XY(i,1)/delt
       a2u_XY(i,2)=(x1**2+x2**2+x3**2)*y1-(x1*y1+x2*y2+x3*y3)*x1
       a2u_XY(i,2)=a2u_XY(i,2)/delt
       a2u_XY(i,3)=(x1**2+x2**2+x3**2)*y2-(x1*y1+x2*y2+x3*y3)*x2
       a2u_XY(i,3)=a2u_XY(i,3)/delt
       a2u_XY(i,4)=(x1**2+x2**2+x3**2)*y3-(x1*y1+x2*y2+x3*y3)*x3
       a2u_XY(i,4)=a2u_XY(i,4)/delt
     end if

     x1=vx_TMP(nv(i,1))-xc_TMP(i)
     x2=vx_TMP(nv(i,2))-xc_TMP(i)
     x3=vx_TMP(nv(i,3))-xc_TMP(i)
     y1=vy_TMP(nv(i,1))-yc_TMP(i)
     y2=vy_TMP(nv(i,2))-yc_TMP(i)
     y3=vy_TMP(nv(i,3))-yc_TMP(i)


     ai1=y2-y3
     ai2=y3-y1
     ai3=y1-y2
     bi1=x3-x2
     bi2=x1-x3
     bi3=x2-x1
     ci1=x2*y3-x3*y2
     ci2=x3*y1-x1*y3
     ci3=x1*y2-x2*y1

     aw0_XY(i,1)=-ci1/2./art(i)
     aw0_XY(i,2)=-ci2/2./art(i)
     aw0_XY(i,3)=-ci3/2./art(i)
     awx_XY(i,1)=-ai1/2./art(i)
     awx_XY(i,2)=-ai2/2./art(i)
     awx_XY(i,3)=-ai3/2./art(i)
     awy_XY(i,1)=-bi1/2./art(i)
     awy_XY(i,2)=-bi2/2./art(i)
     awy_XY(i,3)=-bi3/2./art(i)
   end do

   DEALLOCATE(XC_TMP,YC_TMP,VX_TMP,VY_TMP)
   
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: shape_coef_xy"

   return
   end subroutine shape_coef_xy

!==============================================================================|

!==============================================================================|
!==============================================================================|
   SUBROUTINE ADV_Q_XY(XFLUX,PQPX,PQPY,PQPXD,PQPYD,VISCOFF,Q,UQ,VQ,K,UQ1,VQ1,CETA)               

!==============================================================================|
!   Calculate the Turbulent Kinetic Energy and Mixing Length Based on          |
!   The Mellor-Yamada Level 2.5 Turbulent Closure Model                        |
!==============================================================================|


!------------------------------------------------------------------------------|

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:MT,KB)     :: XFLUX,Q
   REAL(SP), DIMENSION(M)           :: PQPX,PQPY,PQPXD,PQPYD,VISCOFF
   REAL(SP), DIMENSION(3*(NT),KBM1)      :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2 
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF   
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PQPX_TMP,PQPY_TMP,PQPXD_TMP,PQPYD_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2

   REAL(SP) :: XIJE1_TMP,YIJE1_TMP,XIJE2_TMP,YIJE2_TMP
   REAL(SP) :: Q1MIN, Q1MAX, Q2MIN, Q2MAX

   REAL(SP), DIMENSION(0:NT,KB)    :: UQ,VQ

   REAL(SP), DIMENSION(0:,:)    :: UQ1, VQ1
   REAL(SP) :: CETA
!------------------------------------------------------------------------------!
   IF (NODE_NORTHPOLE .EQ. 0) RETURN

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: adv_q_xy(K):",K


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
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA,K) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB,K) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     DTIJ(I,K)=DT1(I1)*DZ1(I1,K)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!
   I = NODE_NORTHPOLE

   IF(I==0)  RETURN

   IF(I /= 0)THEN
   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
      
     VX_TMP = REARTH * COS(VY(I)*DEG2RAD) * COS(VX(I)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I)*DEG2RAD))
     VY_TMP = REARTH * COS(VY(I)*DEG2RAD) * SIN(VX(I)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(I)*DEG2RAD))
		     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*DEG2RAD) * COS(VX(NV(I1,J1))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*DEG2RAD))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*DEG2RAD) * SIN(VX(NV(I1,J1))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*DEG2RAD))
		     
     VX2_TMP= REARTH * COS(YC(I1)*DEG2RAD) * COS(XC(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(YC(I1)*DEG2RAD))
     VY2_TMP= REARTH * COS(YC(I1)*DEG2RAD) * SIN(XC(I1)*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(YC(I1)*DEG2RAD))
		     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*DEG2RAD) * COS(VX(NV(I1,J2))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*DEG2RAD))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*DEG2RAD) * SIN(VX(NV(I1,J2))*DEG2RAD) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*DEG2RAD))
		     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)
     
     U_TMP = -VQ(I1,K)*COS(XC(I1)*DEG2RAD)-UQ(I1,K)*SIN(XC(I1)*DEG2RAD)
     V_TMP = -VQ(I1,K)*SIN(XC(I1)*DEG2RAD)+UQ(I1,K)*COS(XC(I1)*DEG2RAD)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2
   VISCOFF(I)=SQRT(TMP1+TMP2)*ART1(I)

   IF(K == KBM1) THEN
     AH_BOTTOM(I) = (FACT*VISCOFF(I) + FM1)*NN_HVC(I)
   END IF
   endif
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     
     IF((IA <= M .AND. IB <= M) .AND. I1 <= N)THEN
       XIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))
       YIJE1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,1)*DEG2RAD))

       XIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       YIJE2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD) &
                  * 2._SP /(1._SP+sin(YIJE(I,2)*DEG2RAD))
       XI_TMP =0.5_SP*(XIJE1_TMP+XIJE2_TMP)
       YI_TMP =0.5_SP*(YIJE1_TMP+YIJE2_TMP)

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         VXA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * COS(VX(IA)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))
         VYA_TMP = REARTH * COS(VY(IA)*DEG2RAD) * SIN(VX(IA)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IA)*DEG2RAD))

         VXB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * COS(VX(IB)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))
         VYB_TMP = REARTH * COS(VY(IB)*DEG2RAD) * SIN(VX(IB)*DEG2RAD) &
                   * 2._SP /(1._SP+sin(VY(IB)*DEG2RAD))

         DXA=XI_TMP-VXA_TMP
         DYA=YI_TMP-VYA_TMP
         DXB=XI_TMP-VXB_TMP
         DYB=YI_TMP-VYB_TMP

        IF(IA == NODE_NORTHPOLE)THEN
	  PQPX_TMP=-PQPY(IB)*COS(VX(IB)*DEG2RAD)-PQPX(IB)*SIN(VX(IB)*DEG2RAD)
          PQPY_TMP=-PQPY(IB)*SIN(VX(IB)*DEG2RAD)+PQPX(IB)*COS(VX(IB)*DEG2RAD)
   
	  PQPXD_TMP=-PQPYD(IB)*COS(VX(IB)*DEG2RAD)-PQPXD(IB)*SIN(VX(IB)*DEG2RAD)
          PQPYD_TMP=-PQPYD(IB)*SIN(VX(IB)*DEG2RAD)+PQPXD(IB)*COS(VX(IB)*DEG2RAD)
   
          FIJ1=Q(IA,K)+DXA*PQPX(IA)+DYA*PQPY(IA)
          FIJ2=Q(IB,K)+DXB*PQPX_TMP+DYB*PQPY_TMP

          ! VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))/HPRNU

          TXX=0.5_SP*(PQPXD(IA)+PQPXD_TMP)*VISCOF
          TYY=0.5_SP*(PQPYD(IA)+PQPYD_TMP)*VISCOF
        ELSE IF(IB == NODE_NORTHPOLE)THEN
	  PQPX_TMP=-PQPY(IA)*COS(VX(IA)*DEG2RAD)-PQPX(IA)*SIN(VX(IA)*DEG2RAD)
          PQPY_TMP=-PQPY(IA)*SIN(VX(IA)*DEG2RAD)+PQPX(IA)*COS(VX(IA)*DEG2RAD)
   
	  PQPXD_TMP=-PQPYD(IA)*COS(VX(IA)*DEG2RAD)-PQPXD(IA)*SIN(VX(IA)*DEG2RAD)
          PQPYD_TMP=-PQPYD(IA)*SIN(VX(IA)*DEG2RAD)+PQPXD(IA)*COS(VX(IA)*DEG2RAD)
   
          FIJ1=Q(IA,K)+DXA*PQPX_TMP+DYA*PQPY_TMP
          FIJ2=Q(IB,K)+DXB*PQPX(IB)+DYB*PQPY(IB)

          !VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
          ! David moved HPRNU and added VHC
          VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))/HPRNU

          TXX=0.5_SP*(PQPXD_TMP+PQPXD(IB))*VISCOF
          TYY=0.5_SP*(PQPYD_TMP+PQPYD(IB))*VISCOF
        END IF

       Q1MIN=MINVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MIN=MIN(Q1MIN, Q(IA,K))
       Q1MAX=MAXVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MAX=MAX(Q1MAX, Q(IA,K))
       Q2MIN=MINVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MIN=MIN(Q2MIN, Q(IB,K))
       Q2MAX=MAXVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MAX=MAX(Q2MAX, Q(IB,K))
       IF(FIJ1 < Q1MIN) FIJ1=Q1MIN
       IF(FIJ1 > Q1MAX) FIJ1=Q1MAX
       IF(FIJ2 < Q2MIN) FIJ2=Q2MIN
       IF(FIJ2 > Q2MAX) FIJ2=Q2MAX

!        FXX=-DTIJ(I,K)*TXX*DLTYE(I)
!        FYY= DTIJ(I,K)*TYY*DLTXE(I)

!       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -VQ(I1,K)*COS(XC(I1)*DEG2RAD)-UQ(I1,K)*SIN(XC(I1)*DEG2RAD)
         VIJ_TMP = -VQ(I1,K)*SIN(XC(I1)*DEG2RAD)+UQ(I1,K)*COS(XC(I1)*DEG2RAD)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * COS(XIJE(I,1)*DEG2RAD)
         VY1_TMP = REARTH * COS(YIJE(I,1)*DEG2RAD) * SIN(XIJE(I,1)*DEG2RAD)

         VX2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * COS(XIJE(I,2)*DEG2RAD)
         VY2_TMP = REARTH * COS(YIJE(I,2)*DEG2RAD) * SIN(XIJE(I,2)*DEG2RAD)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       
         FXX=-DTIJ(I,K)*TXX*DLTYE_TMP
         FYY= DTIJ(I,K)*TYY*DLTXE_TMP

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*DTIJ(I,K)*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
       
         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX_TMP+FXX+FYY
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX_TMP-FXX-FYY
         END IF
       END IF
     END IF  
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: ADV_Q_XY(K):",K

   RETURN
   END SUBROUTINE ADV_Q_XY
!==============================================================================|

END MODULE MOD_NORTHPOLE
