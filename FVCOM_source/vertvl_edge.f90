










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
!   CALCULATE THE SIGMA COORDINATE VERTICAL VELOCITY FOR THE 3D MODE (omega)   |
!							                       |
!   DETERMINED FROM EQUATION:						       |
!   									       !
!   d/dt(D) + d/dx(uD) + d/dy(uD) = d/sigma(omega)                             !
!==============================================================================|

   SUBROUTINE VERTVL_EDGE         

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_WD
   USE MOD_NORTHPOLE
 



   IMPLICIT NONE 
   REAL(SP) :: XFLUX(MT,KBM1),WBOTTOM(MT)
   REAL(SP) :: DIJ,UIJ,VIJ,UN,EXFLUX,TMP1,DIJ1,UIJ1,VIJ1
   INTEGER  :: I,K,IA,IB,I1,I2,I3,I4,J,JJ,J1,J2
!------------------------------------------------------------------------------|
   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: Vertvl_edge.F"
!----------------------INITIALIZE FLUX-----------------------------------------!

   XFLUX = 0.0_SP

!----------------------ACCUMULATE FLUX-----------------------------------------!

!!#  if !defined (1)
   DO I=1,NCV
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)

     DO K=1,KBM1
       DIJ=DT1(I1)*DZ1(I1,K)
       UIJ=U(I1,K)
       VIJ=V(I1,K)
       EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))


       XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX


     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!     IA=NIEC(I,1)
!!     IB=NIEC(I,2)

!!     DO K=1,KBM1
!!#      if !defined (SEMI_IMPLICIT)
!!       DIJ=DT1(I1)*DZ1(I1,K)
!!       UIJ=US(I1,K)
!!       VIJ=VS(I1,K)
!!       EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
!!#      else
!!       DIJ=DT1(I1)*DZ1(I1,K)
!!       DIJ1=D1(I1)*DZ1(I1,K)
!!       UIJ=US(I1,K)
!!       VIJ=VS(I1,K)
!!       UIJ1=UF(I1,K)
!!       VIJ1=VF(I1,K)
!!       EXFLUX=( (1.0_SP-IFCETA)*DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))+IFCETA*DIJ1*(-UIJ1*DLTYE(I)+VIJ1*DLTXE(I)) )*ISWETCT(I1)*ISWETC(I1)
!!#      endif
!!       XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
!!       XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
!!     END DO
!!   END DO
!!#  endif

   
!-----------------------NULLIFY BOUNDARY FLUX----------------------------------!
! For "tide + meanflow"/"meanflow only" case, this part should be commented out;
! For "tide only" case, this part may be kept.
! However, the effect of this term is small from my experience.

      DO I=1,M
        DO K=1,KBM1
          IF(ISONB(I) == 2) XFLUX(I,K)=0.0_SP  
        ENDDO
      ENDDO
! can be changed to (no IF statements)
!     DO I=1,IOBCN
!        DO K=1,KBM1
!           XFLUX(I_OBC_N(I),K)=0.0_SP
!        ENDDO
!     ENDDO



!-----------------------FRESH WATER INFLOW-------------------------------------!

   IF(NUMQBC >= 1) THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
       DO J=1,NUMQBC
         JJ=INODEQ(J)
         DO K=1,KBM1
           XFLUX(JJ,K)=XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)    !/DZ(JJ,K)
         END DO
       END DO
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
       DO J=1,NUMQBC
         J1=N_ICELLQ(J,1)
         J2=N_ICELLQ(J,2)
         DO K=1,KBM1
           XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)    !/DZ1(J1,K)
           XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)    !/DZ1(J2,K)
         END DO
       END DO
     END IF
   END IF


!---IF NO FRESH WATER INFLOW, OMEGA IS ZERO AT FREE SURFACE AND BOTTOM---------!

   !CLEAR OLD VALUES
   WBOTTOM = 0.0_SP
   WTS = 0.0_SP


! QXU changed sign of evap/precip
   WTS(:,1) = -(QEVAP+QPREC)*ROFVROS 


   ! SET BOTTOM VELOCITY
   WBOTTOM(1:M)=  BFWDIS(1:M)/ART1(1:M)

!--------------------------CALCULATE OMEGA-------------------------------------!

   DO I=1,M
    IF(ISWETNT(I)*ISWETN(I) == 1)THEN
     DO K=1,KBM1
!       WTS(I,K+1)=WTS(I,K)+DZ(I,K)*(XFLUX(I,K)/ART1(I)+(EL(I)-ET(I))/DTI)
       WTS(I,K+1)=WTS(I,K)+XFLUX(I,K)/ART1(I)+DZ(I,K)*(D(I)-DT(I))/DTI
     END DO
    ELSE
     DO K=1,KBM1
       WTS(I,K+1)=0.0_SP
     END DO
    END IF
   END DO


!-------------------------ADJUSTMENT FOR DAM MODULE----------------------------

!# if defined (THIN_DAM)
!   DO I=1,NODE_DAM1_N
!     I1 = I_NODE_DAM1_N(I,1)
!     I2 = I_NODE_DAM1_N(I,2)
!     if(i1==nlid(2944).or.i2==nlid(2944))print*,'orginal 2944:'&
!          &,kdam(nlid(2944)),wts(nlid(2944),3)
!     if(i1==nlid(4851).or.i2==nlid(4851))print*,'orginal 4851:'&
!          &,kdam(nlid(4851)),wts(nlid(4851),3)
!
!#  if defined (1)
!    IF(ISWETNT(I1)*ISWETN(I1) == 1)THEN
!#  endif
!     DO K=1,KDAM(I1)
!       WTS(I1,K+1)=WTS(I1,K)+(XFLUX(I1,K)+XFLUX(I2,K))/   &
!                   (ART1(I1)+ART1(I2))+DZ(I1,K)*(D(I1)-DT(I1))/DTI
!     END DO
!     DO K=KDAM(I1)+1,KBM1
!       WTS(I1,K+1)=WTS(I1,K)+XFLUX(I1,K)/ART1(I1)+DZ(I1,K)*(D(I1)-DT(I1))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I1,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!
!#  if defined (1)
!    IF(ISWETNT(I2)*ISWETN(I2) == 1)THEN
!#  endif
!     DO K=1,KDAM(I2)
!       WTS(I2,K+1)=WTS(I2,K)+(XFLUX(I1,K)+XFLUX(I2,K))/   &
!                   (ART1(I1)+ART1(I2))+DZ(I2,K)*(D(I2)-DT(I2))/DTI
!     END DO
!     DO K=KDAM(I2)+1,KBM1
!       WTS(I2,K+1)=WTS(I2,K)+XFLUX(I2,K)/ART1(I2)+DZ(I2,K)*(D(I2)-DT(I2))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I2,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!     if(i1==nlid(2944).or.i2==nlid(2944))print*,'adjusted 2944:'&
!          &,kdam(nlid(2944)),wts(nlid(2944),3)
!     if(i1==nlid(4851).or.i2==nlid(4851))print*,'adjusted 4851:'&
!          &,kdam(nlid(4851)),wts(nlid(4851),3)
!   END DO


!   DO I=1,NODE_DAM2_N
!     I1 = I_NODE_DAM2_N(I,1)
!     I2 = I_NODE_DAM2_N(I,2)
!     I3 = I_NODE_DAM2_N(I,3)
!#  if defined (1)
!    IF(ISWETNT(I1)*ISWETN(I1) == 1)THEN
!#  endif
!     DO K=1,KDAM(I1)
!       WTS(I1,K+1)=WTS(I1,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3))+DZ(I1,K)*(D(I1)-DT(I1))/DTI
!     END DO
!     DO K=KDAM(I1)+1,KBM1
!       WTS(I1,K+1)=WTS(I1,K)+XFLUX(I1,K)/ART1(I1)+DZ(I1,K)*(D(I1)-DT(I1))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I1,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!
!#  if defined (1)
!    IF(ISWETNT(I2)*ISWETN(I2) == 1)THEN
!#  endif
!     DO K=1,KDAM(I2)
!       WTS(I2,K+1)=WTS(I2,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3))+DZ(I2,K)*(D(I2)-DT(I2))/DTI
!     END DO
!     DO K=KDAM(I2)+1,KBM1
!       WTS(I2,K+1)=WTS(I2,K)+XFLUX(I2,K)/ART1(I2)+DZ(I2,K)*(D(I2)-DT(I2))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I2,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!
!#  if defined (1)
!    IF(ISWETNT(I3)*ISWETN(I3) == 1)THEN
!#  endif
!     DO K=1,KDAM(I3)
!       WTS(I3,K+1)=WTS(I3,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3))+DZ(I3,K)*(D(I3)-DT(I3))/DTI
!     END DO
!     DO K=KDAM(I3)+1,KBM1
!       WTS(I3,K+1)=WTS(I3,K)+XFLUX(I3,K)/ART1(I3)+DZ(I3,K)*(D(I3)-DT(I3))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I3,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!   END DO
!
!   DO I=1,NODE_DAM3_N
!     I1 = I_NODE_DAM3_N(I,1)
!     I2 = I_NODE_DAM3_N(I,2)
!     I3 = I_NODE_DAM3_N(I,3)
!     I4 = I_NODE_DAM3_N(I,4)
!#  if defined (1)
!    IF(ISWETNT(I1)*ISWETN(I1) == 1)THEN
!#  endif
!     DO K=1,KDAM(I1)
!       WTS(I1,K+1)=WTS(I1,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K)+XFLUX(I4,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))+        &
!		   DZ(I1,K)*(D(I1)-DT(I1))/DTI
!     END DO
!     DO K=KDAM(I1)+1,KBM1
!       WTS(I1,K+1)=WTS(I1,K)+XFLUX(I1,K)/ART1(I1)+DZ(I1,K)*(D(I1)-DT(I1))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I1,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!
!#  if defined (1)
!    IF(ISWETNT(I2)*ISWETN(I2) == 1)THEN
!#  endif
!     DO K=1,KDAM(I2)
!       WTS(I2,K+1)=WTS(I2,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K)+XFLUX(I4,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))+        &
!		   DZ(I2,K)*(D(I2)-DT(I2))/DTI
!     END DO
!     DO K=KDAM(I2)+1,KBM1
!       WTS(I2,K+1)=WTS(I2,K)+XFLUX(I2,K)/ART1(I2)+DZ(I2,K)*(D(I2)-DT(I2))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I2,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!
!#  if defined (1)
!    IF(ISWETNT(I3)*ISWETN(I3) == 1)THEN
!#  endif
!     DO K=1,KDAM(I3)
!       WTS(I3,K+1)=WTS(I3,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K)+XFLUX(I4,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))+        &
!		   DZ(I3,K)*(D(I3)-DT(I3))/DTI
!     END DO
!     DO K=KDAM(I3)+1,KBM1
!       WTS(I3,K+1)=WTS(I3,K)+XFLUX(I3,K)/ART1(I3)+DZ(I3,K)*(D(I3)-DT(I3))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I3,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!
!#  if defined (1)
!    IF(ISWETNT(I4)*ISWETN(I4) == 1)THEN
!#  endif
!     DO K=1,KDAM(I4)
!       WTS(I4,K+1)=WTS(I4,K)+(XFLUX(I1,K)+XFLUX(I2,K)+XFLUX(I3,K)+XFLUX(I4,K))/   &
!                   (ART1(I1)+ART1(I2)+ART1(I3)+ART1(I4))+        &
!		   DZ(I4,K)*(D(I4)-DT(I4))/DTI
!     END DO
!     DO K=KDAM(I4)+1,KBM1
!       WTS(I4,K+1)=WTS(I4,K)+XFLUX(I4,K)/ART1(I4)+DZ(I4,K)*(D(I4)-DT(I4))/DTI
!     END DO
!#  if defined (1)
!    ELSE
!     DO K=1,KBM1
!       WTS(I4,K+1)=0.0_SP
!     END DO
!    END IF
!#  endif
!   END DO
!
!#  endif

!--------------------------ADJUST OMEGA----------------------------------------!
! IMPROVES MASS CONSERVATION

   DO I=1,M
     IF(ABS(WTS(I,KB)-WBOTTOM(I)) > 1.0E-7_SP)THEN
       IF(ISONB(I) /= 2)THEN
         TMP1=ELF(I)*FLOAT(KBM1)-(WTS(I,KB)-WBOTTOM(I))*DTI/DZ(I,1)
         TMP1=TMP1/FLOAT(KBM1)


         DTFA(I)=TMP1+H(I)

         DO K=2,KB
           WTS(I,K)=WTS(I,K)-FLOAT(K-1)/FLOAT(KBM1)*(WTS(I,KB)-WBOTTOM(I))
         END DO
       END IF
     END IF
   END DO

   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,WTS)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,WTS)
!
!----TRANSFER OMEGA TO FACE CENTER---------------------------------------------!
!

   DO I=1,N
     DO K=1,KB
       W(I,K) = ONE_THIRD*(WTS(NV(I,1),K)+WTS(NV(I,2),K)+WTS(NV(I,3),K))
     END DO
   END DO

   DO I=1,N
     DO K=1,KB
       W(I,K) = FLOAT(ISWETC(I))*W(I,K)
     END DO
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: Vertvl_edge.F"

   RETURN
   END SUBROUTINE VERTVL_EDGE
!==============================================================================|
