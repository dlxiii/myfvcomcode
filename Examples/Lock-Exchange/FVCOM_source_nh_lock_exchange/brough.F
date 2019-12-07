! $Id: brough.F,v 1.4.2.6 2008/04/24 14:26:06 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.4.2.6 $
!==============================================================================|
!  Calculate Bottom Drag Coefficient based on Bottom Roughness                 !
!   note:                                                                      !
!   when the log function derived from the constant stress log-viscous         !
!   layer is applied to an estuary, if the value of z0 is close to             !
!   (zz(kbm1)-z(kb)*dt1, drag coefficient "cbc" could become a huge            !
!   number due to near-zero value of alog function. In our application         !
!   we simply cutoff at cbc=0.005. One could adjust this cutoff value          !
!   based on observations or his or her experiences.                           !   
!   CALCULATES:   WUBOT(N), WVBOT(N) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|

SUBROUTINE BOTTOM_ROUGHNESS

  !==============================================================================!
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_WD
  USE MOD_PAR

  IMPLICIT NONE
  INTEGER :: I,II
  REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
  REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
  REAL(SP)             :: ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM
  !==============================================================================!


  !
  !  SET CONSTANTS
  !

  SELECT CASE(BOTTOM_ROUGHNESS_TYPE) 
  !==============================================================================|
  CASE(BR_ORIG) !USE ORIGINAL FVCOM FORM FOR BOTTOM FRICTION  |
  !==============================================================================|

     ! SET A EFFECTIVE MAXIMUM FOR CBC USING THE DEPTH
     WHERE (DT1 > 3.0_SP)
        CBC = VK2/(LOG((ZZ1(:,KBM1)-Z1(:,KB))*DT1/CC_Z0B))**2
     ELSEWHERE
        CBC = VK2/(LOG((ZZ1(:,KBM1)-Z1(:,KB))*3.0/CC_Z0B))**2
     END WHERE
     
     ! SET A MINIMUM VALUE FOR CBC
     WHERE (CBC < CBCMIN)
        CBC=CBCMIN
     END WHERE


  !==============================================================================|
  CASE(BR_GOTM) !GOTM FORMULATION FOR BOTTOM FRICTION    |
  !==============================================================================|

     !----Convert Input Z0B to GOTMS H0B
     !     H0B = Z0B/.03  
     ! DAS fixed bug to match gotm's friction.f90
     DO I=1,N
        U_TAUB = 0.0_SP
        DO II=1,40       
           IF (UMOL <= 0.0_SP) THEN
              Z0B_GOTM=CC_Z0B(I)   !0.03*H0B 
           ELSE
              Z0B_GOTM=0.1_SP*UMOL/MAX(UMOL,U_TAUB)+CC_Z0B(I) !0.03*H0B
           END IF
           ztemp=(zz1(I,kbm1)-z1(I,kb))*dt1(i)
           RR=KAPPA/(LOG((Z0B_GOTM+ZTEMP)/Z0B_GOTM))
           U_TAUB = RR*SQRT( U(I,KBM1)*U(I,KBM1) + V(I,KBM1)*V(I,KBM1) )
        END DO
        CBC(I) =   RR*RR
     END DO


  CASE DEFAULT
     CALL FATAL_ERROR ("BROUGH: UNKNOWN BOTTOM_ROUGHNESS_TYPE:"&
          & ,TRIM(BOTTOM_ROUGHNESS_TYPE) )
  END SELECT

  !==============================================================================|
  !  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
  !==============================================================================|
  DO  I = 1, N
#      if !defined (TWO_D_MODEL)     
     BTPS= CBC(I)*SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
     WUBOT(I) = -BTPS * U(I,KBM1)
     WVBOT(I) = -BTPS * V(I,KBM1)
#      else
     BTPS= CBC(I)*SQRT(UA(I)**2+VA(I)**2)
     WUBOT(I) = -BTPS * UA(I)
     WVBOT(I) = -BTPS * VA(I)
#      endif       
  
  END DO


  RETURN
END SUBROUTINE BOTTOM_ROUGHNESS
!==============================================================================|
