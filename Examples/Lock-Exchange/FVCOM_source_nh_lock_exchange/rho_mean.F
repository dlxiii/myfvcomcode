!==============================================================================|
SUBROUTINE RHO_MEAN
  !==============================================================================|
  USE ALL_VARS
  USE MOD_UTILS
  IMPLICIT NONE

  INTEGER, PARAMETER :: KBB=601
  INTEGER, PARAMETER :: KBBM1=KBB-1

  REAL(SP) DELTZ
  REAL(SP), DIMENSION(KBBM1)   :: PHY_Z   !Depth(m) in every standary Z levels 
  REAL(SP), DIMENSION(KBBM1)   :: RHOZTMP !density in standary Z levels 
  !REAL(SP), DIMENSION(M,KBBM1) :: RHOZ    !density in standary Z levals,in Nodes 
  REAL(SP), ALLOCATABLE :: RHOZ(:,:)
  REAL(SP), DIMENSION(KBBM1)   :: RHOA    !density mean in standary Z levals

  REAL(SP), DIMENSION(KBM1)    :: ZM      !Depth (m) in every sigma levels for giving node
  REAL(SP), DIMENSION(KBM1)    :: RHOS    !density for giving node
  INTEGER ICOUNT,I,K,status

  !-- ALLOCATE LARGE ARRAY

  ALLOCATE(RHOZ(M,KBBM1),stat=status)
  IF(STATUS/=0) CALL FATAL_ERROR("RHO_MEAN: CAN NOT ALLOCATE RHOZ ARRAY!")

  !--CALCULATE Z-LEVELS TO MAX DEPTH---------------------------------------------|

  DELTZ=HMAX/FLOAT(KBBM1)

  DO K=1,KBBM1
     PHY_Z(K)=(0.5_SP-FLOAT(K))*DELTZ
  END DO

  !--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT Z LEVELS-------------------|

  DO I=1,M
     DO K=1,KBM1
        ZM(K)=ZZ(I,K)*DT(I)+EL(I)
        RHOS(K)= RHO1(I,K)                   !=DBLE(RHO1(I,K))
     END DO

     CALL SINTER_P(ZM,RHOS,PHY_Z,RHOZTMP,KBM1,KBBM1)

     RHOZ(I,:)=RHOZTMP(:)

  END DO

  !--DO THE AVERAGE OVER Z_levels 

  DO K=1,KBBM1
     RHOA(K)=0.0
     ICOUNT=0
     DO I=1,M
        IF(-H(I).LE.PHY_Z(K)) THEN
           ICOUNT = ICOUNT + 1
           RHOA(K)=RHOA(K)+RHOZ(I,K)
        END IF
     ENDDO

     IF(ICOUNT.GE.1) then
        RHOA(K) = RHOA(K)/float(ICOUNT)
     ELSE 
        IF(K .LT.1) CALL FATAL_ERROR("RHO_MEAN: BAD STATE, AVOIDING SIGSEV...")

        RHOA(K) = RHOA(K-1) 
     END IF
     !        write(*,'("K,RHOA= ",f10.6)') K,RHOA(K)
  END DO

  DEALLOCATE(RHOZ,stat=status)
  IF(STATUS/=0) CALL FATAL_ERROR("RHO_MEAN: CAN NOT DEALLOCATE RHOZ ARRAY!")

  !--LINEARLY INTERPOLATE TO OBTAIN DENSITY VALUES AT SIGMA LEVELS-------------------|


  DO I=1,M
     DO K=1,KBM1
        ZM(K)=ZZ(I,K)*DT(I)+EL(I)
     END DO

     CALL SINTER(PHY_Z,RHOA,ZM,RHOS,KBBM1,KBM1)
     DO K=1,KBM1
        RMEAN1(I,K)=RHOS(K)                  !=REAL(RHOS(K))
     END DO
  END DO

  RMEAN1(:,KB)=0.0_SP


  CALL N2E3D(RMEAN1,RMEAN)

  RETURN
END SUBROUTINE RHO_MEAN
