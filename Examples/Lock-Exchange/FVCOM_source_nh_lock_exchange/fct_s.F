! $Id: fct_s.F,v 1.3.2.4 2008/04/24 22:49:58 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.3.2.4 $
!==============================================================================|
!  FLUX CONTROL FOR SALINITY                                                        |
!==============================================================================|

SUBROUTINE FCT_S
  !#  if defined (WET_DRY)

  !==============================================================================|
  USE ALL_VARS
  USE MOD_UTILS
  USE BCS
  USE MOD_OBCS
  IMPLICIT NONE
  REAL(SP):: SMAX,SMIN
  INTEGER :: I,J,K
  !==============================================================================|

  IF(HEATING_TYPE == 'body') RETURN

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: fct_s"

  nodes: DO I=1,M

     ! SKIP OPEN BOUNDARY NODES
     IF(IOBCN > 0)THEN
        DO J=1,IOBCN
           IF(I == I_OBC_N(J)) CYCLE nodes
        END DO
     END IF

     ! SKIP RIVER INFLOW POINTS
     IF(NUMQBC > 0)THEN
        DO J=1,NUMQBC
           IF(RIVER_INFLOW_LOCATION == 'node')THEN
              IF(I == INODEQ(J)) CYCLE nodes
           END IF
           IF(RIVER_INFLOW_LOCATION == 'edge')THEN
              IF(I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2)) CYCLE nodes
           END IF
        END DO
     END IF

     ! SKIP GROUND WATER INFLOW POINTS
     IF(BFWDIS(I) .GT. 0.0_SP .and. GROUNDWATER_SALT_ON) CYCLE nodes

     DO K=1,KBM1
        SMAX = MAXVAL(S1(NBSN(I,1:NTSN(I)),K))
        SMIN = MINVAL(S1(NBSN(I,1:NTSN(I)),K))

        IF(K == 1)THEN
           SMAX = MAX(SMAX,(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
           SMIN = MIN(SMIN,(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K+1)))
        ELSE IF(K == KBM1)THEN
           SMAX = MAX(SMAX,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
           SMIN = MIN(SMIN,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)))
        ELSE
           SMAX = MAX(SMAX,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
           SMIN = MIN(SMIN,(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/  &
                (DZ(I,K)+DZ(I,K-1)),                             &
                (S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/           &
                (DZ(I,K)+DZ(I,K+1)))
        END IF

        IF(SMIN-SF1(I,K) > 0.0_SP)SF1(I,K) = SMIN
        IF(SF1(I,K)-SMAX > 0.0_SP)SF1(I,K) = SMAX

     END DO

  END DO nodes

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: fct_s"
  !#  endif
END SUBROUTINE FCT_S
!==============================================================================|


