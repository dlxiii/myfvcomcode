! $Id: fct_q2l.F,v 1.2.2.3 2008/04/03 16:27:06 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.2.2.3 $
!==============================================================================|
!  FLUX CONTROL FOR SALINITY                                                        |
!==============================================================================|

   SUBROUTINE FCT_Q2L
!#  if defined (WET_DRY)

!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_UTILS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP):: Q2LMAX,Q2LMIN
   INTEGER :: I,J,K
!==============================================================================|

   IF(HEATING_TYPE == 'body') return

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: fct_q2"

nodes:DO I=1,M

     IF(IOBCN > 0)THEN
       DO J=1,IOBCN
         IF(I == I_OBC_N(J)) CYCLE nodes
       END DO
     END IF  	 
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
     DO K=2,KBM1
       Q2LMAX = MAXVAL(Q2L(NBSN(I,1:NTSN(I)),K))
       Q2LMIN = MINVAL(Q2L(NBSN(I,1:NTSN(I)),K))

       IF(K == 2)THEN
         Q2LMAX = MAX(Q2LMAX,(Q2L(I,K)*DZZ(I,K)+Q2L(I,K+1)*DZZ(I,K-1))/  &
                  (DZZ(I,K)+DZZ(I,K-1)))
         Q2LMIN = MIN(Q2LMIN,(Q2L(I,K)*DZZ(I,K)+Q2L(I,K+1)*DZZ(I,K-1))/  &
                  (DZZ(I,K)+DZZ(I,K-1)))
       ELSE IF(K == KBM1)THEN
         Q2LMAX = MAX(Q2LMAX,(Q2L(I,K)*DZZ(I,K-2)+Q2L(I,K-1)*DZZ(I,K-1))/  &
                  (DZZ(I,K-1)+DZZ(I,K-2)))
         Q2LMIN = MIN(Q2LMIN,(Q2L(I,K)*DZZ(I,K-2)+Q2L(I,K-1)*DZZ(I,K-1))/  &
                  (DZZ(I,K-1)+DZZ(I,K-2)))
       ELSE
         Q2LMAX = MAX(Q2LMAX,(Q2L(I,K)*DZZ(I,K-2)+Q2L(I,K-1)*DZZ(I,K-1))/  &         
                  (DZZ(I,K-1)+DZZ(I,K-2)),                                 &
                  (Q2L(I,K)*DZZ(I,K)+Q2L(I,K+1)*DZZ(I,K-1))/               &
                  (DZZ(I,K)+DZZ(I,K-1)))
         Q2LMIN = MIN(Q2LMIN,(Q2L(I,K)*DZZ(I,K-2)+Q2L(I,K-1)*DZZ(I,K-1))/  &
                  (DZZ(I,K-1)+DZZ(I,K-2)),                                 &
                  (Q2L(I,K)*DZZ(I,K)+Q2L(I,K+1)*DZZ(I,K-1))/               &
                  (DZZ(I,K)+DZZ(I,K-1)))
       END IF

       IF(Q2LMIN-Q2LF(I,K) > 0.0_SP)Q2LF(I,K) = Q2LMIN
       IF(Q2LF(I,K)-Q2LMAX > 0.0_SP)Q2LF(I,K) = Q2LMAX

     END DO

   END DO nodes

   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: fct_q2"

   END SUBROUTINE FCT_Q2L
!==============================================================================|


