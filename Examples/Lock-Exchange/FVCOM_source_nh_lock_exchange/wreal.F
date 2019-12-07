! $Id: wreal.F,v 1.2.2.3 2008/04/03 17:10:48 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.2.2.3 $
!==============================================================================|
!  Compute Cartesian Vertical Velocity                                         |
!==============================================================================|
   SUBROUTINE WREAL               
!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP) :: DDDX,DDDY,DEDX,DEDY,ETF1AA,WW1,WW2 
   INTEGER  :: I,K,J1,J2,J3
!==============================================================================|

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: wreal"

!------------------------------------------------------------------------------!
!  SAVE OMEGA VELOCITY FROM PREVIOUS TIMESTEP (USED FOR LAGRANGIAN TRACKING)   !
!------------------------------------------------------------------------------!

   WTTS = WTS

!------------------------------------------------------------------------------!
!  CALCULATE A NEW OMEGA VELOCITY                                              !
!------------------------------------------------------------------------------!

   DO I=1,N
#  if defined (WET_DRY)
    IF(ISWETC(I) == 1)THEN
#  endif
     J1=NV(I,1)
     J2=NV(I,2)
     J3=NV(I,3)
     DDDX=AWX(I,1) * D(J1)+AWX(I,2) * D(J2)+AWX(I,3)*D(J3)
     DDDY=AWY(I,1) * D(J1)+AWY(I,2) * D(J2)+AWY(I,3)*D(J3)
     DEDX=AWX(I,1)*ELF(J1)+AWX(I,2)*ELF(J2)+AWX(I,3)*ELF(J3)
     DEDY=AWY(I,1)*ELF(J1)+AWY(I,2)*ELF(J2)+AWY(I,3)*ELF(J3)
     ETF1AA=ONE_THIRD*(EL(NV(I,1))+EL(NV(I,2))+EL(NV(I,3)))
     DO K=1,KBM1
      WW1=0.5_SP*(W(I,K)+W(I,K+1))+U(I,K)*(ZZ1(I,K)*DDDX+DEDX)+ &
                                V(I,K)*(ZZ1(I,K)*DDDY+DEDY)
      WW2=(ZZ1(I,K)+1.)*(ETF1AA-ET1(I))/DTI
      WW(I,K)=WW1+WW2
     END DO
#  if defined (WET_DRY)
    ELSE
     DO K=1,KBM1
      WW(I,K)=0.0_SP
     END DO
    END IF
#  endif
   END DO

   IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: wreal"

   END SUBROUTINE WREAL
!==============================================================================|
