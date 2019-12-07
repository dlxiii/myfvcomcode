! $Id: mod_bulk.F,v 1.1.2.3 2008/04/03 16:39:00 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.3 $
MODULE mod_bulk
  USE MOD_PREC
  implicit none

CONTAINS

  SUBROUTINE SIMPLE_DRAG(spdx,spdy,strx,stry)
    IMPLICIT NONE
    REAL(SP), POINTER, INTENT(IN)  :: SPDX(:),SPDY(:)
    REAL(SP), POINTER, INTENT(OUT) :: STRX(:),STRY(:)
    INTEGER :: I, N
    REAL(SP) :: CD, WDS, TX, TY


    N = UBOUND(SPDX,1)

    IF(N /= UBOUND(SPDY,1)) WRITE(5,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"
    IF(N /= UBOUND(STRY,1)) WRITE(5,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"
    IF(N /= UBOUND(STRX,1)) WRITE(5,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"

    DO I=1,N
       TX = SPDX(I)
       TY = SPDY(I)
       WDS=SQRT(TX*TX+TY*TY)
       CD=1.2E-3
       IF (WDS >= 11.0_SP) THEN
          CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
       ELSE IF (WDS >= 25.0_SP)THEN
          CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP
       END IF

       STRX(I) = 1.2_SP*CD*TX*WDS
       STRY(I) = 1.2_SP*CD*TY*WDS

    END DO
    


  END SUBROUTINE SIMPLE_DRAG
  


END MODULE mod_bulk
