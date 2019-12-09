












!****************************************************************
!
  SUBROUTINE ICETT2
!
!****************************************************************
!
  USE M_GENARR
  USE SWCOMM2
  USE SWCOMM3
  USE SWCOMM4
  USE OCPCOMM4
  USE M_PARALL
  USE ALL_VARS
  USE VARS_WAVE

  INTEGER :: I,ISS,ID

  DO I = 1,MT
   DO ISS = 1,MSC
    DO ID = 1,MDC
     if (Sice(ID,ISS,I)>0) then
       AC2(ID,ISS,I)=AC2(ID,ISS,I)-Sice(ID,ISS,I)
       AC2(ID,ISS,I)=MAX(0.0,AC2(ID,ISS,I))
     endif
    END DO
   END DO
  END DO

  RETURN
  END SUBROUTINE ICETT2

