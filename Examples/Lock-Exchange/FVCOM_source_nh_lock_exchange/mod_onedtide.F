! $Id: mod_onedtide.F,v 1.1.2.3 2008/05/09 17:20:26 jqi Exp $
! $Name: New_Input $
! $Revision: 1.1.2.3 $
MODULE MOD_ONEDTIDE
#  if defined (ONE_D_MODEL)
   USE MOD_TYPES
   USE BCS
   
   IMPLICIT NONE
   SAVE
   
   INTEGER, PARAMETER :: NUMTIDES_MAX = 8
   INTEGER  :: NUMTIDES
   REAL(SP) :: UMAX(NUMTIDES_MAX),PMAX(NUMTIDES_MAX)

   NAMELIST /NML_ONEDTIDE/      &
        & NUMTIDES,             &
	& UMAX,                 &
	& PMAX
	   
   CONTAINS
!============================================================================!
!
!============================================================================!
   SUBROUTINE GET_UMAX
   USE ALL_VARS
   USE MOD_UTILS
   
   IMPLICIT NONE
   integer :: ios, i
   Character(Len=120):: FNAME
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "Subroutine Begins: get_umax;"

    ios = 0

    FNAME = "./"//trim(casename)//"_run.nml"

    if(DBG_SET(dbg_io)) &
         & write(IPT,*) "Get_umax: File: ",trim(FNAME)

    CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

    !READ NAME LIST FILE

    ! Read name list nml_onedtide
    READ(UNIT=NMLUNIT, NML=NML_ONEDTIDE,IOSTAT=ios)
    if(ios .NE. 0 ) then
       if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_ONEDTIDE)
       Call Fatal_Error("Can Not Read NameList NML_ONEDTIDE from file: "//trim(FNAME))
    end if

    REWIND(NMLUNIT)

    if(DBG_SET(dbg_scl)) &
         & write(IPT,*) "Read_Name_List:"

    if(DBG_SET(dbg_scl)) &
         & write(UNIT=IPT,NML=NML_ONEDTIDE)


    IF(NUMTIDES > nTideComps)THEN
      CALL FATAL_ERROR("NUMBER OF TIDAL CONSTITUTES IS GREATER THAN nTideComps IN 1D MODULE")
    END IF
      
    CLOSE(NMLUNIT)

   RETURN
   END SUBROUTINE GET_UMAX
!===========================================================================!
!
!===========================================================================!        
   SUBROUTINE ONEDTIDE(PSTX_TM,PSTY_TM)
   USE ALL_VARS
   USE BCS
   IMPLICIT NONE
   TYPE(TIME):: TIME_ELAPSED
   INTEGER :: I,J,K
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB),TIME1
   REAL(DP):: TIME_FLT
!**************     EXECUTABLE     ***************
!   TIME1 = TIME * 86400.0_SP
   TIME_ELAPSED = RKTime - SpecTime
   
   PSTX_TM = 0.0_SP
   PSTY_TM = 0.0_SP
!***************  INITIALIZATION   ***************
   IF (IINT == 1) THEN
     DO K=1,KBM1
       DO I=1,NT
         DO J=1,6
           PSTX_TM(I,K) = PSTX_TM(I,K)+UMAX(J)/SQRT(2._SP)/DTI
         END DO
       END DO
     END DO
   END IF 

!***************  COMPUTATION   ***************
   DO K=1,KBM1
     DO I=1,NT
       DO J=1,nTideComps    !6
         TIME_FLT = SECONDS(TIME_ELAPSED * PI2/PERIOD(J))
         PSTX_TM(I,K) = PSTX_TM(I,K)-UMAX(J)*PI2/(PERIOD(J)*SQRT(2.0_SP))*  &
	                SIN(TIME_FLT-PMAX(J))
!	                SIN(PI2/PERIOD(J)*TIME1-PMAX(J))
       END DO
       PSTX_TM(I,K) = PSTX_TM(I,K)*DT1(I)
       PSTY_TM(I,K) = PSTX_TM(I,K)
     END DO
   END DO
			
   RETURN
   END SUBROUTINE ONEDTIDE
!===========================================================================!
!
!===========================================================================!   
   SUBROUTINE NAME_LIST_INITIALIZE_ONED
   USE CONTROL
   
   IMPLICIT NONE
   
   !--Parameters in NameList NML_ONEDTIDE
   NUMTIDES = 0
   UMAX     = 0.0
   PMAX     = 0.0

   RETURN
   END SUBROUTINE NAME_LIST_INITIALIZE_ONED
!===========================================================================!
!
!===========================================================================!   
   SUBROUTINE NAME_LIST_PRINT_ONED
   USE CONTROL
   
   IMPLICIT NONE
   
   write(UNIT=IPT,NML=NML_ONEDTIDE)
   
   RETURN
   END SUBROUTINE NAME_LIST_PRINT_ONED
!===========================================================================!   
      
#  endif
END MODULE MOD_ONEDTIDE   
   
    
