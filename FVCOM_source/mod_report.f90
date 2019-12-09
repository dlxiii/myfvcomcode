










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

MODULE MOD_REPORT
  IMPLICIT NONE
  
  CONTAINS

  FUNCTION REPORT_NOW(IINT,IREPORT)
    USE MOD_TIME
    IMPLICIT NONE
    LOGICAL :: REPORT_NOW
    INTEGER(ITIME) :: IINT
    INTEGER :: IREPORT
    INTEGER(ITIME) :: IREPORT_DB

    REPORT_NOW = .FALSE.
    IF(ireport == 0) return
    IREPORT_DB = IREPORT
    IF(mod(IINT,IREPORT_DB) /= 0)  return

    REPORT_NOW = .TRUE.

  END FUNCTION REPORT_NOW
  !==============================================================================|
  SUBROUTINE REPORT(INFO_STRING)
    !==============================================================================|
    !     REPORT INITIAL INFORMATION                                               |
    !==============================================================================|
    USE ALL_VARS
    USE MOD_WD
    USE MOD_ICE
    IMPLICIT NONE
    CHARACTER(LEN=*) :: INFO_STRING     !!INFORMATION STRING
    INTEGER :: E3TOT,ESTOT,IERR
    REAL(DP), DIMENSION(22) :: SBUF,RBUF1,RBUF2,RBUF3

    REAL(SP), ALLOCATABLE :: AICE1(:), VICE1(:)
    REAL(SP), ALLOCATABLE :: Q21(:,:),Q2L1(:,:),L1(:,:)
    REAL(SP), ALLOCATABLE :: KH1(:,:),KQ1(:,:)
    INTEGER :: I,J,K
    !==============================================================================|

    ALLOCATE(Q21(1:N,KBM1));   Q21   = 0.0_SP
    ALLOCATE(Q2L1(1:N,KBM1));  Q2L1  = 0.0_SP
    ALLOCATE(L1(1:N,KBM1));    L1    = 0.0_SP
    ALLOCATE(KH1(1:N,KBM1));   KH1   = 0.0_SP
    ALLOCATE(KQ1(1:N,KBM1));   KQ1   = 0.0_SP






    DO K=1,KBM1
       DO I=1,N
          DO J=1,3
             Q21(I,K)  = Q21(I,K)+Q2(NV(I,J),K)
             Q2L1(I,K) = Q2L1(I,K)+Q2L(NV(I,J),K)
             L1(I,K)   = L1(I,K)+L(NV(I,J),K)
             KH1(I,K)  = KH1(I,K)+KH(NV(I,J),K)
             KQ1(I,K)  = KQ1(I,K)+KQ(NV(I,J),K)
          END DO
          Q21(I,K)  = Q21(I,K)/3.0_SP
          Q2L1(I,K) = Q2L1(I,K)/3.0_SP
          L1(I,K)   = L1(I,K)/3.0_SP
          KH1(I,K)  = KH1(I,K)/3.0_SP
          KQ1(I,K)  = KQ1(I,K)/3.0_SP
       END DO
    END DO

    SBUF = 0.0_DP

    SBUF(1)  = SUM(DBLE(UA(1:N)))
    SBUF(2)  = SUM(DBLE(VA(1:N)))
    SBUF(3)  = SUM(DBLE(EL1(1:N)))
    SBUF(4)  = SUM(DBLE(H1(1:N)))
    SBUF(5)  = SUM(DBLE(U(1:N,1:KBM1)))
    SBUF(6)  = SUM(DBLE(V(1:N,1:KBM1)))
    SBUF(7)  = SUM(DBLE(S(1:N,1:KBM1)))
    SBUF(8)  = SUM(DBLE(T(1:N,1:KBM1)))
    SBUF(9)  = SUM(DBLE(Q21(1:N,2:KBM1)))
    SBUF(10) = SUM(DBLE(Q2L1(1:N,2:KBM1)))
    SBUF(11) = SUM(DBLE(L1(1:N,1:KBM1)))
    SBUF(12) = SUM(DBLE(KM1(1:N,1:KBM1)))
    SBUF(13) = SUM(DBLE(KQ1(1:N,1:KBM1)))
    SBUF(14) = SUM(DBLE(KH1(1:N,1:KBM1)))
    SBUF(15) = SUM(DBLE(RHO(1:N,1:KBM1)))
    SBUF(16) = SUM(DBLE(D1(1:N)))
    SBUF(17) = SUM(ISWETC(1:N))


    RBUF1 = SBUF
    IF(PAR)CALL MPI_REDUCE(SBUF,RBUF1,22,MPI_DP,MPI_SUM,MSRID-1,MPI_FVCOM_GROUP,IERR)

    SBUF = 0.0_DP

    SBUF(1)  = MAXVAL(UA(1:N))
    SBUF(2)  = MAXVAL(VA(1:N))
    SBUF(3)  = MAXVAL(EL(1:M))
    SBUF(4)  = MAXVAL(H(1:M))
    SBUF(5)  = MAXVAL(U(1:N,1:KBM1))
    SBUF(6)  = MAXVAL(V(1:N,1:KBM1))
    SBUF(7)  = MAXVAL(S1(1:M,1:KBM1))
    SBUF(8)  = MAXVAL(T1(1:M,1:KBM1))
    SBUF(9)  = MAXVAL(Q2(1:M,1:KBM1))
    SBUF(10) = MAXVAL(Q2L(1:M,1:KBM1))
    SBUF(11) = MAXVAL(L(1:M,1:KBM1))
    SBUF(12) = MAXVAL(KM(1:M,1:KBM1))
    SBUF(13) = MAXVAL(KQ(1:M,1:KBM1))
    SBUF(14) = MAXVAL(KH(1:M,1:KBM1))
    SBUF(15) = MAXVAL(RHO1(1:M,1:KBM1))
    SBUF(16) = MAXVAL(D(1:M))



    RBUF2 = SBUF
    IF(PAR)CALL MPI_REDUCE(SBUF,RBUF2,22,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)

    SBUF = 0.0_DP

    SBUF(1)  = MINVAL(UA(1:N))
    SBUF(2)  = MINVAL(VA(1:N))
    SBUF(3)  = MINVAL(EL(1:M))
    SBUF(4)  = MINVAL(H(1:M))
    SBUF(5)  = MINVAL(U(1:N,1:KBM1))
    SBUF(6)  = MINVAL(V(1:N,1:KBM1))
    SBUF(7)  = MINVAL(S1(1:M,1:KBM1))
    SBUF(8)  = MINVAL(T1(1:M,1:KBM1))
    SBUF(9)  = MINVAL(Q2(1:M,1:KBM1))
    SBUF(10)  = MINVAL(Q2L(1:M,1:KBM1))
    SBUF(11) = MINVAL(L(1:M,1:KBM1))
    SBUF(12) = MINVAL(KM(1:M,1:KBM1))
    SBUF(13) = MINVAL(KQ(1:M,1:KBM1))
    SBUF(14) = MINVAL(KH(1:M,1:KBM1))
    SBUF(15) = MINVAL(RHO1(1:M,1:KBM1))
    SBUF(16) = MINVAL(D(1:M))



    RBUF3 = SBUF
    IF(PAR)CALL MPI_REDUCE(SBUF,RBUF3,22,MPI_DP,MPI_MIN,MSRID-1,MPI_FVCOM_GROUP,IERR)

    IF(MSR)THEN
       IF(LEN_TRIM(INFO_STRING) /= 0)THEN
          WRITE(IPT,*  )'!===================',TRIM(INFO_STRING),'======================'
       END IF
       RBUF1(15) = (RBUF1(15)+NGL*KBM1)*1000.
       RBUF2(15) = (RBUF2(15)+1.)*1000.
       RBUF3(15) = (RBUF3(15)+1.)*1000.
       E3TOT = DBLE(NGL*KBM1)
       ESTOT = DBLE(NGL)
       WRITE(IPT,*  )'!  QUANTITY              :     AVG           MAX         MIN'
       WRITE(IPT,100)'!  EXTERNAL UVEL         :',RBUF1(1)/ESTOT,RBUF2(1),RBUF3(1)
       WRITE(IPT,100)'!  EXTERNAL VVEL         :',RBUF1(2)/ESTOT,RBUF2(2),RBUF3(2)
       WRITE(IPT,100)'!  FREE SURFACE          :',RBUF1(3)/ESTOT,RBUF2(3),RBUF3(3)
       WRITE(IPT,100)'!  BATH DEPTH            :',RBUF1(4)/ESTOT,RBUF2(4),RBUF3(4)
       WRITE(IPT,100)'!  INTERNAL UVEL         :',RBUF1(5)/E3TOT,RBUF2(5),RBUF3(5)
       WRITE(IPT,100)'!  INTERNAL VVEL         :',RBUF1(6)/E3TOT,RBUF2(6),RBUF3(6)
       WRITE(IPT,100)'!  SALINITY              :',RBUF1(7)/E3TOT,RBUF2(7),RBUF3(7)
       WRITE(IPT,100)'!  TEMPERATURE           :',RBUF1(8)/E3TOT,RBUF2(8),RBUF3(8)
       WRITE(IPT,100)'!  TURBULENT KE          :',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
       WRITE(IPT,100)'!  TURB KE*L             :',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)
       WRITE(IPT,100)'!  TURB LENGTH SCALE     :',RBUF1(11)/E3TOT,RBUF2(11),RBUF3(11)
       WRITE(IPT,100)'!  KM                    :',RBUF1(12)/E3TOT,RBUF2(12),RBUF3(12)
       WRITE(IPT,100)'!  KQ                    :',RBUF1(13)/E3TOT,RBUF2(13),RBUF3(13)
       WRITE(IPT,100)'!  KH                    :',RBUF1(14)/E3TOT,RBUF2(14),RBUF3(14)
       WRITE(IPT,100)'!  DENSITY               :',RBUF1(15)/E3TOT,RBUF2(15),RBUF3(15)
       WRITE(IPT,100)'!  DEPTH                 :',RBUF1(16)/ESTOT,RBUF2(16),RBUF3(16)

       WRITE(IPT,*  )'!  WET/DRY INFO          :   #WET       #DRY             %WET'
       IF(RBUF1(17) == FLOAT(NGL))THEN
          WRITE(IPT,*)'!  NO DRY POINTS          '
       ELSE
          WRITE(IPT,101)'!  WET/DRY DATA          :',INT(RBUF1(17)),NGL-INT(RBUF1(17)),100.*RBUF1(17)/FLOAT(NGL)
       END IF

    END IF


    DEALLOCATE(Q21,Q2L1,L1)
    DEALLOCATE(KH1,KQ1)

    RETURN
100 FORMAT(1X,A26,3F12.6)
101 FORMAT(1X,A26,2I12,F12.6)
  END SUBROUTINE REPORT
  !==============================================================================|

END MODULE MOD_REPORT
