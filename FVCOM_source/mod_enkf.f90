










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

MODULE ENKFVAL
END MODULE ENKFVAL

!=============================================
!MODULE MOD_STARTUP_MIMIC
!# if defined (ENKF)
!  USE MOD_UTILS
!  USE MOD_NCTOOLS
!  USE MOD_INPUT
!  USE ALL_VARS
!  USE EQS_OF_STATE
!  USE MOD_WD
!  USE SINTER
!
!
!#  if defined (ICE)
!  USE MOD_ICE
!  USE MOD_ICE2D
!#  endif
!
!# if defined (WATER_QUALITY)
!  USE MOD_WQM
!# endif
!  
!# if defined (BioGen)
!  USE MOD_BIO_3D
!# endif
!
!# if defined (NH)
!  USE NON_HYDRO, ONLY: W4ZT, NHQDRX, NHQDRY, NHQDRZ, NHQ2DX, NHQ2DY
!# endif
!  
!  IMPLICIT NONE
!  
!  PRIVATE
!
!  PUBLIC :: READ_SSH1
!  PUBLIC :: READ_UV1
!  PUBLIC :: READ_TURB1
!  PUBLIC :: READ_TS1
!  PUBLIC :: READ_WETDRY1
!
!  
!CONTAINS
!!==============================================================================!
!  SUBROUTINE READ_SSH1(NC_START)
!#   if defined (SEDIMENT)
!    USE MOD_SED, only : morpho_model,sed_hot_start 
!#   endif
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"
!    
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!    ! LOAD EL
!    VAR => FIND_VAR(NC_START,'zeta',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'zeta'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, EL)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ET
!    VAR => FIND_VAR(NC_START,'et',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'et'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ET)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    !----------------------------------------------------------------
!    ! Read the most recent bathymetry if Morphodynamics is Active
!    !----------------------------------------------------------------
!#   if defined(SEDIMENT)
!    IF(MORPHO_MODEL .and. SED_HOT_START)THEN
!    VAR => FIND_VAR(NC_START,'h',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'h'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, h)
!    CALL NC_READ_VAR(VAR,STKCNT)
!    ENDIF
!#   endif
!
!
!    !----------------------------------------------------------------
!    ! Given SSH and Bathy, Update the Bathymetry 
!    !----------------------------------------------------------------
!    D  = H + EL
!    DT = H + ET
!
!    
!    CALL N2E2D(H,H1)
!    CALL N2E2D(EL,EL1)
!    CALL N2E2D(D,D1)
!    CALL N2E2D(DT,DT1)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_SSH"
!
!  END SUBROUTINE READ_SSH1
!!==============================================================================!
!!==============================================================================!
!  SUBROUTINE READ_WETDRY1(NC_START)
!    USE MOD_WD
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_WETDRY"
!
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!    ! LOAD ISWETN
!    VAR => FIND_VAR(NC_START,'wet_nodes',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETN)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ISWETC
!    VAR => FIND_VAR(NC_START,'wet_cells',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETC)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    ! LOAD ISWETNT
!    VAR => FIND_VAR(NC_START,'wet_nodes_prev_int',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_nodes_prev_int'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETNT)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ISWETCT
!    VAR => FIND_VAR(NC_START,'wet_cells_prev_int',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_int'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETCT)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD ISWETCE
!    VAR => FIND_VAR(NC_START,'wet_cells_prev_ext',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'wet_cells_prev_ext'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, ISWETC)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_WETDRY"
!    
!  END SUBROUTINE READ_WETDRY1
!!==============================================================================!
!  SUBROUTINE READ_TS1(NC_START)
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT, K
!    REAL(SP), DIMENSION(0:MT,KB) :: PRESSURE
!    
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "Start: READ_TS" 
!   
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!
!    ! LOAD TEMPERATURE
!    VAR => FIND_VAR(NC_START,'temp',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'temp'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, T1)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    ! LOAD MEAN INITIAL TEMPERATURE
!    VAR => FIND_VAR(NC_START,'tmean1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tmean1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, tmean1)
!    CALL NC_READ_VAR(VAR)
!
!    ! LOAD SALINITY
!    VAR => FIND_VAR(NC_START,'salinity',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'saltinity'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, S1)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD MEAN INITIAL SALINITY
!    VAR => FIND_VAR(NC_START,'smean1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'smean1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, smean1)
!    CALL NC_READ_VAR(VAR)
!
!
!    ! LOAD DENSITY
!    VAR => FIND_VAR(NC_START,'rho1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rho1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, RHO1)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    ! LOAD MEAN DENSITY
!    VAR => FIND_VAR(NC_START,'rmean1',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'rmean1'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, rmean1)
!    CALL NC_READ_VAR(VAR)
!
!    ! AVERAGE FROM CELLS TO FACE CENTERS
!
!
!!JQI    SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
!!JQI    CASE(SW_DENS1)
!
!!JQI       ! SET MEAN DENSITY
!!JQI       DO K=1,KBM1
!!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
!!JQI       END DO
!!JQI       CALL FOFONOFF_MILLARD(SMEAN1,TMEAN1,PRESSURE,0.0_SP,RMEAN1)
!!JQI       RMEAN1(0,:)=0.0_SP
!!JQI       RMEAN1(:,KB)=0.0_SP
!
!!JQI       ! SET REAL DENSITY
!!JQI       CALL DENS1 ! GENERIC CALL TO FOFONOFF_MILLARD FOR S1,T1...
!       
!!JQI    CASE(SW_DENS2)
!!JQI       ! SET MEAN DENSITY
!!JQI       CALL DENS2G(SMEAN1,TMEAN1,RMEAN1)
!!JQI       RMEAN1(0,:)=0.0_SP
!!JQI       RMEAN1(:,KB)=0.0_SP
!
!!JQI       ! SET REAL DENSITY
!!JQI       CALL DENS2 ! GENERIC CALL TO DENS2G FOR S1,T1...
!
!!JQI    CASE(SW_DENS3)
!
!!JQI       ! SET MEAN DENSITY
!!JQI       DO K=1,KBM1
!!JQI          PRESSURE(:,K) = -GRAV_N*1.025_SP*(ZZ(:,K)*D(:))*0.1_SP
!!JQI       END DO
!!JQI       CALL JACKET_MCDOUGALL(SMEAN1,TMEAN1,PRESSURE,RMEAN1)
!!JQI       RMEAN1(0,:)=0.0_SP
!!JQI       RMEAN1(:,KB)=0.0_SP
!
!!JQI       ! SET REAL DENSITY
!!JQI       CALL DENS3 ! GENERIC CALL TO JACKET_MCDOUGALL FOR S1,T1.
!
!!JQI    CASE DEFAULT
!!JQI       CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
!!JQI            & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
!!JQI    END SELECT
!    
!    CALL N2E3D(T1,T)
!    CALL N2E3D(S1,S)
!    CALL N2E3D(Tmean1,Tmean)
!    CALL N2E3D(Smean1,Smean)
!    CALL N2E3D(Rmean1,Rmean)
!
!   
!   IF(DBG_SET(DBG_SBR)) write(ipt,*) "End: READ_TS" 
!   
! END SUBROUTINE READ_TS1
!!==============================================================================!
!
!  SUBROUTINE READ_UV1(NC_START)
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_UV" 
!
!
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!    VAR => FIND_VAR(NC_START,'u',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'u'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, U)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    VAR => FIND_VAR(NC_START,'v',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'v'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, V)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'omega',FOUND)
!    IF(FOUND) THEN
!       CALL NC_CONNECT_AVAR(VAR, WTS)
!       CALL NC_READ_VAR(VAR,STKCNT)
!       
!       CALL N2E3D(WTS,W)
!    ELSE
!       VAR => FIND_VAR(NC_START,'w',FOUND)
!       IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'w' &
!            & or 'omega' IN THE HOTSTART FILE OBJECT")
!       CALL NC_CONNECT_AVAR(VAR, W)
!       CALL NC_READ_VAR(VAR,STKCNT)
!    END IF
!
!    VAR => FIND_VAR(NC_START,'ua',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'ua'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, UA)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'va',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'va'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, VA)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_UV" 
!
!  END SUBROUTINE READ_UV1
!!==============================================================================!
!  SUBROUTINE READ_TURB1(NC_START)
!    IMPLICIT NONE
!    TYPE(NCFILE),POINTER :: NC_START
!    TYPE(NCVAR),  POINTER :: VAR
!    TYPE(NCDIM),  POINTER :: DIM
!    LOGICAL :: FOUND
!    INTEGER :: STKCNT
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "START: READ_TURB" 
!
!    STKCNT = NC_START%FTIME%PREV_STKCNT
!
!#    if defined (GOTM)
!
!    VAR => FIND_VAR(NC_START,'tke',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'tke'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, TKE)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'teps',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'teps'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, TEPS)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    L = .001 
!    L(1:MT,2:KBM1) = (.5544**3)*TKE(1:MT,2:KBM1)**1.5/TEPS(1:MT,2:KBM1)
!    
!# else
!
!    VAR => FIND_VAR(NC_START,'q2',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, Q2)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'q2l',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'q2l'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, Q2L)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'l',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'l'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, L)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!# endif
!
!    VAR => FIND_VAR(NC_START,'km',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'km'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, km)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'kq',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kq'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, KQ)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!    VAR => FIND_VAR(NC_START,'kh',FOUND)
!    IF(.not. FOUND) CALL FATAL_ERROR("COULD NOT FIND VARIABLE 'kh'&
!         & IN THE HOTSTART FILE OBJECT")
!    CALL NC_CONNECT_AVAR(VAR, KH)
!    CALL NC_READ_VAR(VAR,STKCNT)
!
!
!    CALL N2E3D(KM,KM1)
!
!    IF(DBG_SET(DBG_SBR)) write(ipt,*) "END: READ_TURB" 
!
!  END SUBROUTINE READ_TURB1
!!==============================================================================|
!# endif
!END MODULE MOD_STARTUP_MIMIC

MODULE MOD_ENKF
   USE ENKFVAL
   USE mod_ncdio, only : update_iodata 
   USE MOD_INPUT, only : NC_START
   USE CONTROL
   USE MOD_UTILS
   USE MOD_NCTOOLS
   IMPLICIT NONE
   SAVE
   

   LOGICAL  :: ENKF_ON

END MODULE MOD_ENKF

