MODULE MOD_ICE2D
#  if defined (ICE)
   USE ALL_VARS
   USE MOD_PREC
#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif
#  endif

# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif

!===============================================================================|
!   use CICE      
!===============================================================================|
   USE ICE_KINDS_MOD
   USE ICE_STATE
   use ice_mechred
   USE ICE_CONSTANTS
   use ice_calendar, only: dtice,dyn_dt
 
   IMPLICIT NONE
   SAVE
   
   REAL(SP), PARAMETER   :: ECC = 4.0
   REAL(SP), PARAMETER   :: ECCI = C1I/ECC
   real(SP), ALLOCATABLE :: strocnxE(:),strocnyE(:)! ice-ocean stress
   real(SP), ALLOCATABLE :: STRAIRXE(:),STRAIRYE(:)
   
   REAL(SP), ALLOCATABLE :: AICETMP(:),AIU(:)
   REAL(SP), ALLOCATABLE :: UICE2(:),VICE2(:),UICE2F(:),VICE2F(:)
   REAL(SP), ALLOCATABLE :: UICE2RK(:),VICE2RK(:)
   REAL(SP), ALLOCATABLE :: PSTXICE(:),PSTYICE(:)
   REAL(SP), ALLOCATABLE :: ADVSTRX(:),ADVSTRY(:)
   REAL(SP), ALLOCATABLE :: ADVUICE(:),ADVVICE(:)
   REAL(SP), ALLOCATABLE :: IMASS1(:),IMASS(:),PICE(:)
   REAL(SP), ALLOCATABLE :: TAUXWI(:),TAUYWI(:),VREL(:)
   INTEGER,  ALLOCATABLE :: ISICEN(:),ISICEC(:),ISICEC_OLD(:)
   
   real (kind=dbl_kind), dimension (:,:) ,allocatable,save ::&

!     &   shear        ! strain rate II component (1/s)
         divu       &  ! strain rate I component, velocity divergence (1/s)
      ,  Delta         ! function of strain rates (1/s)
!     &,  rcon         ! for damping criterion (kg/s)
!     &,  prs_sig      ! replacement pressure, for stress calc

!     logical (kind=log_kind) :: &
!        evp_damping  ! if true, use evp damping procedure

!      integer (kind=int_kind) ::  &
!        kdyn        & ! type of dynamics ( 1 = evp )
!      ,  ndte         ! number of subcycles:  ndte=dyn_dt/dte

      !ice_calendar.F:     &,  ndyn_dt   ! reduced timestep for dynamics: ndyn_dt=dt/dyn_dt
      !ice_calendar.F:     &,  dyn_dt             ! dynamics/transport/ridging timestep (s)
      !ice_calendar.F:      dyn_dt = dt/real(ndyn_dt,kind=dbl_kind) ! dynamics et al timestep
      !ice_dyn_evp.F:     &,  ndte         ! number of subcycles:  ndte=dyn_dt/dte
      !ice_dyn_evp.F:      dte = dyn_dt/real(ndte,kind=dbl_kind)        ! s
      !ice_calendar.F:      dyn_dt = dt/real(ndyn_dt,kind=dbl_kind) ! dynamics et al timestep
      !ice_dyn_evp.F:      dtei = c1/dte              ! 1/s
      !ice_dyn_evp.F:      rcon(:,:) = 1230._dbl_kind*eyc*dyn_dt*dtei**2  ! kg/s
      !ice_calendar.F:     &,  dtei               ! 1/dte, where dte is the EVP timestep (1/s)
      !ice_dyn_evp.F:      rcon(:,:) = 1230._dbl_kind*eyc*dyn_dt*dtei**2  ! kg/s
     
     REAL(DP), PARAMETER :: eyc =0.36
!     REAL(DP), PARAMETER :: rcon = 1230.*eyc*(dtice)*(1._SP/(dtice/4._SP))**2
     REAL(DP) :: rcon
  
   CHARACTER(LEN=80) :: RHEOLOGY

   REAL(SP), PARAMETER :: DRAGW = 0.00536_SP*RHOW
                                ! drag coefficient for water on ice *rhow (kg/m^3)
   REAL(SP), PARAMETER :: COSW = C1I  !cos(ocean turning angle) !turning angle = 0
   REAL(SP), PARAMETER :: SINW = C0I  !sin(ocean turning angle) !turning angle = 0
   
   CONTAINS

!===============================================================================|
!===============================================================================|
   SUBROUTINE ALLOC_UVICE

!   USE MOD_INP
   
   IMPLICIT NONE
  
   CHARACTER(LEN=120) :: FNAME
   INTEGER  :: ISCAN
   
   ALLOCATE(UICE2(0:NT));  UICE2 = 0.0_SP
   ALLOCATE(VICE2(0:NT));  VICE2 = 0.0_SP
   ALLOCATE(AIU(0:NT))  ;    AIU = 0.0_SP

   ALLOCATE(strocnxE(0:NT)); strocnxE =0.0_SP
   ALLOCATE(strocnyE(0:NT)); strocnyE =0.0_SP

   ALLOCATE(ISICEC(0:NT));  ISICEC = 0
   ALLOCATE(ISICEC_OLD(0:NT));  ISICEC_OLD = 0
   ALLOCATE(ISICEN(0:MT));  ISICEN = 0
!  for ridge redistribution
   allocate(divu(M,1)); divu  =0.0
   allocate(Delta(M,1)); Delta =0.0

!!$   IF(MSR) THEN 
!!$   FNAME = "./"//trim(casename)//"_run.dat"
!!$   ISCAN = SCAN_FILE(FNAME,"RHEOLOGY",CVAL = RHEOLOGY)
!!$   IF(ISCAN /= 0)THEN
!!$     WRITE(IPT,*)'ERROR READING RHEOLOGY: ',ISCAN
!!$     CALL PSTOP 
!!$   END IF
!!$   END IF

   rcon = 1230.*eyc*(dtice)*(1._SP/(dtice/4._SP))**2

   dyn_dt=dtice
   IF(MSR) write(ipt,*)' ICE dynamice time step  dyn_dt =' , dyn_dt

   RETURN
   END SUBROUTINE ALLOC_UVICE
!===============================================================================|
!
!===============================================================================|
   SUBROUTINE ICE_RUNGE_KUTTA

   IMPLICIT NONE
!   integer (kind=int_kind), intent(in) :: kstrngth
   integer (kind=int_kind) :: kstrngth
   INTEGER :: K   
   integer :: I,J
   real, allocatable :: rsub(:,:)
   real(SP) :: aice_tmp
   
!===============================================================================|
   !  temp use for calculate the air -ice stress
   real :: vmag_tmp      &! surface wind magnitude   (m/s)
      ,      rd_tmp      &! sqrt of exchange coefficient (momentum)
      ,     rdn_tmp      &! sqrt of exchange coefficient (momentum)
      ,     tau_tmp      &! stress at zlvl
      ,   ustar_tmp      &! ustar (m/s)
      ,    rhoa_tmp      &! air density (kg/m^3)
      ,windspeed_tmp

!===============================================================================|
   ! major/minor axis length ratio, squared
      
   ALLOCATE(AICETMP(0:MT));   AICETMP = 0.0_SP
   ALLOCATE(UICE2F(0:NT)) ;   UICE2F  = 0.0_SP
   ALLOCATE(VICE2F(0:NT)) ;   VICE2F  = 0.0_SP
   ALLOCATE(UICE2RK(0:NT));   UICE2RK = 0.0_SP
   ALLOCATE(VICE2RK(0:NT));   VICE2RK = 0.0_SP
   ALLOCATE(PSTXICE(0:NT));   PSTXICE = 0.0_SP
   ALLOCATE(PSTYICE(0:NT));   PSTYICE = 0.0_SP
   ALLOCATE(ADVSTRX(0:NT));   ADVSTRX = 0.0_SP
   ALLOCATE(ADVSTRY(0:NT));   ADVSTRY = 0.0_SP
   ALLOCATE(ADVUICE(0:NT));   ADVUICE = 0.0_SP
   ALLOCATE(ADVVICE(0:NT));   ADVVICE = 0.0_SP
   ALLOCATE(IMASS1(0:NT)) ;   IMASS1  = 0.0_SP
   ALLOCATE(PICE(1:MT))   ;   PICE    = 0.0_SP
   ALLOCATE(IMASS(0:MT))  ;   IMASS   = 0.0_SP
   ALLOCATE(TAUXWI(0:NT)) ;   TAUXWI  = 0.0_SP
   ALLOCATE(TAUYWI(0:NT)) ;   TAUYWI  = 0.0_SP
   ALLOCATE(VREL(0:NT))   ;   VREL  = 0.0_SP
         
   !-----------------------------------------------------------------
   ! pressure and forcing terms; set sigma=0 for no ice;
   ! initialize ice velocity in cells previously empty to ocn current
   !-----------------------------------------------------------------

   CALL ICE_STRENGTH(KSTRNGTH)
   DO I=1,M !T
      PICE(I) = STRENGTH(I,1)
   END DO
   PICE =PICE*RAMP

   !-----------------------------------------------------------------
   ! total mass of ice and snow, centered in T-cell
   ! NOTE: vice and vsno must be up to date in all grid cells,
   !       including ghost cells
   !-----------------------------------------------------------------

    ISICEC_OLD = ISICEC  !  ice in last step
    ISICEC = 0
    ISICEN = 0
    IMASS  =0.0_SP
    IMASS1 =0.0_SP

   DO I=1,M !T
      IF (AICE(I,1) >p001) then
      IMASS(I) = (RHOI*VICE(I,1) + RHOS*VSNO(I,1)) ! kg/m^2
      IF(IMASS(I) >P01)  ISICEN(I) =1
      END IF
   END DO
# if defined (MULTIPROCESSOR)
!       if(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,IMASS)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,IMASS)
# endif

   CALL  N2E2D(IMASS,IMASS1)

   DO I=1,M !T 
     AICETMP(I) = AICE(I,1)
   END DO
   AIU =0.0_SP

# if defined (MULTIPROCESSOR)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,AICETMP)
# endif

   CALL N2E2D(AICETMP,AIU)

     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,IMASS1,AIU)
#    endif
   
   DO I=1,N !T
      IF(sum(ISICEN(NV(I,1:3))) > 1 .and.IMASS1(I)>p01  &
         .and. AIU(I)>P001) &
      ISICEC(i)=1
   END DO

   DO I=1,N !T
     IF(ISICEC(I)==1) THEN
     DO J=1,3
!     IMASS1(I)=IMASS1(I)+IMASS(NV(I,J))*ISICEN(NV(I,J))
     END DO
!     IMASS1(I) =IMASS1(I)/(float(sum(ISICEN(NV(I,1:3)))))
     END IF
   END DO

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,IMASS1)
#    endif

! set new Ice speed --->ocean speed
   DO I=1,N  !T
     IF(ISICEC(I) ==1 .and. ISICEC_OLD(I)==0) THEN
!       UICE2(I) =U(I,1)  
!       VICE2(I) =V(I,1)  
     ENDIF
     UICE2(I) =UICE2(I)*ISICEC(I)
     VICE2(I) =VICE2(I)*ISICEC(I)

   END DO

#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UICE2,VICE2)
#    endif

   UICE2RK  = UICE2
   VICE2RK  = VICE2

   !-----------------------------------------------------------   
   !  calculate the wind stress on ice
   !-----------------------------------------------------------
     allocate(STRAIRXE(0:NT),STRAIRYE(0:NT))
     STRAIRXE = 0.0_SP
     STRAIRYE = 0.0_SP

      !------------------------------------------------------------
      ! Define functions
      !------------------------------------------------------------
      !------------------------------------------------------------
      ! Compute turbulent flux coefficients, wind stress, and 
      !------------------------------------------------------------
      !   rdn_tmp  = vonkar/log(zref/iceruf) ! neutral coefficient
        rdn_tmp = 0.4/log(10./0.0005)
         DO I=1,N
          windspeed_tmp=SQRT(uuwind(i)*uuwind(i)+vvwind(i)*vvwind(i))
!           if ( aicen(i,j,ni) > puny) then
      !------------------------------------------------------------
      ! define some needed variables
      !------------------------------------------------------------
         vmag_tmp = max(1.0, windspeed_tmp)
         ! ustar,tstar,qstar
         ustar_tmp = rdn_tmp * vmag_tmp
      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
      ! tau = rhoa(i,j) * ustar * ustar 
      ! strx = tau * uatm(i,j) / vmag 
      ! stry = tau * vatm(i,j) / vmag 
      !------------------------------------------------------------

         !tau = rhoa(i,1) * ustar_tmp * ustar_tmp ! not the stress at zlvl(i,j)
         tau_tmp = 1.3 * ustar_tmp * ustar_tmp ! not the stress at zlvl(i,j)
         !STRAIRXE(I) = tau_tmp * uatm(i,1)*AIU(I) * RAMP
         !STRAIRYE(I) = tau_tmp * vatm(i,1)*AIU(I) * RAMP    

         !FORCEX=aice*rhoair*dragair*windspeed*uatm
         STRAIRXE(I) =1.3*1.0E-3*windspeed_tmp*UUWIND(i)*AIU(I) * RAMP
         STRAIRYE(I) =1.3*1.0E-3*windspeed_tmp*VVwind(i)*AIU(I) * RAMP
       END DO
   !-----------------------------------------------------------
   !-----------------------------------------------------------
  
!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!

     CALL ICEOCEAN_STRESS
   DO K=1,4
     CALL ICE_STRESS(K)
     !CALL ICEOCEAN_STRESS

!     CALL FLUX_OBN(K)
     IF(trim(RHEOLOGY) == 'VP')THEN
!       CALL ADVAVE_ICE(ADVUICE,ADVVICE)           
     END IF

     CALL ICE_UV(K)

!     CALL BCOND_GCN(2)

     !UPDATE ICE VELOCITY FIELD
     UICE2  = UICE2F
     VICE2  = VICE2F

     !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
#    if defined (MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,UICE2,VICE2)
#    endif

!         !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
!         IF(K == 4.AND.IOBCN > 0) THEN
!           DO I=1,IOBCN
!             J=I_OBC_N(I)
!             TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
!             UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
!            END DO
!         END IF

   END DO     !! END RUNGE-KUTTA LOOP
  
!==========================================================================
!==========================================================================
      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling
      ! scale to full grid cell
      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling

!      strocnXE=0.0_SP
!      strocnYE=0.0_SP

!      DO I=1,NT 
!       IF(ISICEC(I) == 1)THEN
!       vrel(I) = dragw*sqrt((u(i,1) - uice2(i))**2 +   &
!     &                    (v(i,1) - vice2(i))**2)  ! m/s
!        strocnxE(I) = vrel(I)*((u(i,1)*cosw - v(i,1)*sinw)   &
!     &               - (uice2(i)*cosw - vice2(i)*sinw))
!        strocnyE(I) = vrel(I)*((v(i,1)*cosw + u(i,1)*sinw)   &
!     &               - (vice2(i)*cosw + uice2(i)*sinw))
!       END IF
!      END DO

!==========================================================================
!==========================================================================
    
   DEALLOCATE(AICETMP)
   DEALLOCATE(UICE2F,VICE2F)
   DEALLOCATE(UICE2RK,VICE2RK)
   DEALLOCATE(PSTXICE,PSTYICE)
   DEALLOCATE(ADVSTRX,ADVSTRY)
   DEALLOCATE(ADVUICE,ADVVICE)
   DEALLOCATE(IMASS1,IMASS,PICE)
   DEALLOCATE(TAUXWI,TAUYWI,vrel)
   DEALLOCATE(STRAIRXE,STRAIRYE)
   
   RETURN
   END SUBROUTINE ICE_RUNGE_KUTTA
!==========================================================================
!
!==========================================================================   
!  CALCULATE FLUXES OF FREE SURFACE ELEVATION (CONTINUITY) EQUATION            |
!==============================================================================|
   SUBROUTINE ICE_STRESS(K)       
!------------------------------------------------------------------------------|
!   USE BCS
!   USE MOD_OBCS
   IMPLICIT NONE

   REAL(SP), ALLOCATABLE :: PUPX(:),PVPX(:),PUPY(:),PVPY(:)
   REAL(SP), ALLOCATABLE :: STRESSX1(:),STRESSX2(:),STRESSX3(:),STRESSX4(:),STRESSX5(:)
   REAL(SP), ALLOCATABLE :: STRESSY1(:),STRESSY2(:),STRESSY3(:),STRESSY4(:),STRESSY5(:)
   REAL(SP), ALLOCATABLE :: PS1PX(:),PS2PX(:),PS3PX(:),PS4PX(:),PS5PX(:)
   REAL(SP), ALLOCATABLE :: PS1PY(:),PS2PY(:),PS3PY(:),PS4PY(:),PS5PY(:)

   REAL(SP) :: ELIJ,UIJ,VIJ,EXFLUX
   INTEGER  :: I,K,I1,IA,IB,J1,J2
   
   REAL(SP) :: E11,E12,E22,DD,DT,DS,DELT,ZETA,ETA
   REAL(SP) :: STRESSX1IJ,STRESSX2IJ,STRESSX3IJ,STRESSX4IJ,STRESSX5IJ
   REAL(SP) :: STRESSY1IJ,STRESSY2IJ,STRESSY3IJ,STRESSY4IJ,STRESSY5IJ

!#  if defined (SPHERICAL)
   REAL(SP), ALLOCATABLE :: UICE2N(:), VICE2N(:)
   REAL(SP) :: XTMP,XTMP1
!#  endif        
!==============================================================================|

!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   ALLOCATE(PUPX(0:MT))     ;PUPX= 0.0_SP
   ALLOCATE(PVPX(0:MT))     ;PVPX= 0.0_SP
   ALLOCATE(PUPY(0:MT))     ;PUPY= 0.0_SP
   ALLOCATE(PVPY(0:MT))     ;PVPY= 0.0_SP
   ALLOCATE(STRESSX1(0:MT)) ;STRESSX1 = 0.0_SP
   ALLOCATE(STRESSX2(0:MT)) ;STRESSX2 = 0.0_SP
   ALLOCATE(STRESSX3(0:MT)) ;STRESSX3 = 0.0_SP
   ALLOCATE(STRESSX4(0:MT)) ;STRESSX4 = 0.0_SP
   ALLOCATE(STRESSX5(0:MT)) ;STRESSX5 = 0.0_SP

   ALLOCATE(STRESSY1(0:MT)) ;STRESSY1 = 0.0_SP
   ALLOCATE(STRESSY2(0:MT)) ;STRESSY2 = 0.0_SP
   ALLOCATE(STRESSY3(0:MT)) ;STRESSY3 = 0.0_SP
   ALLOCATE(STRESSY4(0:MT)) ;STRESSY4 = 0.0_SP
   ALLOCATE(STRESSY5(0:MT)) ;STRESSY5 = 0.0_SP

   ALLOCATE(PS1PX(0:NT))    ;PS1PX =0.0_SP
   ALLOCATE(PS2PX(0:NT))    ;PS2PX =0.0_SP
   ALLOCATE(PS3PX(0:NT))    ;PS3PX =0.0_SP
   ALLOCATE(PS4PX(0:NT))    ;PS4PX =0.0_SP
   ALLOCATE(PS5PX(0:NT))    ;PS5PX =0.0_SP

   ALLOCATE(PS1PY(0:NT))    ;PS1PY =0.0_SP
   ALLOCATE(PS2PY(0:NT))    ;PS2PY =0.0_SP
   ALLOCATE(PS3PY(0:NT))    ;PS3PY =0.0_SP
   ALLOCATE(PS4PY(0:NT))    ;PS4PY =0.0_SP
   ALLOCATE(PS5PY(0:NT))    ;PS5PY =0.0_SP

   Delta= 0.0_SP
   divu = 0.0_SP
   
!---------ACCUMULATE STRESS BY LOOPING OVER CONTROL VOLUME HALF EDGES----------!
!#  if defined (SPHERICAL)
    ALLOCATE(UICE2N(0:MT)); UICE2N =0.0_SP
    ALLOCATE(VICE2N(0:MT)); VICE2N =0.0_SP

   DO I=1,M
      DO I1=1, NTVE(I)
          UICE2N(I) =UICE2N(I)+UICE2(NBVE(I,I1))
          VICE2N(I) =VICE2N(I)+VICE2(NBVE(I,I1))
      END DO
      UICE2N(I) = UICE2N(I) /NTVE(I)
      VICE2N(I) = VICE2N(I) /NTVE(I)
   END DO
!#  endif
# if defined(MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,UICE2N,VICE2N)
# endif

   DO I=1,NCV_I
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)

     IF(ISICEN(IA) == 1 .OR. ISICEN(IB) == 1)THEN
       UIJ = UICE2(I1)
       VIJ = VICE2(I1)

       EXFLUX = -UIJ*DLTYE(I)                  
       PUPX(IA) = PUPX(IA)-EXFLUX
       PUPX(IB) = PUPX(IB)+EXFLUX

       EXFLUX = -VIJ*DLTYE(I)                  
       PVPX(IA) = PVPX(IA)-EXFLUX
       PVPX(IB) = PVPX(IB)+EXFLUX

       EXFLUX =  UIJ*DLTXE(I)  
       PUPY(IA) = PUPY(IA)-EXFLUX
       PUPY(IB) = PUPY(IB)+EXFLUX

       EXFLUX =  VIJ*DLTXE(I)  
       PVPY(IA) = PVPY(IA)-EXFLUX
       PVPY(IB) = PVPY(IB)+EXFLUX

#  if defined (SPHERICAL)
       PUPY(IA) = PUPY(IA) +UICE2N(IA)*TAN(DEG2RAD*YC(I1))/REARTH
       PVPY(IA) = PVPY(IA) +VICE2N(IA)*TAN(DEG2RAD*YC(I1))/REARTH
       PUPY(IB) = PUPY(IB) +UICE2N(IB)*TAN(DEG2RAD*YC(I1))/REARTH
       PVPY(IB) = PVPY(IB) +VICE2N(IB)*TAN(DEG2RAD*YC(I1))/REARTH
!! ggao change for test 1115/2007

#  endif
     END IF  
   END DO

!#  if defined (SPHERICAL)
    DEALLOCATE(UICE2N,VICE2N) 
!#  endif

# if defined(MULTIPROCESSOR)
!     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,PUPX,PVPX)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,PUPX,PVPX)
!     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,PUPY,PVPY)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,PUPY,PVPY)
# endif

   PUPX = PUPX/ART1
   PVPX = PVPX/ART1
   PUPY = PUPY/ART1
   PVPY = PVPY/ART1

   DO I=1,M
     IF(ISICEN(I) == 1)THEN
       E11 = PUPX(I)
       E12 = 0.5_SP*(PUPY(I)+PVPX(I))
       E22 = PVPY(I)
   
       DD = E11+E22
       DT = E11-E22
       DS = 2.0_SP*E12
   
       DELT = DD*DD+ECCI*(DT*DT+DS*DS)
       DELT = SQRT(DELT)

!--------------------------------------------------  
!       calculate stress for ridging
!--------------------------------------------------  
       IF(K==4) THEN
       Delta (I,1)= DELT 
       divu  (I,1)= DD
       END IF
!--------------------------------------------------  
!--------------------------------------------------  
   
!       ZETA = 0.5_SP*PICE(I)/DELT
!       ETA  = 0.5_SP*PICE(I)/(DELT*ECC)

       !ice_grid.F:     tinyarea(i,j) = puny*tarea(i,j)
       !ice_grid.F:!echmod   tinyarea(i,j) = 1.e-7*tarea(i,j) ! for automatic differenti
       !ice_dyn_evp.F:     &,  rcon         ! for damping criterion (kg/s)
       !ice_dyn_evp.F:      rcon(:,:) = 1230._dbl_kind*eyc*dyn_dt*dtei**2  ! kg/s
       !if (evp_damping) then

       ! enforce damping criterion
       !--------------------------------------------------  
        ZETA = 0.5_SP*min(PICE(I)/max(DELT,1.e-7*art1(i)),rcon)
        !ETA  = 0.5_SP*PICE(I)/(DELT*ECC)

       !--------------------------------------------------  
       !   original version
        !ZETA = 0.5_SP*PICE(i)/max(DELT,1.e-7*art1(i))
        ETA  = ZETA*ECCI
       !--------------------------------------------------  
   
       STRESSX1(I) = (ETA+ZETA)*PUPX(I)
       STRESSX2(I) = ETA*PUPY(I)
       STRESSX3(I) = (ZETA-ETA)*PVPY(I)
       STRESSX4(I) = ETA*PVPX(I)
       STRESSX5(I) = 0.5_SP*PICE(I)
   
       STRESSY1(I) = (ETA+ZETA)*PVPY(I)
       STRESSY2(I) = ETA*PVPX(I)
       STRESSY3(I) = (ZETA-ETA)*PUPX(I)
       STRESSY4(I) = ETA*PUPY(I)
       STRESSY5(I) = 0.5_SP*PICE(I)
     END IF  
   END DO  

# if defined(MULTIPROCESSOR)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,STRESSX1,STRESSX2)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,STRESSX3,STRESSX4)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,STRESSX5)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,STRESSX1,STRESSX2)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,STRESSX3,STRESSX4)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,STRESSX5)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,STRESSY1,STRESSY2)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,STRESSY3,STRESSY4)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,STRESSy5)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,STRESSY1,STRESSY2)
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,STRESSY3,STRESSY4,STRESSY5)
# endif

   PS1PX = 0.0_SP
   PS2PX = 0.0_SP
   PS3PX = 0.0_SP
   PS4PX = 0.0_SP
   PS5PX = 0.0_SP

   PS1PY = 0.0_SP
   PS2PY = 0.0_SP
   PS3PY = 0.0_SP
   PS4PY = 0.0_SP
   PS5PY = 0.0_SP

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     
     IF(ISICEC(IA) == 1 .OR. ISICEC(IB) == 1)THEN
       ELIJ=0.5_SP*(EL(J1)+EL(J2))
       !  ggao 0516
       ELIJ=0.0    
 
       STRESSX1IJ=0.5_SP*(STRESSX1(J1)+STRESSX1(J2))
       STRESSX2IJ=0.5_SP*(STRESSX2(J1)+STRESSX2(J2))
       STRESSX3IJ=0.5_SP*(STRESSX3(J1)+STRESSX3(J2))
       STRESSX4IJ=0.5_SP*(STRESSX4(J1)+STRESSX4(J2))

       STRESSX5IJ=0.5_SP*(STRESSX5(J1)+STRESSX5(J2))

       STRESSY1IJ=0.5_SP*(STRESSY1(J1)+STRESSY1(J2))
       STRESSY2IJ=0.5_SP*(STRESSY2(J1)+STRESSY2(J2))
       STRESSY3IJ=0.5_SP*(STRESSY3(J1)+STRESSY3(J2))
       STRESSY4IJ=0.5_SP*(STRESSY4(J1)+STRESSY4(J2))
       STRESSY5IJ=0.5_SP*(STRESSY5(J1)+STRESSY5(J2))

#      if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o latitude         
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF  

       EXFLUX    = STRESSX1IJ*DLTYC(I)
       PS1PX(IA) = PS1PX(IA) - EXFLUX/ART(IA)
       PS1PX(IB) = PS1PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSX2IJ*XTMP*COS(DEG2RAD*YC(IA))
       PS2PY(IA) = PS2PY(IA) + EXFLUX/ART(IA)
       EXFLUX    = STRESSX2IJ*XTMP*COS(DEG2RAD*YC(IB))
       PS2PY(IB) = PS2PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSX3IJ*DLTYC(I)
       PS3PX(IA) = PS3PX(IA) - EXFLUX/ART(IA)
       PS3PX(IB) = PS3PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSX4IJ*XTMP*COS(DEG2RAD*YC(IA))
       PS4PY(IA) = PS4PY(IA) + EXFLUX/ART(IA)
       EXFLUX    = STRESSX4IJ*XTMP*COS(DEG2RAD*YC(IB))
       PS4PY(IB) = PS4PY(IB) - EXFLUX/ART(IB)

!       EXFLUX    = STRESSX5IJ*DLTYC(I)
!       PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)
!       PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = DLTYC(I)
       PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)*STRESSX5IJ
       PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)*STRESSX5IJ
!-------------------------------------------------------------
!-------------------------------------------------------------
       EXFLUX    = STRESSY1IJ*XTMP*COS(DEG2RAD*YC(IA))
       PS1PY(IA) = PS1PY(IA) + EXFLUX/ART(IA)
       EXFLUX    = STRESSY1IJ*XTMP*COS(DEG2RAD*YC(IB))
       PS1PY(IB) = PS1PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSY2IJ*DLTYC(I)
       PS2PX(IA) = PS2PX(IA) - EXFLUX/ART(IA)
       PS2PX(IB) = PS2PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSY3IJ*XTMP*COS(DEG2RAD*YC(IA))
       PS3PY(IA) = PS3PY(IA) + EXFLUX/ART(IA)
       EXFLUX    = STRESSY3IJ*XTMP*COS(DEG2RAD*YC(IB))
       PS3PY(IB) = PS3PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSY4IJ*DLTYC(I)
       PS4PX(IA) = PS4PX(IA) - EXFLUX/ART(IA)
       PS4PX(IB) = PS4PX(IB) + EXFLUX/ART(IB)

!       EXFLUX    = STRESSY5IJ*XTMP*COS(DEG2RAD*YC(IA))
!       PS5PY(IA) = PS5PY(IA) - EXFLUX/ART(IA)
!       EXFLUX    = STRESSY5IJ*XTMP*COS(DEG2RAD*YC(IB))
!       PS5PY(IB) = PS5PY(IB) + EXFLUX/ART(IB)

       EXFLUX    = XTMP*COS(DEG2RAD*YC(IA))
       PS5PY(IA) = PS5PY(IA) + EXFLUX/ART(IA)*STRESSY5IJ
       EXFLUX    = XTMP*COS(DEG2RAD*YC(IB))
       PS5PY(IB) = PS5PY(IB) - EXFLUX/ART(IB)*STRESSY5IJ

!      ACCUMULATE BAROTROPIC FLUX
       PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*DLTYC(I)/ART(IA)
       PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*XTMP*COS(DEG2RAD*YC(IA))/ART(IA)

       PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*DLTYC(I)/ART(IB)
       PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*XTMP*COS(DEG2RAD*YC(IB))/ART(IB)

#      else
       EXFLUX    = STRESSX1IJ*DLTYC(I)
       PS1PX(IA) = PS1PX(IA) - EXFLUX/ART(IA)
       PS1PX(IB) = PS1PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSX2IJ*DLTXC(I)
       PS2PY(IA) = PS2PY(IA) + EXFLUX/ART(IA)
       PS2PY(IB) = PS2PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSX3IJ*DLTYC(I)
       PS3PX(IA) = PS3PX(IA) - EXFLUX/ART(IA)
       PS3PX(IB) = PS3PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSX4IJ*DLTXC(I)
       PS4PY(IA) = PS4PY(IA) + EXFLUX/ART(IA)
       PS4PY(IB) = PS4PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSX5IJ*DLTYC(I)
       PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)
       PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSY1IJ*DLTXC(I)
       PS1PY(IA) = PS1PY(IA) + EXFLUX/ART(IA)
       PS1PY(IB) = PS1PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSY2IJ*DLTYC(I)
       PS2PX(IA) = PS2PX(IA) - EXFLUX/ART(IA)
       PS2PX(IB) = PS2PX(IB) + EXFLUX/ART(IB)

       EXFLUX    = STRESSY3IJ*DLTXC(I)
       PS3PY(IA) = PS3PY(IA) + EXFLUX/ART(IA)
       PS3PY(IB) = PS3PY(IB) - EXFLUX/ART(IB)

       EXFLUX    = STRESSY4IJ*DLTYC(I)
       PS4PX(IA) = PS4PX(IA) - EXFLUX/ART(IA)
       PS4PX(IB) = PS4PX(IB) + EXFLUX/ART(IB)

!       EXFLUX    = STRESSY5IJ*DLTXC(I)
!       PS5PY(IA) = PS5PY(IA) - EXFLUX/ART(IA)
!       PS5PY(IB) = PS5PY(IB) + EXFLUX/ART(IB)

       EXFLUX    = DLTXC(I)
       PS5PY(IA) = PS5PY(IA) + EXFLUX/ART(IA)*STRESSY5IJ
       PS5PY(IB) = PS5PY(IB) - EXFLUX/ART(IB)*STRESSY5IJ

!      ACCUMULATE BAROTROPIC FLUX
       PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*DLTYC(I)/ART(IA)
       PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*DLTXC(I)/ART(IA)
       PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*DLTYC(I)/ART(IB)
       PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*DLTXC(I)/ART(IB)

#      endif     
     END IF
   END DO
   
#  if defined (SPHERICAL)
#  if defined (NORTHPOLE)
   CALL ICE_STRESS_XY(PS1PX,PS2PX,PS3PX,PS4PX,PS5PX,       &
                      PS1PY,PS2PY,PS3PY,PS4PY,PS5PY,K      &        
      ,STRESSX1,STRESSX2,STRESSX3,STRESSX4,STRESSX5        &
      ,STRESSY1,STRESSY2,STRESSY3,STRESSY4,STRESSY5)

#  endif
#  endif

   ADVSTRX=0.0_SP
   ADVSTRY=0.0_SP


!! this part need to double check  the sign
!!  ICE_STRESS_XY
!!  ggao 2007--1003/
!!  change end
   DO I=1,N
     IF(sum(ISICEN(NV(I,1:3)))==3) THEN
     ADVSTRX(I) = PS1PX(I)+PS2PY(I)+PS3PX(I)+PS4PY(I) -PS5PX(I)
     ADVSTRY(I) = PS1PY(I)+PS2PX(I)+PS3PY(I)+PS4PX(I) -PS5PY(I)
     END IF
   END DO

# if defined(MULTIPROCESSOR)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ADVSTRX,ADVSTRY)
# endif


   DEALLOCATE(PUPX) 
   DEALLOCATE(PVPX) 
   DEALLOCATE(PUPY) 
   DEALLOCATE(PVPY)
   DEALLOCATE(STRESSX1)
   DEALLOCATE(STRESSX2)
   DEALLOCATE(STRESSX3)
   DEALLOCATE(STRESSX4)
   DEALLOCATE(STRESSX5)
   DEALLOCATE(STRESSY1)
   DEALLOCATE(STRESSY2)
   DEALLOCATE(STRESSY3)
   DEALLOCATE(STRESSY4)
   DEALLOCATE(STRESSY5)

   DEALLOCATE(PS1PX)  
   DEALLOCATE(PS2PX)  
   DEALLOCATE(PS3PX)  
   DEALLOCATE(PS4PX)  
   DEALLOCATE(PS5PX)

   DEALLOCATE(PS1PY)    
   DEALLOCATE(PS2PY)   
   DEALLOCATE(PS3PY)   
   DEALLOCATE(PS4PY)   
   DEALLOCATE(PS5PY)   


   RETURN
   END SUBROUTINE ICE_STRESS
!==============================================================================|
!
!==============================================================================|
   SUBROUTINE ICEOCEAN_STRESS
   
   IMPLICIT NONE
   
!   REAL(SP), PARAMETER :: DRAGW = 0.00536_SP*RHOW 
!                                ! drag coefficient for water on ice *rhow (kg/m^3)
!   REAL(SP), PARAMETER :: COSW = C1I  !cos(ocean turning angle) !turning angle = 0
!   REAL(SP), PARAMETER :: SINW = C0I  !sin(ocean turning angle) !turning angle = 0
   REAL(SP) :: WATERX,WATERY
   INTEGER  :: I

   DO I=1,NT
     IF(ISICEC(I) == 1)THEN
       WATERX = U(I,1)*COSW - V(I,1)*SINW
       WATERY = V(I,1)*COSW + U(I,1)*SINW
       ! (magnitude of relative ocean current)*rhow*drag*aice
       VREL(I) = AIU(i)*DRAGW*SQRT((U(I,1)-UICE2(I))**2+(V(I,1)-VICE2(I))**2)  ! m/s
       ! ice/ocean stress
       TAUXWI(I) = VREL(I)*WATERX ! NOTE this is not the entire
       TAUYWI(I) = VREL(I)*WATERY ! ocn stress term

     IF(ISBCE(I) == 1) THEN
!       TAUXWI(I) = 0.0 ! NOTE this is not the entire
!       TAUYWI(I) = 0.0 ! ocn stress term
     END IF


     END IF 
   END DO  

# if defined(MULTIPROCESSOR)
!     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,TAUXWI,TAUYWI)
# endif

   RETURN
   END SUBROUTINE ICEOCEAN_STRESS
!==============================================================================|
!
!==============================================================================|
!     ACCUMLATE FLUXES FOR ICE                                                 |
!==============================================================================|

   SUBROUTINE ICE_UV(K)       

!==============================================================================|

#  if defined (SPHERICAL)
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif   
#  endif   

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESXICE,RESYICE
   REAL(SP) :: UICE2FT,VICE2FT
   REAL(SP) :: UI,VI
   INTEGER  :: I,J
   real(SP) :: WIND_SP
   real(SP) :: cca,ccb,ab2,cc1,cc2
   real(SP) :: TAUX, TAUY

!==============================================================================|
!
!--ACCUMULATE RESIDUALS FOR ICE EQUATIONS--------------------------------------|
!

   UICE2FT = UICE2F(0)
   VICE2FT = VICE2F(0)

   DO I=1,N
     IF(ISICEC(I) == 1)THEN
!       RESXICE(I) = ADVUICE(I)-ADVSTRX(I)+IMASS1(I)*(PSTXICE(I)    &
!                   -COR(I)*VICE2(I)*ART(I))                               &
!                   -(TAUXWI(I)+STRAIRXE(I))*ART(I)
!       RESYICE(I) = ADVVICE(I)-ADVSTRY(I)+IMASS1(I)*(PSTYICE(I)           &
!                   +COR(I)*UICE2(I)*ART(I))                               &
!                   -(TAUYWI(I)+STRAIRYE(I))*ART(I)

!       RESXICE(I) = ADVUICE(I)/ART(I)-ADVSTRX(I)+IMASS1(I)*(PSTXICE(I)    &
!                   -COR(I)*VICE2(I))-(TAUXWI(I)+STRAIRXE(I))
!       RESYICE(I) = ADVVICE(I)/ART(I)-ADVSTRY(I)+IMASS1(I)*(PSTYICE(I)           &
!                   +COR(I)*UICE2(I))-(TAUYWI(I)+STRAIRYE(I))

#      if defined (SPHERICAL)
       IF(trim(RHEOLOGY) == 'VP')THEN
!         RESXICE(I) = RESXICE(I)-UICE2(I)*VICE2(I)/REARTH*TAN(DEG2RAD*YC(I))
!         RESYICE(I) = RESYICE(I)+UICE2(I)*UICE2(I)/REARTH*TAN(DEG2RAD*YC(I))
       END IF
#      endif
!
!--UPDATE----------------------------------------------------------------------|

        ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
!        cca = imass1(i)/dtice + vrel * cosw         ! alpha, kg/m^2 s
        cca = imass1(i)/dyn_dt + vrel(I) * cosw*ALPHA_RK(K) ! alpha, kg/m^2 s
!        ccb = fm(i,j)        + vrel * sinw        ! beta,  kg/m^2 s
        ccb = (cor(i)*imass1(i) + vrel(I) * sinw)*ALPHA_RK(K)  ! beta,  kg/m^2 s
        ab2 = cca**2 + ccb**2

        ! divergence of the internal stress tensor
        !  ADVSTRX(I),ADVSTRY(I)

        ! finally, the velocity components
        cc1 = (ADVSTRX(I) - IMASS1(I)*PSTXICE(I) + TAUXWI(I) +STRAIRXE(I))&
     &       *ALPHA_RK(K) + imass1(I)/dyn_dt*UICE2RK(I)
        cc2 = (ADVSTRY(I) - IMASS1(I)*PSTYICE(I) + TAUYWI(I) +STRAIRYE(I))&
     &        *ALPHA_RK(K)+ imass1(I)/dyn_dt*VICE2RK(I)

       UICE2F(I) = (cca*cc1 + ccb*cc2)/ab2              ! m/s
       VICE2F(I) = (cca*cc2 - ccb*cc1)/ab2
     END IF  !  end isicec   
   END DO

#  if defined (SPHERICAL)
#  if defined (NORTHPOLE)
   CALL ICE_UV_XY(K)
#  endif
#  endif

   DO I=1,N
     IF(ISBCE(I) == 1) THEN
       UI= UICE2F(I)*(SIN(ALPHA(I)))**2-VICE2F(I)*SIN(ALPHA(I))*COS(ALPHA(I))
       VI=-UICE2F(I)*SIN(ALPHA(I))*COS(ALPHA(I))+VICE2F(I)*(COS(ALPHA(I)))**2
!  ggao 0908/2007
!      UICE2F(I)=UI
!      VICE2F(I)=VI
     END IF  
     UICE2F(I) =UICE2F(I)*ISICEC(I)
     VICE2F(I) =VICE2F(I)*ISICEC(I)
   END DO

   VICE2F(0) = VICE2FT
   UICE2F(0) = UICE2FT

   RETURN
   END SUBROUTINE ICE_UV
!==============================================================================|
!
!==============================================================================|
!   CALCULATE CONVECTION FLUXES FOR ICE                                        !
!==============================================================================|
   SUBROUTINE ADVAVE_ICE(XFLUX,YFLUX)
!==============================================================================|

#  if defined (SPHERICAL)
   USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
   USE MOD_NORTHPOLE
#  endif   
#  endif
!   USE BCS
!   USE MOD_OBCS

   IMPLICIT NONE
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3
   REAL(SP) :: XIJ,YIJ,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,UN_TMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)

#  if defined (SPHERICAL)
   REAL(DP) :: XTMP,XTMP1
#  endif      

!
!-------------------------INITIALIZE FLUXES------------------------------------!
!
   XFLUX = 0.0_SP
   YFLUX = 0.0_SP
!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

!     IF(ISICEC(IA) == 1 .OR. ISICEC(IB) == 1)THEN
!      FLUX FROM LEFT
       K1=NBE(IA,1)
       K2=NBE(IA,2)
       K3=NBE(IA,3)
         
       COFA1=A1U(IA,1)*UICE2(IA)+A1U(IA,2)*UICE2(K1)+A1U(IA,3)*UICE2(K2)+A1U(IA,4)*UICE2(K3)
       COFA2=A2U(IA,1)*UICE2(IA)+A2U(IA,2)*UICE2(K1)+A2U(IA,3)*UICE2(K2)+A2U(IA,4)*UICE2(K3)
       COFA5=A1U(IA,1)*VICE2(IA)+A1U(IA,2)*VICE2(K1)+A1U(IA,3)*VICE2(K2)+A1U(IA,4)*VICE2(K3)
       COFA6=A2U(IA,1)*VICE2(IA)+A2U(IA,2)*VICE2(K1)+A2U(IA,3)*VICE2(K2)+A2U(IA,4)*VICE2(K3)
     
#      if defined (SPHERICAL)
       UIJ1=UICE2(IA)+COFA1*DLTXNE(I,1)+COFA2*DLTYNE(I,1)
       VIJ1=VICE2(IA)+COFA5*DLTXNE(I,1)+COFA6*DLTYNE(I,1)
#      else
       XIJ=XIJC(I)-XC(IA)
       YIJ=YIJC(I)-YC(IA)
       UIJ1=UICE2(IA)+COFA1*XIJ+COFA2*YIJ
       VIJ1=VICE2(IA)+COFA5*XIJ+COFA6*YIJ
#      endif

!      FLUX FROM RIGHT
       K1=NBE(IB,1)
       K2=NBE(IB,2)
       K3=NBE(IB,3)
          
       COFA3=A1U(IB,1)*UICE2(IB)+A1U(IB,2)*UICE2(K1)+A1U(IB,3)*UICE2(K2)+A1U(IB,4)*UICE2(K3)
       COFA4=A2U(IB,1)*UICE2(IB)+A2U(IB,2)*UICE2(K1)+A2U(IB,3)*UICE2(K2)+A2U(IB,4)*UICE2(K3)
       COFA7=A1U(IB,1)*VICE2(IB)+A1U(IB,2)*VICE2(K1)+A1U(IB,3)*VICE2(K2)+A1U(IB,4)*VICE2(K3)
       COFA8=A2U(IB,1)*VICE2(IB)+A2U(IB,2)*VICE2(K1)+A2U(IB,3)*VICE2(K2)+A2U(IB,4)*VICE2(K3)
     
#      if defined (SPHERICAL)
       UIJ2=UICE2(IB)+COFA3*DLTXNE(I,2)+COFA4*DLTYNE(I,2)
       VIJ2=VICE2(IB)+COFA7*DLTXNE(I,2)+COFA8*DLTYNE(I,2)
#      else
       XIJ=XIJC(I)-XC(IB)
       YIJ=YIJC(I)-YC(IB)
       UIJ2=UICE2(IB)+COFA3*XIJ+COFA4*YIJ
       VIJ2=VICE2(IB)+COFA7*XIJ+COFA8*YIJ
#      endif

!      NORMAL VELOCITY
       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       UN_TMP=-UIJ*DLTYC(I) + VIJ*DLTXC(I)

!      ADD CONVECTIVE AND VISCOUS FLUXES
       XADV=UN_TMP*&
            ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
       YADV=UN_TMP* &
            ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP

!      ACCUMULATE FLUX
       XFLUX(IA)=XFLUX(IA)+XADV*(1.0_SP-ISBC(I))
       YFLUX(IA)=YFLUX(IA)+YADV*(1.0_SP-ISBC(I))
       XFLUX(IB)=XFLUX(IB)-XADV*(1.0_SP-ISBC(I))
       YFLUX(IB)=YFLUX(IB)-YADV*(1.0_SP-ISBC(I))
!     END IF  
   END DO

#  if defined (SPHERICAL)
#  if defined (NORTHPOLE)
   CALL ADVAVE_ICE_XY(XFLUX,YFLUX)
#  endif 
#  endif  

!
!-------------------------SET BOUNDARY VALUES----------------------------------!
!

!  MODIFY BOUNDARY FLUX
!      DO I=1,N
!        IF(ISBCE(I) == 2) THEN
!          XFLUX(I)=(XFLUX(I)+Fluxobn(I)*UICE2(I))*IUCP(I)
!          YFLUX(I)=(YFLUX(I)+Fluxobn(I)*VICE2(I))*IUCP(I)
!        ENDIF
!      END DO

   RETURN
   END SUBROUTINE ADVAVE_ICE
!==============================================================================|


#  if defined (NORTHPOLE)
!
!==============================================================================|
   SUBROUTINE ICE_STRESS_XY(PS1PX,PS2PX,PS3PX,PS4PX,PS5PX,            &
                            PS1PY,PS2PY,PS3PY,PS4PY,PS5PY,K           &
                    ,STRESSX1,STRESSX2,STRESSX3,STRESSX4,STRESSX5     &
                    ,STRESSY1,STRESSY2,STRESSY3,STRESSY4,STRESSY5)


!------------------------------------------------------------------------------|
!   USE BCS
!   USE MOD_OBCS

   IMPLICIT NONE

   REAL(SP), DIMENSION(0:MT) :: PUPX,PVPX,PUPY,PVPY
   REAL(SP), DIMENSION(0:MT) :: STRESSX1,STRESSX2,STRESSX3,STRESSX4,STRESSX5
   REAL(SP), DIMENSION(0:MT) :: STRESSY1,STRESSY2,STRESSY3,STRESSY4,STRESSY5
   REAL(SP), DIMENSION(0:NT) :: PS1PX,PS2PX,PS3PX,PS4PX,PS5PX
   REAL(SP), DIMENSION(0:NT) :: PS1PY,PS2PY,PS3PY,PS4PY,PS5PY

   REAL(DP) :: ELIJ,UIJ,VIJ,EXFLUX
   INTEGER  :: I,II,K,I1,IA,IB,J1,J2
   
   REAL(SP) :: E11,E12,E22,DD,DT,DS,DELT,ZETA,ETA
   REAL(SP) :: STRESSX1IJ,STRESSX2IJ,STRESSX3IJ,STRESSX4IJ,STRESSX5IJ
   REAL(SP) :: STRESSY1IJ,STRESSY2IJ,STRESSY3IJ,STRESSY4IJ,STRESSY5IJ
   
   REAL(DP) :: UIJ_TMP,VIJ_TMP
   REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
   REAL(DP) :: DLTXE_TMP,DLTYE_TMP,DLTXC_TMP,DLTYC_TMP
   REAL(DP) :: STXJ1,STYJ1,STXJ2,STYJ2

!==============================================================================|

!----------INITIALIZE STRESS ARRAY ----------------------------------------------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       PUPX(IA) = 0.0_SP
       PVPX(IA) = 0.0_SP
       PUPY(IA) = 0.0_SP
       PVPY(IA) = 0.0_SP
     END IF
     IF(IB == NODE_NORTHPOLE)THEN  
       PUPX(IB) = 0.0_SP
       PVPX(IB) = 0.0_SP
       PUPY(IB) = 0.0_SP
       PVPY(IB) = 0.0_SP
     END IF  
   END DO  

!---------ACCUMULATE STRESS BY LOOPING OVER CONTROL VOLUME HALF EDGES----------!

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)

     IF(ISICEN(IA) == 1 .OR. ISICEN(IB) == 1)THEN
       UIJ = UICE2(I1)
       VIJ = VICE2(I1)

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -VIJ*COS(XC(I1)*PI/180.)-UIJ*SIN(XC(I1)*PI/180.)
         VIJ_TMP = -VIJ*SIN(XC(I1)*PI/180.)+UIJ*COS(XC(I1)*PI/180.)
       
         VX1_TMP = REARTH * DCOS(YIJE(I,1)*PI/180.0_DP) * DCOS(XIJE(I,1)*PI/180.0_DP)&
                   * 2._DP /(1._DP+Dsin(YIJE(I,1)*PI/180.0_DP))
         VY1_TMP = REARTH * DCOS(YIJE(I,1)*PI/180.0_DP) * DSIN(XIJE(I,1)*PI/180.0_DP)&
                   * 2._DP /(1._DP+Dsin(YIJE(I,1)*PI/180.0_DP))

         VX2_TMP = REARTH * DCOS(YIJE(I,2)*PI/180.0_DP) * DCOS(XIJE(I,2)*PI/180.0_DP)&
                   * 2._DP /(1._DP+Dsin(YIJE(I,2)*PI/180.0_DP))
         VY2_TMP = REARTH * DCOS(YIJE(I,2)*PI/180.0_DP) * DSIN(XIJE(I,2)*PI/180.0_DP)&
                   * 2._DP /(1._DP+Dsin(YIJE(I,2)*PI/180.0_DP))

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP

         EXFLUX = -UIJ_TMP*DLTYE_TMP                  
         IF(IA == NODE_NORTHPOLE) THEN    
           PUPX(IA) = PUPX(IA)-EXFLUX/ART1(IA)
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           PUPX(IB) = PUPX(IB)+EXFLUX/ART1(IB)
         END IF    

         EXFLUX = -VIJ_TMP*DLTYE_TMP                  
         IF(IA == NODE_NORTHPOLE) THEN    
           PVPX(IA) = PVPX(IA)-EXFLUX/ART1(IA)
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           PVPX(IB) = PVPX(IB)+EXFLUX/ART1(IB)
         END IF    

         EXFLUX =  UIJ_TMP*DLTXE_TMP  
         IF(IA == NODE_NORTHPOLE) THEN    
           PUPY(IA) = PUPY(IA)-EXFLUX/ART1(IA)
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           PUPY(IB) = PUPY(IB)+EXFLUX/ART1(IB)
         END IF    

         EXFLUX =  VIJ_TMP*DLTXE_TMP  
         IF(IA == NODE_NORTHPOLE) THEN    
           PVPY(IA) = PVPY(IA)-EXFLUX/ART1(IA)
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           PVPY(IB) = PVPY(IB)+EXFLUX/ART1(IB)
         END IF 
       END IF  
     END IF 
   END DO

   I = NODE_NORTHPOLE
   IF(I == 0) RETURN
   IF(ISICEN(I) == 0)RETURN
   
   E11 = PUPX(I)
   E12 = 0.5_SP*(PUPY(I)+PVPX(I))
   E22 = PVPY(I)
   
   DD = E11+E22
   DT = E11-E22
   DS = 2.0_SP*E12
   
   DELT = DD*DD+ECCI*(DT*DT+DS*DS)
   DELT = SQRT(DELT)

!--------------------------------------------------
!       calculate stress for ridging
!--------------------------------------------------
       IF(K==4) THEN
       Delta (I,1)= DELT
       divu  (I,1)= DD
       END IF
!--------------------------------------------------
!--------------------------------------------------

!       ZETA = 0.5_SP*PICE(I)/DELT
!       ETA  = 0.5_SP*PICE(I)/(DELT*ECC)

       !ice_grid.F:     tinyarea(i,j) = puny*tarea(i,j)
       !ice_grid.F:!echmod   tinyarea(i,j) = 1.e-7*tarea(i,j) ! for automatic differenti
       !ice_dyn_evp.F:     &,  rcon         ! for damping criterion (kg/s)
       !ice_dyn_evp.F:      rcon(:,:) = 1230._dbl_kind*eyc*dyn_dt*dtei**2  ! kg/s
!      if (evp_damping) then
!        ! enforce damping criterion
!       ZETA = 0.5_SP*min(PICE(I)/max(DELT,1.e-7*art1(i)),rcon)
       ZETA = 0.5_SP*min(PICE(I)/max(DELT,1.e-7*art1(i)),rcon)
!       ETA  = 0.5_SP*PICE(I)/(DELT*ECC)
!      else
!        ! original version
!        ZETA = 0.5_SP*PICE(i)/max(DELT,1.e-7*art1(i))
        ETA  = ZETA*ECCI
!      end if
  
!   ZETA = 0.5_SP*PICE(I)/DELT
!   ETA  = 0.5_SP*PICE(I)/(DELT*ECC)
   
   STRESSX1(I) = (ETA+ZETA)*PUPX(I)
   STRESSX2(I) = ETA*PUPY(I)
   STRESSX3(I) = (ZETA-ETA)*PVPY(I)
   STRESSX4(I) = ETA*PVPX(I)
   STRESSX5(I) = 0.5_SP*PICE(I)
   
   STRESSY1(I) = (ETA+ZETA)*PVPY(I)
   STRESSY2(I) = ETA*PVPX(I)
   STRESSY3(I) = (ZETA-ETA)*PUPX(I)
   STRESSY4(I) = ETA*PUPY(I)
   STRESSY5(I) = 0.5_SP*PICE(I)

   
   DO II=1,NPE
     I = NPEDGE_LST(II)
     IA = IEC(I,1)
     IB = IEC(I,2)
     
     IF(CELL_NORTHAREA(IA) == 1)THEN
       PS1PX(IA) = 0.0_SP
       PS2PX(IA) = 0.0_SP
       PS3PX(IA) = 0.0_SP
       PS4PX(IA) = 0.0_SP
       PS5PX(IA) = 0.0_SP

       PS1PY(IA) = 0.0_SP
       PS2PY(IA) = 0.0_SP
       PS3PY(IA) = 0.0_SP
       PS4PY(IA) = 0.0_SP
       PS5PY(IA) = 0.0_SP
       
       PSTXICE(IA) = 0.0_SP
       PSTYICE(IA) = 0.0_SP
     END IF  
     IF(CELL_NORTHAREA(IB) == 1)THEN  
       PS1PX(IB) = 0.0_SP
       PS2PX(IB) = 0.0_SP
       PS3PX(IB) = 0.0_SP
       PS4PX(IB) = 0.0_SP
       PS5PX(IB) = 0.0_SP

       PS1PY(IB) = 0.0_SP
       PS2PY(IB) = 0.0_SP
       PS3PY(IB) = 0.0_SP
       PS4PY(IB) = 0.0_SP
       PS5PY(IB) = 0.0_SP
        
       PSTXICE(IB) = 0.0_SP
       PSTYICE(IB) = 0.0_SP
    END IF  
   END DO  
     
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)
     
!     IF(ISICEC(IA) == 1 .OR. ISICEC(IB) == 1)THEN
       ELIJ=0.5_SP*(EL(J1)+EL(J2))

       ! ggao 0516
         ELIJ=0.0

       STRESSX1IJ=0.5_SP*(STRESSX1(J1)+STRESSX1(J2))
       STRESSX2IJ=0.5_SP*(STRESSX2(J1)+STRESSX2(J2))
       STRESSX3IJ=0.5_SP*(STRESSX3(J1)+STRESSX3(J2))
       STRESSX4IJ=0.5_SP*(STRESSX4(J1)+STRESSX4(J2))
       STRESSX5IJ=0.5_SP*(STRESSX5(J1)+STRESSX5(J2))

       STRESSY1IJ=0.5_SP*(STRESSY1(J1)+STRESSY1(J2))
       STRESSY2IJ=0.5_SP*(STRESSY2(J1)+STRESSY2(J2))
       STRESSY3IJ=0.5_SP*(STRESSY3(J1)+STRESSY3(J2))
       STRESSY4IJ=0.5_SP*(STRESSY4(J1)+STRESSY4(J2))
       STRESSY5IJ=0.5_SP*(STRESSY5(J1)+STRESSY5(J2))
   

       VX1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_DP) * COS(VX(IENODE(I,1))*PI/180.0_DP) &
                 * 2._DP /(1._DP+sin(VY(IENODE(I,1))*PI/180.0_DP))
       VY1_TMP = REARTH * COS(VY(IENODE(I,1))*PI/180.0_DP) * SIN(VX(IENODE(I,1))*PI/180.0_DP) &
                 * 2._DP /(1._DP+sin(VY(IENODE(I,1))*PI/180.0_DP))

       VX2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_DP) * COS(VX(IENODE(I,2))*PI/180.0_DP) &
                 * 2._DP /(1._DP+sin(VY(IENODE(I,2))*PI/180.0_DP))
       VY2_TMP = REARTH * COS(VY(IENODE(I,2))*PI/180.0_DP) * SIN(VX(IENODE(I,2))*PI/180.0_DP) &
                 * 2._DP /(1._DP+sin(VY(IENODE(I,2))*PI/180.0_DP))

       DLTXC_TMP = VX2_TMP-VX1_TMP
       DLTYC_TMP = VY2_TMP-VY1_TMP

       EXFLUX    = STRESSX1IJ*DLTYC_TMP

       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS1PX(IA) = PS1PX(IA) - EXFLUX/ART(IA)
         PS1PX(IB) = PS1PX(IB) + EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS1PX(IA) = PS1PX(IA) - EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS1PX(IB) = PS1PX(IB) + EXFLUX/ART(IB)
       END IF
     
       EXFLUX    = STRESSX2IJ*DLTXC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS2PY(IA) = PS2PY(IA) + EXFLUX/ART(IA)
         PS2PY(IB) = PS2PY(IB) - EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS2PY(IA) = PS2PY(IA) + EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS2PY(IB) = PS2PY(IB) - EXFLUX/ART(IB)
       END IF

       EXFLUX    = STRESSX3IJ*DLTYC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS3PX(IA) = PS3PX(IA) - EXFLUX/ART(IA)
         PS3PX(IB) = PS3PX(IB) + EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS3PX(IA) = PS3PX(IA) - EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS3PX(IB) = PS3PX(IB) + EXFLUX/ART(IB)
       END IF

       EXFLUX    = STRESSX4IJ*DLTXC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS4PY(IA) = PS4PY(IA) + EXFLUX/ART(IA)
         PS4PY(IB) = PS4PY(IB) - EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS4PY(IA) = PS4PY(IA) + EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS4PY(IB) = PS4PY(IB) - EXFLUX/ART(IB)
       END IF

!       EXFLUX    = -STRESSX5IJ*DLTYC_TMP
       EXFLUX    = DLTYC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
!         PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)
!         PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)
         PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)*STRESSX5IJ
         PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)*STRESSX5IJ

       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
!         PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)
         PS5PX(IA) = PS5PX(IA) - EXFLUX/ART(IA)*STRESSX5IJ
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
!         PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)
         PS5PX(IB) = PS5PX(IB) + EXFLUX/ART(IB)*STRESSX5IJ
       END IF

       EXFLUX    = STRESSY1IJ*DLTXC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS1PY(IA) = PS1PY(IA) + EXFLUX/ART(IA)
         PS1PY(IB) = PS1PY(IB) - EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS1PY(IA) = PS1PY(IA) + EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS1PY(IB) = PS1PY(IB) - EXFLUX/ART(IB)
       END IF

       EXFLUX    = STRESSY2IJ*DLTYC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS2PX(IA) = PS2PX(IA) - EXFLUX/ART(IA)
         PS2PX(IB) = PS2PX(IB) + EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS2PX(IA) = PS2PX(IA) - EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS2PX(IB) = PS2PX(IB) + EXFLUX/ART(IB)
       END IF

       EXFLUX    = STRESSY3IJ*DLTXC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS3PY(IA) = PS3PY(IA) + EXFLUX/ART(IA)
         PS3PY(IB) = PS3PY(IB) - EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS3PY(IA) = PS3PY(IA) + EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS3PY(IB) = PS3PY(IB) - EXFLUX/ART(IB)
       END IF

       EXFLUX    = STRESSY4IJ*DLTYC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PS4PX(IA) = PS4PX(IA) - EXFLUX/ART(IA)
         PS4PX(IB) = PS4PX(IB) + EXFLUX/ART(IB)
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PS4PX(IA) = PS4PX(IA) - EXFLUX/ART(IA)
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PS4PX(IB) = PS4PX(IB) + EXFLUX/ART(IB)
       END IF

!       EXFLUX    = STRESSY5IJ*DLTXC_TMP
       EXFLUX    = DLTXC_TMP
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
!         PS5PY(IA) = PS5PY(IA) - EXFLUX/ART(IA)
!         PS5PY(IB) = PS5PY(IB) + EXFLUX/ART(IB)
         PS5PY(IA) = PS5PY(IA) + EXFLUX/ART(IA)*STRESSY5IJ
         PS5PY(IB) = PS5PY(IB) - EXFLUX/ART(IB)*STRESSY5IJ

       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
!         PS5PY(IA) = PS5PY(IA) - EXFLUX/ART(IA)
         PS5PY(IA) = PS5PY(IA) + EXFLUX/ART(IA)*STRESSY5IJ
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
!         PS5PY(IB) = PS5PY(IB) - EXFLUX/ART(IB)
         PS5PY(IB) = PS5PY(IB) - EXFLUX/ART(IB)*STRESSY5IJ
       END IF

!      ACCUMULATE BAROTROPIC FLUX
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*(DLTYC_TMP/ART(IA))/          &
                     (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
         PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*(DLTXC_TMP/ART(IA))/           &
                     (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
         PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*(DLTYC_TMP/ART(IB))/           &
                     (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
         PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*(DLTXC_TMP/ART(IB))/           &
                     (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         PSTXICE(IA)=PSTXICE(IA)-GRAV_E(IA)*ELIJ*(DLTYC_TMP/ART(IA))/           &
                     (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
         PSTYICE(IA)=PSTYICE(IA)+GRAV_E(IA)*ELIJ*(DLTXC_TMP/ART(IA))/           &
                     (2._SP /(1._SP+sin(YC(IA)*PI/180.0_SP)))
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN  
         PSTXICE(IB)=PSTXICE(IB)+GRAV_E(IB)*ELIJ*(DLTYC_TMP/ART(IB))/           &
                     (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
         PSTYICE(IB)=PSTYICE(IB)-GRAV_E(IB)*ELIJ*(DLTXC_TMP/ART(IB))/           &
                     (2._SP /(1._SP+sin(YC(IB)*PI/180.0_SP)))
       END IF
!     END IF  !  is ice

   END DO
   
   RETURN
   END SUBROUTINE ICE_STRESS_XY
!==============================================================================|
!
!==============================================================================|
!     ACCUMLATE FLUXES FOR ICE                                                 |
!==============================================================================|

   SUBROUTINE ICE_UV_XY(K)       

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESXICE,RESYICE
   REAL(SP) :: UICE2FT,VICE2FT
   INTEGER  :: I,II
   
   REAL(SP) :: UICE_TMP,VICE_TMP,STRAIRX_TMP,STRAIRY_TMP
   REAL(SP) :: PSTXICE_TMP,PSTYICE_TMP
   REAL(SP) :: UICE2RK_TMP,VICE2RK_TMP,UICE2F_TMP,VICE2F_TMP
   REAL(SP) :: DLTYC_TMP ,DLTXC_TMP

   real(SP) :: cca,ccb,ab2,cc1,cc2
   real(SP) :: UW_TMP,VW_TMP
   real(SP) :: TAUXWI_tmp,TAUYWI_tmp
!==============================================================================|
!
!--ACCUMULATE RESIDUALS FOR ICE EQUATIONS--------------------------------------|
!

   DO II=1,NP
     I=NP_LST(II)
!     IF(ISICEC(I) == 1)THEN
       UICE_TMP = -VICE2(I)*COS(XC(I)*PI/180.)-UICE2(I)*SIN(XC(I)*PI/180.)
       VICE_TMP = -VICE2(I)*SIN(XC(I)*PI/180.)+UICE2(I)*COS(XC(I)*PI/180.)

       UW_TMP = -V(I,1)*COS(XC(I)*PI/180.)-U(I,1)*SIN(XC(I)*PI/180.)
       VW_TMP = -V(I,1)*SIN(XC(I)*PI/180.)+U(I,1)*COS(XC(I)*PI/180.)

       TAUXWI_TMP = -TAUYWI(I)*COS(XC(I)*PI/180.)-TAUXWI(I)*SIN(XC(I)*PI/180.)
       TAUYWI_TMP = -TAUYWI(I)*SIN(XC(I)*PI/180.)+TAUXWI(I)*COS(XC(I)*PI/180.)

       STRAIRX_TMP = -STRAIRYE(I)*COS(XC(I)*PI/180.)-STRAIRXE(I)*SIN(XC(I)*PI/180.)
       STRAIRY_TMP = -STRAIRYE(I)*SIN(XC(I)*PI/180.)+STRAIRXE(I)*COS(XC(I)*PI/180.)

       PSTXICE_TMP = -PSTYICE(I)*COS(XC(I)*PI/180.)-PSTXICE(I)*SIN(XC(I)*PI/180.)
       PSTYICE_TMP = -PSTYICE(I)*SIN(XC(I)*PI/180.)+PSTXICE(I)*COS(XC(I)*PI/180.)

!       RESXICE(I) = ADVUICE(I)-ADVSTRX(I)+IMASS1(I)*(PSTXICE(I)           &
!                   -COR(I)*VICE_TMP*ART(I))                               &
!                   -(TAUXWI_TMP+STRAIRX_TMP)*ART(I)
!       RESYICE(I) = ADVVICE(I)-ADVSTRY(I)+IMASS1(I)*(PSTYICE(I)           &
!                   +COR(I)*UICE_TMP*ART(I))                               &
!                   -(TAUYWI_TMP+STRAIRY_TMP)*ART(I)

!       RESXICE(I) = ADVUICE(I)/ART(I)-ADVSTRX(I)+IMASS1(I)*(PSTXICE(I)           &
!                   -COR(I)*VICE_TMP) -(TAUXWI_TMP+STRAIRX_TMP)
!       RESYICE(I) = ADVVICE(I)/ART(I)-ADVSTRY(I)+IMASS1(I)*(PSTYICE(I)           &
!                   +COR(I)*UICE_TMP) -(TAUYWI_TMP+STRAIRY_TMP)

!
!--UPDATE----------------------------------------------------------------------|
!

        ! alpha, beta are defined in Hunke and Dukowicz (1997), section 3.2
!        cca = imass1(i)/dtice + vrel * cosw         ! alpha, kg/m^2 s
        cca = imass1(i)/DYN_DT + vrel(I) * cosw*ALPHA_RK(K) ! alpha, kg/m^2 s
!        ccb = fm(i,j)        + vrel * sinw        ! beta,  kg/m^2 s
        ccb = (cor(i)*imass1(i) + vrel(I) * sinw)*ALPHA_RK(K)  ! beta,  kg/m^2 s

        ab2 = cca**2 + ccb**2

        ! divergence of the internal stress tensor
        !  ADVSTRX(I),ADVSTRY(I)

         UICE2RK_TMP = -VICE2RK(I)*COS(XC(I)*PI/180.)-UICE2RK(I)*SIN(XC(I)*PI/180.)
         VICE2RK_TMP = -VICE2RK(I)*SIN(XC(I)*PI/180.)+UICE2RK(I)*COS(XC(I)*PI/180.)

        ! finally, the velocity components
        cc1 = (ADVSTRX(I) - IMASS1(I)*PSTXICE_TMP + TAUXWI_TMP+STRAIRX_TMP)&
     &       *ALPHA_RK(K) + imass1(I)/DYN_DT*UICE2RK_TMP
        cc2 = (ADVSTRY(I) - IMASS1(I)*PSTYICE_TMP + TAUYWI_TMP+STRAIRY_TMP)&
     &        *ALPHA_RK(K)+ imass1(I)/DYN_DT*VICE2RK_TMP

       UICE2F_TMP = (cca*cc1 + ccb*cc2)/ab2              ! m/s
       VICE2F_TMP = (cca*cc2 - ccb*cc1)/ab2

       UICE2F(I)  = VICE2F_TMP*COS(XC(I)*PI/180.)-UICE2F_TMP*SIN(XC(I)*PI/180.)
       VICE2F(I)  = UICE2F_TMP*COS(XC(I)*PI/180.)+VICE2F_TMP*SIN(XC(I)*PI/180.)
       VICE2F(I)  = -VICE2F(I)


!    endif
   END DO

   RETURN
   END SUBROUTINE ICE_UV_XY
!==============================================================================|
!
!==============================================================================|
!   CALCULATE CONVECTION FLUXES FOR ICE                                        !
!==============================================================================|
   SUBROUTINE ADVAVE_ICE_XY(XFLUX,YFLUX)
!==============================================================================|

!   USE BCS
!   USE MOD_OBCS

   IMPLICIT NONE
   INTEGER  :: I,II,IA,IB,J1,J2,K1,K2,K3
   REAL(SP) :: XIJ,YIJ,UIJ,VIJ,UIJ1,VIJ1,UIJ2,VIJ2
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,UN_TMP
   REAL(SP) :: XFLUX(0:NT),YFLUX(0:NT)
   
   REAL(SP) :: UICEIA,VICEIA,UICEIB,VICEIB
   REAL(SP) :: UICEK1,VICEK1,UICEK2,VICEK2,UICEK3,VICEK3
   REAL(SP) :: XIJC_TMP,YIJC_TMP,XIJ_TMP,YIJ_TMP
   REAL(SP) :: XCIA_TMP,YCIA_TMP,XCIB_TMP,YCIB_TMP
   REAL(SP) :: DLTYC_TMP , DLTXC_TMP

!
!--Initialize Variables--------------------------------------------------------|
!
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     IF(CELL_NORTHAREA(IA) == 1)THEN
       XFLUX(IA) = 0.0_SP
       YFLUX(IA) = 0.0_SP
     END IF  
     IF(CELL_NORTHAREA(IB) == 1)THEN  
       XFLUX(IB) = 0.0_SP
       YFLUX(IB) = 0.0_SP
     END IF  
   END DO  
!
!-------------------------ACCUMULATE FLUX OVER ELEMENT EDGES-------------------!
!
   DO II=1,NPE
     I=NPEDGE_LST(II)
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

     IF(ISICEC(IA) == 1 .OR. ISICEC(IB) == 1)THEN

!      FLUX FROM LEFT
       K1=NBE(IA,1)
       K2=NBE(IA,2)
       K3=NBE(IA,3)
         
       UICEIA = -VICE2(IA)*COS(XC(IA)*PI/180.)-UICE2(IA)*SIN(XC(IA)*PI/180.)
       VICEIA = -VICE2(IA)*SIN(XC(IA)*PI/180.)+UICE2(IA)*COS(XC(IA)*PI/180.)
       UICEK1 = -VICE2(K1)*COS(XC(K1)*PI/180.)-UICE2(K1)*SIN(XC(K1)*PI/180.)
       VICEK1 = -VICE2(K1)*SIN(XC(K1)*PI/180.)+UICE2(K1)*COS(XC(K1)*PI/180.)
       UICEK2 = -VICE2(K2)*COS(XC(K2)*PI/180.)-UICE2(K2)*SIN(XC(K2)*PI/180.)
       VICEK2 = -VICE2(K2)*SIN(XC(K2)*PI/180.)+UICE2(K2)*COS(XC(K2)*PI/180.)
       UICEK3 = -VICE2(K3)*COS(XC(K3)*PI/180.)-UICE2(K3)*SIN(XC(K3)*PI/180.)
       VICEK3 = -VICE2(K3)*SIN(XC(K3)*PI/180.)+UICE2(K3)*COS(XC(K3)*PI/180.)
     
       COFA1=A1U(IA,1)*UICEIA+A1U(IA,2)*UICEK1+A1U(IA,3)*UICEK2+A1U(IA,4)*UICEK3
       COFA2=A2U(IA,1)*UICEIA+A2U(IA,2)*UICEK1+A2U(IA,3)*UICEK2+A2U(IA,4)*UICEK3
       COFA5=A1U(IA,1)*VICEIA+A1U(IA,2)*VICEK1+A1U(IA,3)*VICEK2+A1U(IA,4)*VICEK3
       COFA6=A2U(IA,1)*VICEIA+A2U(IA,2)*VICEK1+A2U(IA,3)*VICEK2+A2U(IA,4)*VICEK3
     
       XIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * COS(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
       YIJC_TMP = REARTH * COS(YIJC(I)*PI/180.0_SP) * SIN(XIJC(I)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YIJC(I)*PI/180.0_SP))
       XCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * COS(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
       YCIA_TMP = REARTH * COS(YC(IA)*PI/180.0_SP) * SIN(XC(IA)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IA)*PI/180.0_SP))
  
       XIJ_TMP = XIJC_TMP-XCIA_TMP
       YIJ_TMP = YIJC_TMP-YCIA_TMP
       UIJ1=UICEIA+COFA1*XIJ_TMP+COFA2*YIJ_TMP
       VIJ1=VICEIA+COFA5*XIJ_TMP+COFA6*YIJ_TMP

!      FLUX FROM RIGHT
       K1=NBE(IB,1)
       K2=NBE(IB,2)
       K3=NBE(IB,3)
          
       UICEIB = -VICE2(IB)*COS(XC(IB)*PI/180.)-UICE2(IB)*SIN(XC(IB)*PI/180.)
       VICEIB = -VICE2(IB)*SIN(XC(IB)*PI/180.)+UICE2(IB)*COS(XC(IB)*PI/180.)
       UICEK1 = -VICE2(K1)*COS(XC(K1)*PI/180.)-UICE2(K1)*SIN(XC(K1)*PI/180.)
       VICEK1 = -VICE2(K1)*SIN(XC(K1)*PI/180.)+UICE2(K1)*COS(XC(K1)*PI/180.)
       UICEK2 = -VICE2(K2)*COS(XC(K2)*PI/180.)-UICE2(K2)*SIN(XC(K2)*PI/180.)
       VICEK2 = -VICE2(K2)*SIN(XC(K2)*PI/180.)+UICE2(K2)*COS(XC(K2)*PI/180.)
       UICEK3 = -VICE2(K3)*COS(XC(K3)*PI/180.)-UICE2(K3)*SIN(XC(K3)*PI/180.)
       VICEK3 = -VICE2(K3)*SIN(XC(K3)*PI/180.)+UICE2(K3)*COS(XC(K3)*PI/180.)
     
       COFA3=A1U(IB,1)*UICEIB+A1U(IB,2)*UICEK1+A1U(IB,3)*UICEK2+A1U(IB,4)*UICEK3
       COFA4=A2U(IB,1)*UICEIB+A2U(IB,2)*UICEK1+A2U(IB,3)*UICEK2+A2U(IB,4)*UICEK3
       COFA7=A1U(IB,1)*VICEIB+A1U(IB,2)*VICEK1+A1U(IB,3)*VICEK2+A1U(IB,4)*VICEK3
       COFA8=A2U(IB,1)*VICEIB+A2U(IB,2)*VICEK1+A2U(IB,3)*VICEK2+A2U(IB,4)*VICEK3
     
       XCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * COS(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
       YCIB_TMP = REARTH * COS(YC(IB)*PI/180.0_SP) * SIN(XC(IB)*PI/180.0_SP) &
                  * 2._SP /(1._SP+SIN(YC(IB)*PI/180.0_SP))
  
       XIJ_TMP=XIJC_TMP-XCIB_TMP
       YIJ_TMP=YIJC_TMP-YCIB_TMP
       UIJ2=UICEIB+COFA3*XIJ_TMP+COFA4*YIJ_TMP
       VIJ2=VICEIB+COFA7*XIJ_TMP+COFA8*YIJ_TMP

!      NORMAL VELOCITY
       UIJ=0.5_SP*(UIJ1+UIJ2)
       VIJ=0.5_SP*(VIJ1+VIJ2)
       UN_TMP=-UIJ*DLTYC_TMP + VIJ*DLTXC_TMP

!      ADD CONVECTIVE AND VISCOUS FLUXES
       XADV=UN_TMP*&
            ((1.0_SP-SIGN(1.0_SP,UN_TMP))*UIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*UIJ1)*0.5_SP
       YADV=UN_TMP* &
            ((1.0_SP-SIGN(1.0_SP,UN_TMP))*VIJ2+(1.0_SP+SIGN(1.0_SP,UN_TMP))*VIJ1)*0.5_SP

!      ACCUMULATE FLUX
       IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) == 1)THEN
         XFLUX(IA)=XFLUX(IA)+XADV*(1.0_SP-ISBC(I))
         YFLUX(IA)=YFLUX(IA)+YADV*(1.0_SP-ISBC(I))
         XFLUX(IB)=XFLUX(IB)-XADV*(1.0_SP-ISBC(I))
         YFLUX(IB)=YFLUX(IB)-YADV*(1.0_SP-ISBC(I))
       ELSE IF(CELL_NORTHAREA(IA) == 1 .AND. CELL_NORTHAREA(IB) /= 1)THEN
         XFLUX(IA)=XFLUX(IA)+XADV*(1.0_SP-ISBC(I))
         YFLUX(IA)=YFLUX(IA)+YADV*(1.0_SP-ISBC(I))
       ELSE IF(CELL_NORTHAREA(IB) == 1 .AND. CELL_NORTHAREA(IA) /= 1)THEN
         XFLUX(IB)=XFLUX(IB)-XADV*(1.0_SP-ISBC(I))
         YFLUX(IB)=YFLUX(IB)-YADV*(1.0_SP-ISBC(I))
       END IF 
     END IF
   END DO

   RETURN
   END SUBROUTINE ADVAVE_ICE_XY

!==============================================================================|
#  endif
!end if defined (NORTHPOLE)


!==============================================================================|
!   Calculate Advection  for ICE variables                                     |
!==============================================================================|

   SUBROUTINE ICE_ADV_fvcom(var_ice,cuice,cvice)       

!------------------------------------------------------------------------------|
   IMPLICIT NONE

   REAL(SP), DIMENSION(0:MT)    :: XFLUX
   REAL(SP), DIMENSION(0:MT)    :: PUPX,PUPY,PVPX,PVPY  
   REAL(DP), DIMENSION(0:MT)    :: PIPX,PIPY
   REAL(SP), DIMENSION(3*(NT))  :: DTIJ 
   REAL(SP), DIMENSION(3*(NT))  :: UVN

   REAL(SP), DIMENSION(0:NT)  :: cuice,cvice
   REAL(SP), DIMENSION(0:MT)  :: var_ice,var_ice1

   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
#  if defined (SPHERICAL)
   REAL(DP) :: TY,TXPI,TYPI
   REAL(DP) :: XTMP1,XTMP
   REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
   REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
# if defined (NORTHPOLE)
   REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
   REAL(DP) :: TXPI_TMP,TYPI_TMP
# endif
#  endif

   REAL(SP) :: VI1MIN, VI1MAX, VI2MIN, VI2MAX
   REAL(SP) :: XMIN, XMAX


!------------------------------------------------------------------------------!
!
!--Initialize Fluxes-----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   UVN       =0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO I=1,NCV
     I1=NTRG(I)
     UVN(I) = cvice(I1)*DLTXE(I) - cuice(I1)*DLTYE(I)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

      PIPX  = 0.0_SP
      PIPY  = 0.0_SP
      DO I=1,M
         DO J=1,NTSN(I)-1
            I1=NBSN(I,J)
            I2=NBSN(I,J+1)
            FFD=0.5_SP*(var_ice(I1)+var_ice(I2))
            FF1=0.5_SP*(var_ice(I1)+var_ice(I2))


!!$#        if defined (SPHERICAL)
!!$         XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$         XTMP1 = VX(I2)-VX(I1)
!!$         IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$         ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$         END IF  
!!$!====ggao 2007 1120====================================================
!!$!         TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$        TXPI=XTMP*VAL_COS_VY(I)
!!$!======================================================================
!!$
!!$         TYPI=(VY(I1)-VY(I2))*tpi
!!$#    if defined (NORTHPOLE)
!!$         IF(NODE_NORTHAREA(I) == 1)THEN
!!$           VX1_TMP = REARTH * COS(VY(I1)*PI/180.0_SP) * COS(VX(I1)*PI/180.0_SP) &
!!$                     * 2._SP /(1._SP+SIN(VY(I1)*PI/180.0_SP))
!!$           VY1_TMP = REARTH * COS(VY(I1)*PI/180.0_SP) * SIN(VX(I1)*PI/180.0_SP) &
!!$                     * 2._SP /(1._SP+SIN(VY(I1)*PI/180.0_SP))
!!$
!!$           VX2_TMP = REARTH * COS(VY(I2)*PI/180.0_SP) * COS(VX(I2)*PI/180.0_SP) &
!!$                     * 2._SP /(1._SP+SIN(VY(I2)*PI/180.0_SP))
!!$           VY2_TMP = REARTH * COS(VY(I2)*PI/180.0_SP) * SIN(VX(I2)*PI/180.0_SP) &
!!$                     * 2._SP /(1._SP+SIN(VY(I2)*PI/180.0_SP))
!!$
!!$           TXPI = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP)))
!!$           TYPI = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP)))
!!$           IF(I /= NODE_NORTHPOLE)THEN
!!$             TXPI_TMP = TYPI*COS(VX(I)*PI/180.)-TXPI*SIN(VX(I)*PI/180.)
!!$             TYPI_TMP = TXPI*COS(VX(I)*PI/180.)+TYPI*SIN(VX(I)*PI/180.)
!!$             TYPI_TMP = -TYPI_TMP
!!$       
!!$             TXPI = TXPI_TMP
!!$             TYPI = TYPI_TMP
!!$           END IF  
!!$          END IF 
!!$#    endif	 
!!$         PIPX(I)=PIPX(I)+FF1*TYPI
!!$         PIPY(I)=PIPY(I)+FF1*TXPI
!!$#        else
!!$         PIPX(I)=PIPX(I)+FF1*(VY(I1)-VY(I2))
!!$         PIPY(I)=PIPY(I)+FF1*(VX(I2)-VX(I1))
!!$#        endif
!!$       END DO
!!$       PIPX(I)=PIPX(I)/ART2(I)
!!$       PIPY(I)=PIPY(I)/ART2(I)
            
            PIPX(I)=PIPX(I)+FF1*DLTYTRIE(i,j)
            PIPY(I)=PIPY(I)+FF1*DLTXTRIE(i,j)
            
         END DO
         PIPX(I)=PIPX(I)/ART2(I)
         PIPY(I)=PIPY(I)/ART2(I)
         
      END DO
      
      DO I=1,NCV_I
         
         I1=NTRG(I)
         
         !      IF(ISICEC(I1)==1) THEN
         UVN(I) = cvice(I1)*DLTXE(I) - cuice(I1)*DLTYE(I)
         
         IA=NIEC(I,1)
         IB=NIEC(I,2)


!!$       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!!$       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!$#      if defined (SPHERICAL)
!!$       X1_DP=XIJE(I,1)
!!$       Y1_DP=YIJE(I,1)
!!$       X2_DP=XIJE(I,2)
!!$       Y2_DP=YIJE(I,2)
!!$!======ggao 1120/2007
!!$!       CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!!$         XII=XCG2(I)
!!$         YII=YCG2(I)
!!$!======ggao 1120/2007
!!$
!!$
!!$       XI=XII
!!$       XTMP  = XI*TPI-VX(IA)*TPI
!!$       XTMP1 = XI-VX(IA)
!!$       IF(XTMP1 >  180.0_SP)THEN
!!$         XTMP = -360.0_SP*TPI+XTMP
!!$       ELSE IF(XTMP1 < -180.0_SP)THEN
!!$         XTMP =  360.0_SP*TPI+XTMP
!!$       END IF 
!!$       DXA=XTMP*COS(DEG2RAD*VY(IA))  
!!$       DYA=(YI-VY(IA))*TPI
!!$
!!$       XTMP  = XI*TPI-VX(IB)*TPI
!!$       XTMP1 = XI-VX(IB)
!!$       IF(XTMP1 >  180.0_SP)THEN
!!$         XTMP = -360.0_SP*TPI+XTMP
!!$       ELSE IF(XTMP1 < -180.0_SP)THEN
!!$         XTMP =  360.0_SP*TPI+XTMP
!!$       END IF 
!!$
!!$       DXB=XTMP*COS(DEG2RAD*VY(IB))
!!$       DYB=(YI-VY(IB))*TPI
!!$#      else
!!$       DXA=XI-VX(IA)
!!$       DYA=YI-VY(IA)
!!$       DXB=XI-VX(IB)
!!$       DYB=YI-VY(IB)
!!$#      endif
!!$       FIJ1=var_ice(IA)+DXA*PIPX(IA)+DYA*PIPY(IA)
!!$       FIJ2=var_ice(IB)+DXB*PIPX(IB)+DYB*PIPY(IB)

         FIJ1=var_ice(IA)+DLTXNCVE(I,1)*PIPX(IA)+DLTYNCVE(I,1)*PIPY(IA)
         FIJ2=var_ice(IB)+DLTXNCVE(I,2)*PIPX(IB)+DLTYNCVE(I,2)*PIPY(IB)
         


         VI1MIN=MINVAL(var_ice(NBSN(IA,1:NTSN(IA)-1)))
         VI1MIN=MIN(VI1MIN, var_ice(IA))
         VI1MAX=MAXVAL(var_ice(NBSN(IA,1:NTSN(IA)-1)))
         VI1MAX=MAX(VI1MAX, var_ice(IA))
         VI2MIN=MINVAL(var_ice(NBSN(IB,1:NTSN(IB)-1)))
         VI2MIN=MIN(VI2MIN, var_ice(IB))
         VI2MAX=MAXVAL(var_ice(NBSN(IB,1:NTSN(IB)-1)))
         VI2MAX=MAX(VI2MAX, var_ice(IB))
         IF(FIJ1 < VI1MIN) FIJ1=VI1MIN
         IF(FIJ1 > VI1MAX) FIJ1=VI1MAX
         IF(FIJ2 < VI2MIN) FIJ2=VI2MIN
         IF(FIJ2 > VI2MAX) FIJ2=VI2MAX
         
         UN=UVN(I)
         
         EXFLUX=-UN* &
              ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP 
         
         XFLUX(IA)=XFLUX(IA)+EXFLUX
         XFLUX(IB)=XFLUX(IB)-EXFLUX
         
         !      END IF  ! just calculate the ice
         
      END DO
      
#    if defined (SPHERICAL)
#    if defined (NORTHPOLE)
     CALL ICE_ADV_XY(XFLUX,PIPX,PIPY,var_ice,cuice,cvice)
#    endif
#    endif  

# if defined (MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,XFLUX)
# endif

!--Update ice variables----------------------------------------------------------!
!
    DO I=1,M
       var_ICE1(I)=var_ice(I)-XFLUX(I)*dtice/ART1(I)

       XMAX = MAX(MAXVAL(VAR_ICE(NBSN(I,1:NTSN(I)))),VAR_ICE(I))
       XMIN = MIN(MINVAL(VAR_ICE(NBSN(I,1:NTSN(I)))),VAR_ICE(I))

       IF(XMIN-VAR_ICE1(I) > 0.0_SP)VAR_ICE1(I) = XMIN
       IF(VAR_ICE1(I)-XMAX > 0.0_SP)VAR_ICE1(I) = XMAX

       !IF(ISICEN(I)==1)THEN
       var_ICE(I)=var_ICE1(I)
       !ELSE
       !VAR_ICE(I)=0.0_SP
       !END IF
       !if (VAR_ICE(i) < puny)VAR_ICE(I) = c0i
       if (VAR_ICE(I) < 0.0)VAR_ICE(I) = c0i
    END DO


   RETURN
   END SUBROUTINE ICE_ADV_fvcom
!==============================================================================|

#  if defined (NORTHPOLE)

!==============================================================================|
   SUBROUTINE ICE_ADV_XY(XFLUX,PIPX,PIPY,var_ice,cuice,cvice)               

!------------------------------------------------------------------------------|

   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT)     :: XFLUX,XFLUX_ADV
   REAL(DP), DIMENSION(0:MT)     :: PIPX,PIPY
   REAL(SP), DIMENSION(3*(NT))   :: DTIJ 
   REAL(SP) :: XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2 
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF   
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: TXPI,TYPI

   REAL(SP) :: VX_TMP,VY_TMP,VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP,VX3_TMP,VY3_TMP
   REAL(SP) :: XI_TMP,YI_TMP,VXA_TMP,VYA_TMP,VXB_TMP,VYB_TMP
   REAL(SP) :: UIJ_TMP,VIJ_TMP,DLTXE_TMP,DLTYE_TMP,UVN_TMP,EXFLUX_TMP
   REAL(SP) :: PUPX_TMP,PUPY_TMP,PVPX_TMP,PVPY_TMP
   REAL(SP) :: PIPX_TMP,PIPY_TMP
   REAL(SP) :: U_TMP,V_TMP
   REAL(SP) :: X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2
   integer :: KK

   REAL(SP), DIMENSION(0:NT)  :: cuice,cvice
   REAL(SP), DIMENSION(0:MT)  :: var_ice,var_ice1

!------------------------------------------------------------------------------!
!
!--Initialize Fluxes-----------------------------------------------------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     IA = NIEC(I,1)
     IB = NIEC(I,2)
     IF(IA == NODE_NORTHPOLE)THEN
       XFLUX(IA) = 0.0_SP
       XFLUX_ADV(IA) = 0.0_SP
     ELSE IF(IB == NODE_NORTHPOLE)THEN  
       XFLUX(IB) = 0.0_SP
       XFLUX_ADV(IB) = 0.0_SP
     END IF  
   END DO  
     
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!
   I = NODE_NORTHPOLE

   IF(I==0)  RETURN

   PUPX_TMP=0.0_SP
   PUPY_TMP=0.0_SP
   PVPX_TMP=0.0_SP
   PVPY_TMP=0.0_SP

   DO J=1,NTVE(I)
     I1=NBVE(I,J)
     JTMP=NBVT(I,J)
     J1=JTMP+1-(JTMP+1)/4*3
     J2=JTMP+2-(JTMP+2)/4*3
       
     VX_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * COS(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
     VY_TMP = REARTH * COS(VY(I)*PI/180.0_SP) * SIN(VX(I)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(I)*PI/180.0_SP))
     
     VX1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * COS(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
     VY1_TMP= REARTH * COS(VY(NV(I1,J1))*PI/180.0_SP) * SIN(VX(NV(I1,J1))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J1))*PI/180.0_SP))
     
     VX2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * COS(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
     VY2_TMP= REARTH * COS(YC(I1)*PI/180.0_SP) * SIN(XC(I1)*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(YC(I1)*PI/180.0_SP))
     
     VX3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * COS(VX(NV(I1,J2))*PI/180.0_SP) &
                     * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
     VY3_TMP= REARTH * COS(VY(NV(I1,J2))*PI/180.0_SP) * SIN(VX(NV(I1,J2))*PI/180.0_SP) &
                    * 2._SP /(1._SP+SIN(VY(NV(I1,J2))*PI/180.0_SP))
     
     X11=0.5_SP*(VX_TMP+VX1_TMP)
     Y11=0.5_SP*(VY_TMP+VY1_TMP)
     X22=VX2_TMP
     Y22=VX2_TMP
     X33=0.5_SP*(VX_TMP+VX3_TMP)
     Y33=0.5_SP*(VY_TMP+VY3_TMP)

     U_TMP = -cvice(I1)*COS(XC(I1)*PI/180.)-cuice(I1)*SIN(XC(I1)*PI/180.)
     V_TMP = -cvice(I1)*SIN(XC(I1)*PI/180.)+cuice(I1)*COS(XC(I1)*PI/180.)

     PUPX_TMP=PUPX_TMP+U_TMP*(Y11-Y33)
     PUPY_TMP=PUPY_TMP+U_TMP*(X33-X11)
     PVPX_TMP=PVPX_TMP+V_TMP*(Y11-Y33)
     PVPY_TMP=PVPY_TMP+V_TMP*(X33-X11)
   END DO

   PUPX_TMP=PUPX_TMP/ART1(I)
   PUPY_TMP=PUPY_TMP/ART1(I)
   PVPX_TMP=PVPX_TMP/ART1(I)
   PVPY_TMP=PVPY_TMP/ART1(I)
   TMP1=PUPX_TMP**2+PVPY_TMP**2
   TMP2=0.5_SP*(PUPY_TMP+PVPX_TMP)**2

   DO II=1,NPCV
     I = NCEDGE_LST(II)
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)
     
     IF((IA <= M .AND. IB <= M) .AND. I1 <= N)THEN
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))

       IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         XI_TMP = REARTH * COS(YI*PI/180.0_SP) * COS(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))
         YI_TMP = REARTH * COS(YI*PI/180.0_SP) * SIN(XI*PI/180.0_SP) &
                  * 2._SP /(1._SP+sin(YI*PI/180.0_SP))

         VXA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * COS(VX(IA)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))
         VYA_TMP = REARTH * COS(VY(IA)*PI/180.0_SP) * SIN(VX(IA)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IA)*PI/180.0_SP))

         VXB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * COS(VX(IB)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))
         VYB_TMP = REARTH * COS(VY(IB)*PI/180.0_SP) * SIN(VX(IB)*PI/180.0_SP) &
                   * 2._SP /(1._SP+sin(VY(IB)*PI/180.0_SP))

!         IF(IA == NODE_NORTHPOLE)THEN
         DXA=XI_TMP-VXA_TMP
         DYA=YI_TMP-VYA_TMP
!         ELSE IF(IB == NODE_NORTHPOLE)THEN
         DXB=XI_TMP-VXB_TMP
         DYB=YI_TMP-VYB_TMP
!	 END IF
!       END IF

        IF(IA == NODE_NORTHPOLE)THEN
          PIPX_TMP=-PIPY(IB)*COS(VX(IB)*PI/180.)-PIPX(IB)*SIN(VX(IB)*PI/180.)
          PIPY_TMP=-PIPY(IB)*SIN(VX(IB)*PI/180.)+PIPX(IB)*COS(VX(IB)*PI/180.)
   
   
          FIJ1=var_ice(IA)+DXA*PIPX(IA)+DYA*PIPY(IA)
          FIJ2=var_ice(IB)+DXB*PIPX_TMP+DYB*PIPY_TMP

        ELSE IF(IB == NODE_NORTHPOLE)THEN
          PIPX_TMP=-PIPY(IA)*COS(VX(IA)*PI/180.)-PIPX(IA)*SIN(VX(IA)*PI/180.)
          PIPY_TMP=-PIPY(IA)*SIN(VX(IA)*PI/180.)+PIPX(IA)*COS(VX(IA)*PI/180.)
   
          FIJ1=var_ice(IA)+DXA*PIPX_TMP+DYA*PIPY_TMP
          FIJ2=var_ice(IB)+DXB*PIPX(IB)+DYB*PIPY(IB)
        END IF

!      IF(IA == NODE_NORTHPOLE .OR. IB == NODE_NORTHPOLE)THEN
         UIJ_TMP = -cvice(I1)*COS(XC(I1)*PI/180.)-cuice(I1)*SIN(XC(I1)*PI/180.)
         VIJ_TMP = -cvice(I1)*SIN(XC(I1)*PI/180.)+cuice(I1)*COS(XC(I1)*PI/180.)
       
         VX1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * COS(XIJE(I,1)*PI/180.0_SP)
         VY1_TMP = REARTH * COS(YIJE(I,1)*PI/180.0_SP) * SIN(XIJE(I,1)*PI/180.0_SP)

         VX2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * COS(XIJE(I,2)*PI/180.0_SP)
         VY2_TMP = REARTH * COS(YIJE(I,2)*PI/180.0_SP) * SIN(XIJE(I,2)*PI/180.0_SP)

         DLTXE_TMP = VX2_TMP-VX1_TMP
         DLTYE_TMP = VY2_TMP-VY1_TMP
       

         UVN_TMP = VIJ_TMP*DLTXE_TMP - UIJ_TMP*DLTYE_TMP
         EXFLUX_TMP = -UVN_TMP*((1.0_SP+SIGN(1.0_SP,UVN_TMP))*FIJ2+   &
                      (1.0_SP-SIGN(1.0_SP,UVN_TMP))*FIJ1)*0.5_SP
       
         IF(IA == NODE_NORTHPOLE)THEN
           XFLUX(IA)=XFLUX(IA)+EXFLUX_TMP 
         ELSE IF(IB == NODE_NORTHPOLE)THEN
           XFLUX(IB)=XFLUX(IB)-EXFLUX_TMP 
         END IF
       END IF
     END IF  
   END DO
# if defined (MULTIPROCESSOR)
   !IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,XFLUX)
# endif




   RETURN
   END SUBROUTINE ICE_ADV_XY
!==============================================================================|
#  endif
!end if defined (NORTHPOLE)

#  endif
END MODULE MOD_ICE2D
 
