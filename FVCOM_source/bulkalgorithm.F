
!/===========================================================================/
! Copyright (c) 2020, The University of Tokyo 
! Produced at Jun sasaki lab, Yulong Wang
! All rights reserved.
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!------------------------------------------------------------------------------
# if defined (HEATING_CALCULATED)
# if defined (DLW_CALCULATED)
  subroutine bulkalgorithm(wind,tair,humidity,pair,TC,cloud,sql,sqh,sqe,rain,evap)

! Based on MEC model, bulk.f
! Yulong Wang
! Oct 27, 2020

    implicit none

    real :: wind,tair,humidity,pair,TC,cloud,sql,sqh,sqe
    real :: Ca,Sss,sigma,Ccc,Lll,Ce,Ced,Ch,Chd
    real :: windabs,dt,TaK,es,ea,rhoa,qs,qa,dq
    real :: tbifu
    real :: rain,evap

!===========================================================================================|
!   Inp/Otp Meteorological Boundary Conditions for Calculating Heat Flux                    |
!===========================================================================================|
!    [I]  wind speed at height 2m:             m/s          "wind"                          |
!    [I]  bulk air temperature at height 2m:   degree(C)    "tair"                          |
!    [I]  relative humidity at height 2m:      (%)          "humidity"                      |
!    [I]  surface pressure:                    mb[0.1kPa]   "pair"                          |
!    [I]  calculated water temperature:        degree(C)    "TC"                            |
!    [I]  cloud cover:                         (0-1)        "cloud"                         |
!    [O]  Longwave radiation from sea surface: w/m^2        "sql"                           |
!    [O]  Sensible heat flux due to con&trans: w/m^2        "sqh"                           |
!    [O]  Latent heat flux due to evaporation: w/m^2        "sqe"                           |
!    [I]  precipitation:                       kg/m2/s      "rain",ocean lose water is negative|   
!    [O]  evaporation:                         kg/m2/s      "evap",ocean lose water is negative|  
!===========================================================================================|

!***********  meteorological constants ***********************************
    Ca=1010.         ! Specific heat of air at constant pressure, [J/kg/K]
    Sss=0.96         ! Emissivity, [-]
    sigma=5.67e-8    ! Stefan-Boltzmann constant, [W/m2/K4]
    Ccc=0.65         ! Constant due to latitude, [-]
    Lll=2.45e6      ! Latent heat by evaporation, [J/kg]
    Ce=1.1e-3       ! Bulk coefficient of latent heat transfer, [-]
    Ced=1.1e-3       ! Bulk coefficient of latent heat transfer, [-]
    Ch=1.1e-3       ! Bulk coefficient of sensible heat transfer, [-]
    Chd=1.1e-3       ! Bulk coefficient of sensible heat transfer, [-]

!***********  basic values *********************************************
    windabs=ABS(wind)
    dt=TC-tair
    TaK=tair+273.15
    es=6.1078*10**(7.5*tair/TaK)                                    ! [hPa]
    ea=es*humidity/100.                                             ! [hPa]
    rhoa=1.293*(273.15/TaK)*(pair/1013.25)*(1.0-0.378*(ea/pair))    ! [kg/m3]
    qs=0.622*es/(pair-0.378*es)
    qa=0.622*ea/(pair-0.378*ea)
    dq=qs-qa

!*********** long wave radiation ***************
    sql=Sss*sigma*TaK**4*(0.39-0.058*sqrt(ea))* &
    (1.-Ccc*cloud**2)+4.*Sss*sigma*TaK**(3)*dt

!*********** sensible heat & latent heat***************
    IF (windabs>=1.) THEN
        sqh=Ca*rhoa*Ch*dt*windabs
        sqe=max(Lll*rhoa*Ce*dq*windabs,0.)
    ELSE
        tbifu=dt+0.61*TaK*dq
        IF (tbifu< 0.) THEN
            sqh=-Ca*rhoa*Chd*-tbifu**(1./3.)*dt
            sqe=0.
        ELSE    
            sqh=Ca*rhoa*Chd*tbifu**(1./3.)*dt
            sqe=Lll*rhoa*Ced*tbifu**(1./3.)*dq
        END IF
    END IF

!*********** adjust them to upward value negitive ***************
    sql = sql * -1.
    sqh = sqh * -1.
    sqe = sqe * -1.

!*********** rainfall and evaporation ***************
    IF (windabs>=1.) THEN
        evap = rhoa*Ce*dt*windabs * -1.
    ELSE
        tbifu=dt+0.61*TaK*dq
        IF (tbifu< 0.) THEN
            evap = rain * -1.
        ELSE    
            evap = rhoa*Ced*tbifu**(1./3.)*dq * -1.
        END IF
    END IF

  end subroutine bulkalgorithm
# endif
# endif
