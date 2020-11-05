
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

!------------------------------------------------------------------------------
# if defined (HEATING_CALCULATED)
# if defined (DLW_CALCULATED)
  subroutine bulkalgorithm(wind,tair,humidity,pair,TC,cloud,sql,sqh,sqe)

! Based on MEC model
! Yulong Wang
! Oct 27, 2020

    implicit none

    real :: wind,tair,humidity,pair,TC,cloud,sql,sqh,sqe
    real :: Ca,Sss,sigma,Ccc,Lll,Ce,Ced,Ch,Chd
    real :: windabs,dt,TaK,es,ea,rhoa,qs,qa,dq
    real :: tbifu


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
    es=6.1078*10**(7.5*tair/TaK)
    ea=es*humidity/100.
    rhoa=1.293*(273.15/TaK)*(pair/1013.25)*(1.0-0.378*(ea/pair))
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

  end subroutine bulkalgorithm
# endif
# endif
