










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

!=======================================================================
! FVCOM Sediment Module 
!
! Copyright:    2005(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Contact:      G. Cowles 
!               School for Marine Science and Technology, Umass-Dartmouth
!
! Based on the Community Sediment Transport Model (CSTM) as implemented
!     in ROMS by J. Warner (USGS)
!
! Comments:     Sediment Dynamics Module 
!
! Current FVCOM dependency
!
!   init_sed.F:  - user defined sediment model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration/bottom/bed fields
!   fvcom.F:     - main calls sediment setup 
!   internal.F:  - calls sediment advance
!
! History
!   Feb 7, 2008: added initialization of bottom(:,:) to 0 (w/ T. Hamada)
!              : fixed loop bounds in hot start and archive for conc (w/ T. Hamada)
!              : added comments describing theoretical bases of dynamics
!   Feb 14,2008: added non-constant settling velocity for cohesive sediments (w/ T. Hamada) 
!              : updated vertical flux routine to handle non-constant vertical velocity (w/ T. Hamada)
!              : added a user-defined routine to calculate settling velocity based on concentration (w/ T. Hamada)
!              : added a user-defined routine to calculate erosion for a general case (w/ T. Hamada)
!
!  PLEASE NOTE!!!!!!!!!!! 
!  Do NOT USE INTEL FORTRAN COMPILER VERSION 11.0 IT HAS KNOWN BUGS WHEN DEALING WITH TYPES WITH ALLOCATABLE
!    COMPONENTS.  YOU WILL SEE WEIRD BEHAVIOR.  VERSION 11.1 IS OK.
!   
!
!  Later
!   1.) Modify vertical flux routines to work with general vertical coordinate
!   2.) Add divergence term for bedload transport calc 
!   3.) Add ripple roughness calculation
!   4.) Add morphological change (bathymetry + vertical velocity condition) 
!   5.) Eliminate excess divisions and recalcs
!
!=======================================================================
Module Mod_Sed  
   

End Module Mod_Sed 
