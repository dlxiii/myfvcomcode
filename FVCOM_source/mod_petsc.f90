










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

!  PETSC 2.3.2 / PETSC 2.3.3 COMPATIBILITY ISSUE FOR VecScatter:
! 2.3.2 : CALL VecScatterBegin(BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,L2G_EL,IERR);CHKERRQ(IERR)
!
! 2.3.3 : CALL VecScatterBegin(L2G_EL,BL_EL,B_EL,INSERT_VALUES,SCATTER_FORWARD,IERR);CHKERRQ(IERR)
!


MODULE MOD_PETSc

!=========================================================================================
!
!=========================================================================================
END MODULE MOD_PETSc

