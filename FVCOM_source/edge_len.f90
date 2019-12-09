










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

SUBROUTINE edge_len
  USE all_vars
  use mod_spherical
  use mod_northpole
  use mod_utils

  implicit none
  INTEGER  :: I,IP,I1,I2,IA,IB,J,J1
  REAL(SP) :: XI, YI, X11,X33,Y11,Y33

  INTEGER :: JTMP,J2



  ! Distance between control volue edge and the node

  DO I=1,NCV_I
     IA=NIEC(I,1)
     IB=NIEC(I,2)     
     XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
     YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
     DLTXNCVE(I,1)=XI-VX(IA)
     DLTYNCVE(I,1)=YI-VY(IA)
     DLTXNCVE(I,2)=XI-VX(IB)
     DLTYNCVE(I,2)=YI-VY(IB)

  END DO


  ! Set the distance between Nodes

  IF(MAXVAL(NTSN) > 13) CALL FATAL_ERROR &
       & ("THERE ARE MORE THAN 12 NODES AROUND ONE NODE:",&
       "PLEASE INCREASE THE SIZE OF DLTXPI AND DLTYPI IN MOD_MAIN",&
       "BUT REALLY, WHAT IS WRONG WITH YOUR MESH?")

  DO I=1,M
     DO J=1,NTSN(I)-1
        I1=NBSN(I,J)
        I2=NBSN(I,J+1)

        DLTYTRIE(i,j) = VY(I1)-VY(I2)
        DLTXTRIE(i,j) = VX(I2)-VX(I1)

     END DO
  END DO

  ! Set the distance between Nodes for the North Pole region

  ! Set the distance between triangle edge centers
  IF(MAXVAL(NTVE) > 13) CALL FATAL_ERROR &
       & ("THERE ARE MORE THAN 12 CELLS AROUND ONE NODE:",&
       "PLEASE INCREASE THE SIZE OF DLVISCXPI AND DLVISCYPI IN MOD_MAIN",&
       "BUT REALLY, WHAT IS WRONG WITH YOUR MESH?")


  DO I=1,M
     DO J=1,NTVE(I)
        I1=NBVE(I,J)
        JTMP=NBVT(I,J)
        J1=JTMP+1-(JTMP+1)/4*3
        J2=JTMP+2-(JTMP+2)/4*3
        X11=0.5_SP*(VX(I)+VX(NV(I1,J1)))
        Y11=0.5_SP*(VY(I)+VY(NV(I1,J1)))
!        X22=XC(I1)
!        Y22=YC(I1)
        X33=0.5_SP*(VX(I)+VX(NV(I1,J2)))
        Y33=0.5_SP*(VY(I)+VY(NV(I1,J2)))


        DLTYECEC(I,J)=(Y11-Y33)
        DLTXECEC(I,J)=(X33-X11)



! Set the distance between the node and the edge Center
! NOTE: THE SIGN MATTERS!
         DLTYNEC(I,J)=(VY(I)-Y11)
         DLTXNEC(I,J)=(X11-VX(I))


     END DO
  END DO


   END SUBROUTINE edge_len
