










!
!     SWAN - SERVICE ROUTINES
!
!  Contents of this file
!
!     READXY
!     REFIXY
!     INFRAM
!     DISTR
!     KSCIP1
!     AC2TST
!     CVCHEK                                                              30.60
!     CVMESH                                                              30.60
!     NEWTON                                                              30.60
!     EVALF                                                               30.60
!     SWOBST                                                              30.60
!     TCROSS                                                              40.04
!     SWTRCF
!     SSHAPE                                                              40.00
!     SINTRP                                                              40.00
!     HSOBND                                                              32.01
!     CHGBAS                                                              40.00
!     GAMMA                                                               40.00
!     WRSPEC                                                              40.00
!TIMG!     SWTSTA                                                              40.23
!TIMG!     SWTSTO                                                              40.23
!TIMG!     SWPRTI                                                              40.23
!     TXPBLA                                                              40.23
!     INTSTR                                                              40.23
!     NUMSTR                                                              40.23
!     SWCOPI                                                              40.23
!     SWCOPR                                                              40.23
!     SWI2B                                                               40.30
!     SWR2B                                                               40.30
!
!  functions:
!  ----------
!  DEGCNV  (converts from cartesian convention to nautical and            32.01
!           vice versa)                                                   32.01
!  ANGRAD  (converts radians to degrees)                                  32.01
!  ANGDEG  (converts degrees to radians)                                  32.01
!
!  subroutines:
!  ------------
!  HSOBND  (Hs is calculated after a SWAN computation at all sides.       32.01
!           The calculated wave height from SWAN is then compared with    32.01
!           the wave heigth as provided by the user                       32.01
!
!***********************************************************************
!                                                                      *
   SUBROUTINE READXY (NAMX, NAMY, XX, YY, KONT, XSTA, YSTA)
!                                                                      *
!***********************************************************************
!
   USE OCPCOMM1                                                        
   USE SWCOMM2                                                         
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!     40.22: John Cazes and Tim Campbell
!     40.13: Nico Booij
!     40.51: Marcel Zijlema
!
!  1. UPDATE
!
!       Nov. 1996               offset values are added to standard values
!                               because they will be subtracted later
!     40.13, Nov. 01: a valid value for YY is required if a valid value
!                     for XX has been given; ocpcomm1.inc reactivated
!     40.51, Feb. 05: correction to location points equal to offset values
!
!  2. PURPOSE
!
!       Read x and y, initialize offset values XOFFS and YOFFS
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NAMX, NAMY   inp char    names of the two coordinates as given in
!                                the user manual
!       XX, YY       out real    values of x and y taking into account offset
!       KONT         inp char    what to be done if values are missing
!                                see doc. of INDBLE (Ocean Pack doc.)
!       XSTA, YSTA   inp real    standard values of x and y
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
!       INDBLE (Ocean Pack)
   LOGICAL EQREAL
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Read x and y in double prec.
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
   DOUBLE PRECISION XTMP, YTMP
   CHARACTER  NAMX *(*), NAMY *(*), KONT *(*)
   SAVE  IENT
   DATA  IENT /0/
   CALL  STRACE (IENT,'READXY')
!
!JQI   CALL INDBLE (NAMX, XTMP, KONT, DBLE(XSTA)+DBLE(XOFFS))
   IF(CHGVAL)THEN                                                    
!       a valid value was given for XX                                    
!JQI     CALL INDBLE (NAMY, YTMP, 'REQ', DBLE(YSTA)+DBLE(YOFFS))           
   ELSE                                                                
!JQI     CALL INDBLE (NAMY, YTMP, KONT, DBLE(YSTA)+DBLE(YOFFS))
   ENDIF                                                               
   IF(.NOT.LXOFFS)THEN
     XOFFS = REAL(XTMP)
     YOFFS = REAL(YTMP)
     LXOFFS = .TRUE.
   ENDIF
!JQI   IF(.NOT.EQREAL(XOFFS,REAL(XTMP)))THEN                             
!JQI     XX = REAL(XTMP-DBLE(XOFFS))
!JQI   ELSE IF(OPTG == 3)THEN                                            
!JQI     XX = 1.E-5                                                       
!JQI   ELSE                                                                
!JQI     XX = 0.                                                          
!JQI   END IF                                                              
!JQI   IF(.NOT.EQREAL(YOFFS,REAL(YTMP)))THEN                             
!JQI     YY = REAL(YTMP-DBLE(YOFFS))
!JQI   ELSE IF(OPTG == 3)THEN                                            
!JQI     YY = 1.E-5                                                       
!JQI   ELSE                                                                
!JQI     YY = 0.                                                          
!JQI   END IF                                                              
!
   RETURN
   END SUBROUTINE READXY
 

!****************************************************************
!
   SUBROUTINE TXPBLA(TEXT,IF,IL)
!
!****************************************************************
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     determines the position of the first and the last non-blank
!     (or non-tabulator) character in the text-string
!
!  4. Argument variables
!
!     IF          position of the first non-blank character in TEXT
!     IL          position of the last non-blank character in TEXT
!     TEXT        text string
!
   INTEGER IF, IL
   CHARACTER*(*) TEXT
!
!  6. Local variables
!
!     FOUND :     TEXT is found or not
!     ITABVL:     integer value of tabulator character
!     LENTXT:     length of TEXT
!
   INTEGER LENTXT, ITABVL
   LOGICAL FOUND
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
!DOS      ITABVL = ICHAR('	')
!UNIX      ITABVL = ICHAR('	')
   LENTXT = LEN (TEXT)
   IF = 1
   FOUND = .FALSE.
100 IF(IF <= LENTXT .AND. .NOT. FOUND)THEN
     IF(.NOT. (TEXT(IF:IF) ==  ' ' .OR.                   &
       ICHAR(TEXT(IF:IF)) == ITABVL))THEN
       FOUND = .TRUE.
     ELSE
       IF = IF + 1
     ENDIF
     GOTO 100
   ENDIF
   IL = LENTXT + 1
   FOUND = .FALSE.
200 IF(IL > 1 .AND. .NOT. FOUND)THEN
    IL = IL - 1
     IF(.NOT. (TEXT(IL:IL) ==  ' ' .OR.                   &
       ICHAR(TEXT(IL:IL)) == ITABVL))THEN
       FOUND = .TRUE.
     ENDIF
     GOTO 200
   ENDIF

   RETURN
   END SUBROUTINE TXPBLA
 
!****************************************************************
!
    CHARACTER*20 FUNCTION NUMSTR ( IVAL, RVAL, FORM )
!
!****************************************************************
!
    USE OCPCOMM4                                                        
!
    IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Convert integer or real to string with given format
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!     FORM        given format
!     RVAL        real to be converted
!
    INTEGER   IVAL
    REAL      RVAL
    CHARACTER FORM*20
!
!  6. Local variables
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
    IF(IVAL /= INAN)THEN
      WRITE(NUMSTR,FORM) IVAL
    ELSE IF( RVAL /= RNAN)THEN
      WRITE (NUMSTR,FORM) RVAL
    ELSE
      NUMSTR = ''
    END IF

    RETURN
    END FUNCTION NUMSTR
!****************************************************************
!
   SUBROUTINE SWCOPI ( IARR1, IARR2, LENGTH )
!
!****************************************************************
!
   USE OCPCOMM4                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copies integer array IARR1 to IARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IARR1       source array
!     IARR2       target array
!     LENGTH      array length
!
   INTEGER LENGTH
   INTEGER IARR1(LENGTH), IARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
   INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
   SAVE IENT
   DATA IENT/0/
   IF(LTRACE) CALL STRACE (IENT,'SWCOPI')

!  --- check array length

   IF(LENGTH <= 0)THEN
     CALL MSGERR( 3, 'Array length should be positive' )
   END IF

!  --- copy elements of array IARR1 to IARR2

   DO I = 1, LENGTH
     IARR2(I) = IARR1(I)
   END DO

   RETURN
   END SUBROUTINE SWCOPI
 

!***********************************************************************
!                                                                      *
    SUBROUTINE KSCIP1(MMT,SIG,D,K,CG,N,ND)
!                                                                      *
!***********************************************************************
!
    USE OCPCOMM4                                                        
    USE SWCOMM3                                                         
!
    IMPLICIT NONE
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!  0. Authors
!
!  1. Updates
!
!  2. Purpose
!
!     Calculation of the wave number, group velocity, group number N
!     and the derivative of N w.r.t. depth (=ND)
!
!  3. Method
!
!     --
!
!  4. Argument variables
!
!     MMT     input    number of frequency points
!
      INTEGER   MMT
!
!     CG      output   group velocity
!     D       input    local depth
!     K       output   wave number
!     N       output   ratio of group and phase velocity
!     ND      output   derivative of N with respect to D
!     SIG     input    rel. frequency for which wave parameters
!                      must be determined
!
      REAL      CG(MMT), D,                                    &
                K(MMT), N(MMT), ND(MMT), SIG(MMT)
!
!  6. Local variables
!
!     C         phase velocity
!     FAC1      auxiliary factor
!     FAC2      auxiliary factor
!     FAC3      auxiliary factor
!     IENT      number of entries
!     IS        counter in frequency (sigma-space)
!     KND       dimensionless wave number
!     ROOTDG    square root of D/GRAV
!     WGD       square root of GRAV*D
!     SND       dimensionless frequency
!     SND2      = SND*SND
!
      INTEGER   IENT, IS
      REAL      KND, ROOTDG, SND, WGD, SND2, C, FAC1, FAC2, FAC3
!
!  8. Subroutines used
!
!     --
!
!  9. Subroutines calling
!
!     SWOEXA, SWOEXF (Swan/Output)
!
! 10. Error messages
!
!     --
!
! 11. Remarks
!
!     --
!
! 12. Structure
!
!     -----------------------------------------------------------------
!      Compute non-dimensional frequency SND
!      IF SND >= 2.5, then
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND according to
!        deep water theory
!      ELSE IF SND =< 1.e-6
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND
!        according to extremely shallow water
!      ELSE
!        Compute wave number K, group velocity CGO and the ratio of
!        group and phase velocity N by Pade and other simple formulas.
!        Compute the derivative of N w.r.t. D = ND.
!     -----------------------------------------------------------------
!
! 13. Source text
!
!    SAVE IENT
!    DATA IENT /0/
!    IF (LTRACE) CALL STRACE (IENT, 'KSCIP1')
!
    ROOTDG = SQRT(D/GRAV_W)                                               
    WGD    = ROOTDG*GRAV_W                                                
    DO IS = 1, MMT
!     SND is dimensionless frequency
      SND = SIG(IS) * ROOTDG
      IF(SND >= 2.5)THEN
!     ******* deep water *******
        K(IS)  = SIG(IS) * SIG(IS) / GRAV_W                               
        CG(IS) = 0.5 * GRAV_W / SIG(IS)                                   
        N(IS)  = 0.5
        ND(IS) = 0.
      ELSE IF(SND < 1.E-6)THEN
!     *** very shallow water ***  
!        print*,'IN VERY SHALLOW WATER'
!	stop                                      
        K(IS)  = SND/D                                                  
        CG(IS) = WGD
        N(IS)  = 1.
        ND(IS) = 0.
      ELSE
        SND2  = SND*SND                                                 
        C     = SQRT(GRAV_W*D/(SND2+1./(1.+0.666*SND2+0.445*SND2**2      &
	        -0.105*SND2**3+0.272*SND2**4)))   
        K(IS) = SIG(IS)/C                                               
        KND   = K(IS)*D                                                 
        FAC1  = 2.*KND/SINH(2.*KND)                                     
        N(IS) = 0.5*(1.+FAC1)                                           
        CG(IS)= N(IS)*C                                                 
        FAC2  = SND2/KND                                                
        FAC3  = 2.*FAC2/(1.+FAC2*FAC2)                                  
        ND(IS)= FAC1*(0.5/D - K(IS)/FAC3)                               
      ENDIF
    END DO  

    RETURN
    END SUBROUTINE KSCIP1
!***********************************************************************
!                                                                      *
    SUBROUTINE KSCIP2(S,D,CG)
!                                                                      *
!***********************************************************************
!
    USE OCPCOMM4
    USE SWCOMM3
!
    IMPLICIT NONE
      REAL      CG(MSC), D, K,S(MSC)
      INTEGER   IENT, IS
      REAL      KND, ROOTDG, SND, WGD, SND2, C, FAC1,N
    ROOTDG = SQRT(D/GRAV_W)
    WGD    = ROOTDG*GRAV_W
!     SND is dimensionless frequency
     DO IS=1,MSC
      SND = S(IS) * ROOTDG
      IF(SND >= 2.5)THEN
!     ******* deep water *******
        K = S(IS) * S(IS) / GRAV_W
        CG(IS) = 0.5 * GRAV_W / S(IS)
      ELSE IF(SND < 1.E-6)THEN
!     *** very shallow water ***  
!        print*,'IN VERY SHALLOW WATER'
!       stop                                      
        K  = SND/D
        CG(IS) = WGD
      ELSE
        SND2  = SND*SND
        C     = SQRT(GRAV_W*D/(SND2+1./(1.+0.666*SND2+0.445*SND2**2      &
                -0.105*SND2**3+0.272*SND2**4)))
        K = S(IS)/C
        KND   = K*D
        FAC1  = 2.*KND/SINH(2.*KND)
        N = 0.5*(1.+FAC1)
        CG(IS)= N*C
      ENDIF
    END DO
    RETURN
    END SUBROUTINE KSCIP2



!****************************************************************

!****************************************************************
!
      CHARACTER*20 FUNCTION INTSTR ( IVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Convert integer to string
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!
      INTEGER IVAL
!
!  6. Local variables
!
!     CVAL  :     character represented an integer of mantisse
!     I     :     counter
!     IPOS  :     position in mantisse
!     IQUO  :     whole quotient
!
      INTEGER I, IPOS, IQUO
      CHARACTER*1, ALLOCATABLE :: CVAL(:)
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IPOS = 1
 100  CONTINUE
      IF (IVAL/10**IPOS.GE.1.) THEN
         IPOS = IPOS + 1
         GO TO 100
      END IF
      ALLOCATE(CVAL(IPOS))

      DO I=IPOS,1,-1
         IQUO=IVAL/10**(I-1)
         CVAL(IPOS-I+1)=CHAR(INT(IQUO)+48)
         IVAL=IVAL-IQUO*10**(I-1)
      END DO

      WRITE (INTSTR,*) (CVAL(I), I=1,IPOS)

      RETURN
      END FUNCTION INTSTR

!
!********************************************************************
!                                                                   *
      REAL FUNCTION GAMMA(XX)
!                                                                   *
!********************************************************************
!
!   Updates
!     ver 30.70, Oct 1997 by N.Booij: new subroutine
!
!   Purpose
!     Compute the transcendental function Gamma
!
!   Subroutines used
!     GAMMLN  (Numerical Recipes)
!
      REAL XX, YY, ABIG                                                   
      SAVE IENT, ABIG
      DATA IENT /0/, ABIG /30./
      CALL STRACE (IENT, 'GAMMA')
      YY = GAMMLN(XX)
      IF (YY.GT.ABIG) YY = ABIG
      IF (YY.LT.-ABIG) YY = -ABIG
      GAMMA = EXP(YY)
      RETURN
      END FUNCTION GAMMA

!********************************************************************
!                                                                   *
      FUNCTION GAMMLN(XX)
!                                                                   *
!********************************************************************
!
!   Method:
!     function is copied from: Press et al., "Numerical Recipes"
!
      DOUBLE PRECISION  COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,     &
          -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END FUNCTION GAMMLN

!
!***********************************************************************
!                                                                      *
      SUBROUTINE SSHAPE (ACLOC, SPCSIG, SPCDIR, FSHAPL, DSHAPL)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        
      USE SWCOMM3                                                         
!
!  0. Authors
!
!  1. Updates
!
!  2. Purpose
!
!     Calculating of energy density at boundary point (x,y,sigma,theta)
!
!  3. Method (updated...)
!
!     see: M. Yamaguchi: Approximate expressions for integral properties
!          of the JONSWAP spectrum; Proc. JSCE, No. 345/II-1, pp. 149-152,
!          1984.
!
!     computation of mean period: see Swan system documentation
!
!  4. Argument variables
!
!   o ACLOC : Energy density at a point in space
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL    ACLOC(MDC,MSC)
      REAL    SPCDIR(MDC,6)                                               
      REAL    SPCSIG(MSC)                                                 
!
! i   DSHAPL: Directional distribution
! i   FSHAPL: Shape of spectrum:
!             =1; Pierson-Moskowitz spectrum
!             =2; Jonswap spectrum
!             =3; bin
!             =4; Gauss curve
!             (if >0: period is interpreted as peak per.
!              if <0: period is interpreted as mean per.)
!
      INTEGER FSHAPL, DSHAPL                                              
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!     LSHAPE   absolute value of FSHAPL
!
      INTEGER  ID, IS, LSHAPE
!
!     PKPER    peak period                                                30.80
!     APSHAP   aux. var. used in computation of spectrum
!     AUX1     auxiliary variable
!     AUX2     auxiliary variable
!     AUX3     auxiliary variable
!     COEFF    coefficient for behaviour around the peak (Jonswap)
!     CPSHAP   aux. var. used in computation of spectrum
!     CTOT     total energy
!     CTOTT    total energy (used for comparison)
!     DIFPER   auxiliary variable used to select bin closest
!              to given frequency
!     MPER
!     MS       power in directional distribution
!     RA       action density
!     SALPHA
!     SF       frequency (Hz)
!     SF4      SF**4
!     SF5      SF**5
!     FPK      frequency corresponding to peak period (1/PKPER)           30.80
!     FPK4     FPK**4
!     SYF      peakedness parameter
!
      REAL     APSHAP, AUX1, AUX2, AUX3
      REAL     COEFF ,SYF   ,MPER  ,CTOT  ,CTOTT,PKPER  ,DIFPER
      REAL     MS
      REAL     RA    ,SALPHA,SF   ,SF4   ,SF5   ,FPK   ,FPK4, FAC
!
!     LOGPM    indicates whether peak or mean frequency is used
!     DVERIF   logical used in verification of incident direction
!
      LOGICAL  LOGPM, DVERIF                                              
!
!     PSHAPE   coefficients of spectral distribution (see remarks)
!     SPPARM   array containing integral wave parameters (see remarks)
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
! 10. Error messages
!
! 11. Remarks
!
!     PSHAPE(1): SY0, peak enhancement factor (gamma) in Jonswap spectrum
!     PSHAPE(2): spectral width in case of Gauss spectrum in rad/s
!
!     SPPARM    real     input    incident wave parameters (Hs, Period,
!                                 direction, Ms (dir. spread))
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by the user (either peak or mean)      30.80
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!
!     ---------------------------------------------------------------------
!
!     In the case of a JONSWAP spectrum the initial conditions are given by
!                   _               _       _       _       _
!                  |       _   _ -4  |     |       | S - S   |
!             2    |      |  S  |    |     |       |      p  |
!          a g     |      |  _  |    |  exp|-1/2 * |________ |* 2/pi COS(T-T  )
! E(S,D )= ___  exp|-5/4 *|  S  |    | G   |       | e * S   |              wi
!      wa    5     |      |   p |    |     |_      |_     p _|
!           S      |      |_   _|    |
!                  |_               _|
!
!   where
!         S   : rel. frequency
!
!         D   : Dir. of wave component
!          wa
!
!         a   : equili. range const. (Phillips' constant)
!         g   : gravity acceleration
!
!         S   : Peak frequency
!          p
!
!         G   : Peak enhancement factor
!         e   : Peak width
!
!         T   : local wind direction
!          wi
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       case shape
!       =1:   calculate value of Pierson-Moskowitz spectrum
!       =2:   calculate value of Jonswap spectrum
!       =3:   calculate value of bin spectrum
!       =4:   calculate value of Gauss spectrum
!       else: Give error message because of wrong shape
!       ----------------------------------------------------------------
!       if LOGPM is True
!       then calculate average period
!            if it differs from given average period
!            then recalculate peak period
!                 restart procedure to compute spectral shape
!       ----------------------------------------------------------------
!       for all spectral bins do
!            multiply all action densities by directional distribution
!       ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE(IENT,'SSHAPE')
!
      IF (ITEST >= 80) WRITE (PRTEST, 8) FSHAPL, DSHAPL,          &
                                         (SPPARM(JJ), JJ = 1,4)
   8  FORMAT (' entry SSHAPE ', 2I3, 4E12.4)
      IF(FSHAPL < 0)THEN
        LSHAPE = - FSHAPL
        LOGPM  = .FALSE.
      ELSE
        LSHAPE = FSHAPL
        LOGPM  = .TRUE.
      ENDIF

!
      IF(SPPARM(1) <= 0.)                                         &
         CALL MSGERR(1,' sign. wave height at boundary is not positive')  
!
      PKPER = SPPARM(2)
      ITPER = 0
      IF(LSHAPE == 3)THEN
!       select bin closest to given period
        DIFPER = 1.E10
        DO IS = 1, MSC
          IF(ABS(PKPER - PI2_W/SPCSIG(IS)) < DIFPER)THEN
            ISP = IS
            DIFPER = ABS(PKPER - PI2_W/SPCSIG(IS))
          ENDIF
        ENDDO
      ENDIF
!
!     compute spectral shape using peak period PKPER                      
!
      FAC  = 1.
 100  FPK  = (1./PKPER)                                                   
      FPK4 = FPK**4
      IF(LSHAPE == 1)THEN
        SALPHA = ((SPPARM(1) ** 2) * (FPK4)) * 5. / 16.
      ELSE IF(LSHAPE == 2)THEN
!       *** SALPHA = alpha*(grav**2)/(2.*pi)**4)
        SALPHA = (SPPARM(1)**2 * FPK4) /                            &
	         ((0.06533*(PSHAPE(1)**0.8015)+0.13467)*16.)
      ELSE IF(LSHAPE == 4)THEN
        AUX1 = SPPARM(1)**2 / ( 16.* SQRT (PI2_W) * PSHAPE(2))
        AUX3 = 2. * PSHAPE(2)**2
      ENDIF
!
      CTOTT = 0.
      DO IS = 1, MSC                                                  
!
        IF(LSHAPE == 1)THEN
!         *** LSHAPE = 1 : Pierson and Moskowitz ***
          SF = SPCSIG(IS) / PI2_W
          SF4 = SF**4
          SF5 = SF**5
          RA = (SALPHA/SF5)*EXP(-(5.*FPK4)/(4.*SF4))/(PI2_W*SPCSIG(IS))
          ACLOC(MDC,IS) = RA
        ELSE IF(LSHAPE == 2)THEN
!         *** LSHAPE = 2 : JONSWAP ***
          SF = SPCSIG(IS)/(PI2_W)
          SF4 = SF**4
          SF5 = SF**5
          CPSHAP = 1.25 * FPK4 / SF4
          IF(CPSHAP > 10.)THEN                                         
            RA = 0.
          ELSE
            RA = (SALPHA/SF5) * EXP(-CPSHAP)
          ENDIF
          IF(SF < FPK)THEN
            COEFF = 0.07
          ELSE
            COEFF = 0.09
          ENDIF
          APSHAP =  0.5 * ((SF-FPK) / (COEFF*FPK)) **2
          IF(APSHAP > 10.)THEN                                         
            SYF = 1.
          ELSE
            PPSHAP = EXP(-APSHAP)
            SYF = PSHAPE(1)**PPSHAP
          ENDIF
          RA = SYF*RA/(SPCSIG(IS)*PI2_W)
          ACLOC(MDC,IS) = RA
          IF(ITEST >= 120) WRITE (PRTEST, 112)                  &
	            SF, SALPHA, CPSHAP, APSHAP, SYF, RA
 112      FORMAT (' SSHAPE freq. ', 8E12.4)
        ELSE IF(LSHAPE == 3)THEN
!
!         *** all energy concentrated in one BIN ***
!
          IF(IS == ISP)THEN
            ACLOC(MDC,IS) = ( SPPARM(1)**2 ) /                   &
	                    ( 16. * SPCSIG(IS)**2 * FRINTF )
          ELSE
            ACLOC(MDC,IS) = 0.
          ENDIF
        ELSE IF(LSHAPE == 4)THEN
!
!         *** energy Gaussian distributed (wave-current tests) ***
!
          AUX2 = ( SPCSIG(IS) - ( PI2_W / PKPER ) )**2
          RA = AUX1 * EXP ( -1. * AUX2 / AUX3 ) / SPCSIG(IS)
          ACLOC(MDC,IS) = RA
        ELSE
          IF(IS == 1)THEN
            CALL MSGERR (2,'Wrong type for frequency shape')
            WRITE (PRINTF, *) ' -> ', FSHAPL, LSHAPE
          ENDIF
        ENDIF
        IF (ITEST >= 10)                                         &
	        CTOTT = CTOTT + FRINTF * ACLOC(MDC,IS) * SPCSIG(IS)**2
      END DO   !END DO IS
      IF(ITEST >= 10)THEN
        IF(SPPARM(1) > 0.01)THEN
          HSTMP = 4. * SQRT(CTOTT)
          IF(ABS(HSTMP-SPPARM(1)) > 0.1*SPPARM(1))           &
	      WRITE (PRINTF, 303) SPPARM(1), HSTMP
 303      FORMAT (' SSHAPE, deviation in Hs, should be ', F8.3,      &
                  ', calculated ', F8.3)
        ENDIF
      ENDIF
!
!     if mean frequency was given recalculate PKPER and restart
!
      IF (.NOT.LOGPM .AND. ITPER < 10)THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.
        AM1 = 0.
        DO IS = 1, MSC
          AS2 = ACLOC(MDC,IS) * (SPCSIG(IS))**2
          AS3 = AS2 * SPCSIG(IS)
          AM0 = AM0 + AS2
          AM1 = AM1 + AS3
        ENDDO
!       contribution of tail to total energy density
        PPTAIL = PWTAIL(1) - 1.                                           
        APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))              
        AM0 = AM0 * FRINTF + APTAIL * AS2                                 
        PPTAIL = PWTAIL(1) - 2.                                           
        EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))              
        AM1 = AM1 * FRINTF + EPTAIL * AS3                                 
!       Mean period:
        IF( AM1 /= 0.)THEN                                             
           MPER = PI2_W * AM0 / AM1
        ELSE                                                              
           CALL MSGERR(3, ' first moment is zero in calculating the')     
           CALL MSGERR(3, ' spectrum at boundary using param. bc.')       
        END IF                                                            
        IF(ITEST >= 80) WRITE (PRTEST, 72) ITPER, SPPARM(2), MPER,     &
	                                    PKPER
  72    FORMAT (' SSHAPE iter=', I2, '  period values:', 3F7.2)
        IF(ABS(MPER-SPPARM(2)) > 0.01*SPPARM(2))THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPARM(2) / MPER) * PKPER                              
          GOTO 100
        ENDIF
      ENDIF
!
      IF (ITPER >= 10)THEN
        CALL MSGERR(3, 'No convergence calculating the spectrum')         
        CALL MSGERR(3, 'at the boundary using parametric bound. cond.')   
      ENDIF
!
!     now introduce distribution over directions
!
      ADIR = PI_W * DEGCNV(SPPARM(3)) / 180.                                
      IF(DSHAPL == 1)THEN
        DSPR = PI_W * SPPARM(4) / 180.
        MS = MAX (DSPR**(-2) - 2., 1.)
      ELSE
        MS = SPPARM(4)
      ENDIF
      IF(MS < 12.)THEN
        CTOT = (2.**MS) * (GAMMA(0.5*MS+1.))**2 / (PI_W * GAMMA(MS+1.))
      ELSE
        CTOT =  SQRT (0.5*MS/PI_W) / (1. - 0.25/MS)
      ENDIF
      IF(ITEST >= 100)THEN
        ESOM = 0.
        DO IS = 1, MSC
          ESOM = ESOM + FRINTF * SPCSIG(IS)**2 * ACLOC(MDC,IS)
        ENDDO
        WRITE (PRTEST, *) ' SSHAPE dir ', 4.*SQRT(ABS(ESOM)),         &
	       SPPARM(1), CTOT, MS, GAMMA(0.5*MS+1.), GAMMA(MS+1.),   &
	       CTOT                                                     
      ENDIF
      DVERIF = .FALSE.
      CTOTT = 0.
      DO ID = 1, MDC
        ACOS = COS(SPCDIR(ID,1) - ADIR)
        IF(ACOS > 0.)THEN
          CDIR = CTOT * MAX (ACOS**MS, 1.E-10)
          IF(.NOT.FULCIR)THEN
            IF(ACOS >= COS(DDIR)) DVERIF = .TRUE.
          ENDIF
        ELSE
          CDIR = 0.
        ENDIF
        IF(ITEST >= 10) CTOTT = CTOTT + CDIR * DDIR
        IF(ITEST >= 100) WRITE (PRTEST, 360) ID,SPCDIR(ID,1),CDIR
 360    FORMAT (' ID Spcdir Cdir: ',I3,3(1X,E11.4))
        DO IS = 1, MSC
          ACLOC(ID,IS) = CDIR * ACLOC(MDC,IS)
        ENDDO
      ENDDO
      IF(ITEST >= 10)THEN
        IF(ABS(CTOTT-1.) > 0.1) WRITE (PRINTF, 363) CTOTT
 363    FORMAT (' SSHAPE, integral of Cdir is not 1, but:', F6.3)
      ENDIF
      IF (.NOT.FULCIR .AND. .NOT.DVERIF)                   &
         CALL MSGERR (1, 'incident direction is outside sector')
!
      RETURN

      END SUBROUTINE SSHAPE

!***********************************************************************
!                                                                      *
      REAL FUNCTION DEGCNV (DEGREE)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM3                                                         

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees from nautical to cartesian or vice versa.
!
!  3. METHOD
!
!       DEGCNV = 180 + dnorth - degree
!
!  4. PARAMETERLIST
!
!       DEGCNV      direction in cartesian or nautical degrees.
!       DEGREE      direction in nautical or cartesian degrees.
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!           Nautical convention           Cartesian convention
!
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
!  9. STRUCTURE
!
!     ---------------------------------
!     IF (NAUTICAL DEGREES) THEN
!       CONVERT DEGREES
!     IF (DEGREES > 360 OR < 0) THEN
!       CORRECT DEGREES WITHIN 0 - 360
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      INTEGER :: IENT
      REAL    :: DEGREE
      
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'DEGCNV')
!
      IF ( BNAUT ) THEN
        DEGCNV = 180. + DNORTH - DEGREE
      ELSE
        DEGCNV = DEGREE
      ENDIF
!
      IF (DEGCNV >= 360.) THEN
        DEGCNV = MOD (DEGCNV, 360.)
      ELSE IF (DEGCNV < 0.) THEN
        DEGCNV = MOD (DEGCNV, 360.) + 360.
      ELSE
!       DEGCNV between 0 and 360; do nothing
      ENDIF
!
      RETURN
      END FUNCTION DEGCNV

!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWANOUT(DEP2)    !AC2   ,SPCSIG,HSIBC,KGRPNT)                     
      RETURN
      END SUBROUTINE SWANOUT
