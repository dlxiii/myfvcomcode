MODULE Mod_Spherical
!-------------------------------Constants---------------------------------------------------!
INTEGER, PARAMETER                    :: DP    = SELECTED_REAL_KIND(12,300)
INTEGER, PARAMETER                    :: SP    = SELECTED_REAL_KIND(6,30)
REAL(SP), PARAMETER, DIMENSION(4) :: ALPHA_RK = (/0.2500_SP,0.333333_SP,0.5000_SP,1.0_SP/)
REAL(SP), PARAMETER :: GRAV      = 9.81_SP
REAL(SP), PARAMETER :: PI        = 3.141592653_SP
REAL(SP), PARAMETER :: PI2       = 6.283185307_SP
REAL(SP), PARAMETER :: ZERO      = 0.0_SP 
REAL(SP), PARAMETER :: ONE_THIRD = 1.0_SP/3.0_SP 
REAL(SP), PARAMETER :: REARTH    = 6371.0E03_SP   !!Earth Radius in Meters
REAL(SP), PARAMETER :: DEG2RAD   = PI2/360.0_SP   !!Radians/Degree
REAL(SP), PARAMETER :: TPI       = DEG2RAD*REARTH !TPI=pi*rearth/180.=3.14159265/180.0*6371.*1000.
REAL(SP), PARAMETER :: ROFVROS   = 0.9775171065_SP!!RATIO OF THE DENSITY OF FRESH AND SEA WATER 1000./1023.   

 
CONTAINS

   SUBROUTINE ARC(XX1,YY1,XX2,YY2,ARCL)
!----------------------------------------------------------------------------
!      function:
!           calculate the arc lenth for given two point on the spherical plane
!      input:
!           xx1,yy1,xx2,yy2 :are longitude and latitude of two points
!      output:
!           arcl :  arc lenth of two points in spherical plane
!-----------------------------------------------------------------------------       
       
!  solve the arc length through the earth center
   IMPLICIT NONE
   REAL(DP) :: X1,Y1,X2,Y2,XA,YA,ZA,XB,YB,ZB,AB,AOB,ARCL
   REAL(DP) :: XX1,YY1,XX2,YY2

   X1=XX1*DEG2RAD
   Y1=YY1*DEG2RAD

   X2=XX2*DEG2RAD
   Y2=YY2*DEG2RAD

   XA=COS(Y1)*COS(X1)
   YA=COS(Y1)*SIN(X1)
   ZA=SIN(Y1)

   XB=COS(Y2)*COS(X2)
   YB=COS(Y2)*SIN(X2)
   ZB=SIN(Y2)

   AB=SQRT((XB-XA)**2+(YB-YA)**2+(ZB-ZA)**2)
   AOB=(2.-AB*AB)/2.
   AOB=ACOS(AOB)
   ARCL=REARTH*AOB

   RETURN
   END SUBROUTINE ARC



   SUBROUTINE AREA(SIDE1,SIDE2,SIDE3,AREA1)
!--------------------------------------------------------------------
!      function:
!           calculate the area of a triangle on a spherical plane
!      input:
!           side1,side2 and side3: are 3 arc lenth for one triangle
!      output:
!           areal: is area of a triangle on a spherical plane
!--------------------------------------------------------------------
   IMPLICIT NONE
   REAL(DP) :: SIDE1,SIDE2,SIDE3,AREA1
   REAL(DP) :: PSUM,PM,QMJC

   SIDE1=SIDE1/REARTH
   SIDE2=SIDE2/REARTH
   SIDE3=SIDE3/REARTH
   IF(SIDE1 == 0. .OR. SIDE2 == 0. .OR. SIDE3 == 0.)THEN
     AREA1=0.
   ELSE
     PSUM=0.5*(SIDE1+SIDE2+SIDE3)
     PM=SIN(PSUM)*SIN(PSUM-SIDE1)*SIN(PSUM-SIDE2)*SIN(PSUM-SIDE3)
     PM=SQRT(PM)/(2.0*COS(SIDE1*0.5)*COS(SIDE2*0.5)*COS(SIDE3*0.5))
     QMJC = 2.0*ASIN(PM)

     AREA1=REARTH*REARTH*QMJC

   END IF

   RETURN
   END SUBROUTINE AREA

   SUBROUTINE ARCC(XX1,YY1,XX2,YY2,XXC,YYC)
   IMPLICIT NONE
   REAL(DP) :: XXC,YYC,XX1,YY1,XX2,YY2
   REAL(DP) :: X1,Y1,X2,Y2

   X1=XX1*DEG2RAD
   Y1=YY1*DEG2RAD

   X2=XX2*DEG2RAD
   Y2=YY2*DEG2RAD

   XXC=DCOS(Y1)*DSIN(X1)+DCOS(Y2)*DSIN(X2)
!   XXC=XXC/(COS(Y1)*COS(X1)+COS(Y2)*COS(X2))
!   XXC=ATAN(XXC)
   XXC=DATAN2(XXC,(DCOS(Y1)*DCOS(X1)+DCOS(Y2)*DCOS(X2)))
   XXC=XXC/DEG2RAD
 
!   IF(XXC .LT. 0.0) XXC=180.0+XXC
   IF(XXC < 0.0) XXC=360.0+XXC
   
   YYC=DCOS(Y1)*DCOS(Y1)+DCOS(Y2)*DCOS(Y2)+2.*DCOS(Y1)*DCOS(Y2)*DCOS(X1-X2)
!   YYC=SQRT(YYC)/(SIN(Y1)+SIN(Y2))
   YYC=DATAN2(DSQRT(YYC),(DSIN(Y1)+DSIN(Y2)))
!   YYC=ATAN(YYC)
   YYC=90.-YYC/DEG2RAD

   RETURN
   END SUBROUTINE ARCC


   SUBROUTINE ARCX(XX1,YY1,XX2,YY2,ARCX1)

   IMPLICIT NONE
   INTEGER I,NX
   PARAMETER(NX=500)
   REAL(DP) :: XX1,YY1,XX2,YY2,ARCX1
   REAL(DP) :: X1,Y1,X2,Y2,TY
   REAL(DP) :: XTMP	      

   IF(XX1 == XX2)THEN
     ARCX1=0.
   ELSE
     X1=XX1*DEG2RAD
     Y1=YY1*DEG2RAD

     X2=XX2*DEG2RAD
     Y2=YY2*DEG2RAD

     XTMP  = X2-X1
     IF(XTMP >  PI)THEN
       XTMP = -2*PI+XTMP
     ELSE IF(XTMP < -PI)THEN
       XTMP =  2*PI+XTMP
     END IF  

     TY=0.5*(Y2+Y1)
     ARCX1=REARTH*COS(TY)*XTMP
   END IF
   
   RETURN
   END SUBROUTINE ARCX

   SUBROUTINE ARCX_BACK(XX1,YY1,XX2,YY2,ARCX1)

   IMPLICIT NONE
   INTEGER I,NX
   PARAMETER(NX=500)
   REAL(DP) :: XX1,YY1,XX2,YY2,ARCX1
   REAL(DP) :: X1,Y1,X2,Y2,TY,A1,A2,B1,B2,C1,C2,A,B,C,X(NX+1),Y(NX+1)
   REAL(DP) :: XTMP	      

   IF(XX1 == XX2)THEN
     ARCX1=0.
   ELSE
     X1=XX1*DEG2RAD
     Y1=YY1*DEG2RAD

     X2=XX2*DEG2RAD
     Y2=YY2*DEG2RAD

     X(1)=X1
     Y(1)=Y1
     X(NX+1)=X2
     Y(NX+1)=Y2

     XTMP=X(NX+1)-X(1)
     IF(XTMP >  PI)THEN
       XTMP = -2*PI+XTMP
     ELSE IF(XTMP < -PI)THEN
       XTMP =  2*PI+XTMP
     END IF  

     DO I=2,NX
       X(I)=X(I-1)+XTMP/FLOAT(NX)
!       x(i)=x(i-1)+(x(nx+1)-x(1))/float(nx)
     END DO

     A1=COS(Y(1))*COS(X(1))
     A2=COS(Y(NX+1))*COS(X(NX+1))

     B1=COS(Y(1))*SIN(X(1))
     B2=COS(Y(NX+1))*SIN(X(NX+1))

     C1=SIN(Y(1))
     C2=SIN(Y(NX+1))

     A=A1*B2-A2*B1
     B=B1*C2-B2*C1
     C=A2*C1-A1*C2

     DO I=2,NX
       Y(I)=-B*COS(X(I))-C*SIN(X(I))
       Y(I)=Y(I)/A
       Y(I)=ATAN(Y(I))
     END DO

     ARCX1=0.
     DO I=1,NX
       TY=0.5*(Y(I)+Y(I+1))
       XTMP=X(I+1)-X(I)
       IF(XTMP >  PI)THEN
         XTMP = -2*PI+XTMP
       ELSE IF(XTMP < -PI)THEN
         XTMP =  2*PI+XTMP
       END IF  
       ARCX1=ARCX1+REARTH*COS(TY)*XTMP
!       arcx1=arcx1+rearth*cos(ty)*(x(i+1)-x(i))
     END DO
   END IF

   RETURN
   END SUBROUTINE ARCX_BACK

END MODULE Mod_Spherical