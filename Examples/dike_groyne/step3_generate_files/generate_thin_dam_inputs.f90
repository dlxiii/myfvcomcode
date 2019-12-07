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


module mod_sub
  use Mod_Spherical
  implicit none

  INTEGER, PARAMETER :: IPT = 6

  INTEGER NGL                !!GLOBAL NUMBER OF ELEMENTS
  INTEGER MGL                !!GLOBAL NUMBER OF NODES
  INTEGER IBFW_GL            !!GLOBAL NUMBER OF GROUNDWATER NODES
  INTEGER NUMQBC_GL          !!GLOBAL NUMBER OF FRESHWATER INFLOW NODES
  INTEGER NOBCGEO_GL         !!GLOBAL NUMBER OF OPEN BOUNDARY 
  INTEGER NDRFT_GL           !!GLOBAL NUMBER OF LAGRANGIAN TRACKING PARTICLES 

  INTEGER N                  !!LOCAL NUMBER OF ELEMENTS 
  INTEGER M                  !!LOCAL NUMBER OF NODES
  INTEGER IBFW               !!LOCAL NUMBER OF GROUNDWATER NODES
  INTEGER NUMQBC             !!LOCAL NUMBER OF FRESHWATER INFLOW NODES
  INTEGER NOBCGEO            !!LOCAL NUMBER OF OPEN BOUNDARY 
  INTEGER NDRFT              !!LOCAL NUMBER OF LAGRANGIAN TRACKING PARTICLES
  INTEGER NISBCE_1           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 1
  INTEGER NISBCE_2           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 2
  INTEGER NISBCE_3           !!LOCAL NUMBER OF ELEMENTS WITH ISBCE = 3

  INTEGER KB                 !!NUMBER OF SIGMA LEVELS
  INTEGER KBM1               !!NUMBER OF SIGMA LEVELS-1
  INTEGER KBM2               !!NUMBER OF SIGMA LEVELS-2
  INTEGER MYID               !!UNIQUE PROCESSOR ID (1 => NPROCS)
  INTEGER KSL                !!NUMBER OF STANDARD SEA LEVELS 
  INTEGER NPROCS             !!NUMBER OF PROCESSORS
  INTEGER NE                 !!NUMBER OF UNIQUE EDGES (LOCAL DOMAIN ONLY)
  INTEGER NCV                !!NUMBER OF INTERNAL CONTROL VOLUMES (EXTENDED LOCAL ONLY)

  INTEGER NCV_I              !!NUMBER OF INTERNAL CONTROL VOLUMES (LOCAL ONLY)
  INTEGER NT                 !!TOTAL OF LOCAL INTERNAL + HALO ELEMENTS
  INTEGER MT                 !!TOTAL OF LOCAL INTERNAL + HALO NODES
  INTEGER MX_NBR_ELEM        !!MAX NUMBER OF ELEMENTS SURROUNDING A NODE
  !----------------Node, Boundary Condition, and Control Volume-----------------------!

  INTEGER, ALLOCATABLE :: NV(:,:),NV2(:,:)    !!NODE NUMBERING FOR ELEMENTS
  INTEGER, ALLOCATABLE :: NBE(:,:)            !!INDICES OF ELMNT NEIGHBORS
  INTEGER, ALLOCATABLE :: NTVE(:)      
  INTEGER, ALLOCATABLE :: NTSN(:)      
  INTEGER, ALLOCATABLE :: ISONB(:)            !!NODE MARKER = 0,1,2
  INTEGER, ALLOCATABLE :: ISBC(:)     
  INTEGER, ALLOCATABLE :: ISBCE(:)     
  INTEGER, ALLOCATABLE :: IEC(:,:)
  INTEGER, ALLOCATABLE :: IENODE(:,:)
  INTEGER, ALLOCATABLE :: NBSN(:,:)
  INTEGER, ALLOCATABLE :: NIEC(:,:)
  INTEGER, ALLOCATABLE :: NTRG(:)
  INTEGER, ALLOCATABLE :: NBVE(:,:)
  INTEGER, ALLOCATABLE :: NBVT(:,:)
  INTEGER, ALLOCATABLE :: LISBCE_1(:)          !!LIST OF ELEMENTS WITH ISBCE=1
  INTEGER, ALLOCATABLE :: LISBCE_2(:)          !!LIST OF ELEMENTS WITH ISBCE=2
  INTEGER, ALLOCATABLE :: LISBCE_3(:)          !!LIST OF ELEMENTS WITH ISBCE=3
  REAL(SP),ALLOCATABLE :: DLTXC(:)
  REAL(SP),ALLOCATABLE :: DLTYC(:)
  REAL(SP),ALLOCATABLE :: DLTXYC(:)
  REAL(SP),ALLOCATABLE :: DLTXE(:)
  REAL(SP),ALLOCATABLE :: DLTYE(:)
  REAL(SP),ALLOCATABLE :: DLTXYE(:)
  REAL(SP),ALLOCATABLE :: SITAC(:) 
  REAL(SP),ALLOCATABLE :: SITAE(:) 
  REAL(SP),ALLOCATABLE :: XIJC(:) 
  REAL(SP),ALLOCATABLE :: YIJC(:)
  REAL(SP),ALLOCATABLE :: XIJE(:,:) 
  REAL(SP),ALLOCATABLE :: YIJE(:,:) 
  REAL(SP),ALLOCATABLE :: EPOR(:)            !!ELEMENT FLUX POROSITY (=0. IF ISBCE = 2)
  INTEGER, ALLOCATABLE :: IBCGEO(:)        !!LOCAL GEOSTROPHIC FRICTION CORRECTION NODES
  INTEGER, ALLOCATABLE :: N_ICELLQ(:,:)    !!FLUX ANGLE 

  !--------------------------Global Grid Variables------------------------------------!
  REAL(SP), ALLOCATABLE :: XG(:),XG2(:)              !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE :: YG(:),YG2(:)               !!GLOBAL X-COORD AT NODE 
  REAL(SP), ALLOCATABLE :: HG(:),HG2(:)               !!GLOBAL DEPTH AT NODE 
  REAL(SP), ALLOCATABLE :: XCG(:)              !!GLOBAL X-COORD AT FACE CENTER 
  REAL(SP), ALLOCATABLE :: YCG(:)              !!GLOBAL X-COORD AT FACE CENTER 
  !--------------------------Grid Metrics---------------------------------------------!

  REAL(SP)              :: VXMIN,VYMIN,VXMAX,VYMAX
  REAL(SP), ALLOCATABLE :: XC(:)               !!X-COORD AT FACE CENTER 
  REAL(SP), ALLOCATABLE :: YC(:)               !!Y-COORD AT FACE CENTER
  REAL(SP), ALLOCATABLE :: VX(:)               !!X-COORD AT GRID POINT
  REAL(SP), ALLOCATABLE :: VY(:)               !!Y-COORD AT GRID POINT
  REAL(SP), ALLOCATABLE :: ART(:)              !!AREA OF ELEMENT
  REAL(SP), ALLOCATABLE :: ART1(:)             !!AREA OF NODE-BASE CONTROl VOLUME
  REAL(SP), ALLOCATABLE :: ART2(:)             !!AREA OF ELEMENTS AROUND NODE

  CHARACTER(LEN=255):: FNAME1,FNAME2,CASENAME

  REAL                 :: above_msl ! dike top above the mean sea level 
  INTEGER              :: cell_count

  INTEGER              :: INTDM,INTDAMS,INTNODE
  INTEGER              :: NTDCELL_GL,NTDCELL
  INTEGER              :: NTDNODE_GL,NTDNODE_I
  INTEGER              :: NCROSS_GL,NCROSS_I
  INTEGER              :: NEND_GL,NEND_I
  INTEGER, DIMENSION(200,3000):: STRDAM  
  INTEGER,ALLOCATABLE,DIMENSION(:):: DAMNODE
  INTEGER,ALLOCATABLE,DIMENSION(:):: dyke
  !which is the node number of individual dam
  INTEGER, ALLOCATABLE :: I_TDCELL_GL(:),I_TDCELL_N(:)
  !which are neighboring cells along dams
  INTEGER, ALLOCATABLE :: I_TDNODE_GL(:),I_TDNODE_N(:)
  !which are all nodes of dams including crossing nodes and ending nodes
  INTEGER, ALLOCATABLE :: I_NCROSS_GL(:,:),I_NCROSS_N(:,:)
  !which is the crossing node of several dams, global and local
  ! I_NCROSS_GL(I,1:10)
  ! I_NCROSS_GL(I,1) GLOBAL ID OF CROSSING NODE
  ! I_NCROSS_GL(I,2) HOW MANY DAM STRINGS CROSS THIS NODE
  ! I_NCROSS_GL(I,3) HOW NANY NEIBORING DAM NODE
  ! I_NCROSS_GL(I,4:15) ID OF NEIBORING DAM NODE
  INTEGER, ALLOCATABLE :: I_NEND_GL(:),  I_NEND_N(:) 
  !which is the ending node of dam, global and local   

  !Declare the mirror matrix upon dyke node(non-crossing node)
  !scalar variables including temp,salinity,elvation etc.
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::DK_TEMP,DK_SAL,DK_EL
  !Declare the mirror matrix upon dyke node(crossing node)
  !scalar variables including temp,salinity,elvation etc.
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::DKC_TEMP,DKC_SAL,DKC_EL  

  !Declare the tge arrays,see the info in the heading of tge.F
  !non-crossing dyke node
  INTEGER, ALLOCATABLE :: DK_NBVE(:,:,:)
  INTEGER, ALLOCATABLE :: DK_NBVT(:,:,:)
  INTEGER, ALLOCATABLE :: DK_NBSN(:,:,:)
  INTEGER, ALLOCATABLE :: DK_NTVE(:,:)      
  INTEGER, ALLOCATABLE :: DK_NTSN(:,:)
  !crossing dyke node
  INTEGER, ALLOCATABLE :: DKC_NBVE(:,:,:)
  INTEGER, ALLOCATABLE :: DKC_NBVT(:,:,:)
  INTEGER, ALLOCATABLE :: DKC_NBSN(:,:,:)
  INTEGER, ALLOCATABLE :: DKC_NTVE(:,:)      
  INTEGER, ALLOCATABLE :: DKC_NTSN(:,:)

  INTEGER, DIMENSION(200,3000)    :: add_nodes
  INTEGER,DIMENSION(200) :: add_NODE_num
  INTEGER                 :: add_num

  INTEGER, DIMENSION(200,3000)    :: submerged_nodes
  INTEGER,DIMENSION(200) :: submerged_NODE_num
  INTEGER                 :: submerged_num

  INTEGER, DIMENSION(200,3000)    :: linear_nodes
  INTEGER,DIMENSION(200) :: linear_NODE_num
  INTEGER                 :: linear_num
contains


  SUBROUTINE READ_MESH
    !------------------------------------------------------------------------!
    !                       Original mesh reading                            !
    !------------------------------------------------------------------------!
    INTEGER:: I,J,TMP
    CHARACTER(LEN=3)::E3T
    CHARACTER(LEN=2)::ND
    ALLOCATE(NV(0:NGL,4))          ;NV       = 0  !!NODE NUMBERING FOR ELEMENTS
    ALLOCATE(NV2(0:NGL,4))         ;NV2      = 0  !!NODE NUMBERING FOR ELEMENTS
    ALLOCATE(XG(MGL),YG(MGL),HG(MGL))
    ALLOCATE(XG2(MGL*2),YG2(MGL*2),HG2(MGL*2))
    ALLOCATE(XCG(NGL),YCG(NGL))
    OPEN(1,FILE=TRIM(FNAME1))
    READ(1,*)
    DO I=1,NGL
       READ(1,*)E3T,TMP,NV(I,1:3)
    END DO
    DO I=1,MGL
       READ(1,*)ND,TMP,XG(I),YG(I),HG(I)
    END DO
    CLOSE(1)
    WRITE(IPT,*)'!MESH READING COMPLETED!'
    WRITE(IPT,*)'          '
    WRITE(IPT,*)'!--------------Global mesh infomation-----------------!'
    WRITE(IPT,*)'! NGL = ',NGL
    WRITE(IPT,*)'! MGL = ',MGL
    WRITE(IPT,*)'          '

    MT=MGL
    NT=NGL
    M=MGL
    N=NGL
    NV2=NV
    XG2=0;XG2(1:MGL)=XG(1:MGL)
    YG2=0;YG2(1:MGL)=YG(1:MGL)
    HG2=0;HG2(1:MGL)=HG(1:MGL)
    ALLOCATE(NBE(0:NT,3))         ;NBE      = 0  !!INDICES OF ELMNT NEIGHBORS
    ALLOCATE(ISONB(0:MT))         ;ISONB    = 0  !!NODE MARKER = 0,1,2
    ALLOCATE(ISBCE(0:NT))         ;ISBCE    = 0 
    ALLOCATE(NTVE(0:MT))          ;NTVE     = 0 
    ALLOCATE(NTSN(MT))            ;NTSN     = 0
    DO I=1,NGL
       XCG(I)=(XG(NV(I,1))+XG(NV(I,2))+XG(NV(I,3)))/3.0
       YCG(I)=(YG(NV(I,1))+YG(NV(I,2))+YG(NV(I,3)))/3.0
    END DO
  END SUBROUTINE READ_MESH




  !==========================================================================
  SUBROUTINE READ_NODESTRING
    !--------------------------------------------------------------------------|
    !  READ THE INFORMATION ABOUT THIN DAM                                     |
    !--------------------------------------------------------------------------|
    INTEGER :: I,J,NUM,processor
    CHARACTER(LEN=2)::NS
    INTEGER,DIMENSION(10000,10)::N_STRING,add_string,submerged_string,linear_string
    INTEGER, DIMENSION(200,3000)    :: TEMP,domain
    INTEGER,DIMENSION(200) :: N_NODE

    integer :: ntype1,ntype2,ntype3

    integer :: istart,iend,count,k,icount(2)
    logical :: is_start,is_end
    real(dp)::distance

    add_string=0
    add_node_num=0
    add_num=1



    OPEN(1,FILE='tst_ns_adding.dat',status='old')
    DO I=1,10000
       READ(1,*,END=998)NS,add_STRING(I,1:10)
       IF(MINVAL(add_STRING(I,:))>=0)THEN
          add_NODE_num(add_NUM)=add_NODE_num(add_NUM)+10
          add_nodes(add_NUM,add_NODE_num(add_NUM)-10+1:add_NODE_num(add_NUM))=add_STRING(I,1:10)
       ELSE
          DO J=1,10
             add_NODE_num(add_NUM)=add_NODE_num(add_NUM)+1
             add_nodes(add_NUM,add_NODE_num(add_NUM))=add_STRING(I,J)
             IF(add_STRING(I,J)<0)THEN
                add_num=add_num+1
                EXIT
             END IF
          END DO
       END IF
    END DO
998 CONTINUE
    CLOSE(1)

    submerged_string=0
    submerged_node_num=0
    submerged_num=1
    OPEN(1,FILE='tst_ns_submerged.dat',status='old')
    DO I=1,10000
       READ(1,*,END=997)NS,submerged_STRING(I,1:10)
       IF(MINVAL(submerged_STRING(I,:))>=0)THEN
          submerged_NODE_num(submerged_NUM)=submerged_NODE_num(submerged_NUM)+10
          submerged_nodes(submerged_NUM,submerged_NODE_num(submerged_NUM)-10+1:submerged_NODE_num(submerged_NUM))=submerged_STRING(I,1:10)
       ELSE
          DO J=1,10
             submerged_NODE_num(submerged_NUM)=submerged_NODE_num(submerged_NUM)+1
             submerged_nodes(submerged_NUM,submerged_NODE_num(submerged_NUM))=submerged_STRING(I,J)
             IF(submerged_STRING(I,J)<0)THEN
                submerged_num=submerged_num+1
                EXIT
             END IF
          END DO
       END IF
    END DO
997 CONTINUE
    CLOSE(1)
    do n=1,submerged_num
       do m=1,submerged_node_num(n)
          if(submerged_nodes(n,m)<0)submerged_nodes(n,m)=-submerged_nodes(n,m)
       end do
    end do

    linear_string=0
    linear_node_num=0
    linear_num=1
    OPEN(1,FILE='tst_ns_linear_groyne.dat',status='old')
    DO I=1,10000
       READ(1,*,END=996)NS,linear_STRING(I,1:10)
       IF(MINVAL(linear_STRING(I,:))>=0)THEN
          linear_NODE_num(linear_NUM)=linear_NODE_num(linear_NUM)+10
          linear_nodes(linear_NUM,linear_NODE_num(linear_NUM)-10+1:linear_NODE_num(linear_NUM))=linear_STRING(I,1:10)
       ELSE
          DO J=1,10
             linear_NODE_num(linear_NUM)=linear_NODE_num(linear_NUM)+1
             linear_nodes(linear_NUM,linear_NODE_num(linear_NUM))=linear_STRING(I,J)
             IF(linear_STRING(I,J)<0)THEN
                linear_num=linear_num+1
                EXIT
             END IF
          END DO
       END IF
    END DO
996 CONTINUE
    CLOSE(1)
    do n=1,linear_num
       do m=1,linear_node_num(n)
          if(linear_nodes(n,m)<0)linear_nodes(n,m)=-linear_nodes(n,m)
       end do
    end do

    N_STRING = 0
    N_NODE=0
    NUM=1
    OPEN(1,FILE='tst_ns_after_renumber.dat',status='old')
    DO I=1,10000
       READ(1,*,END=999)NS,processor,N_STRING(I,1:10)
       IF(MINVAL(N_STRING(I,:))>=0)THEN
          N_NODE(NUM)=N_NODE(NUM)+10
          TEMP(NUM,N_NODE(NUM)-10+1:N_NODE(NUM))=N_STRING(I,1:10)
          domain(NUM,N_NODE(NUM)-10+1:N_NODE(NUM))=processor
       ELSE
          DO J=1,10
             N_NODE(NUM)=N_NODE(NUM)+1
             TEMP(NUM,N_NODE(NUM))=N_STRING(I,J)
             domain(NUM,N_NODE(NUM))=processor
             IF(N_STRING(I,J)<0)THEN
                NUM=NUM+1
                EXIT
             END IF
          END DO
       END IF
    END DO
999 CONTINUE
    CLOSE(1)
    NUM=NUM-1
    OPEN(1,FILE=TRIM(CASENAME)//'_tdm.dat',status='replace')
    WRITE(1,*)'! THIN DAM'
    WRITE(1,*)NUM
    DO I=1,NUM
       DO J=1,N_NODE(I)
          WRITE(1,*)'TD',TEMP(I,J),domain(i,j)
       END DO
    END DO
    CLOSE(1)

    open(1,file=TRIM(CASENAME)//'_dam_cell.dat',status='replace')
    cell_count = 0
    do i=1,num
       do j=1,n_node(i)-1
          istart=TEMP(I,J)
          if(istart<0)istart= -istart
          iend  =TEMP(I,J+1)
          if(iend<0)iend= -iend
          count = 0
          do n=1,ngl
             is_start = .false.
             is_end   = .false.
             do k=1,3
                call arc(dble(xg(nv(n,k))),dble(yg(nv(n,k))),dble(xg(istart)),dble(yg(istart)),distance)
                if(distance<20.0)then
                   is_start = .true.
                end if
                call arc(dble(xg(nv(n,k))),dble(yg(nv(n,k))),dble(xg(iend)),dble(yg(iend)),distance)
                if(distance<20.0)then
                   is_end = .true.
                end if
             end do
             if(is_start.and.is_end)then
                count = count +1
                icount(count)=n
             end if
          end do
          cell_count = cell_count + 1
       end do
    end do

    write(1,"(i8)")cell_count
    do i=1,num
       do j=1,n_node(i)-1
          istart=TEMP(I,J)
          if(istart<0)istart= -istart
          iend  =TEMP(I,J+1)
          if(iend<0)iend= -iend
          count = 0
          do n=1,ngl
             is_start = .false.
             is_end   = .false.
             do k=1,3
                call arc(dble(xg(nv(n,k))),dble(yg(nv(n,k))),dble(xg(istart)),dble(yg(istart)),distance)
                if(distance<20.0)then
                   is_start = .true.
                end if
                call arc(dble(xg(nv(n,k))),dble(yg(nv(n,k))),dble(xg(iend)),dble(yg(iend)),distance)
                if(distance<20.0)then
                   is_end = .true.
                end if
             end do
             if(is_start.and.is_end)then
                count = count +1
                icount(count)=n
             end if
          end do
          write(1,"(3i8)")icount(1:2),domain(i,j)
       end do
    end do
    close(1)


  END SUBROUTINE READ_NODESTRING

  !==========================================================================
  SUBROUTINE READ_THINDAM
    !--------------------------------------------------------------------------|
    !  READ THE INFORMATION ABOUT THIN DAM                                     |
    !--------------------------------------------------------------------------|
    IMPLICIT NONE
    CHARACTER(LEN=2)             :: TD
    INTEGER                      :: INTDM 
    INTEGER                      :: I,J,K,N,NN
    INTEGER                      :: TOTAL
    LOGICAL                      :: ISSAME,ISBEGIN,ISEND
    !--------ASSUMING THE NUMBER OF DAMS IS LESS THAN 200,
    !--------THE NODES OF STRING IS LESS THAN 3000
    INTEGER, DIMENSION(200,3000)    :: TEMP1,TMP1  
    REAL(SP),DIMENSION(200,3000)    :: TEMP2,TMP2
    INTEGER, DIMENSION(200*3000,2)  :: TEMP3,TMP3 

    integer  , allocatable   :: node(:),node1(:,:),node2(:,:),node3(:,:)
    real(sp) , allocatable   :: dep1(:,:),dep2(:,:),dep3(:,:)
    integer  , dimension(10) :: tmp
    real(sp) , dimension(10) :: htmp
    integer                  :: ntype1,ntype2,ntype3,num_node,count
    logical                  :: exist
    real(dp)                 :: distance,relax_factor,relax_distance
    real(dp)                 :: x0,y0,x1,y1,x2,y2,dis1,dis2,dis0,height
    logical                  :: is_add
    logical                  :: is_submerged
    logical                  :: is_linear

    !   REWIND(INTDM)
    INTDM = 47
    OPEN(INTDM,FILE=TRIM(CASENAME)//'_tdm.dat',status='old')
    READ(INTDM,*)
    READ(INTDM,*)INTDAMS
    WRITE(IPT,*)'!-----------------DAM INFORMATION---------------------!'
    WRITE(IPT,*)'!Global IDs of nodes along Dams:',INTDAMS
    ALLOCATE(DAMNODE(INTDAMS))

    TOTAL = 0
    IF(INTDAMS>0)THEN
       !----------------read all nodes of dams-----------------------
       N=1
       DO I=1,INTDAMS
          DO J=1,3000
             READ(INTDM,*)TD,TEMP1(I,J)
             print*,TEMP1(I,J)
             !          IF(MSR)WRITE(IPT,*)TD,TEMP1(I,J),TEMP2(I,J)
             ISSAME = .FALSE.
             DO K=1,N-1
                IF(TEMP1(I,J)==TEMP3(K,1))THEN
                   TEMP3(K,2)=TEMP3(K,2)+1
                   ISSAME=.TRUE.
                   EXIT
                ENDIF
             END DO
             IF(.NOT.ISSAME)THEN
                TEMP3(N,1)=INT(ABS(TEMP1(I,J)))
                N=N+1
             ENDIF
             IF(TEMP1(I,J)<0)EXIT
          END DO
          IF(TEMP1(I,J)>0.AND.J>3000)THEN
             WRITE(IPT,*)'CHANGE FORMAT STATEMENT BELOW TO ACCOMODATE'
             WRITE(IPT,*)'MAXIMUM DAM NODES IN mod_thindam.F'
             STOP
          ELSEIF(TEMP1(I,J)<0)THEN
             TEMP1(I,J)=-1*TEMP1(I,J)
             DAMNODE(I)=J
             TOTAL = TOTAL +J
             TMP1(I,:)=TEMP1(I,:)
          ENDIF
       END DO
       STRDAM=TEMP1
       !-------------------------------------------------------------------------
       !   separate the crossed nodes from all nodes
       !--------------------------------------------------------------------------
       NTDNODE_GL=N-1
       ALLOCATE(I_TDNODE_GL (NTDNODE_GL))
       I_TDNODE_GL(1:NTDNODE_GL)=TEMP3(1:NTDNODE_GL,1)
       WRITE(IPT,'(5I7)')I_TDNODE_GL(1:NTDNODE_GL)
       WRITE(IPT,*)'!TOTAL NODES OF DAMS:',NTDNODE_GL
       TEMP3(:,2)=TEMP3(:,2)+1
       DO K=1,NTDNODE_GL
          IF(TEMP3(K,2)>1)NCROSS_GL=NCROSS_GL+1
       END DO
       WRITE(IPT,*)'!NCROSS_GL :',NCROSS_GL
       ALLOCATE(I_NCROSS_GL(NCROSS_GL,15))
       I_NCROSS_GL=0
       N=1
       DO K=1,NTDNODE_GL
          IF(TEMP3(K,2)>1)THEN
             I_NCROSS_GL(N,1)=TEMP3(K,1)
             I_NCROSS_GL(N,2)=TEMP3(K,2)
             !           IF(MSR)WRITE(IPT,*)' !CROSS NODE OF DAM: ',I_NCROSS_GL(N,1),I_NCROSS_GL(N,2)
             N=N+1           
          END IF
       END DO
       !--------FIND NEBORING DYKE NODE AROUND CROSSING NODES--------!
       DO I=1,INTDAMS 
          DO J=1,DAMNODE(I)
             DO N=1,NCROSS_GL
                IF(J>1.AND.J<DAMNODE(I))THEN
                   IF(TMP1(I,J+1)==I_NCROSS_GL(N,1))THEN
                      I_NCROSS_GL(N,3)=I_NCROSS_GL(N,3)+1
                      I_NCROSS_GL(N,3+I_NCROSS_GL(N,3))=TMP1(I,J)
                   END IF
                   IF(TMP1(I,J-1)==I_NCROSS_GL(N,1))THEN
                      I_NCROSS_GL(N,3)=I_NCROSS_GL(N,3)+1
                      I_NCROSS_GL(N,3+I_NCROSS_GL(N,3))=TMP1(I,J)
                   END IF
                END IF
             END DO
          END DO
       END DO

       DO N=1,NCROSS_GL
          WRITE(IPT,'(A20,15I6)')' !CROSSING NODE INFO:'&
               &,(I_NCROSS_GL(N,K),K=1,3+I_NCROSS_GL(N,3)) 
       END DO

       NEND_GL=INTDAMS*2
       DO I=1,INTDAMS       
          DO N=1,NCROSS_GL
             IF(TEMP1(I,1)==I_NCROSS_GL(N,1))NEND_GL=NEND_GL -1
             IF(TEMP1(I,DAMNODE(I))==I_NCROSS_GL(N,1))NEND_GL=NEND_GL -1
          END DO
       END DO
       WRITE(IPT,*)'!NEND_GL :',NEND_GL
       ALLOCATE(I_NEND_GL(NEND_GL));I_NEND_GL=0
       NN=1

       DO I=1,INTDAMS 
          ISBEGIN = .TRUE.
          ISEND   = .TRUE.      
          DO N=1,NCROSS_GL
             IF(TEMP1(I,1)==I_NCROSS_GL(N,1))ISBEGIN = .FALSE.
             IF(TEMP1(I,DAMNODE(I))==I_NCROSS_GL(N,1)) ISEND = .FALSE.
          END DO
          IF(ISBEGIN)THEN
             I_NEND_GL(NN)=TEMP1(I,1)
             NN=NN+1
          END IF
          IF(ISEND)THEN
             I_NEND_GL(NN)=TEMP1(I,DAMNODE(I))
             NN=NN+1
          END IF
       END DO
       !      DBNODE=NTDNODE_GL-NEND_GL-NCROSS_GL
    END IF

    !--------FIND NEBORING DYKE NODE AROUND CROSSING NODES--------!
    allocate(node(3000))
    allocate(node1(3000,2),node2(3000,3),node3(3000,4))
    allocate(dep1(3000,2),dep2(3000,3),dep3(3000,4))
    node1  = 0
    node2  = 0
    node3  = 0
    ntype1 = 0
    ntype2 = 0
    ntype3 = 0

    num_node = 0
    DO I=1,INTDAMS 
       DO J=1,DAMNODE(I)
          exist = .false.
          do k=1,num_node
             if(node(k)==TEMP1(i,j))exist = .true.
          end do
          if(exist)cycle
          num_node = num_node + 1
          node(num_node) = TEMP1(i,j)
       end do
    end do
    print*,'Total Dam Node = ',num_node
    do i=1,num_node
       count = 0
       tmp   = 0
       htmp  = 0
       do n=1,mgl
          call arc(dble(xg(n)),dble(yg(n)),dble(xg(node(i))),dble(yg(node(i))),distance)
          if(distance<20.0)then
             count = count + 1
             tmp(count) = n
             htmp(count) = hg(n)
          end if
       end do
       if(count == 2)then
          ntype1 = ntype1 + 1
          node1(ntype1,1:count)=tmp(1:count)
          dep1(ntype1,1:count)=htmp(1:count)
       end if
       if(count == 3)then
          ntype2 = ntype2 + 1
          node2(ntype2,1:count)=tmp(1:count)
          dep2(ntype2,1:count)=htmp(1:count)
       end if
       if(count == 4)then
          ntype3 = ntype3 + 1
          node3(ntype3,1:count)=tmp(1:count)
          dep3(ntype3,1:count)=htmp(1:count)
       end if
    end do
    relax_factor = 0.005
    relax_distance = 1000
    open(1,file=TRIM(CASENAME)//'_dam_node.dat',status='replace')
    open(2,file=TRIM(CASENAME)//'_dam_spg.dat',status='replace')
    write(1,*)' DAM TYPE 1'
    write(1,*)ntype1
    do i=1,ntype1
       is_add=.false.
       do n=1,add_num
          do m=1,add_node_num(n)
             if(node1(i,1)==add_nodes(n,m))is_add=.true.
             if(node1(i,2)==add_nodes(n,m))is_add=.true.
          end do
       end do

       is_submerged=.false.
       do n=1,submerged_num
          do m=1,submerged_node_num(n)
             if(node1(i,1)==submerged_nodes(n,m))is_submerged=.true.
             if(node1(i,2)==submerged_nodes(n,m))is_submerged=.true.
          end do
       end do

       is_linear=.false.
       do n=1,linear_num
          do m=1,linear_node_num(n)
             if(node1(i,1)==linear_nodes(n,m).and.node1(i,1)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node1(i,1))
                y0=yg(node1(i,1))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep1(i,:)-height)<0)then
                   write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),0.0,0.0
                else
                   write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),dep1(i,:)-height
                end if
             end if
             if(node1(i,2)==linear_nodes(n,m).and.node1(i,2)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node1(i,2))
                y0=yg(node1(i,2))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep1(i,:)-height)<0)then
                   write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),0.0,0.0
                else
                   write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),dep1(i,:)-height
                end if
             end if
          end do
       end do

       if(is_add)then
          write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),max(dep1(i,:)+4.37,0.0)
       elseif(is_submerged)then
          if(sum(dep1(i,:))>0.0)then
             write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),dep1(i,:)-5.0
          else
             write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),0.0,0.0
          end if
       elseif(is_linear==.false.)then
          if(sum(dep1(i,:))>0.0)then
             write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),dep1(i,:)+above_msl
          else
             write(1,"(2i6,2f8.3,f7.0,f7.4)")node1(i,:),max(dep1(i,:)+above_msl,0.0)
          end if
       end if
       write(2,"(i6,f7.0,f7.4)")node1(i,1),relax_distance,relax_factor
       write(2,"(i6,f7.0,f7.4)")node1(i,2),relax_distance,relax_factor
    end do
    write(1,*)' DAM TYPE 2'
    write(1,*)ntype2
    do i=1,ntype2
       is_add=.false.
       do n=1,add_num
          do m=1,add_node_num(n)
             if(node2(i,1)==add_nodes(n,m))is_add=.true.
             if(node2(i,2)==add_nodes(n,m))is_add=.true.
             if(node2(i,3)==add_nodes(n,m))is_add=.true.
          end do
       end do

       is_submerged=.false.
       do n=1,submerged_num
          do m=1,submerged_node_num(n)
             if(node2(i,1)==submerged_nodes(n,m))is_submerged=.true.
             if(node2(i,2)==submerged_nodes(n,m))is_submerged=.true.
             if(node2(i,3)==submerged_nodes(n,m))is_submerged=.true.
          end do
       end do


       is_linear=.false.
       do n=1,linear_num
          do m=1,linear_node_num(n)
             if(node2(i,1)==linear_nodes(n,m).and.node2(i,1)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node2(i,1))
                y0=yg(node2(i,1))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep2(i,:)-height)<0)then
                   write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),0.0,0.0,0.0
                else
                   write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),dep2(i,:)-height
                end if
             end if
             if(node2(i,2)==linear_nodes(n,m).and.node2(i,2)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node2(i,2))
                y0=yg(node2(i,2))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,m))
                y2=yg(linear_nodes(n,m))
                dis1=sqrt((x1-x0)**2+(y1-y0)**2)
                dis2=sqrt((x2-x0)**2+(y2-y0)**2)
                dis0=sqrt((x2-x1)**2+(y2-y1)**2)
                height=2.0-dis1/dis0*2.0
                if(sum(dep2(i,:)-height)<0)then
                   write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),0.0,0.0,0.0
                else
                   write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),dep2(i,:)-height
                end if
             end if
             if(node2(i,3)==linear_nodes(n,m).and.node2(i,3)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node2(i,3))
                y0=yg(node2(i,3))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep2(i,:)-height)<0)then
                   write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),0.0,0.0,0.0
                else
                   write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),dep2(i,:)-height
                end if
             end if
          end do
       end do

       if(is_add)then
          write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),max(dep2(i,:)+4.37,0.0)
       elseif(is_submerged)then
          if(sum(dep2(i,:))>0.0)then
             write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),dep2(i,:)-5.0
          else
             write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),0.0,0.0,0.0
          end if
       elseif(is_linear==.false.)then
          if(sum(dep2(i,:))>0.0)then
             write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),dep2(i,:)+above_msl
          else
             write(1,"(3i6,3f8.3,f7.0,f7.4)")node2(i,:),max(dep2(i,:)+above_msl,0.0)
          end if
       end if
       write(2,"(i6,f7.0,f7.4)")node2(i,1),relax_distance,relax_factor
       write(2,"(i6,f7.0,f7.4)")node2(i,2),relax_distance,relax_factor
       write(2,"(i6,f7.0,f7.4)")node2(i,3),relax_distance,relax_factor
    end do
    write(1,*)' DAM TYPE 3'
    write(1,*)ntype3
    do i=1,ntype3
       is_add=.false.
       do n=1,add_num
          do m=1,add_node_num(n)
             if(node3(i,1)==add_nodes(n,m))is_add=.true.
             if(node3(i,2)==add_nodes(n,m))is_add=.true.
             if(node3(i,3)==add_nodes(n,m))is_add=.true.
             if(node3(i,4)==add_nodes(n,m))is_add=.true.
          end do
       end do

       is_submerged=.false.
       do n=1,submerged_num
          do m=1,submerged_node_num(n)
             if(node3(i,1)==submerged_nodes(n,m))is_submerged=.true.
             if(node3(i,2)==submerged_nodes(n,m))is_submerged=.true.
             if(node3(i,3)==submerged_nodes(n,m))is_submerged=.true.
             if(node3(i,4)==submerged_nodes(n,m))is_submerged=.true.
          end do
       end do



       is_linear=.false.
       do n=1,linear_num
          do m=1,linear_node_num(n)
             if(node3(i,1)==linear_nodes(n,m).and.node3(i,1)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node3(i,1))
                y0=yg(node3(i,1))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep3(i,:)-height)<0)then
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),0.0,0.0,0.0,0.0
                else
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),dep3(i,:)-height
                end if
             end if
             if(node3(i,2)==linear_nodes(n,m).and.node3(i,2)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node3(i,2))
                y0=yg(node3(i,2))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep3(i,:)-height)<0)then
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),0.0,0.0,0.0,0.0
                else
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),dep3(i,:)-height
                end if
             end if
             if(node3(i,3)==linear_nodes(n,m).and.node3(i,3)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node3(i,3))
                y0=yg(node3(i,3))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep3(i,:)-height)<0)then
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),0.0,0.0,0.0,0.0
                else
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),dep3(i,:)-height
                end if
             end if
             if(node3(i,4)==linear_nodes(n,m).and.node3(i,4)/=linear_nodes(n,linear_node_num(n)))then
                is_linear=.true.
                x0=xg(node3(i,4))
                y0=yg(node3(i,4))
                x1=xg(linear_nodes(n,2))
                y1=yg(linear_nodes(n,2))
                x2=xg(linear_nodes(n,linear_node_num(n)))
                y2=yg(linear_nodes(n,linear_node_num(n)))
                call arc(x0,y0,x1,y1,dis1)
                call arc(x0,y0,x2,y2,dis2)
                call arc(x1,y1,x2,y2,dis0)
                height=2.0-dis1/dis0*1.5
                if(sum(dep3(i,:)-height)<0)then
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),0.0,0.0,0.0,0.0
                else
                   write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),dep3(i,:)-height
                end if
             end if
          end do
       end do


       if(is_add)then
          write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),max(dep3(i,:)+4.37,0.0)
       elseif(is_submerged)then
          if(sum(dep3(i,:))>0.0)then
             write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),dep3(i,:)-5.0
          else
             write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),0.0,0.0,0.0,0.0
          end if
       elseif(is_linear==.false.)then
          if(sum(dep3(i,:))>0.0)then
             write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),dep3(i,:)+above_msl
          else
             write(1,"(4i6,4f8.3,f7.0,f7.4)")node3(i,:),max(dep3(i,:)+above_msl,0.0)
          end if
       end if
       write(2,"(i6,f7.0,f7.4)")node3(i,1),relax_distance,relax_factor
       write(2,"(i6,f7.0,f7.4)")node3(i,2),relax_distance,relax_factor
       write(2,"(i6,f7.0,f7.4)")node3(i,3),relax_distance,relax_factor
       write(2,"(i6,f7.0,f7.4)")node3(i,4),relax_distance,relax_factor
    end do
    close(1)
    close(2)
  END SUBROUTINE READ_THINDAM


end module mod_sub





!========================================================================!
! Renumber_Nodestring: to check the dam nodes and cells after mesh       !
! renumbering in SMS.
!========================================================================!
program renumber_nodestring

use mod_sub

implicit none

CASENAME='tst'
FNAME1='./tst_dike_renumbered.2dm'
NGL= 9250
MGL= 4890

above_msl= -5.0 ! dike top is 5.0m deep below the mean sea level

CALL READ_MESH
CALL READ_NODESTRING
CALL READ_THINDAM

end program renumber_nodestring

