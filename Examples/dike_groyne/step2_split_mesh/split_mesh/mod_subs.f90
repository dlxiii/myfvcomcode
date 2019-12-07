MODULE MOD_SUBS
USE ALL_VARS
USE MOD_THINDAM
IMPLICIT NONE
SAVE

CONTAINS

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
OPEN(1,FILE=TRIM(FNAME1),STATUS='old')
READ(1,*)
DO I=1,NGL
   READ(1,*)E3T,TMP,NV(I,1:3)
END DO
DO I=1,MGL
!   READ(1,"(a2,i7,1x,e15.8,1x,e15.8,1x,e15.8)")ND,TMP,XG(I),YG(I),HG(I)
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

SUBROUTINE DIKE_MESH_GENERATE
!------------------------------------------------------------------------!
!                       Insert dike boundary nodes                       !
!------------------------------------------------------------------------!
INTEGER :: I,J,K,NODE,NN,II,JJ,KK,SUM,NUM,NCROSS
INTEGER :: PREV ! The previous node in current dike string
INTEGER :: NEXT ! The next node in current dike string
REAL(SP) :: ANGLE_K,ANGLE_P
! dike(n) 
!         =3 refers crossing node 
!         =2 refers dike node
!         =1 refers ending dike node
!         =0 MEAN NOT.  
NODE=0   ! nodes added to the mesh     

open(1,file=TRIM(CASENAME)//'_dam_node.dat',status='replace')

write(1,*)'DAM TYPE 1'
DO I=1,INTDAMS 
    DO J=1,DAMNODE(I)
        K=STRDAM(I,J)
        PRINT*,K,DIKE(K),ISONB(K)

        ! ending nodes and it's also a solid boundary node
        IF(DIKE(K)==1.AND.ISONB(K)==1)THEN
           NODE=NODE+1
           IF(J==1)THEN
             PREV=STRDAM(I,J)
             NEXT=STRDAM(I,J+1)
           ELSE
             PREV=STRDAM(I,J-1)
             NEXT=STRDAM(I,J)
           END IF
           DO NN=1,NTVE(K)
               ANGLE_K=ATAN((YG(NEXT)-YG(PREV))/(XG(NEXT)-XG(PREV)))   ! DIKE LINE
               IF((XG(NEXT)-XG(PREV))==0.AND.(YG(NEXT)-YG(PREV))>0)ANGLE_K=90
               IF((XG(NEXT)-XG(PREV))==0.AND.(YG(NEXT)-YG(PREV))<0)ANGLE_K=270
               IF((XG(NEXT)-XG(PREV))>0 .AND.(YG(NEXT)-YG(PREV))==0)ANGLE_K=360
               IF((XG(NEXT)-XG(PREV))<0 .AND.(YG(NEXT)-YG(PREV))==0)ANGLE_K=180
               IF(YG(NEXT)-YG(PREV)>0.AND.XG(NEXT)-XG(PREV)>0)ANGLE_K=ANGLE_K/PI*180.0      !in the first quadrant
               IF(YG(NEXT)-YG(PREV)<0.AND.XG(NEXT)-XG(PREV)>0)ANGLE_K=360+ANGLE_K/PI*180.0  !in the forth quadrant
               IF(YG(NEXT)-YG(PREV)>0.AND.XG(NEXT)-XG(PREV)<0)ANGLE_K=180+ANGLE_K/PI*180.0  !in the second quadrant
               IF(YG(NEXT)-YG(PREV)<0.AND.XG(NEXT)-XG(PREV)<0)ANGLE_K=180+ANGLE_K/PI*180.0  ! in the third quadrant
               ANGLE_P=ATAN((YCG(NBVE(K,NN))-YG(PREV))/(XCG(NBVE(K,NN))-XG(PREV)))   ! Face center of the cell which attached to the dike line
               IF((XCG(NBVE(K,NN))-XG(PREV))==0.AND.(YCG(NBVE(K,NN))-YG(PREV))>0)ANGLE_P=90
               IF((XCG(NBVE(K,NN))-XG(PREV))==0.AND.(YCG(NBVE(K,NN))-YG(PREV))<0)ANGLE_P=270
               IF((XCG(NBVE(K,NN))-XG(PREV))>0 .AND.(YCG(NBVE(K,NN))-YG(PREV))==0)ANGLE_P=360
               IF((XCG(NBVE(K,NN))-XG(PREV))<0 .AND.(YCG(NBVE(K,NN))-YG(PREV))==0)ANGLE_P=180
               IF(YCG(NBVE(K,NN))-YG(PREV)>0.AND.XCG(NBVE(K,NN))-XG(PREV)>0)ANGLE_P=ANGLE_P/PI*180.0      !in the first quadrant
               IF(YCG(NBVE(K,NN))-YG(PREV)<0.AND.XCG(NBVE(K,NN))-XG(PREV)>0)ANGLE_P=360+ANGLE_P/PI*180.0  !in the forth quadrant
               IF(YCG(NBVE(K,NN))-YG(PREV)>0.AND.XCG(NBVE(K,NN))-XG(PREV)<0)ANGLE_P=180+ANGLE_P/PI*180.0  !in the second quadrant
               IF(YCG(NBVE(K,NN))-YG(PREV)<0.AND.XCG(NBVE(K,NN))-XG(PREV)<0)ANGLE_P=180+ANGLE_P/PI*180.0  ! in the third quadrant
               IF(ABS(ANGLE_P-ANGLE_K)>180)THEN
                      IF(ANGLE_P-ANGLE_K>0)THEN
                         ANGLE_K=ANGLE_K+360
                      ELSE
                         ANGLE_P=ANGLE_P+360
                      END IF
               END IF
!               PRINT*,ANGLE_K,ANGLE_P
               IF(ANGLE_K>ANGLE_P)THEN
                 DO II=1,3
                   IF(K==NV(NBVE(K,NN),II))THEN
                     PRINT*,'ADD',K,MGL+NODE
                     NV2(NBVE(K,NN),II)=MGL+NODE
                     XG2(MGL+NODE)=XG(NV(NBVE(K,NN),II))
                     YG2(MGL+NODE)=YG(NV(NBVE(K,NN),II))
                     HG2(MGL+NODE)=HG(NV(NBVE(K,NN),II))
                   END IF
                 END DO
               END IF
          
           END DO
!           PAUSE
           write(1,"(2i6,2f6.2,f8.0,f8.4)") K,MGL+NODE,DAM_HEIGHT,DAM_HEIGHT,SPG_R,SPG_C  !HG2(K)+0.3,HG2(K)+0.3
        END IF
        ! dike nodes but non-crossing and non-ending
        IF(DIKE(K)==2)THEN
           PREV=STRDAM(I,J-1)
           NEXT=STRDAM(I,J+1)
!           PRINT*,K,NTVE(K)
           DO NN=1,NTVE(K)
              SUM=0
              DO II=1,3
                 IF(NEXT==NV(NBVE(K,NN),II))SUM=SUM+1
                 IF(   K==NV(NBVE(K,NN),II))SUM=SUM+1
              END DO

              IF(SUM==2)THEN
!                 PRINT*,'NTVE:',NBVE(K,NN)
                 ANGLE_K=ATAN((YG(NEXT)-YG(K))/(XG(NEXT)-XG(K)))   ! DIKE LINE
                 IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))>0)ANGLE_K=90
                 IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))<0)ANGLE_K=270
                 IF((XG(NEXT)-XG(K))>0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=360
                 IF((XG(NEXT)-XG(K))<0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=180
                 IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=ANGLE_K/PI*180.0      !in the first quadrant
                 IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=360+ANGLE_K/PI*180.0  !in the forth quadrant
                 IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  !in the second quadrant
                 IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  ! in the third quadrant
                 ANGLE_P=ATAN((YCG(NBVE(K,NN))-YG(K))/(XCG(NBVE(K,NN))-XG(K)))   ! Face center of the cell which attached to the dike line
                 IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))>0)ANGLE_P=90
                 IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))<0)ANGLE_P=270
                 IF((XCG(NBVE(K,NN))-XG(K))>0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=360
                 IF((XCG(NBVE(K,NN))-XG(K))<0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=180
                 IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=ANGLE_P/PI*180.0      !in the first quadrant
                 IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=360+ANGLE_P/PI*180.0  !in the forth quadrant
                 IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  !in the second quadrant
                 IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  ! in the third quadrant 
                 IF(ABS(ANGLE_P-ANGLE_K)>180)THEN
                      IF(ANGLE_P-ANGLE_K>0)THEN
                         ANGLE_K=ANGLE_K+360
                      ELSE
                         ANGLE_P=ANGLE_P+360
                      END IF
                 END IF          
                 IF(ANGLE_K>ANGLE_P)EXIT
              END IF              
           END DO
!           PRINT*,ANGLE_K,ANGLE_P,NN,NBVE(K,NN)
           NODE=NODE+1
           DO WHILE(.TRUE.)
             IF(NEXT==PREV)EXIT
!             PRINT*,NEXT,PREV,K,NV(NBVE(K,NN),1:3)
             DO II=1,3
                IF(K==NV(NBVE(K,NN),II))THEN
                   PRINT*,'ADD',K,MGL+NODE
                   NV2(NBVE(K,NN),II)=MGL+NODE
                   XG2(MGL+NODE)=XG(NV(NBVE(K,NN),II))
                   YG2(MGL+NODE)=YG(NV(NBVE(K,NN),II))
                   HG2(MGL+NODE)=HG(NV(NBVE(K,NN),II))
                END IF
             END DO
             DO II=1,3
               IF( (K/=NV(NBVE(K,NN),II)) .AND. (NEXT/=NV(NBVE(K,NN),II)) )THEN
                  NEXT=NV(NBVE(K,NN),II)
                  EXIT
               END IF
             END DO
             DO NN=1,NTVE(K)
                SUM=0
                DO II=1,3
                  IF(NEXT==NV(NBVE(K,NN),II))SUM=SUM+1
                  IF(   K==NV(NBVE(K,NN),II))SUM=SUM+1
                END DO
                IF(SUM==2)THEN
                   ANGLE_K=ATAN((YG(NEXT)-YG(K))/(XG(NEXT)-XG(K)))   ! DIKE LINE
                   IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))>0)ANGLE_K=90
                   IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))<0)ANGLE_K=270
                   IF((XG(NEXT)-XG(K))>0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=360
                   IF((XG(NEXT)-XG(K))<0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=180
                   IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=ANGLE_K/PI*180.0      !in the first quadrant
                   IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=360+ANGLE_K/PI*180.0  !in the forth quadrant
                   IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  !in the second quadrant
                   IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  ! in the third quadrant
                   ANGLE_P=ATAN((YCG(NBVE(K,NN))-YG(K))/(XCG(NBVE(K,NN))-XG(K)))   ! Face center of the cell which attached to the dike line
                   IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))>0)ANGLE_P=90
                   IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))<0)ANGLE_P=270
                   IF((XCG(NBVE(K,NN))-XG(K))>0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=360
                   IF((XCG(NBVE(K,NN))-XG(K))<0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=180
                   IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=ANGLE_P/PI*180.0      !in the first quadrant
                   IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=360+ANGLE_P/PI*180.0  !in the forth quadrant
                   IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  !in the second quadrant
                   IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  ! in the third quadrant
                   IF(ABS(ANGLE_P-ANGLE_K)>180)THEN
                      IF(ANGLE_P-ANGLE_K>0)THEN
                         ANGLE_K=ANGLE_K+360
                      ELSE
                         ANGLE_P=ANGLE_P+360
                      END IF
                   END IF
!                   PRINT*,NEXT,NV(NBVE(K,NN),1:3),ANGLE_K,ANGLE_P
!                   PAUSE
                   IF(ANGLE_K>ANGLE_P)EXIT
                END IF 
             END DO
           END DO
           write(1,"(2i6,2f6.2,f8.0,f8.4)") K,MGL+NODE,DAM_HEIGHT,DAM_HEIGHT,SPG_R,SPG_C  !HG2(K)+0.3,HG2(K)+0.3
        END IF

!    write(1,"(2i6,2f6.2)") K,MGL+NODE,HG2(K)+0.3,HG2(K)+0.3
        
    END DO

END DO


write(1,*)'DAM TYPE 2'
! crossing dike nodes
PRINT*,'!crossing dike nodes'
 DO KK=1,NCROSS_GL
    NCROSS=I_NCROSS_GL(KK,3)
    NUM=KK
    K=I_NCROSS_GL(KK,1)
    print*,I_NCROSS_GL(NUM,4:4+NCROSS-1)
    NEXT=I_NCROSS_GL(NUM,1+3)

    IF(ISONB(K)>0)NCROSS=NCROSS+1  ! If the crossing node is on the boundary, the loop count should be added once more.

    DO JJ=1,NCROSS-1
        NODE=NODE+1
        DO NN=1,NTVE(K)
            SUM=0
            DO II=1,3
              IF(NEXT==NV(NBVE(K,NN),II))SUM=SUM+1
              IF(   K==NV(NBVE(K,NN),II))SUM=SUM+1
            END DO
            IF(SUM==2)THEN
              ANGLE_K=ATAN((YG(NEXT)-YG(K))/(XG(NEXT)-XG(K)))   ! DIKE LINE
              IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))>0)ANGLE_K=90
              IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))<0)ANGLE_K=270
              IF((XG(NEXT)-XG(K))>0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=360
              IF((XG(NEXT)-XG(K))<0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=180
              IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=ANGLE_K/PI*180.0      !in the first quadrant
              IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=360+ANGLE_K/PI*180.0  !in the forth quadrant
              IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  !in the second quadrant
              IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  ! in the third quadrant
              ANGLE_P=ATAN((YCG(NBVE(K,NN))-YG(K))/(XCG(NBVE(K,NN))-XG(K)))   ! Face center of the cell which attached to the dike line
              IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))>0)ANGLE_P=90
              IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))<0)ANGLE_P=270
              IF((XCG(NBVE(K,NN))-XG(K))>0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=360
              IF((XCG(NBVE(K,NN))-XG(K))<0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=180
              IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=ANGLE_P/PI*180.0      !in the first quadrant
              IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=360+ANGLE_P/PI*180.0  !in the forth quadrant
              IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  !in the second quadrant
              IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  ! in the third quadrant
              IF(ABS(ANGLE_P-ANGLE_K)>180)THEN
                 IF(ANGLE_P-ANGLE_K>0)THEN
                    ANGLE_K=ANGLE_K+360
                 ELSE
                    ANGLE_P=ANGLE_P+360
                 END IF
              END IF
              IF(ANGLE_K>ANGLE_P)EXIT
            END IF              
        END DO

        DO WHILE(.TRUE.)
          DO II=1,3
            IF(K==NV(NBVE(K,NN),II))THEN
                
                NV2(NBVE(K,NN),II)=MGL+NODE
                XG2(MGL+NODE)=XG(NV(NBVE(K,NN),II))
                YG2(MGL+NODE)=YG(NV(NBVE(K,NN),II))
                HG2(MGL+NODE)=HG(NV(NBVE(K,NN),II))
            END IF
          END DO
          DO II=1,3
             IF( (K/=NV(NBVE(K,NN),II)) .AND. (NEXT/=NV(NBVE(K,NN),II)) )THEN
                NEXT=NV(NBVE(K,NN),II)
                EXIT
             END IF
          END DO
          PRINT*,'ADD',K,MGL+NODE,NBVE(K,NN),NEXT,JJ
          IF(NEXT==I_NCROSS_GL(NUM,5))THEN
             NEXT=I_NCROSS_GL(NUM,5)
             EXIT
          END IF
          IF(NEXT==I_NCROSS_GL(NUM,6))THEN
             NEXT=I_NCROSS_GL(NUM,6)
             EXIT
          END IF
          IF(NEXT==I_NCROSS_GL(NUM,4))THEN
             NEXT=I_NCROSS_GL(NUM,4)
             EXIT
          END IF
          IF(NEXT==I_NCROSS_GL(NUM,7))THEN
             NEXT=I_NCROSS_GL(NUM,7)
             EXIT
          END IF
          DO NN=1,NTVE(K)
            SUM=0
            DO II=1,3
              IF(NEXT==NV(NBVE(K,NN),II))SUM=SUM+1
              IF(   K==NV(NBVE(K,NN),II))SUM=SUM+1
            END DO
            IF(SUM==2)THEN
                ANGLE_K=ATAN((YG(NEXT)-YG(K))/(XG(NEXT)-XG(K)))   ! DIKE LINE
                IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))>0)ANGLE_K=90
                IF((XG(NEXT)-XG(K))==0.AND.(YG(NEXT)-YG(K))<0)ANGLE_K=270
                IF((XG(NEXT)-XG(K))>0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=360
                IF((XG(NEXT)-XG(K))<0 .AND.(YG(NEXT)-YG(K))==0)ANGLE_K=180
                IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=ANGLE_K/PI*180.0      !in the first quadrant
                IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)>0)ANGLE_K=360+ANGLE_K/PI*180.0  !in the forth quadrant
                IF(YG(NEXT)-YG(K)>0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  !in the second quadrant
                IF(YG(NEXT)-YG(K)<0.AND.XG(NEXT)-XG(K)<0)ANGLE_K=180+ANGLE_K/PI*180.0  ! in the third quadrant
                ANGLE_P=ATAN((YCG(NBVE(K,NN))-YG(K))/(XCG(NBVE(K,NN))-XG(K)))   ! Face center of the cell which attached to the dike line
                IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))>0)ANGLE_P=90
                IF((XCG(NBVE(K,NN))-XG(K))==0.AND.(YCG(NBVE(K,NN))-YG(K))<0)ANGLE_P=270
                IF((XCG(NBVE(K,NN))-XG(K))>0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=360
                IF((XCG(NBVE(K,NN))-XG(K))<0 .AND.(YCG(NBVE(K,NN))-YG(K))==0)ANGLE_P=180
                IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=ANGLE_P/PI*180.0      !in the first quadrant
                IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)>0)ANGLE_P=360+ANGLE_P/PI*180.0  !in the forth quadrant
                IF(YCG(NBVE(K,NN))-YG(K)>0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  !in the second quadrant
                IF(YCG(NBVE(K,NN))-YG(K)<0.AND.XCG(NBVE(K,NN))-XG(K)<0)ANGLE_P=180+ANGLE_P/PI*180.0  ! in the third quadrant
                IF(ABS(ANGLE_P-ANGLE_K)>180)THEN
                IF(ANGLE_P-ANGLE_K>0)THEN
                    ANGLE_K=ANGLE_K+360
                ELSE
                    ANGLE_P=ANGLE_P+360
                END IF
                END IF
                IF(ANGLE_K>ANGLE_P)EXIT
            END IF 
         END DO
        END DO
        
    END DO
    if(jj==3)write(1,"(3i6,3f6.2,f8.0,f8.4)") K,MGL+NODE-1,MGL+NODE,DAM_HEIGHT,DAM_HEIGHT,DAM_HEIGHT,SPG_R,SPG_C !HG2(K)+0.3,HG2(K)+0.3,HG2(K)+0.3
    if(jj==4)write(1,"(4i6,4f6.2,f8.0,f8.4)") K,MGL+NODE-2,MGL+NODE-1,MGL+NODE,DAM_HEIGHT,DAM_HEIGHT,DAM_HEIGHT,DAM_HEIGHT,SPG_R,SPG_C !HG2(K)+0.3,HG2(K)+0.3,HG2(K)+0.3,HG2(K)+0.3
 END DO

close(1)

PRINT*,NODE
OPEN(1,FILE=TRIM(CASENAME)//'_dike.2dm')
WRITE(1,*)'MESH2D'
DO I=1,NGL
  WRITE(1,*)'E3T',I,NV2(I,1:3),1
END DO
DO I=1,MGL+NODE
  WRITE(1,"(a2,i7,1x,e15.9,1x,e15.9,1x,e15.9)")'ND',I,XG2(I),YG2(I),HG2(I)
END DO
CLOSE(1)



END SUBROUTINE DIKE_MESH_GENERATE

















END MODULE MOD_SUBS
