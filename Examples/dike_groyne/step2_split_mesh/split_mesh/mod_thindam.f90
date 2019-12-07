MODULE MOD_THINDAM
   USE ALL_VARS
   USE Mod_Spherical
   IMPLICIT NONE
   SAVE 
   INTEGER              :: INTDM,INTDAMS,INTNODE
   INTEGER              :: NTDCELL_GL,NTDCELL
   INTEGER              :: NTDNODE_GL,NTDNODE_I
   INTEGER              :: NCROSS_GL,NCROSS_I
   INTEGER              :: NEND_GL,NEND_I
   INTEGER, DIMENSION(200,3000):: STRDAM  
   INTEGER,ALLOCATABLE,DIMENSION(:):: DAMNODE
   INTEGER,ALLOCATABLE,DIMENSION(:):: dike
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

  !Declare the mirror matrix upon dike node(non-crossing node)
  !scalar variables including temp,salinity,elvation etc.
   REAL,ALLOCATABLE,DIMENSION(:,:,:)::DK_TEMP,DK_SAL,DK_EL
  !Declare the mirror matrix upon dike node(crossing node)
  !scalar variables including temp,salinity,elvation etc.
   REAL,ALLOCATABLE,DIMENSION(:,:,:)::DKC_TEMP,DKC_SAL,DKC_EL  

  !Declare the tge arrays,see the info in the heading of tge.F
  !non-crossing dike node
   INTEGER, ALLOCATABLE :: DK_NBVE(:,:,:)
   INTEGER, ALLOCATABLE :: DK_NBVT(:,:,:)
   INTEGER, ALLOCATABLE :: DK_NBSN(:,:,:)
   INTEGER, ALLOCATABLE :: DK_NTVE(:,:)      
   INTEGER, ALLOCATABLE :: DK_NTSN(:,:)
  !crossing dike node
   INTEGER, ALLOCATABLE :: DKC_NBVE(:,:,:)
   INTEGER, ALLOCATABLE :: DKC_NBVT(:,:,:)
   INTEGER, ALLOCATABLE :: DKC_NBSN(:,:,:)
   INTEGER, ALLOCATABLE :: DKC_NTVE(:,:)      
   INTEGER, ALLOCATABLE :: DKC_NTSN(:,:)
   CONTAINS

!==========================================================================
   SUBROUTINE READ_NODESTRING
!--------------------------------------------------------------------------|
!  READ THE INFORMATION ABOUT THIN DAM                                     |
!--------------------------------------------------------------------------|
   INTEGER :: I,J,NUM
   CHARACTER(LEN=2)::NS
   INTEGER,DIMENSION(10000,10)::N_STRING
   INTEGER, DIMENSION(200,3000)    :: TEMP
   INTEGER,DIMENSION(200) :: N_NODE 
   
   integer :: istart,iend,count,k,icount(2)
   logical :: is_start,is_end
   
   N_STRING = 0
   N_NODE=0
   NUM=1
   OPEN(1,FILE=TRIM(CASENAME)//'_ns.dat',status='old')
   DO I=1,10000
     READ(1,*,END=999)NS,N_STRING(I,1:10)
     IF(MINVAL(N_STRING(I,:))>=0)THEN
        N_NODE(NUM)=N_NODE(NUM)+10
        TEMP(NUM,N_NODE(NUM)-10+1:N_NODE(NUM))=N_STRING(I,1:10)
     ELSE
        DO J=1,10
           N_NODE(NUM)=N_NODE(NUM)+1
           TEMP(NUM,N_NODE(NUM))=N_STRING(I,J)
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
   open(2,file=trim(casename)//'_xy.dat',status='replace')
   WRITE(1,*)'! THIN DAM'
   WRITE(1,*)NUM
   DO I=1,NUM
     DO J=1,N_NODE(I)
        write(2,*)xg(int(abs(temp(i,j)))),yg(int(abs(temp(i,j))))
        WRITE(1,*)'TD',TEMP(I,J)
     END DO
     write(2,*)
   END DO
   CLOSE(1)
   
   open(1,file=TRIM(CASENAME)//'_dam_cell.dat',status='replace')
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
              if(nv(n,k)==istart)then
                 is_start = .true.
              end if
              if(nv(n,k)==iend)then
                 is_end = .true.
              end if
           end do
           if(is_start.and.is_end)then
              count = count +1
              icount(count)=n
           end if
        end do
        write(1,"(2i6)")icount(1:2)
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
             IF(int(abs(TEMP1(I,J)))==TEMP3(K,1))THEN
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
           WRITE(IPT,*)' !CROSS NODE OF DAM: ',I_NCROSS_GL(N,1),I_NCROSS_GL(N,2)
!           IF(MSR)WRITE(IPT,*)' !CROSS NODE OF DAM: ',I_NCROSS_GL(N,1),I_NCROSS_GL(N,2)
           N=N+1           
        END IF
      END DO

!--------FIND NEBORING DIKE NODE AROUND CROSSING NODES--------!
      DO I=1,INTDAMS 
        DO J=1,DAMNODE(I)
           DO N=1,NCROSS_GL
              IF(J>=1.AND.J<=DAMNODE(I))THEN
                IF(J==1)THEN
                  IF(TMP1(I,J+1)==I_NCROSS_GL(N,1))THEN
                    I_NCROSS_GL(N,3)=I_NCROSS_GL(N,3)+1
                    I_NCROSS_GL(N,3+I_NCROSS_GL(N,3))=TMP1(I,J)
                  END IF
                ELSEIF(J==DAMNODE(I))THEN
                  IF(TMP1(I,J-1)==I_NCROSS_GL(N,1))THEN
                    I_NCROSS_GL(N,3)=I_NCROSS_GL(N,3)+1
                    I_NCROSS_GL(N,3+I_NCROSS_GL(N,3))=TMP1(I,J)
                  END IF
                ELSE
                  IF(TMP1(I,J+1)==I_NCROSS_GL(N,1))THEN
                    I_NCROSS_GL(N,3)=I_NCROSS_GL(N,3)+1
                    I_NCROSS_GL(N,3+I_NCROSS_GL(N,3))=TMP1(I,J)
                  END IF
                  IF(TMP1(I,J-1)==I_NCROSS_GL(N,1))THEN
                    I_NCROSS_GL(N,3)=I_NCROSS_GL(N,3)+1
                    I_NCROSS_GL(N,3+I_NCROSS_GL(N,3))=TMP1(I,J)
                  END IF
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
   
   END SUBROUTINE READ_THINDAM




!==========================================================================
   SUBROUTINE INITIALIZE_THINDAM
!--------------------------------------------------------------------------|
!  INITIALIZE THE INFORMATION ABOUT THIN DAM MODULE                        |
!--------------------------------------------------------------------------|
   LOGICAL :: ISEND,ISCROSS,LTEND,LTCROSS
   INTEGER :: I,J,K,TMP

! divide the dam node into two mirror nodes of different solid
! boundaries except ending nodes and crossing nodes. so do the edges.
   
   NTDNODE_I = 0 
   NCROSS_I  = 0
   ALLOCATE(DIKE(MGL)) !  =1 MEANS THIS IS A
                        !  =3 refers crossing node 
                        !  =2 refers dike node
                        !  =1 refers ending dike node
                        !  =0 MEAN NOT.
   DIKE = 0
   DO I=1,MGL
      DO J=1,NTDNODE_GL
         IF(I==I_TDNODE_GL(J))DIKE(I)=2
      END DO
      DO J=1,NCROSS_GL
         IF(I==I_NCROSS_GL(J,1))DIKE(I)=3
      END DO
      DO J=1,NEND_GL
         IF(I==I_NEND_GL(J))DIKE(I)=1
      END DO
   END DO
   PRINT*,''
   PRINT*,''
!      NTDNODE_I=NTDNODE_GL
!      ALLOCATE(I_TDNODE_N(NTDNODE_I))
!      I_TDNODE_N=I_TDNODE_GL
!      NCROSS_I = NCROSS_GL
!      ALLOCATE(I_NCROSS_N(NCROSS_I,15))
!      I_NCROSS_N = I_NCROSS_GL
!      DO I=1,MGL
!        DO J=1,NTDNODE_I
!           IF(I==I_TDNODE_N(J))DIKE(I)=1
!        END DO
!      END DO

!----------MAPPING TO LOCAL DOMAIN----------------------!

   END SUBROUTINE INITIALIZE_THINDAM


END MODULE MOD_THINDAM
