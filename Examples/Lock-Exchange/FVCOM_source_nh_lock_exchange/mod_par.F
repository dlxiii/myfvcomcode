! $Id: mod_par.F,v 1.5.2.13 2008/05/16 02:06:23 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.5.2.13 $
MODULE MOD_PAR  
   USE MOD_TYPES
   USE MOD_UTILS
   USE LIMS, ONLY : NGL, MGL
   IMPLICIT NONE
   SAVE

!
!--Global Information
!
   INTEGER, ALLOCATABLE :: EL_PID(:)         !!PROCESSOR OWNER OF GLOBAL ELEMENT
   INTEGER, ALLOCATABLE,TARGET :: ELID(:)    !!LOCAL VALUE OF GLOBAL ELEMENT
   INTEGER, ALLOCATABLE,TARGET :: NLID(:)    !!LOCAL VALUE OF GLOBAL NODE 
   INTEGER, POINTER :: ELID_X(:)             !!LOCAL VALUE OF GLOBAL ELEMENT INCLUDING HALOS
   INTEGER, POINTER :: NLID_X(:)             !!LOCAL VALUE OF GLOBAL NODE INCLUDING HALOS 
!
!--Global to Local Boundary Condition Transforms
!
   INTEGER, ALLOCATABLE :: OBN_GL2LOC(:)     !!GIVES GLOBAL IDENTITY OF LOCAL BC ARRAY
!
!--Internal Information (Local)
!
   INTEGER, ALLOCATABLE,TARGET :: EGID(:)    !!GLOBAL ID OF LOCAL ELEMENT
   INTEGER, ALLOCATABLE,TARGET :: NGID(:)    !!GLOBAL ID OF LOCAL NODE 
   INTEGER, POINTER :: EGID_X(:)             !!GLOBAL ID OF LOCAL ELEMENT
   INTEGER, POINTER :: NGID_X(:)             !!GLOBAL ID OF LOCAL NODE 
!
!--Boundary Information: Halo Elements          
!
   INTEGER              :: NHE               !!NUMBER OF HALO ELEMENTS
   INTEGER, ALLOCATABLE :: HE_LST(:)         !!GLOBAL IDENTITIES OF HALO ELEMENTS
   INTEGER, ALLOCATABLE :: HE_OWN(:)         !!OWNER OF HALO ELEMENTS

!
!--Boundary Information: Internal Boundary Nodes
!
   INTEGER :: NBN                            !!NUMBER OF BOUNDARY NODES
   INTEGER :: MX_MLT                         !!MAX MULTIPLICITY OF BOUNDARY NODES
   INTEGER, ALLOCATABLE :: BN_LST(:)         !!GLOBAL IDENTITY OF BOUNDARY NODES
   INTEGER, ALLOCATABLE :: BN_LOC(:)         !!LOCAL IDENTITY OF BOUNDARY NODES
   INTEGER, ALLOCATABLE :: BN_MLT(:)         !!MULTIPLICITY OF BOUNDARY NODES
   INTEGER, ALLOCATABLE :: BN_NEY(:,:)       !!NODE OWNER LIST
   INTEGER, ALLOCATABLE :: NDE_ID(:)         !! = 0 IF INTERNAL, 1 IF ON INTERNAL BOUNDARY

!
!--Boundary Information: Halo Nodes
!
   INTEGER :: NHN  !!NUMBER OF HALO NODES
   INTEGER, ALLOCATABLE :: HN_LST(:)         !!LIST OF HALO NODES 
   INTEGER, ALLOCATABLE :: HN_OWN(:)         !!PRIMARY OWNER OF HALO NODES

# if defined (MULTIPROCESSOR)

!
!--Communication Objects    [SIZE: NPROCS] 
!
   TYPE(COMM), ALLOCATABLE, DIMENSION(:) :: EC,NC,BNC 

!
!--MPI TYPE Objects    [SIZE: NComponents] 
!
   TYPE:: TYPE_DEF
      INTEGER :: NCOMPS
      INTEGER :: OFFSET(100)
      INTEGER :: BLOCKTYPE(100)
      INTEGER :: BLOCKCOUNT(100)
   END TYPE TYPE_DEF

!
!--Maps for Global Array Reconstruction   [SIZE: NPROCS]
!
   TYPE(MAP), ALLOCATABLE, DIMENSION(:) :: EMAP,NMAP, LSFMAP
   TYPE(MAP), ALLOCATABLE, DIMENSION(:) :: EXMAP,NXMAP,BCMAP

   TYPE MAPLINK
      TYPE(MAP), POINTER, DIMENSION(:) :: MAP
      INTEGER :: GLOBAL_SIZE
      TYPE(MAPLINK), POINTER :: NEXT
   END TYPE MAPLINK

   TYPE(MAPLINK) :: HALO_MAPS

   TYPE(MAPLINK) :: INTERNAL_MAPS

!
!--Statistics Calculation   [SIZE: NPROCS]
!

   INTEGER, ALLOCATABLE :: PNE(:)        !!NUMBER OF ELEMENTS IN EACH PROC
   INTEGER, ALLOCATABLE :: PNN(:)        !!NUMBER OF NODES IN EACH PROC
   INTEGER, ALLOCATABLE :: PNHE(:)       !!NUMBER OF HALO ELEMENTS IN EACH PROC
   INTEGER, ALLOCATABLE :: PNBN(:)       !!NUMBER OF INTERNAL BOUNDARY NODES IN EACH PROC
   INTEGER, ALLOCATABLE :: PMBM(:)       !!MAX MULTIPLICITY OF INTERNAL BOUNDARY NODES
   INTEGER, ALLOCATABLE :: PNHN(:)       !!NUMBER OF HALO NODES IN EACH PROC

   ! DEAL FOR POINTERS
   INTERFACE PDEAL
      MODULE PROCEDURE VEC_INT_PDEAL
      MODULE PROCEDURE ARR_INT_PDEAL
      MODULE PROCEDURE CUB_INT_PDEAL

      MODULE PROCEDURE VEC_FLT_PDEAL
      MODULE PROCEDURE ARR_FLT_PDEAL
      MODULE PROCEDURE CUB_FLT_PDEAL

      MODULE PROCEDURE VEC_DBL_PDEAL
      MODULE PROCEDURE ARR_DBL_PDEAL
      MODULE PROCEDURE CUB_DBL_PDEAL
   END INTERFACE
   ! DEAL FOR ALLOCATABLES
   INTERFACE ADEAL
      MODULE PROCEDURE VEC_INT_ADEAL
      MODULE PROCEDURE ARR_INT_ADEAL
      MODULE PROCEDURE CUB_INT_ADEAL

      MODULE PROCEDURE VEC_FLT_ADEAL
      MODULE PROCEDURE ARR_FLT_ADEAL
      MODULE PROCEDURE CUB_FLT_ADEAL

      MODULE PROCEDURE VEC_DBL_ADEAL
      MODULE PROCEDURE ARR_DBL_ADEAL
      MODULE PROCEDURE CUB_DBL_ADEAL
   END INTERFACE

   ! COLLECT FOR POINTER
   INTERFACE PCOLLECT
      MODULE PROCEDURE VEC_INT_PCOLLECT
      MODULE PROCEDURE ARR_INT_PCOLLECT
      MODULE PROCEDURE CUB_INT_PCOLLECT

      MODULE PROCEDURE VEC_FLT_PCOLLECT
      MODULE PROCEDURE ARR_FLT_PCOLLECT
      MODULE PROCEDURE CUB_FLT_PCOLLECT

      MODULE PROCEDURE VEC_DBL_PCOLLECT
      MODULE PROCEDURE ARR_DBL_PCOLLECT
      MODULE PROCEDURE CUB_DBL_PCOLLECT
   END INTERFACE

   ! COLLECT FOR ALLOCATABLES
   INTERFACE ACOLLECT
      MODULE PROCEDURE VEC_INT_ACOLLECT
      MODULE PROCEDURE ARR_INT_ACOLLECT
      MODULE PROCEDURE CUB_INT_ACOLLECT

      MODULE PROCEDURE VEC_FLT_ACOLLECT
      MODULE PROCEDURE ARR_FLT_ACOLLECT
      MODULE PROCEDURE CUB_FLT_ACOLLECT

      MODULE PROCEDURE VEC_DBL_ACOLLECT
      MODULE PROCEDURE ARR_DBL_ACOLLECT
      MODULE PROCEDURE CUB_DBL_ACOLLECT
   END INTERFACE

   INTERFACE PEXCHANGE
      MODULE PROCEDURE VEC_FLT_PEXCHANGE
      MODULE PROCEDURE VEC_INT_PEXCHANGE
      MODULE PROCEDURE VEC_DBL_PEXCHANGE

      MODULE PROCEDURE ARR_FLT_PEXCHANGE
      MODULE PROCEDURE ARR_INT_PEXCHANGE
      MODULE PROCEDURE ARR_DBL_PEXCHANGE
   END INTERFACE

   INTERFACE AEXCHANGE
      MODULE PROCEDURE VEC_FLT_AEXCHANGE
      MODULE PROCEDURE VEC_INT_AEXCHANGE
      MODULE PROCEDURE VEC_DBL_AEXCHANGE

      MODULE PROCEDURE ARR_FLT_AEXCHANGE
      MODULE PROCEDURE ARR_INT_AEXCHANGE
      MODULE PROCEDURE ARR_DBL_AEXCHANGE
   END INTERFACE

   INTERFACE PPRINT
      MODULE PROCEDURE PPRINT_ARR
      MODULE PROCEDURE PPRINT_VEC
   END INTERFACE

   INTERFACE APRINT
      MODULE PROCEDURE APRINT_ARR
      MODULE PROCEDURE APRINT_VEC
   END INTERFACE
!===================================================================================|
   CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
!===================================================================================|

   SUBROUTINE INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID) 
!===================================================================================|
!  INITIALIZE MPI ENVIRONMENT                                                       |
!===================================================================================|
     INTEGER, INTENT(OUT) :: MYID,NPROCS,MSRID
     LOGICAL, INTENT(OUT) :: SERIAL,PAR,MSR
     INTEGER IERR
     
     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "STARTING INIT_MPI_ENV"
     
     IERR=0
     
     CALL MPI_INIT(IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_INIT"
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_COMM_RANK"
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERR)
     IF(IERR/=0) WRITE(*,*) "BAD MPI_COMM_SIZE"
     
     MYID = MYID + 1
     MSRID = 1
     IF(NPROCS > 1) SERIAL=.FALSE.
     IF(NPROCS > 1) PAR   =.TRUE.
     IF(MYID /=  1) MSR   =.FALSE.
     
     
     ! INITIALIZE THE LIST OF MAPS
     nullify(halo_maps%next)
     nullify(halo_maps%map)
     
     nullify(internal_maps%next)
     nullify(internal_maps%map)
     
     
     if(DBG_SET(dbg_sbr)) &
          & write(IPT,*) "END INIT_MPI_ENV"
     
     RETURN
   END SUBROUTINE INIT_MPI_ENV
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! THE FOLLOWING THREE SUBROUTINES ARE USED TO DEFINE MPI TYPES WHICH
! CAN TRANSMIT FORTRAN TYPE DATA BETWEEN PROCESSORS. USE
! 'ADD_TO_MPI_TYPE' TO ADD COMPONENTS TO AN MPI TYPE. THE ORDER
! MATTERS AND MUST MATCH THE ORDER IN WHICH THE DATA IS DECLARED
! IN THE FORTRAN TYPE. SEE EXAMPLE IN particle.F
!
!
! MPI AND FORTRAN DATA TYPE NAMES
!
!
!$MPI datatype          $FORTRAN datatype
!MPI_INTEGER            INTEGER
!MPI_REAL               REAL
!MPI_REAL8              REAL*8
!MPI_DOUBLE_PRECISION   DOUBLE PRECISION
!MPI_COMPLEX            COMPLEX
!MPI_LOGICAL            LOGICAL
!MPI_CHARACTER          CHARACTER
!MPI_BYTE               -	 
!MPI_PACKED             -
!
! mod_prec.F DECLARES THE FOLLOWING EQUIVELENCE
!
! MPI_F                REAL(SP) (FOR DOUBLE AND SINGLE PRECISION MODELS)
! MPI_DP               REAL(DP) 
!
!==============================================================================|
   FUNCTION INIT_TYPE_DEF() RESULT(DEF)
     IMPLICIT NONE
     TYPE(TYPE_DEF) :: DEF
     DEF%NCOMPS     = 0
     DEF%BLOCKTYPE  = 0
     DEF%BLOCKCOUNT = 0
     DEF%OFFSET     = 0
   END FUNCTION INIT_TYPE_DEF
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   SUBROUTINE ADD_TO_MPI_TYPE(DEF,MYTYPE,COUNT)
     IMPLICIT NONE

     TYPE(TYPE_DEF) :: DEF
     INTEGER, INTENT(IN) :: MYTYPE
     INTEGER, INTENT(IN) :: COUNT

     INTEGER :: N, EXTENT,IERR


     ! INCRIMENT THE NUMBER OF COMPONENTS
     DEF%NCOMPS = DEF%NCOMPS+1

     IF (DEF%NCOMPS .GT. 99) CALL FATAL_ERROR&
          &("ADD_TO_MPI_TYPE: More than 100 components!",&
          & "You must edit the TYPE_DEF and recompile!" )

     N = DEF%NCOMPS

     ! SET THE NEW BLOCKTYPE
     DEF%BLOCKTYPE(N) = MYTYPE 

     ! SET THE NEW BLOCKCOUNT
     DEF%BLOCKCOUNT(N) = COUNT 

     call MPI_TYPE_EXTENT(MYTYPE, EXTENT, IERR)
     IF(IERR /= 0) CALL FATAL_ERROR&
          & ("ADD_TO_MPI_TYPE: COULD NOT GET EXTENT FOR THE TYPE?")
        
     ! SET THE NEXT OFFSET (LAST ONE IS NOT USED)
     DEF%OFFSET(N+1) = DEF%OFFSET(N)+ COUNT * EXTENT


   END SUBROUTINE ADD_TO_MPI_TYPE
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   FUNCTION CREATE_MPI_TYPE(DEF) RESULT(MYMPITYPE)
     IMPLICIT NONE

     TYPE(TYPE_DEF),INTENT(IN) :: DEF
     INTEGER :: MYMPITYPE
     
     INTEGER :: IERR

     CALL MPI_TYPE_STRUCT(DEF%NCOMPS, DEF%blockcount, DEF%offset, DEF%blocktype, mympitype, ierr) 
     IF(IERR /= 0) CALL FATAL_ERROR&
          & ("CREATE_MPI_TYPE: COULD NOT CREATE MPI_TYPE_STRUCT?")


      CALL MPI_TYPE_COMMIT(MYMPITYPE, ierr)
      IF(IERR /= 0) CALL FATAL_ERROR&
          & ("CREATE_MPI_TYPE: COULD NOT COMMIT MPI TYPE")

   END FUNCTION CREATE_MPI_TYPE
!==============================================================================|
   
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
  SUBROUTINE DOMDEC(NVG,NPROCS,EL_PID,MSR)              
!==============================================================================|
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVG(0:NGL,4)
   INTEGER, INTENT(IN)  :: NPROCS
   INTEGER, INTENT(OUT) :: EL_PID(NGL)
   LOGICAL, INTENT(IN)  :: MSR
   INTEGER, ALLOCATABLE :: NVT(:)
   INTEGER :: I,J,COUNT,NTEMP,IERR,ii
!==============================================================================|
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING DOMDEC"

!
!----------------CONVERT NVG to 1D-----------------------!
! 

   IF(MSR) THEN
      ALLOCATE(NVT(3*NGL))
      DO I=1,NGL
         DO J = 1,3
            COUNT = (I-1)*3 + J
            ! FVCOM SWITCHES NODE ORDER
            NVT(COUNT) = NVG(I,5-J)
         END DO
      END DO
      
!
!-------------DECOMPOSE ELEMENTS USING METIS GRAPH PARTITIONING ---------------!
!
      
      
      CALL PARTITION(NPROCS,NGL,MAXVAL(NVT),loc(NVT),loc(EL_PID))
      EL_PID = EL_PID + 1
      DEALLOCATE(NVT)
   END IF

   !first make sure it is an 2D mesh that can be decomposed exactly
   if(mod(ngl,8) /= 0)then
     write(*,*)'mesh size: ',ngl
     write(*,*)'this mesh size is not a factor of 8, not a 2D mesh'
     call pstop
   endif
   if(mod(ngl/8,nprocs) /= 0)then
     write(*,*)'mesh size: ',ngl
     write(*,*)'wont decompose in 2D fashion exactly into #procs',nprocs
     call pstop
   endif

   !decompose
   ntemp = ngl/nprocs
   do i=1,nprocs
     el_pid( (i-1)*ntemp+1 : i*ntemp ) = i
   end do

!---------------------BROADCAST RESULT TO ALL PROCESSORS-----------------------!

   
   CALL MPI_BCAST(EL_PID,NGL,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END DOMDEC"

   END SUBROUTINE DOMDEC
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_FLT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:)

   REAL(SPA), POINTER :: AP(:)
   REAL(SPA), POINTER :: BP(:)
   REAL(SPA), POINTER :: CP(:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE VEC_FLT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_FLT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), POINTER           :: A(:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: B(:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: C(:)
   INTEGER :: NT  
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(SPA), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_FLT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1)               = ",Size(A,1)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
                & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("VEC_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS + 1
       LBUF  = NVARS * CM(I)%NRCV
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_REAL,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 

       N1 = N2+1  ; N2 = N1 + LBUF -1
       SBUF(N1:N2) = A(CM(I)%SNDP(:))
       IF(BYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = B(CM(I)%SNDP(:))
       END IF
       IF(CYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = C(CM(I)%SNDP(:))
       END IF

       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_REAL,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP

     N1 = N2+1 ; N2 = N1 + LBUF -1
     A(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     IF(BYES)THEN
        N1 = N2+1; N2 = N1 + LBUF -1
        B(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
     IF(CYES)THEN
        N1 = N2+1 ; N2 = N1 + LBUF -1
        C(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
  END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_FLT_EXCHANGE"

   RETURN
 END SUBROUTINE VEC_FLT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_INT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:)

   INTEGER, POINTER :: AP(:)
   INTEGER, POINTER :: BP(:)
   INTEGER, POINTER :: CP(:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL VEC_INT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE VEC_INT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_INT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), POINTER           :: A(:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: B(:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: C(:)
   INTEGER :: NT  
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   INTEGER, ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_INT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1)               = ",Size(A,1)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
                & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("VEC_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS + 1
       LBUF  = NVARS * CM(I)%NRCV
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_INTEGER,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 

       N1 = N2+1  ; N2 = N1 + LBUF -1
       SBUF(N1:N2) = A(CM(I)%SNDP(:))
       IF(BYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = B(CM(I)%SNDP(:))
       END IF
       IF(CYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = C(CM(I)%SNDP(:))
       END IF

       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_INTEGER,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP

     N1 = N2+1 ; N2 = N1 + LBUF -1
     A(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     IF(BYES)THEN
        N1 = N2+1; N2 = N1 + LBUF -1
        B(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
     IF(CYES)THEN
        N1 = N2+1 ; N2 = N1 + LBUF -1
        C(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
  END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_INT_EXCHANGE"

   RETURN
 END SUBROUTINE VEC_INT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_DBL_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:)

   REAL(DP), POINTER :: AP(:)
   REAL(DP), POINTER :: BP(:)
   REAL(DP), POINTER :: CP(:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE VEC_DBL_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE VEC_DBL_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), POINTER           :: A(:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: B(:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: C(:)
   INTEGER :: NT  
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(DP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_DBL_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1)               = ",Size(A,1)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
                & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("VEC_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("VEC_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1)) CALL FATAL_ERROR &
             &("VEC_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS + 1
       LBUF  = NVARS * CM(I)%NRCV
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_DP,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 

       N1 = N2+1  ; N2 = N1 + LBUF -1
       SBUF(N1:N2) = A(CM(I)%SNDP(:))
       IF(BYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = B(CM(I)%SNDP(:))
       END IF
       IF(CYES)THEN
          N1 = N2+1 ; N2 = N1 + LBUF -1
          SBUF(N1:N2) = C(CM(I)%SNDP(:))
       END IF

       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_DP,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP

     N1 = N2+1 ; N2 = N1 + LBUF -1
     A(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     IF(BYES)THEN
        N1 = N2+1; N2 = N1 + LBUF -1
        B(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
     IF(CYES)THEN
        N1 = N2+1 ; N2 = N1 + LBUF -1
        C(CM(LPROC)%RCVP(:)) = RBUF(N1:N2) 
     END IF
  END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_DBL_EXCHANGE"

   RETURN
 END SUBROUTINE VEC_DBL_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_FLT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:,:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:,:)
   REAL(SPA), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:,:)

   REAL(SPA), POINTER :: AP(:,:)
   REAL(SPA), POINTER :: BP(:,:)
   REAL(SPA), POINTER :: CP(:,:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE ARR_FLT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_FLT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SPA), INTENT(INOUT), POINTER           :: A(:,:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: B(:,:)
   REAL(SPA), INTENT(INOUT), POINTER, OPTIONAL :: C(:,:)
   INTEGER :: NT,KT
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(SPA), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_FLT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                  = ",MYID
        write(IPT,*) "NPROCS                = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1),Size(A,2)   = ",Size(A,1),Size(A,2)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
             &write(IPT,*) "SIZE(B,1),Size(B,2)   = ",Size(B,1),Size(B,2)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                &write(IPT,*) "SIZE(C,1),Size(C,2)   = ",Size(C,1),Size(C,2)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("ARR_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     KT = Size(A,2)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1) .OR. KT .NE. Size(B,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1) .OR. KT .NE. Size(C,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_REAL,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_REAL,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_FLT_EXCHANGE"

   RETURN
 END SUBROUTINE ARR_FLT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_INT_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:,:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:,:)
   INTEGER, INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:,:)

   INTEGER, POINTER :: AP(:,:)
   INTEGER, POINTER :: BP(:,:)
   INTEGER, POINTER :: CP(:,:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL ARR_INT_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE ARR_INT_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_INT_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   INTEGER, INTENT(INOUT), POINTER           :: A(:,:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: B(:,:)
   INTEGER, INTENT(INOUT), POINTER, OPTIONAL :: C(:,:)
   INTEGER :: NT,KT
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   INTEGER, ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_INT_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                  = ",MYID
        write(IPT,*) "NPROCS                = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1),Size(A,2)   = ",Size(A,1),Size(A,2)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
             &write(IPT,*) "SIZE(B,1),Size(B,2)   = ",Size(B,1),Size(B,2)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                &write(IPT,*) "SIZE(C,1),Size(C,2)   = ",Size(C,1),Size(C,2)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("ARR_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     KT = Size(A,2)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1) .OR. KT .NE. Size(B,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1) .OR. KT .NE. Size(C,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_INTEGER,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_INTEGER,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_INT_EXCHANGE"

   RETURN
 END SUBROUTINE ARR_INT_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_DBL_AEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET           :: A(:,:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: B(:,:)
   REAL(DP), INTENT(INOUT), ALLOCATABLE, TARGET, OPTIONAL :: C(:,:)

   REAL(DP), POINTER :: AP(:,:)
   REAL(DP), POINTER :: BP(:,:)
   REAL(DP), POINTER :: CP(:,:)
 

   IF(ALLOCATED(A))  AP => A

   IF( PRESENT(B) .AND.PRESENT(C)) THEN 
      IF(ALLOCATED(B))  BP => B
      IF(ALLOCATED(C))  CP => C
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,BP,CP) 
   ELSE IF(PRESENT(B)) THEN 
      IF(ALLOCATED(B))  BP => B
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,B=BP) 
   ELSE IF(PRESENT(C)) THEN 
      IF(ALLOCATED(C))  CP => C
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP,C=CP)
   ELSE
      CALL ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,AP)
   END IF

 END SUBROUTINE ARR_DBL_AEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE ARR_DBL_PEXCHANGE(CM,MYID,NPROCS,A,B,C) 
!===================================================================================|
!    PASS ELEMENT/NODE INFORMATION AMONG PROCESSORS                                 |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(DP), INTENT(INOUT), POINTER           :: A(:,:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: B(:,:)
   REAL(DP), INTENT(INOUT), POINTER, OPTIONAL :: C(:,:)
   INTEGER :: NT,KT
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(DP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DBL_PEXCHANGE"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                  = ",MYID
        write(IPT,*) "NPROCS                = ",NPROCS
        IF(ASSOCIATED(A))&
             &write(IPT,*) "SIZE(A,1),Size(A,2)   = ",Size(A,1),Size(A,2)
        IF(PRESENT(B)) THEN
           IF(ASSOCIATED(B))&
             &write(IPT,*) "SIZE(B,1),Size(B,2)   = ",Size(B,1),Size(B,2)
        END IF

        IF(PRESENT(C)) THEN
           IF(ASSOCIATED(C))&
                &write(IPT,*) "SIZE(C,1),Size(C,2)   = ",Size(C,1),Size(C,2)
        END IF

        DO I = 1,NPROCS
           write(IPT,*)"=============MAP INFO=============="
           write(IPT,*) "ID                      = ",I
           write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
           write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="
        
     END if

     IF(.NOT. ASSOCIATED(A)) CALL FATAL_ERROR&
          &("ARR_FLT_PEXCHANGE: PRIMARY ARGUMENT NOT ASSOCIATED")
     NT = Size(A,1)
     KT = Size(A,2)
     IF(PRESENT(B)) THEN
        IF(.NOT. ASSOCIATED(B)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: SECONDARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(B,1) .OR. KT .NE. Size(B,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION SIZES DO NOT MATCH")
     END IF

     IF(PRESENT(C)) THEN
        IF(.NOT. ASSOCIATED(C)) CALL FATAL_ERROR&
             &("ARR_FLT_PEXCHANGE: TERTIARY ARGUMENT NOT ASSOCIATED")
        IF (NT .NE. Size(C,1) .OR. KT .NE. Size(C,2)) CALL FATAL_ERROR &
             &("ARR_FLT_PEXCHANGE: DIMENSION ARUGEMENTS DO NOT MATCH")
     END IF



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_DP,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_DP,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_DBL_EXCHANGE"

   RETURN
 END SUBROUTINE ARR_DBL_PEXCHANGE
!===================================================================================|

!-----------------------------------------------------------------------------------!

!===================================================================================|
   SUBROUTINE NODE_MATCH(IMATCH,NBN,BN_MLT,BN_LOC,CM,NT,KT,MYID,NPROCS,A,B,C) 
!===================================================================================|
! IMATCH=1:   ENFORCE AGREEMENT OF A,B,C ON BOUNDARY NODES                          |
! IMATCH=0:   ACCUMULATE VALUES OF A,B,C AT BOUNDARY NODES                          |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER, INTENT(IN)             :: IMATCH
   INTEGER, INTENT(IN)             :: NBN
   INTEGER, INTENT(IN)             :: BN_MLT(NBN)
   INTEGER, INTENT(IN)             :: BN_LOC(NBN)
   INTEGER, INTENT(IN)             :: NT,KT,MYID,NPROCS
   TYPE(COMM), INTENT(IN)          :: CM(NPROCS)
   REAL(SP), INTENT(INOUT)           :: A(0:NT,KT)
   REAL(SP), INTENT(INOUT), OPTIONAL :: B(0:NT,KT)
   REAL(SP), INTENT(INOUT), OPTIONAL :: C(0:NT,KT)
!------------------------------------------------------------------------------
   LOGICAL             :: BYES,CYES
   INTEGER               ::IREQR(NPROCS),IREQS(NPROCS)
   REAL(SP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),ISTATR(MPI_STATUS_SIZE,NPROCS),IERR,J,N1,N2,NCNT
   INTEGER   I,IFROM,ITO,ISTAG,IRTAG,TRCV,TSND,NVARS,LBUF,LP,NMSG,INDX,LPROC,NSZE
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING NODE_MATCH"

   if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "KT                      = ",KT
        write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        IF(PRESENT(B)) &
             & write(IPT,*) "SIZE(B,1)           = ",Size(B,1)
        IF(PRESENT(C)) &
             & write(IPT,*) "SIZE(C,1)           = ",Size(C,1)

        DO I = 1,NPROCS
        write(IPT,*)"=============MAP INFO=============="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "CM(I)%NSND              = ",CM(I)%NSND
        write(IPT,*) "CM(I)%NRCV              = ",CM(I)%NRCV
        END DO
        write(IPT,*)"==================================="

           
     END if



   NVARS = 1 ; BYES = .FALSE. ; CYES = .FALSE.
   IF(PRESENT(B)) THEN
     NVARS = NVARS + 1
     BYES  = .TRUE.
   END IF
      
   IF(PRESENT(C)) THEN
     NVARS = NVARS + 1
     CYES = .TRUE. 
   END IF

   ALLOCATE(RBUF(NVARS*SUM(CM(1:NPROCS)%NRCV*KT)))
   ALLOCATE(SBUF(NVARS*SUM(CM(1:NPROCS)%NSND*KT)))

!===================================================================================|
!    POST NON-BLOCKING RECEIVES FROM NEIGHBORS                                      |
!===================================================================================|
   TRCV = 0
   DO I=1,NPROCS

     IF(CM(I)%NRCV > 0)THEN
       IFROM = I-1
       IRTAG = I*1000
       TRCV  = TRCV + 1
       LP    = CM(I)%RCPT*NVARS*KT + 1
       LBUF  = NVARS * CM(I)%NRCV *KT
       CALL MPI_IRECV(RBUF(LP),LBUF,MPI_F,IFROM,IRTAG,MPI_COMM_WORLD,IREQR(TRCV),IERR)
     END IF

   END DO

!===================================================================================|
!    SEND DATA TO NEIGHBORS                                                         |
!===================================================================================|
   TSND = 0
   NCNT = 0
   DO I=1,NPROCS
     LBUF = CM(I)%NSND
     IF(LBUF > 0)THEN
       NSZE = LBUF*KT*NVARS
!       ALLOCATE(SBUF(NSZE))
       N2 = NCNT 
       DO J=1,KT
         N1 = N2+1  ; N2 = N1 + LBUF -1
         SBUF(N1:N2) = A(CM(I)%SNDP(:),J)
         IF(BYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = B(CM(I)%SNDP(:),J)
         END IF
         IF(CYES)THEN
           N1 = N2+1 ; N2 = N1 + LBUF -1
           SBUF(N1:N2) = C(CM(I)%SNDP(:),J)
         END IF
       END DO
       TSND  = TSND + 1
       ITO   = I-1
       ISTAG = MYID*1000
       CALL MPI_ISEND(SBUF(NCNT+1),NSZE,MPI_F,ITO,ISTAG,MPI_COMM_WORLD,IREQS(TSND),IERR)
       NCNT = NCNT + LBUF*KT*NVARS 
!       DEALLOCATE(SBUF)
     END IF
   END DO


!===================================================================================|
!    LOOP OVER PROCS UNTIL A MESSAGE IS RECEIVED AND UNPACK                         |
!===================================================================================|
   DO NMSG = 1,TRCV 
     CALL MPI_WAITANY(TRCV,IREQR,INDX,STAT,IERR)
     LPROC = STAT(MPI_SOURCE) +1 
     LP    = CM(LPROC)%RCPT*NVARS*KT 
     LBUF  = CM(LPROC)%NRCV
     N2 = LP
     DO J=1,KT
       N1 = N2+1 ; N2 = N1 + LBUF -1
       A(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) + A(CM(LPROC)%RCVP(:),J) 
       IF(BYES)THEN
         N1 = N2+1; N2 = N1 + LBUF -1
         B(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) + B(CM(LPROC)%RCVP(:),J)
       END IF
       IF(CYES)THEN
         N1 = N2+1 ; N2 = N1 + LBUF -1
         C(CM(LPROC)%RCVP(:),J) = RBUF(N1:N2) + C(CM(LPROC)%RCVP(:),J) 
       END IF
     END DO
   END DO

!===================================================================================|
!    WAIT FOR COMPLETION OF NON-BLOCKING SENDS                                      |
!===================================================================================|

   CALL MPI_WAITALL(TSND,IREQS,ISTATR,IERR)
   DEALLOCATE(RBUF,SBUF)

!===================================================================================|
!  USE MULTIPLICITY OF NODES TO COMPUTE TRUE AVERAGE VALUE                          |
!===================================================================================|
   IF(IMATCH /=1)RETURN

   DO J=1,KT
     DO I=1,NBN
       A( BN_LOC(I),J) = A( BN_LOC(I),J)/FLOAT(BN_MLT(I))
       IF(BYES)B( BN_LOC(I),J) = B( BN_LOC(I),J)/FLOAT(BN_MLT(I))
       IF(CYES)C( BN_LOC(I),J) = C( BN_LOC(I),J)/FLOAT(BN_MLT(I))
     END DO
   END DO

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END NODE_MATCH"

   RETURN
   END SUBROUTINE NODE_MATCH
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of INTEGERS FROM A GLOBAL VECTOR                                 |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                        |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:) :: AP
   INTEGER,   POINTER,DIMENSION(:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: AG
   IF(ALLOCATED(A))  AP => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE VEC_INT_ADEAL

!===================================================================================|
   SUBROUTINE VEC_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of INTEGERS FROM A GLOBAL VECTOR                                 |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                        |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:) :: A
   INTEGER,   POINTER,DIMENSION(:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:),RBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE
   INTEGER, PARAMETER :: TAG = 30111 ! UNIQUE TAG FOR VEC_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_INT_PDEAL"


   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")

      

   if(DBG_SET(dbg_sbrio)) then
      
   
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID
        if(associated(A)) then
           write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1)           = ",Size(AG,1)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"=============MAP INFO=============="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

           
     END if

   if (MYID .EQ. SENDID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED BUT NEED NOT BE GLOBAL SIZE FOR ALL PROCS")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) CYCLE

         ALLOCATE(SBUF(NSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_INT_PDEAL")
         
         DO I=1,NSZE
            SBUF(I) = AG(GM(IP)%LOC_2_GL(I))
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE) = SBUF
            ! THIS MAY BE A PROBLEM, MAKE ASSIGNMENT DIMENSION MORE SPECIFIC!
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("Send Error in VEC_INT_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PDEAL")

      END DO
   else ! IF I AM ONE OF THE RECEIVERS

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0 ) RETURN

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(RBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_INT_PDEAL")
      
      RBUF=0
      CALL MPI_RECV(RBUF,NSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD , STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN VEC_INT_PDEAL")
        
      A(1:NSZE)=RBUF
      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PDEAL")
      
   end if

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_INT_PDEAL"
 END SUBROUTINE VEC_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF INTEGERS FROM LOCAL VECTORS                                |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:)   :: AG
   INTEGER,   POINTER,DIMENSION(:) :: AP
   INTEGER,   POINTER,DIMENSION(:) :: AGP
   IF(ALLOCATED(A))   AP  => A
   IF(ALLOCATED(AG))  AGP => AG
   CALL VEC_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE VEC_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE VEC_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF INTEGERS FROM LOCAL VECTORS                                |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE
   INTEGER, PARAMETER :: TAG =  30112 ! UNIQUE TAG FOR VEC_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL COLLECT AS THE&
        & RECEIVER")


     if(DBG_SET(dbg_sbrio)) then
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "RECVID                  = ",RECVID
        if(associated(A)) then
           write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1)           = ",Size(AG,1)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="
     END if

   if (MYID .EQ. RECVID) then

        if(.not. associated(AG))&
             & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = A(I)
            END DO
         else
            ALLOCATE(RBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in VEC_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = RBUF(I)
            END DO
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PCOLLECT")
         end if
       End DO
   else

      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(SBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_INT_PCOLLECT")
      SBUF = A(1:NSZE)

      DEST = RECVID - 1
      CALL MPI_SEND(SBUF,NSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN VEC_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_INT_PCOLLECT")

   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_INT_PCOLLECT"
 END SUBROUTINE VEC_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:) :: AGP
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:) :: A
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE ARR_INT_ADEAL
!===================================================================================|
   SUBROUTINE ARR_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:) :: A
   INTEGER,   POINTER,DIMENSION(:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:,:),RBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE
   INTEGER, PARAMETER :: TAG = 30121 ! UNIQUE TAG FOR ARR_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_INT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")


   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),  SIZE(A,2)   = ",Size(A,1),Size(A,2)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),  SIZE(AG,2)   = ",Size(AG,1),Size(AG,2)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if


   if (MYID .EQ. SENDID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED BUT NEED NOT BE GLOBAL SIZE FOR ALL PROCS")
      
      PSZE=size(AG,2)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_INT_PDEAL")
         
         DO I=1,NSZE
            SBUF(I,:) = AG(GM(IP)%LOC_2_GL(I),:)
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")
           
            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE,:) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_INT_PDEAL")
         end if
         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PDEAL")
      END DO
   else ! I am a receiver
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0 ) RETURN
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")
      
      PSZE=size(A,2)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
      
      ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_INT_PDEAL")

      RBUF=0
      CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN ARR_INT_PDEAL")
      
      A(1:NSZE,1:PSZE)=RBUF
      
      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PDEAL")
      
     end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_INT_PDEAL"
 END SUBROUTINE ARR_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE ARR_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE ARR_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:,:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:,:),SBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE
   INTEGER, PARAMETER :: TAG =  30122 ! UNIQUE TAG FOR ARR_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & RECEIVER")

     if(DBG_SET(dbg_sbrio)) then
   
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "RECVID                  = ",RECVID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),  SIZE(A,2)   = ",Size(A,1),Size(A,2)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),  SIZE(AG,2)   = ",Size(AG,1),Size(AG,2)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if

   if (MYID .EQ. RECVID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")

     PSZE = size(AG,2)
      
      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:) = A(I,:)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:) = RBUF(I,:)
            END DO
           
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PCOLLECT")
         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")
      PSZE = size(A,2)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_INT_PCOLLECT")

      SBUF= A(1:NSZE,1:PSZE)

      CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN ARR_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_INT_PCOLLECT")
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_INT_PCOLLECT"
 END SUBROUTINE ARR_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_INT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AGP
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:,:) :: A
   INTEGER,   ALLOCATABLE, TARGET,DIMENSION(:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_INT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE CUB_INT_ADEAL
!===================================================================================|
   SUBROUTINE CUB_INT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of INTEGERS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:) :: A
   INTEGER,   POINTER,DIMENSION(:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: SBUF(:,:,:),RBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE
   INTEGER, PARAMETER :: TAG = 30131 ! UNIQUE TAG FOR ARR_INT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_INT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")


   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1), SIZE(A,2), SIZE(A,3)   = ",Size(A,1),Size(A,2),Size(A,3)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),SIZE(AG,2),SIZE(AG,3)   = ",Size(AG,1),Size(AG,2),Size(AG,3)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if


   if (MYID .EQ. SENDID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED BUT NEED NOT BE GLOBAL SIZE FOR ALL PROCS")
      
      PSZE=size(AG,2)
      QSZE=size(AG,3)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_INT_PDEAL")
         
         DO I=1,NSZE
            SBUF(I,:,:) = AG(GM(IP)%LOC_2_GL(I),:,:)
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE,:,:) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_INT_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PDEAL")

      END DO
   else ! I am a receiver
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0 ) RETURN
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")
      
      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
      
      PSZE=size(A,2)
      QSZE=size(A,3)
      
      ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_INT_PDEAL")
      RBUF=0
      
      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN CUB_INT_PDEAL")
      
      A(1:NSZE,1:PSZE,1:QSZE)=RBUF
      
      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_INT_PDEAL"
 END SUBROUTINE CUB_INT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_INT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AP
   INTEGER,   POINTER,DIMENSION(:,:,:) :: AGP
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN) :: A
   INTEGER,   ALLOCATABLE,TARGET,DIMENSION(:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE CUB_INT_ACOLLECT
!===================================================================================|
   SUBROUTINE CUB_INT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF INTEGERS FROM A LOCAL ARRAYS                               |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   INTEGER,   POINTER,DIMENSION(:,:,:),INTENT(IN) :: A
   INTEGER,   POINTER,DIMENSION(:,:,:)   :: AG
!------------------------------------------------------------------------------
   INTEGER, ALLOCATABLE :: RBUF(:,:,:),SBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE
   INTEGER, PARAMETER :: TAG =  30132 ! UNIQUE TAG FOR CUB_INT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_INT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & RECEIVER")

     if(DBG_SET(dbg_sbrio)) then
   
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "RECVID                  = ",RECVID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),SIZE(A,2),SIZE(A,3)   = ",Size(A,1),Size(A,2),Size(A,3)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),SIZE(AG,2),SIZE(AG,3)   = ",Size(AG,1),Size(AG,2),Size(AG,3)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if

   if (MYID .EQ. RECVID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")

     PSZE = size(AG,2)
     QSZE = size(AG,3)
      
      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:,:) = A(I,:,:)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_INT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_INTEGER,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_INT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:,:) = RBUF(I,:,:)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PCOLLECT")

         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")
      PSZE = size(A,2)
      QSZE = size(A,3)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
 

      ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_INT_PCOLLECT")
      SBUF = A(1:NSZE,1:PSZE,1:QSZE)
      
      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_INTEGER,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN ARR_INT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_INT_PCOLLECT")

   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_INT_PCOLLECT"
 END SUBROUTINE CUB_INT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of REALS FROM A GLOBAL VECTOR                                    |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE VEC_FLT_ADEAL
!===================================================================================|
   SUBROUTINE VEC_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of REALS FROM A GLOBAL VECTOR                                    |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:) :: A
   REAL(SPA),   POINTER,DIMENSION(:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:),RBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE
   INTEGER, PARAMETER :: TAG = 30211 ! UNIQUE TAG FOR VEC_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_FLT_PDEAL"


   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")

      

   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID
        if(associated(A)) then
           write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1)           = ",Size(AG,1)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"=============MAP INFO=============="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="
     END if

   if (MYID .EQ. SENDID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED BUT NEED NOT BE GLOBAL SIZE FOR ALL PROCS")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         ALLOCATE(SBUF(NSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PDEAL")
         
         DO I=1,NSZE
            SBUF(I) = AG(GM(IP)%LOC_2_GL(I))
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("Send Error in VEC_FLT_PDEAL")
         end if
        
         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PDEAL")

      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(RBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PDEAL")
      RBUF=0.0_SPA

      CALL MPI_RECV(RBUF,NSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN VEC_FLT_PDEAL")
      
      A(1:NSZE)=RBUF

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_FLT_PDEAL"
 END SUBROUTINE VEC_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF REALS FROM LOCAL VECTORS                                   |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:) :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE VEC_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE VEC_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF REALS FROM LOCAL VECTORS                                   |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE
   INTEGER, PARAMETER :: TAG =  30212 ! UNIQUE TAG FOR VEC_FLT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_FLT_PCOLLECT"
   
  if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL COLLECT AS THE&
        & RECEIVER")


     if(DBG_SET(dbg_sbrio)) then
   
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "RECVID                  = ",RECVID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1)           = ",Size(AG,1)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if

   if (MYID .EQ. RECVID) then
        if(.not. associated(AG))&
             & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,GM(IP)%NSIZE
               AG(GM(IP)%LOC_2_GL(I)) = A(I)
            END DO
         else
            ALLOCATE(RBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in VEC_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = RBUF(I)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PCOLLECT")
            
         end if
       End DO
   else
      
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(SBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_FLT_PCOLLECT")

      SBUF =A(1:NSZE)

      CALL MPI_SEND(SBUF,NSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN VEC_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_FLT_PCOLLECT")
      
   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_FLT_PCOLLECT"
 END SUBROUTINE VEC_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE ARR_FLT_ADEAL
!===================================================================================|
   SUBROUTINE ARR_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:,:),RBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE
   INTEGER, PARAMETER :: TAG = 30221 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_FLT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")


   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),  SIZE(A,2)   = ",Size(A,1),Size(A,2)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),  SIZE(AG,2)   = ",Size(AG,1),Size(AG,2)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if


   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED ON THE SENDING PROC")

      PSZE=size(AG,2)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PDEAL")
         
         DO I=1,NSZE
            SBUF(I,:) = AG(GM(IP)%LOC_2_GL(I),:)
         END DO

         if (IP == MYID) then
            if(.not. associated(A)) CALL FATAL_ERROR&
                 &("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE,:) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_FLT_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PDEAL")

      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")
      PSZE=size(A,2)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
      
      ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PDEAL")
      RBUF=0.0_SPA
     
      CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN ARR_FLT_PDEAL")
      
      A(1:NSZE,:)=RBUF

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_FLT_PDEAL"
 END SUBROUTINE ARR_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE ARR_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE ARR_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:,:),SBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE
   INTEGER, PARAMETER :: TAG =  30222 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_FLT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & RECEIVER")

   
   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID

      if(associated(A)) then
         write(IPT,*) "SIZE(A,1),  SIZE(A,2)   = ",Size(A,1),Size(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "SIZE(AG,1),  SIZE(AG,2)   = ",Size(AG,1),Size(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      PSZE = size(AG,2)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:) = A(I,:)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PCOLLECT")
            RBUF = 0.0_SPA

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:) = RBUF(I,:)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PCOLLECT")

         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
       if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")
      PSZE = size(A,2)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_FLT_PCOLLECT")

      SBUF = A(1:NSZE,1:PSZE)

      CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN ARR_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_FLT_PCOLLECT")


   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_FLT_PCOLLECT"
 END SUBROUTINE ARR_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_FLT_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE CUB_FLT_ADEAL
!===================================================================================|
   SUBROUTINE CUB_FLT_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of REALS FROM A GLOBAL ARRAY                                   |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: SBUF(:,:,:),RBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE
   INTEGER, PARAMETER :: TAG = 30231 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_FLT_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")


   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),SIZE(A,2),SIZE(A,3)= ",Size(A,1),Size(A,2),Size(A,3)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),SIZE(AG,2),SIZE(AG,3)= ",Size(AG,1),Size(AG,2),Size(AG,3)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if


   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED ON THE SENDING PROC")

      PSZE=size(AG,2)
      QSZE=size(AG,3)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PDEAL")
         
         DO I=1,NSZE
            SBUF(I,:,:) = AG(GM(IP)%LOC_2_GL(I),:,:)
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE,:,:) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_FLT_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PDEAL")

      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      PSZE=size(A,2)
      QSZE=size(A,3)

      ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PDEAL")
      RBUF=0.0_SPA

      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN CUB_FLT_PDEAL")

      A(1:NSZE,1:PSZE,1:QSZE)=RBUF

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PDEAL")

   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_FLT_PDEAL"
 END SUBROUTINE CUB_FLT_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_FLT_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(SPA),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(SPA),   ALLOCATABLE,TARGET,DIMENSION(:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE CUB_FLT_ACOLLECT
!===================================================================================|
   SUBROUTINE CUB_FLT_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF REALS FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(SPA),   POINTER,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(SPA),   POINTER,DIMENSION(:,:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(SPA), ALLOCATABLE :: RBUF(:,:,:),SBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE
   INTEGER, PARAMETER :: TAG =  30232 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_FLT_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & RECEIVER")

   
   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID

      if(associated(A)) then
         write(IPT,*) "SIZE(A,1),SIZE(A,2),SIZE(A,3)   = ",Size(A,1),Size(A,2),Size(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "SIZE(AG,1),SIZE(AG,2),SIZE(AG,3)   = ",Size(AG,1),Size(AG,2),Size(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      PSZE = size(AG,2)
      QSZE = size(AG,3)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:,:) = A(I,:,:)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_REAL,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in CUB_FLT_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:,:) = RBUF(I,:,:)
            END DO
            
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PCOLLECT")
                        
         end if
       End DO
   else
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
       if (NSZE == 0)  return
      
      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      PSZE = size(A,2)
      QSZE = size(A,3)

      ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_FLT_PCOLLECT")

      SBUF=A(1:NSZE,1:PSZE,1:QSZE)

      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_REAL,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN CUB_FLT_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_FLT_PCOLLECT")
      
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_FLT_PCOLLECT"
 END SUBROUTINE CUB_FLT_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of DOUBLES FROM A GLOBAL VECTOR                                  |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:) :: AP
   REAL(DP),   POINTER,DIMENSION(:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE VEC_DBL_ADEAL
!===================================================================================|
   SUBROUTINE VEC_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A VECTOR of DOUBLES FROM A GLOBAL VECTOR                                  |
!    INTO LOCAL VECTORS AG(0:NTG) --> A(0:NT) BY MAPPING GM                         |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE

!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:) :: A
   REAL(DP),   POINTER,DIMENSION(:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:),RBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE
   INTEGER, PARAMETER :: TAG = 30311 ! UNIQUE TAG FOR VEC_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_DBL_PDEAL"


   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")

      

   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID

        write(IPT,*) "SIZE(A)                 = ",Size(A)
        write(IPT,*) "SIZE(AG)                = ",SIZE(AG)

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1)           = ",Size(AG,1)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"=============MAP INFO=============="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

           
     END if

   

   if (MYID .EQ. SENDID) then
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED BUT NEED NOT BE GLOBAL SIZE FOR ALL PROCS")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle
         
         ALLOCATE(SBUF(NSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PDEAL")
         
         DO I=1,NSZE
            SBUF(I) = AG(GM(IP)%LOC_2_GL(I))
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("Send Error in VEC_DBL_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PDEAL")

      END DO
   else
      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
      
      ALLOCATE(RBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PDEAL")
         

      CALL MPI_RECV(RBUF,NSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN VEC_DBL_PDEAL")
      
      A(1:NSZE)=RBUF

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PDEAL")
      
   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_DBL_PDEAL"
 END SUBROUTINE VEC_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE VEC_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF DOUBLES FROM LOCAL VECTORS                                 |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:) :: AP
   REAL(DP),   POINTER,DIMENSION(:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL VEC_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE VEC_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE VEC_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT A VECTOR OF DOUBLES FROM LOCAL VECTORS                                 |
!    INTO A GLOBAL VECTOR A(0:NT) --> AG(0:NTG) BY MAPPING GM                       |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:),SBUF(:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE
   INTEGER, PARAMETER :: TAG =  30312 ! UNIQUE TAG FOR VEC_FLT_PCOLLECT
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING VEC_DBL_PCOLLECT"
   
  if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL COLLECT AS THE&
        & RECEIVER")


     if(DBG_SET(dbg_sbrio)) then
   
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "RECVID                  = ",RECVID

        if(associated(A)) then
           write(IPT,*) "SIZE(A,1)            = ",Size(A,1)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1)           = ",Size(AG,1)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if

   if (MYID .EQ. RECVID) then
        if(.not. associated(AG))&
             & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = A(I)
            END DO
         else
            ALLOCATE(RBUF(NSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in VEC_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I)) = RBUF(I)
            END DO
            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
         end if
       End DO
   else
      
      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0)  return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(SBUF(NSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
      SBUF = A(1:NSZE)

      CALL MPI_SEND(SBUF,NSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN VEC_DBL_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN VEC_DBL_PCOLLECT")
   end if
   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END VEC_DBL_PCOLLECT"
 END SUBROUTINE VEC_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
 SUBROUTINE ARR_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE ARR_DBL_ADEAL
!===================================================================================|
 SUBROUTINE ARR_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:) :: A
   REAL(DP),   POINTER,DIMENSION(:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:,:),RBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE
   INTEGER, PARAMETER :: TAG = 30321 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DP_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")


   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID


        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),  SIZE(A,2)   = ",Size(A,1),Size(A,2)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),  SIZE(AG,2)   = ",Size(AG,1),Size(AG,2)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if


   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED ON THE SENDING PROC")

      PSZE=size(AG,2)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0)  cycle
         
         ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PDEAL")
         
         DO I=1,NSZE
            SBUF(I,:) = AG(GM(IP)%LOC_2_GL(I),:)
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE,:) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_DBL_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PDEAL")
      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")
      PSZE=size(A,2)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PDEAL")
      RBUF=0.0_DP

      CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN ARR_DBL_PDEAL")

      A(1:NSZE,1:PSZE)=RBUF

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PDEAL")
   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_DBL_PDEAL"
 END SUBROUTINE ARR_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE ARR_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL ARR_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE ARR_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE ARR_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:,:), SBUF(:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE
   INTEGER, PARAMETER :: TAG =  30322 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DBL_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & RECEIVER")

   
   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "SIZE(A,1),  SIZE(A,2)   = ",Size(A,1),Size(A,2)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "SIZE(AG,1),  SIZE(AG,2)   = ",Size(AG,1),Size(AG,2)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      PSZE = size(AG,2)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:) = A(I,:)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:) = RBUF(I,:)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
         end if
       End DO
   else

      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")
      PSZE = size(A,2)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
      
      ALLOCATE(SBUF(NSZE,PSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
      SBUF=A(1:NSZE,1:PSZE)

      CALL MPI_SEND(SBUF,NSZE*PSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND&
           & IN ARR_DBL_PCOLLECT")

      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PCOLLECT")

   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END ARR_DBL_PCOLLECT"
 END SUBROUTINE ARR_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
 SUBROUTINE CUB_DBL_ADEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AGP,AP)
 END SUBROUTINE CUB_DBL_ADEAL
!===================================================================================|
 SUBROUTINE CUB_DBL_PDEAL(MYID,SENDID,NPROCS,GM,AG,A)
!===================================================================================|
!    DEAL A ARRAY of DOUBLES FROM A GLOBAL ARRAY                                    |
!    INTO LOCAL ARRAYS AG(0:NTG) --> A(0:NT) BY MAPPING GM                          |
!    UPON COMPLETION EACH PROCESSOR HAS A                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,SENDID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: A
   REAL(DP),   POINTER,DIMENSION(:,:,:),INTENT(IN)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: SBUF(:,:,:),RBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE
   INTEGER, PARAMETER :: TAG = 30331 ! UNIQUE TAG FOR ARR_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING ARR_DP_PDEAL"

   if((MYID .GT. NPROCS) .AND. (MYID .NE. SENDID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & SENDER")


   if(DBG_SET(dbg_sbrio)) then
         
        write(IPT,*) "MYID                    = ",MYID
        write(IPT,*) "NPROCS                  = ",NPROCS
        write(IPT,*) "SENDID                  = ",SENDID


        if(associated(A)) then
           write(IPT,*) "SIZE(A,1),SIZE(A,2),SIZE(A,3)   = ",Size(A,1),Size(A,2),Size(A,3)
        else
           write(IPT,*) "A is not associated"
        end if
        
        if(associated(AG)) then
           write(IPT,*) "SIZE(AG,1),SIZE(AG,2),SIZE(AG,3)   = ",Size(AG,1),Size(AG,2),Size(AG,3)
        else
           write(IPT,*) "AG is not associated"
        end if

        DO I = 1,NPROCS
        write(IPT,*)"==================================="
        write(IPT,*) "ID                      = ",I
        write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
        write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
        END DO
        write(IPT,*)"==================================="

     END if


   if (MYID .EQ. SENDID) then
      
      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO DEAL MUST ALREAD&
           &Y BE ASSOCIATED ON THE SENDING PROC")

      PSZE=size(AG,2)
      QSZE=size(AG,3)

      DO IP = 1 , NPROCS

         NSZE = GM(IP)%NSIZE
         if (NSZE == 0)  cycle
         
         ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN CUB_DBL_PDEAL")
         
         DO I=1,NSZE
            SBUF(I,:,:) = AG(GM(IP)%LOC_2_GL(I),:,:)
         END DO

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

            A(1:NSZE,:,:) = SBUF
         else
            DEST = IP - 1
            CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_Send Error in ARR_DBL_PDEAL")
         end if

         DEALLOCATE(SBUF,STAT=IERR)
         if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN CUB_DBL_PDEAL")
      END DO
   else

      SOURCE = SENDID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM DEAL MUST ALREADY BE ASSOCIATED")
      PSZE=size(A,2)
      QSZE=size(A,3)

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
      
      ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN CUB_DBL_PDEAL")
      RBUF=0.0_DP

      CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD, STAT, IERR)
      if(IERR /=0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_RECV IN CUB_DBL_PDEAL")

      A(1:NSZE,1:PSZE,1:QSZE)=RBUF

      DEALLOCATE(RBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN CUB_DBL_PDEAL")
   end if
   

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_DBL_PDEAL"
 END SUBROUTINE CUB_DBL_PDEAL
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE CUB_DBL_ACOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AP
   REAL(DP),   POINTER,DIMENSION(:,:,:) :: AGP
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(DP),   ALLOCATABLE,TARGET,DIMENSION(:,:,:)   :: AG
   IF(ALLOCATED(A))  AP  => A
   IF(ALLOCATED(AG)) AGP => AG
   CALL CUB_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,AP,AGP)
 END SUBROUTINE CUB_DBL_ACOLLECT
!===================================================================================|
   SUBROUTINE CUB_DBL_PCOLLECT(MYID,RECVID,NPROCS,GM,A,AG)
!===================================================================================|
!    COLLECT AN ARRAY OF DOUBLES FROM A LOCAL ARRAYS                                  |
!    INTO A GLOBAL ARRAY A(0:NT) --> AG(0:NTG) BY MAPPING GM                        |
!    UPON COMPLETION ONE PROCESSOR HAS AG                                           |
!===================================================================================|
   IMPLICIT NONE
!------------------------------------------------------------------------------
   INTEGER,   INTENT(IN)           :: MYID,RECVID,NPROCS
   TYPE(MAP), INTENT(IN)           :: GM(NPROCS)
   REAL(DP),   POINTER,DIMENSION(:,:,:),INTENT(IN) :: A
   REAL(DP),   POINTER,DIMENSION(:,:,:)   :: AG
!------------------------------------------------------------------------------
   REAL(DP), ALLOCATABLE :: RBUF(:,:,:),SBUF(:,:,:)
   INTEGER   STAT(MPI_STATUS_SIZE),IERR,I,IP,DEST,SOURCE,NSZE,PSZE,QSZE
   INTEGER, PARAMETER :: TAG =  30332 ! UNIQUE TAG FOR CUB_FLT_PDEAL
!------------------------------------------------------------------------------

   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "STARTING CUB_DBL_PCOLLECT"
   
   if((MYID .GT. NPROCS) .AND. (MYID .NE. RECVID) )&
        &  CALL FATAL_ERROR("IOPROC CAN ONLY CALL DEAL AS THE&
        & RECEIVER")

   
   if(DBG_SET(dbg_sbrio)) then
      
      write(IPT,*) "MYID                    = ",MYID
      write(IPT,*) "NPROCS                  = ",NPROCS
      write(IPT,*) "RECVID                  = ",RECVID
      
      if(associated(A)) then
         write(IPT,*) "SIZE(A,1),SIZE(A,2),SIZE(A,3)   = ",Size(A,1),Size(A,2),Size(A,3)
      else
         write(IPT,*) "A is not associated"
      end if
      
      if(associated(AG)) then
         write(IPT,*) "SIZE(AG,1),SIZE(AG,2),SIZE(AG,3)   = ",Size(AG,1),Size(AG,2),Size(AG,3)
      else
         write(IPT,*) "AG is not associated"
      end if
      
      DO I = 1,NPROCS
         write(IPT,*)"==================================="
         write(IPT,*) "ID                      = ",I
         write(IPT,*) "GM(I)%NSIZE             = ",GM(I)%NSIZE
         write(IPT,*) "SIZE(GM(I)%LOC_2_GL)    = ",SIZE(GM(I)%LOC_2_GL)
      END DO
      write(IPT,*)"==================================="
      
   END if
   
   
   
   if (MYID .EQ. RECVID) then

      if(.not. associated(AG))&
           & CALL FATAL_ERROR("POINTER (AG) PASSED TO COLELCT MUST ALREADY BE ASSOCIATED")
      
      PSZE = size(AG,2)
      QSZE = size(AG,3)

      DO IP = 1 , NPROCS
         NSZE = GM(IP)%NSIZE
         if (NSZE == 0) cycle

         if (IP == MYID) then
            if(.not. associated(A))&
                 & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

            if(Ubound(A,1)<NSZE .or. Size(A,2) .NE. PSZE .or. Size(A,3) .NE. QSZE) &
                 & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")
                        
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:,:) = A(I,:,:)
            END DO
         else
            ALLOCATE(RBUF(NSZE,PSZE,QSZE),STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
            RBUF = 0

            SOURCE = IP - 1
            CALL MPI_RECV(RBUF,NSZE*PSZE*QSZE,MPI_DP,SOURCE,TAG,MPI_COMM_WORLD,STAT,IERR)
            IF (IERR /= 0) CALL FATAL_ERROR("MPI_RECV Error in ARR_DBL_PCOLLECT")
            
            DO I=1,NSZE
               AG(GM(IP)%LOC_2_GL(I),:,:) = RBUF(I,:,:)
            END DO

            DEALLOCATE(RBUF,STAT=IERR)
            if (IERR /= 0) CALL FATAL_ERROR("RECEIVER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
         end if
       End DO
   else

      DEST = RECVID - 1
      NSZE = GM(MYID)%NSIZE
      if (NSZE == 0) return

      if(.not. associated(A))&
           & CALL FATAL_ERROR("POINTER (A) RETURNED FROM COLLECT MUST ALREADY BE ASSOCIATED")

      if(Ubound(A,1)<NSZE) &
           & CALL FATAL_ERROR("POINTER (A) UBOUND DOES NOT MATCH MAP")

      PSZE = size(A,2)
      QSZE = size(A,3)
      ALLOCATE(SBUF(NSZE,PSZE,QSZE),STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT ALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
      SBUF = A(1:NSZE,1:PSZE,1:QSZE)

      CALL MPI_SEND(SBUF,NSZE*PSZE*QSZE,MPI_DP,DEST,TAG,MPI_COMM_WORLD, IERR)
      if(IERR /= 0) CALL FATAL_ERROR("PROCESSOR HIT AN ERROR DURING MPI_SEND IN CUB_DBL_PCOLLECT")
      
      DEALLOCATE(SBUF,STAT=IERR)
      if (IERR /= 0) CALL FATAL_ERROR("SENDER CAN NOT DEALLOCATE MEMORY IN ARR_DBL_PCOLLECT")
   end if

   
   if(DBG_SET(dbg_sbr)) &
        & write(IPT,*) "END CUB_DBL_PCOLLECT"
 END SUBROUTINE CUB_DBL_PCOLLECT
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!===================================================================================|
   SUBROUTINE SORT(N,M,SZE)
!===================================================================================|
!    SORT ELEMENTS IN N AND RETURN ORDER IN M                                       |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN)             :: SZE
   INTEGER, INTENT(INOUT)          :: N(SZE)
   INTEGER, INTENT(OUT)            :: M(SZE)
!------------------------------------------------------------------------------
   INTEGER I,LAST,J(1)
   REAL(SP)  VALM
   INTRINSIC MINLOC
!------------------------------------------------------------------------------

   M = 0
   LAST = 0
   DO I=1,SZE
     J = MINLOC(N,MASK=N > LAST)
     M(I) = J(1)
     LAST = N(J(1))
   END DO

   RETURN
   END SUBROUTINE SORT
!===================================================================================|


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   INTEGER FUNCTION GETLOC(GLOC,MAP,N)
!===================================================================================|
!    DETERMINE LOCAL IDENTITY OF ELEMENT/NODE I USING MAP                           |
!===================================================================================|

!------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: GLOC,N
   INTEGER, INTENT(IN) :: MAP(N)
   INTEGER  I,FOUND
!------------------------------------------------------------------------------
   FOUND = 0
   DO I=1,N
     IF(MAP(I)==GLOC) FOUND = I
   END DO

   GETLOC = FOUND
   RETURN 
   END FUNCTION GETLOC
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   SUBROUTINE ADD_MAP2LIST(HEAD,MP,Gsize)
     IMPLICIT NONE

     TYPE(MAPLINK), target :: HEAD
     TYPE(MAP),TARGET, ALLOCATABLE ::MP(:)
     integer, intent(in) :: Gsize
     TYPE(MAPLINK), POINTER :: current, previous
     integer CNT
     
     CNT = 1
     previous => HEAD
     current => previous%next
     
     DO
        IF(.NOT. ASSOCIATED(CURRENT)) EXIT

        CURRENT  => CURRENT%NEXT
        PREVIOUS => PREVIOUS%NEXT
        CNT = CNT +1
        IF(CNT > 100) CALL FATAL_ERROR&
             &("ADD_MAP_TO_LIST: LOOP COUNT EXCEEDED 100; STOP!")

     END DO
     
     allocate(previous%next)
     previous%next%next =>current
     previous%next%GLOBAL_SIZE = Gsize
     previous%next%MAP => MP

   END SUBROUTINE ADD_MAP2LIST
!===================================================================================|

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!==============================================================================|
   FUNCTION FIND_MAP(HEAD,Gsize,FOUND) RESULT(MP)
     IMPLICIT NONE
     TYPE(MAPLINK),target :: HEAD
     TYPE(MAP), POINTER :: MP(:)
     LOGICAL :: FOUND 
     integer, intent(IN) :: Gsize
     TYPE(MAPLINK), POINTER :: current, previous
     INTEGER CNT

     previous => HEAD
     current => previous%next
     CNT = 0
     FOUND = .FALSE.
     nullify(mp)

     DO
        IF(.NOT. ASSOCIATED(CURRENT)) RETURN

        IF(.NOT. ASSOCIATED(CURRENT%MAP))&
             CALL FATAL_ERROR("FIND_MAP: MAP LIST LINK HAS UNASSOCIATE&
             &D MAP!")

        IF(DBG_SET(DBG_MPI)) WRITE(IPT,*) "SEARCHING MAPS: GLOBAL_SIZE= ",CURRENT&
             &%GLOBAL_SIZE,"; Gsize= ",Gsize

        if(Gsize == CURRENT%GLOBAL_SIZE) then
           MP => CURRENT%MAP
           FOUND = .TRUE.
           RETURN
        end if

        CURRENT  => CURRENT%NEXT
        
        IF (CNT > 100)CALL FATAL_ERROR&
             &("FIND_MAP: LOOP COUNT EXCEEDED 100; STOP!")

     END DO

   END FUNCTION FIND_MAP

# else

   INTERFACE PPRINT
      MODULE PROCEDURE PPRINT_ARR
      MODULE PROCEDURE PPRINT_VEC
   END INTERFACE

   INTERFACE APRINT
      MODULE PROCEDURE APRINT_ARR
      MODULE PROCEDURE APRINT_VEC
   END INTERFACE
!===================================================================================|
   CONTAINS   !!INCLUDED SUBROUTINES FOLLOW
!===================================================================================|
# endif

!==============================================================================|
!  WRITE OUT VARIABLE INFORMATION TO LOCAL FILES                               |
!									       |
!  USAGE EXAMPLES                                                              |
!									       |
!  write u velocity at surface in triangle 256 to file fort.306 with iteration |
!  I1 = LBOUND(U,1) ; I2 = UBOUND(U,1)                                         |
!  CALL PPRINT(306,I1,I2,KB,U,"element",256,1,1,FLOAT(IINT))                   |
!									       |
!  I1 = LBOUND(EL,1) ; I2 = UBOUND(EL,1)                                       |
!  write surface elevation at node 233 to file fort.409 with time in hours     |
!  CALL PPRINT(406,I1,I2,1,EL,"node",233,1,1,THOUR)                            |
!									       |
!  I1 = LBOUND(T1,1) ; I2 = UBOUND(T1,1)                                       |
!  write vertical distribution of salinity at node 422 to file fort.433        |
!  CALL PPRINT(433,I1,I2,KB,T1,"node",422,1,KBM1,THOUR)                        |
!									       |
!  ARGUMENT LIST                                                               |
!       PPRINT(IUNIT,LB1,UB1,UB2,VARP,VART,ILOC,K1,K2,REF)                     |
!    1.) IUNIT - UNIT NUMBER FOR OUTPUT FILE (MUST BE >= 300 .and. <7000)      |
!    2.) LB1   - LOWER BOUND OF 1ST ARGUMENT OF ARRAY TO PRINT (USUAlLY 0)     |
!    3.) LB2   - UPPER BOUND OF 1ST ARRAY DIMENSION (USUALLY NT OR MT)         |
!        NOTE: LB1/LB2 CAN BE DETERMINE AUTOMATICALLY WITH LBOUND/UBOUND       |
!    4.) UB2   - UPPER BOUND OF SECOND ARRAY DIMENSION                         | 
!        UB2   = 1 FOR SURFACE ARRAYS LIKE EL,UA                               |
!        UB2   = KB FOR 3D ARRAYS LIKE U/V                                     |
!    5.) VARP  = VARIABLE TO PRINT (ARRAY NAME = U,V,WW,EL,T1,RHO1, etc)       |
!    6.) VART  = VARIABLE LOCATION ("element" or "node")                       |
!    7.) ILOC  = INDEX OF ELEMENT/NODE TO PRINT                                | 
!    8.) K1    = LOWER RANGE OF SIGMA LEVEL TO PRINT                           |
!    9.) K2    = UPPER RANGE OF SIGMA LEVEL TO PRINT                           |
!        K1 = 1,K2 = 1 FOR SURFACE VALUES ONLY                                 |
!        K1 = 1,K2 = KBM1 FOR ALL LEVELS                                       |
!   10.) REF   = REFERENCE VALUE FOR DATA (MUST BE FLOAT)                      |
!        REF = THOURS FOR CALCULATION TIME IN HOURS                            |
!        REF = FLOAT(IINT) FOR ITERATION NUMBER                                |
!   11.) IPT = UNIT TO WRITE ERRORS TO (USE IPT)                               |
!==============================================================================|

!==============================================================================|
   SUBROUTINE APRINT_VEC(IUNIT,VARP,VART,NOW,ILOC,MSG)
     USE MOD_TIME
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC
     REAL(SP), ALLOCATABLE, INTENT(IN),TARGET :: VARP(:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG
     REAL(SP), POINTER:: VARP_O(:)
     
     VARP_O =>  VARP
     
     IF(PRESENT(MSG))THEN
        CALL PPRINT_VEC(IUNIT,VARP_O,VART,NOW,ILOC,MSG)
     ELSE
        CALL PPRINT_VEC(IUNIT,VARP_O,VART,NOW,ILOC)
     END IF

   END SUBROUTINE APRINT_VEC
   
!==============================================================================|
   SUBROUTINE PPRINT_VEC(IUNIT,VARP,VART,NOW,ILOC,MSG)
     USE CONTROL
     USE LIMS

     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC
     REAL(SP), POINTER, INTENT(IN) :: VARP(:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG

     CHARACTER(LEN=80),parameter :: VAR_E = "element"
     CHARACTER(LEN=80),parameter :: VAR_N = "node"

     CHARACTER(LEN=100) :: STRNG
     CHARACTER(LEN=20) :: short

     INTEGER :: I,J,K,PROCMAX,II,IBND,Kopt,IERR
     LOGICAL :: PRINT_PROC

     !==============================================================================|

     !------------------------------------------------------------------------------|
     !  Process Iunit for Errors                                                    |
     !------------------------------------------------------------------------------|
     IF(IUNIT /= IPT .and. (IUNIT < 300 .or. IUNIT > 7000) )THEN
        CALL FATAL_ERROR('ERROR IN PPRINT',&
             & 'FILE UNIT < 300 AND UNIT > 7000 ARE RESERVED FOR FVCOM I/O',&
             & 'PLEASE INCREASE IUNIT TO 300+')
     END IF

     !------------------------------------------------------------------------------|
     !  Process Vartype for Errors                                                  |
     !------------------------------------------------------------------------------|
     IF(VART /= VAR_E .AND. VART /= VAR_N)THEN
        CALL FATAL_ERROR('VART IN PPRINT NOT CORRECT :'//TRIM(VART),&
             & 'SHOULD BE "'//trim(var_e)//'" or "'//trim(var_n)//'"')
     END IF


     !------------------------------------------------------------------------------|
     !  Process string output                                                       |
     !------------------------------------------------------------------------------|
     IF(PRESENT(MSG)) STRNG=TRIM(MSG)//"; IINT"


     IF(abs(IINT) .lt. 1000) THEN
        WRITE(SHORT,'(I5)') IINT
     ELSE IF(abs(IINT) .lt. 1000000) THEN
        WRITE(SHORT,'(I8)') IINT
     ELSE
        WRITE(SHORT,*) IINT
     END IF

     IF(USE_REAL_WORLD_TIME) THEN
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Date/Time:"&
             &//TRIM(WRITE_DATETIME(NOW,3,TIMEZONE))//"; ILOC= "
     ELSE
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Time(s):"
        WRITE(SHORT,'(f16.8)') SECONDS(NOW)
        STRNG = TRIM(STRNG)//TRIM(SHORT)//"; ILOC="
     END IF

     WRITE(SHORT,'(I8)') ILOC
     STRNG = TRIM(STRNG)//TRIM(SHORT)//"; VALUE="

     !------------------------------------------------------------------------------|
     !  Single Processor Case                                                       |
     !------------------------------------------------------------------------------|
     IF(NPROCS == 1)THEN
        WRITE(IUNIT,*) TRIM(STRNG),VARP(ILOC)
     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Element Based Variable (u,v,ww, etc)              |
     !      Transform to Local Element ID with "ELID"                               |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_e .AND. ELID(ILOC) /= 0)THEN

        WRITE(IUNIT,*) TRIM(STRNG),VARP(ELID(ILOC))

     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Node Based Variable (s1,t1,rho1,e1, etc)          |
     !      Transform to Local Node ID with "NLID"                                  |
     !      If Node is Interprocessor Boundary Node, Choose Processor with Highest  |
     !      ID Number to Write Values to File                                       |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_n .AND. NLID(ILOC) > 0)THEN

        PRINT_PROC = .TRUE.
        IF(NDE_ID(NLID(ILOC)) == 1)THEN   !!BOUNDARY NODE

           DO II=1,NBN
              IF(BN_LST(II) == ILOC) IBND = II
           END DO

           PROCMAX = 10000
           DO J=1,NPROCS
              IF(BN_NEY(IBND,J)==1) THEN
                 IF(J < PROCMAX) PROCMAX = J
              END IF
           END DO

           IF(PROCMAX /=  MYID) PRINT_PROC = .FALSE.  !!NOT RESPONSIBLE FOR OUTPUT 
        END IF

        IF(PRINT_PROC)THEN
           WRITE(IUNIT,*) TRIM(STRNG),VARP(NLID(ILOC))       
        END IF

     END IF


# if defined (MULTIPROCESSOR)
     IF(.NOT. IOPROC) CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
# endif

     RETURN
   END SUBROUTINE PPRINT_VEC
!==============================================================================|
!==============================================================================|
   SUBROUTINE APRINT_ARR(IUNIT,VARP,VART,NOW,ILOC,K1,K2,MSG)
     USE MOD_TIME
     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC,K1
     INTEGER,  INTENT(IN),OPTIONAL :: K2
     REAL(SP), ALLOCATABLE, INTENT(IN),TARGET :: VARP(:,:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG
     REAL(SP), POINTER:: VARP_O(:,:)
     
     VARP_O =>  VARP
     
     IF (PRESENT(K2)) THEN
        IF(PRESENT(MSG))THEN
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1,K2,MSG)
        ELSE
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1,K2)
        END IF
     ELSE
        IF(PRESENT(MSG))THEN
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1,K1,MSG)
        ELSE
           CALL PPRINT_ARR(IUNIT,VARP_O,VART,NOW,ILOC,K1)
        END IF
     END IF
   END SUBROUTINE APRINT_ARR
   
!==============================================================================|
   SUBROUTINE PPRINT_ARR(IUNIT,VARP,VART,NOW,ILOC,K1,K2,MSG)
     USE CONTROL
     USE LIMS

     IMPLICIT NONE
     TYPE(TIME), INTENT(IN) :: NOW
     INTEGER,  INTENT(IN) :: IUNIT,ILOC,K1
     INTEGER,  INTENT(IN),OPTIONAL :: K2
     REAL(SP), POINTER, INTENT(IN) :: VARP(:,:)
     CHARACTER(LEN=*), INTENT(IN) :: VART
     CHARACTER(LEN=*), INTENT(IN),OPTIONAL :: MSG

     CHARACTER(LEN=80),parameter :: VAR_E = "element"
     CHARACTER(LEN=80),parameter :: VAR_N = "node"

     CHARACTER(LEN=100) :: STRNG
     CHARACTER(LEN=20) :: short

     INTEGER :: I,J,K,PROCMAX,II,IBND,Kopt,IERR
     LOGICAL :: PRINT_PROC

     !==============================================================================|

     !------------------------------------------------------------------------------|
     !  Process Iunit for Errors                                                    |
     !------------------------------------------------------------------------------|
     IF(IUNIT /= IPT .and. (IUNIT < 300 .or. IUNIT > 7000) )THEN
        CALL FATAL_ERROR('ERROR IN PPRINT',&
             & 'FILE UNIT < 300 AND UNIT > 7000 ARE RESERVED FOR FVCOM I/O',&
             & 'PLEASE INCREASE IUNIT TO 300+')
     END IF

     !------------------------------------------------------------------------------|
     !  Process Vartype for Errors                                                  |
     !------------------------------------------------------------------------------|
     IF(VART /= VAR_E .AND. VART /= VAR_N)THEN
        CALL FATAL_ERROR('VART IN PPRINT NOT CORRECT :'//TRIM(VART),&
             & 'SHOULD BE "'//trim(var_e)//'" or "'//trim(var_n)//'"')
     END IF

     !------------------------------------------------------------------------------|
     !  Process optional sigma level range                                          |
     !------------------------------------------------------------------------------|

     IF(PRESENT(K2)) THEN
        KOPT=K2
     ELSE
        KOPT=K1
     END IF

     !------------------------------------------------------------------------------|
     !  Process string output                                                       |
     !------------------------------------------------------------------------------|
     IF(PRESENT(MSG)) STRNG=TRIM(MSG)//"; IINT"


     IF(abs(IINT) .lt. 1000) THEN
        WRITE(SHORT,'(I5)') IINT
     ELSE IF(abs(IINT) .lt. 1000000) THEN
        WRITE(SHORT,'(I8)') IINT
     ELSE
        WRITE(SHORT,*) IINT
     END IF

     IF(USE_REAL_WORLD_TIME) THEN
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Date/Time:"&
             &//TRIM(WRITE_DATETIME(NOW,3,TIMEZONE))//"; ILOC "
     ELSE
        STRNG = TRIM(STRNG)//TRIM(SHORT)//", Time(s):"
        WRITE(SHORT,'(f16.8)') SECONDS(NOW)
        STRNG = TRIM(STRNG)//TRIM(SHORT)//"; ILOC="
     END IF

     WRITE(SHORT,'(I8)') ILOC
     STRNG = TRIM(STRNG)//TRIM(SHORT)//"; VALUES="
     !------------------------------------------------------------------------------|
     !  Single Processor Case                                                       |
     !------------------------------------------------------------------------------|
     IF(NPROCS == 1)THEN
        WRITE(IUNIT,*) TRIM(STRNG)
        WRITE(IUNIT,*) (VARP(ILOC,K),K=K1,KOPT)
     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Element Based Variable (u,v,ww, etc)              |
     !      Transform to Local Element ID with "ELID"                               |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_e .AND. ELID(ILOC) /= 0)THEN

        WRITE(IUNIT,*) TRIM(STRNG)
        WRITE(IUNIT,*) (VARP(ELID(ILOC),K),K=K1,KOPT)

     END IF

     !------------------------------------------------------------------------------|
     !  Multi Processor Case with Node Based Variable (s1,t1,rho1,e1, etc)          |
     !      Transform to Local Node ID with "NLID"                                  |
     !      If Node is Interprocessor Boundary Node, Choose Processor with Highest  |
     !      ID Number to Write Values to File                                       |
     !------------------------------------------------------------------------------|

     IF(NPROCS /= 1 .AND. VART == var_n .AND. NLID(ILOC) > 0)THEN

        PRINT_PROC = .TRUE.
        IF(NDE_ID(NLID(ILOC)) == 1)THEN   !!BOUNDARY NODE

           DO II=1,NBN
              IF(BN_LST(II) == ILOC) IBND = II
           END DO

           PROCMAX = 10000
           DO J=1,NPROCS
              IF(BN_NEY(IBND,J)==1) THEN
                 IF(J < PROCMAX) PROCMAX = J
              END IF
           END DO

           IF(PROCMAX /=  MYID) PRINT_PROC = .FALSE.  !!NOT RESPONSIBLE FOR OUTPUT 
        END IF

        IF(PRINT_PROC)THEN
           WRITE(IUNIT,*) TRIM(STRNG)
           WRITE(IUNIT,*) (VARP(NLID(ILOC),K),K=K1,KOPT)        
        END IF

     END IF

# if defined (MULTIPROCESSOR)
     IF(.NOT. IOPROC) CALL MPI_BARRIER(MPI_FVCOM_GROUP,IERR)
# endif

     RETURN
   END SUBROUTINE PPRINT_ARR
!==============================================================================|

END MODULE MOD_PAR  
