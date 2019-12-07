! $Id: genmap_nest.F,v 1.1.2.1 2008/06/05 14:58:38 jqi Exp $
! $Name: New_Input $
! $Revision: 1.1.2.1 $
!==============================================================================|
SUBROUTINE GENMAP_NEST
!==============================================================================|
!
! CREATE A GLOBAL TO LOCAL MAP FOR DOMAIN NESTING
! USES DATA READ INTO: MGL_NEST - The number of nesting nodes
!                      NGL_NEST - The number of nesting elements
!                      MID_NEST_GL - The index to global nesting node  boundary
!                      NID_NEST_GL - The index to global nesting element
!                     
! Creates:             MAP LINK LIST ENTRY FOR IO
!                      NID_NEST,MID_NEST - The transform from the local nesting list to the local domain
!
!==============================================================================|
  USE MOD_NESTING
  USE MOD_PAR
  USE LIMS
  USE CONTROL
  IMPLICIT NONE
  integer :: SENDER,RECVER, ierr, I, NCNT, NSZE, I1, status

  INTEGER, POINTER :: TEMP1(:),TEMP2(:)
  INTEGER, POINTER :: E_G2L(:), N_G2L(:)

  TYPE(MAP), TARGET, SAVE,ALLOCATABLE, DIMENSION(:) :: E_MAP,N_MAP

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "START: GENMAP_NEST"

  IF(SERIAL) THEN

     MT_NEST = MGL_NEST
     NT_NEST = NGL_NEST

     MID_NEST        => MID_NEST_GL
     MID_NEST_OUTPUT => MID_NEST_GL

     NID_NEST        => NID_NEST_GL
     NID_NEST_OUTPUT => NID_NEST_GL

     if (dbg_set(dbg_sbr)) &
          & write(IPT,*) "END: GENMAP_NEST - serial"
     return
  END IF

# if defined (MULTIPROCESSOR)


  IF(.not. IOPROC) THEN


!============================================
! Make a list of the local Nesting nodes
!============================================

     ALLOCATE(TEMP1(MGL_NEST));      TEMP1=0
     ALLOCATE(TEMP2(MGL_NEST));      TEMP2=0
     
     NCNT = 0
     !!SET UP LOCAL OPEN BOUNDARY NODES
     DO I=1,MGL_NEST
        I1 = NLID_X( MID_NEST_GL(I) )
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     MT_NEST = NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(MID_NEST(MT_NEST),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NEST: can not allocate:MID_NEST")

     ALLOCATE(MID_NEST_OUTPUT(MT_NEST),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NEST: can not allocate:MID_NEST_OUTPUT")

     ALLOCATE(N_G2L(MT_NEST),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NEST: can not allocate:N_G2L")

     MID_NEST(1:NCNT)        = TEMP1(1:NCNT)     !MID_NEST_GL(1:NCNT)
     MID_NEST_OUTPUT(1:NCNT) = NGID_X(MID_NEST)
     
     N_G2L(1:NCNT) = TEMP2(1:NCNT)

     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)

!============================================
! Make a list of the local Nesting elements
!============================================

     ALLOCATE(TEMP1(NGL_NEST));      TEMP1=0
     ALLOCATE(TEMP2(NGL_NEST));      TEMP2=0

     NCNT = 0
     !!SET UP LOCAL OPEN BOUNDARY ELEMENTS
     DO I=1,NGL_NEST
        I1 = ELID_X( NID_NEST_GL(I) )
        IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
           TEMP2(NCNT) = I
        END IF
     END DO
     
     ! SET LOCAL NUMBER OF BOUNDARY NODES
     NT_NEST = NCNT

     ! SET GLOBAL TO LOCAL MAP FOR THIS DOMAIN
     ALLOCATE(NID_NEST(NT_NEST),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NEST: can not allocate:NID_NEST")

     ALLOCATE(NID_NEST_OUTPUT(NT_NEST),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NEST: can not allocate:NID_NEST_OUTPUT")

     ALLOCATE(E_G2L(NT_NEST),stat=status)
     if(status /= 0) call fatal_error("GENMAP_NEST: can not allocate:E_G2L")

     NID_NEST(1:NCNT)= TEMP1(1:NCNT)
     NID_NEST_OUTPUT(1:NCNT) = EGID_X(NID_NEST)

     E_G2L(1:NCNT) = TEMP2(1:NCNT)

     DEALLOCATE(TEMP1)
     DEALLOCATE(TEMP2)

  END IF


  !==============================================================================|
  !   SET UP ELEMENT MAPPING FOR GLOBAL 2 LOCAL TRANSFER OF BC'S                 | 
  !   BOUNDARY MAP :: BCMAP(NPROCS)                                              |
  !     BCMAP(1-->NPROCS)%NSIZE  :: NUMBER OF BOUNDARY NODES IN EACH DOM         |
  !     BCMAP(1-->NPROCS)%LOC_2_GL(NSIZE) :: LOCAL TO GLOBAL MAPPING IN EACH DOM |
  !==============================================================================|

  ! ELEMENTS:
  ALLOCATE(E_MAP(NPROCS))

  !--Determine Number of Elements for Each Processor
  DO I=1,NPROCS    
     IF(MYID == I) E_MAP(I)%NSIZE = NT_NEST
     SENDER = I - 1
     CALL MPI_BCAST(E_MAP(I)%NSIZE,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
  END DO

  !--Construct Mapping Array for Each Processor 
  DO I=1,NPROCS
     NSZE = E_MAP(I)%NSIZE
     ALLOCATE(E_MAP(I)%LOC_2_GL(0:NSZE))
     E_MAP(I)%LOC_2_GL=0
     
     if(myid == I) E_MAP(I)%LOC_2_GL(1:NSZE) =  E_G2L(1:NSZE)
     SENDER = I - 1
     CALL MPI_BCAST(E_MAP(I)%LOC_2_GL(1:NSZE),NSZE,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)

  END DO

  DEALLOCATE(E_G2L)

  ! NODES
  ALLOCATE(N_MAP(NPROCS))

  !--Determine Number of Elements for Each Processor
  DO I=1,NPROCS    
     IF(MYID == I) N_MAP(I)%NSIZE = MT_NEST
     SENDER = I - 1
     CALL MPI_BCAST(N_MAP(I)%NSIZE,1,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)
  END DO

  !--Construct Mapping Array for Each Processor 
  DO I=1,NPROCS
     NSZE = N_MAP(I)%NSIZE
     ALLOCATE(N_MAP(I)%LOC_2_GL(0:NSZE))
     N_MAP(I)%LOC_2_GL = 0

     if(myid == I) N_MAP(I)%LOC_2_GL(1:NSZE) =  N_G2L(1:NSZE)
     SENDER = I - 1
     CALL MPI_BCAST(N_MAP(I)%LOC_2_GL(1:NSZE),NSZE,MPI_INTEGER,SENDER,MPI_COMM_WORLD,IERR)

  END DO

  DEALLOCATE(N_G2L)


  ! NOW ADD BOTH MAPS TO THE HALO MAP LIST
!  CALL ADD_MAP2LIST(HALO_MAPS,E_MAP,NGL_NEST)

!  CALL ADD_MAP2LIST(HALO_MAPS,N_MAP,MGL_NEST)



  CALL ADD_MAP2LIST(INTERNAL_MAPS,E_MAP,NGL_NEST)

  CALL ADD_MAP2LIST(INTERNAL_MAPS,N_MAP,MGL_NEST)

! DO NOT DEALLOCATE MAP!!!

# endif

  if (dbg_set(dbg_sbr)) &
       & write(IPT,*) "END: GENMAP_NEST - parallel"   
  RETURN
END SUBROUTINE GENMAP_NEST

