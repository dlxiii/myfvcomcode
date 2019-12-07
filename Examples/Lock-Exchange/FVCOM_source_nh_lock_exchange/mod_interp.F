! $Id: mod_interp.F,v 1.1.2.4 2008/04/03 16:45:14 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.1.2.4 $
MODULE MOD_INTERP
  USE MOD_UTILS
  USE MOD_PREC
  IMPLICIT NONE  
  
  ! SEARCH RADIUS
  REAL(SP) :: search=80000.0_SP


  ! NUMBER OF POINTS PER QUADRENT
  INTEGER, parameter :: nppq=2
  INTEGER, parameter :: quad=4

  TYPE R2PTS
     REAL(SP) :: WEIGHTS(NPPQ,QUAD)
     INTEGER :: QUAD_NUM(NPPQ,QUAD)
  END TYPE R2PTS

  TYPE INTERP_WEIGHTS
     INTEGER :: Nin
     INTEGER :: Nout
     INTEGER, ALLOCATABLE :: INDEX(:)
     TYPE(R2PTS), ALLOCATABLE:: PTW(:)
  END TYPE INTERP_WEIGHTS


!!$  INTERFACE INTERP
!!$     MODULE PROCEDURE INTERP_A
!!$     MODULE PROCEDURE INTERP_P
!!$  END INTERFACE


CONTAINS
  
  SUBROUTINE SETUP_INTERP_A(Xin,Yin,Xout,Yout,WEIGHTS, rzero)
    IMPLICIT NONE 
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Xin(:,:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Yin(:,:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Xout(:)
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Yout(:)
    REAL(SP), OPTIONAL :: rzero

    TYPE(INTERP_WEIGHTS), INTENT(OUT) :: WEIGHTS

    REAL(SP), POINTER :: XinP(:,:)
    REAL(SP), POINTER :: YinP(:,:)
    REAL(SP), POINTER :: XoutP(:)
    REAL(SP), POINTER :: YoutP(:)

    NULLIFY(XinP,YinP,XoutP,YoutP)


    IF(ALLOCATED(Xin)) XinP => Xin

    IF(ALLOCATED(Yin)) YinP => Yin

    IF(ALLOCATED(Xout)) XoutP => Xout

    IF(ALLOCATED(Yout)) YoutP => Yout

    IF(PRESENT(rzero)) THEN
       CALL SETUP_INTERP_P(XinP,YinP,XoutP,YoutP,WEIGHTS,rzero)
    ELSE
       CALL SETUP_INTERP_P(XinP,YinP,XoutP,YoutP,WEIGHTS)
    END IF
  END SUBROUTINE SETUP_INTERP_A




  SUBROUTINE SETUP_INTERP_P(Xin,Yin,Xout,Yout,WEIGHTS,rzero)
    IMPLICIT NONE 
    REAL(SP), POINTER, INTENT(IN) :: Xin(:,:)
    REAL(SP), POINTER, INTENT(IN) :: Yin(:,:)
    REAL(SP), POINTER, INTENT(IN) :: Xout(:)
    REAL(SP), POINTER, INTENT(IN) :: Yout(:)
    REAL(SP), OPTIONAL :: rzero

    TYPE(INTERP_WEIGHTS), INTENT(OUT) :: WEIGHTS


    REAL(SP), ALLOCATABLE :: Xvec(:)
    REAL(SP), ALLOCATABLE :: Yvec(:)  
    REAL(SP), ALLOCATABLE :: XTEMP(:)
    REAL(SP), ALLOCATABLE :: YTEMP(:)
    INTEGER :: I, status, lb, ub, msze,nsze
    
    INTEGER :: UBV(1)

    ! CHECK INPUT ALLOCATION STATUS
    if(.not.associated(Xin)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    if(.not.associated(Yin)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    if(.not.associated(Xout)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    if(.not.associated(Yout)) CALL FATAL_ERROR("SETUP_INTERP: INPUT ARGU&
         &MENTS MUST BE ALLOCATED!")
    

    ! GET THE DIMENSIONS
    UBV = ubound(Xout)
 
    WEIGHTS%Nout= UBV(1) 
    msze=size(Xin,1)
    nsze=size(Xin,2)

    WEIGHTS%Nin  = msze*nsze

    ! ALLOCATE SPACE FOR VECTORIZED POSITION DATA
    ALLOCATE(Xvec(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")

    ALLOCATE(Yvec(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")

    ! ALLOCATE A TEMPORARY ARRAY TO SHIFT DATA AROUND
    ALLOCATE(XTEMP(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")

    ALLOCATE(YTEMP(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")


    ! ALLOCATE THE SPACE FOR THE WEIGHTS OUTPUT AND INDEX
    ALLOCATE(WEIGHTS%INDEX(WEIGHTS%Nin), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")

    ALLOCATE(WEIGHTS%PTW(WEIGHTS%Nout), stat=status)
    IF(STATUS /= 0) CALL FATAL_ERROR("SETUP_INTERP: COULD NOT ALLOCATE SPACE")

    ! THE ORDER HERE MUST MATCH THAT IN INTERP
    DO i = 1,nsze
       lb = (i-1)*msze+1
       ub = i*msze
       XTEMP(lb:ub)=Xin(1:msze,i)
       YTEMP(lb:ub)=Yin(1:msze,i)
    END DO

    CALL SORTRX(WEIGHTS%Nin,XTEMP,WEIGHTS%INDEX)

    DO I = 1,WEIGHTS%Nin
       Xvec(I) = XTEMP(WEIGHTS%INDEX(I))
       Yvec(I) = YTEMP(WEIGHTS%INDEX(I))
    END DO
    

    IF(PRESENT(RZERO)) THEN
       DO I = 1, WEIGHTS%Nout
          CALL GEN_WTS(Xvec,Yvec,Xout(I),Yout(I),WEIGHTS%INDEX,WEIGHTS&
               &%PTW(I),rzero)
       END DO
    ELSE
       DO I = 1, WEIGHTS%Nout
          CALL GEN_WTS(Xvec,Yvec,Xout(I),Yout(I),WEIGHTS%INDEX,WEIGHTS&
               &%PTW(I))
       END DO
    END IF
       

  END SUBROUTINE SETUP_INTERP_P


  ! THIS CODE IS SETUP TO WORK FOR REGULAR GRIDDED DATA
  ! TO USE IT ON OBSERVATION DATA IT WOULD NEED TO BE MODIFIED TO
  ! HANDLE STRANGE CASES WHERE THERE ARE CLUSTERS OF DATA POINTS
  SUBROUTINE GEN_WTS(Xvec,Yvec,Xout,Yout,INDX,PTW,rzero)
    IMPLICIT NONE
    real(SP), INTENT(in) :: Xout,Yout
    real(SP), ALLOCATABLE, INTENT(IN) :: Xvec(:),Yvec(:)
    INTEGER, ALLOCATABLE, INTENT(IN) :: INDX(:)
    TYPE(R2PTS), INTENT(OUT) :: PTW
    REAL(SP), OPTIONAL :: rzero

    INTEGER :: I,Q,K,NB,J
    
    Real(SP) :: DELX, DELY, D, SUMWGHT

    REAL(SP) :: D_MAX(QUAD)
    REAL(SP) :: DIST(nppq,QUAD)
    INTEGER :: NDM
    INTEGER :: UBV(1)

    UBV = ubound(Xvec)
    NB = UBV(1)

    PTW%WEIGHTS=0
    PTW%QUAD_NUM=0

    D_MAX = HUGE(D)
    DIST  = HUGE(D)

    NDM=1

    k = 0
    DO I=1,NB

       DELX = Xvec(I) - Xout

       ! SINCE XVEC IS IN INCREASING ORDER, ONCE DELX IS GREATER THAN
       ! SEARCH WE ARE DONE
       If (DELX.GT.SEARCH) EXIT

       If (ABS(DELX).LE.SEARCH) Then
          DELY = Yvec(I) - Yout
          D = SQRT(DELX**2+DELY**2)
          
          ! WE FOUND A POINT IN THE SEARCH RADIUS
          If (D.LE.SEARCH) Then
             
             k = k +1
             ! DEAL WITH D == 0.0 here?
             Q=QUADRANT(DELX,DELY)

             ! REPLACE THAT VALUE IF IT IS CLOSER THAN THE LAST ONE
             IF(D .LT. D_MAX(Q))THEN
                ! RECORD THE DISTANCE AN UPDATE THE MAX FOR THAT QUADRANT
                UBV=MAXLOC(DIST(:,Q))
                NDM= UBV(1)
                DIST(NDM,Q) = D
                D_MAX(Q) = MAXVAL(DIST(:,Q))
                
                ! RECORD THE INDEX IN QUAD_NUM
                PTW%QUAD_NUM(NDM,Q) = I
             END IF
          End If
       End If
    END DO

    if (k == 0) call warning("FOUND NO Values in interpolation serach radius!")

    IF (PRESENT(RZERO)) THEN
       DIST = Dist/rzero
    END IF


    DO I = 1,quad
       DO K =1 , nppq
          
          IF (DIST(k,I) < HUGE(D)) THEN
!             PTW%WEIGHTS(K,I) = 1/(DIST(K,I)**2)
             PTW%WEIGHTS(K,I) = 1/(.5_SP+DIST(K,I)**2)
          END IF
       END DO
    END DO
    
    ! NORMALIZE THE TOTAL WEIGHT TO ONE
    SUMWGHT = SUM(PTW%WEIGHTS)
    PTW%WEIGHTS = PTW%WEIGHTS / SUMWGHT

  END SUBROUTINE GEN_WTS


  INTEGER FUNCTION QUADRANT(DX,DY)
    IMPLICIT NONE
    REAL(SP), INTENT(IN):: DX,DY
    INTEGER :: IBIT1, IBIT2
    
    IBIT1 = (INT(SIGN(1.,DX))+1) / 2
    IBIT2 = (INT(SIGN(1.,DY))+1) / 2
    !    quadrant number:
    QUADRANT = 2 * IBIT1 + IBIT2 + 1
  END FUNCTION QUADRANT
  

!*****************************************************************
!C Here's a hybrid QuickSort I wrote a number of years ago.  It's
!C based on suggestions in Knuth, Volume 3, and performs much better
!C than a pure QuickSort on short or partially ordered input arrays.  

  SUBROUTINE SORTRX(N,DATA,INDX)
!C===================================================================
!C
!C     SORTRX -- SORT, Real input, indx output
!C
!C
!C     Input:  N     INTEGER
!C             DATA  REAL
!C
!C     Output: INDX INTEGER (DIMENSION N)
!C
!C This routine performs an in-memory sort of the first N elements of
!C array DATA, returning into array INDX the indices of elements of
!C DATA arranged in ascending order.  Thus,
!C
!C    DATA(INDX(1)) will be the smallest number in array DATA;
!C    DATA(INDX(N)) will be the largest number in DATA.
!C
!C The original data is not physically rearranged.  The original order
!C of equal input values is not necessarily preserved.
!C
!C===================================================================
!C
!C SORTRX uses a hybrid QuickSort algorithm, based on several
!C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
!C "pivot key" [my term] for dividing each subsequence is chosen to be
!C the median of the first, last, and middle values of the subsequence;
!C and the QuickSort is cut off when a subsequence has 9 or fewer
!C elements, and a straight insertion sort of the entire array is done
!C at the end.  The result is comparable to a pure insertion sort for
!C very short arrays, and very fast for very large arrays (of order 12
!C micro-sec/element on the 3081K for arrays of 10K elements).  It is
!C also not subject to the poor performance of the pure QuickSort on
!C partially ordered data.
!C
!C Created:  15 Jul 1986  Len Moss
!C
!C===================================================================
    
    INTEGER, INTENT(IN)::  N
    REAL(SP)    :: DATA(N)
    INTEGER :: INDX(N)
    
    INTEGER  :: LSTK(31),RSTK(31),ISTK
    INTEGER  :: L,R,I,J,P,INDXP,INDXT
    REAL     :: DATAP
    
!C     QuickSort Cutoff
!C
!C     Quit QuickSort-ing when a subsequence contains M or fewer
!C     elements and finish off at end with straight insertion sort.
!C     According to Knuth, V.3, the optimum value of M is around 9.
 
    INTEGER, PARAMETER ::   M = 9
    
!C===================================================================
!C
!C     Make initial guess for INDX
 
    DO 50 I=1,N
       INDX(I)=I
50     CONTINUE
 
!C     If array is short, skip QuickSort and go directly to
!C     the straight insertion sort.
 
       IF (N.LE.M) GOTO 900
       
!C===================================================================
!C
!C     QuickSort
!C
!C     The "Qn:"s correspond roughly to steps in Algorithm Q,
!C     Knuth, V.3, PP.116-117, modified to select the median
!C     of the first, last, and middle elements as the "pivot
!C     key" (in Knuth's notation, "K").  Also modified to leave
!C     data in place and produce an INDX array.  To simplify
!C     comments, let DATA[I]=DATA(INDX(I)).
 
!C Q1: Initialize
       ISTK=0
       L=1
       R=N
       
200    CONTINUE
 
!C Q2: Sort the subsequence DATA[L]..DATA[R].
!C
!C     At this point, DATA[l] &lt;= DATA[m] &lt;= DATA[r] for all l &lt; L,
!C     r &gt; R, and L &lt;= m &lt;= R.  (First time through, there is no
!C     DATA for l &lt; L or r &gt; R.)
 
       I=L
       J=R
 
!C Q2.5: Select pivot key
!C
!C     Let the pivot, P, be the midpoint of this subsequence,
!C     P=(L+R)/2; then rearrange INDX(L), INDX(P), and INDX(R)
!C     so the corresponding DATA values are in increasing order.
!C     The pivot key, DATAP, is then DATA[P].
 
       P=(L+R)/2
       INDXP=INDX(P)
       DATAP=DATA(INDXP)
       
       IF (DATA(INDX(L)) .GT. DATAP) THEN
          INDX(P)=INDX(L)
          INDX(L)=INDXP
          INDXP=INDX(P)
          DATAP=DATA(INDXP)
       ENDIF
       
       IF (DATAP .GT. DATA(INDX(R))) THEN
          IF (DATA(INDX(L)) .GT. DATA(INDX(R))) THEN
             INDX(P)=INDX(L)
             INDX(L)=INDX(R)
          ELSE
             INDX(P)=INDX(R)
          ENDIF
          INDX(R)=INDXP
          INDXP=INDX(P)
          DATAP=DATA(INDXP)
       ENDIF
 
!C     Now we swap values between the right and left sides and/or
!C     move DATAP until all smaller values are on the left and all
!C     larger values are on the right.  Neither the left or right
!C     side will be internally ordered yet; however, DATAP will be
!C     in its final position.
 
300    CONTINUE
 
!C Q3: Search for datum on left &gt;= DATAP
!C
!C     At this point, DATA[L] &lt;= DATAP.  We can therefore start scanning
!C     up from L, looking for a value &gt;= DATAP (this scan is guaranteed
!C     to terminate since we initially placed DATAP near the middle of
!C     the subsequence).
 
       I=I+1
       IF (DATA(INDX(I)).LT.DATAP) GOTO 300
       
400    CONTINUE
       
!C Q4: Search for datum on right &lt;= DATAP
!C
!C     At this point, DATA[R] &gt;= DATAP.  We can therefore start scanning
!C     down from R, looking for a value &lt;= DATAP (this scan is guaranteed
!C     to terminate since we initially placed DATAP near the middle of
!C     the subsequence).
 
       J=J-1
       IF (DATA(INDX(J)).GT.DATAP) GOTO 400
 
!C Q5: Have the two scans collided?
 
       IF (I.LT.J) THEN
 
!C Q6: No, interchange DATA[I] &lt;--&gt; DATA[J] and continue
          
          INDXT=INDX(I)
          INDX(I)=INDX(J)
          INDX(J)=INDXT
          GOTO 300
       ELSE
          
!C Q7: Yes, select next subsequence to sort
!C
!C    At this point, I &gt;= J and DATA[l] &lt;= DATA[I] == DATAP &lt;= DATA[r],
!C     for all L &lt;= l &lt; I and J &lt; r &lt;= R.  If both subsequences are
!C     more than M elements long, push the longer one on the stack and
!C     go back to QuickSort the shorter; if only one is more than M
!C     elements long, go back and QuickSort it; otherwise, pop a
!C     subsequence off the stack and QuickSort it.
 
          IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
             ISTK=ISTK+1
             LSTK(ISTK)=J+1
             RSTK(ISTK)=R
             R=I-1
          ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
             ISTK=ISTK+1
             LSTK(ISTK)=L
             RSTK(ISTK)=I-1
             L=J+1
          ELSE IF (R-J .GT. M) THEN
             L=J+1
          ELSE IF (I-L .GT. M) THEN
             R=I-1
          ELSE
!C Q8: Pop the stack, or terminate QuickSort if empty
             IF (ISTK.LT.1) GOTO 900
             L=LSTK(ISTK)
             R=RSTK(ISTK)
             ISTK=ISTK-1
          ENDIF
          GOTO 200
       ENDIF
 
900    CONTINUE
 
!C===================================================================
!C
!C Q9: Straight Insertion sort
 
       DO 950 I=2,N
          !         if(mod(I,10000).eq.0) print*,'i=',i
          IF (DATA(INDX(I-1)) .GT. DATA(INDX(I))) THEN
             INDXP=INDX(I)
             DATAP=DATA(INDXP)
             P=I-1
920          CONTINUE
             INDX(P+1) = INDX(P)
             P=P-1
             IF (P.GT.0) THEN
                IF (DATA(INDX(P)).GT.DATAP) GOTO 920
             ENDIF
             INDX(P+1) = INDXP
          ENDIF
950       CONTINUE
          
!C===================================================================
!C
!C     All done
          RETURN
!       END DO
!    END DO
  END SUBROUTINE SORTRX
!******************************


  SUBROUTINE INTERP_A(zin,Weights,zout)
    IMPLICIT NONE
    REAL(SP), ALLOCATABLE, TARGET, INTENT(IN) :: Zin(:,:)
    TYPE(INTERP_WEIGHTS), INTENT(in) :: WEIGHTS
    REAL(SP), ALLOCATABLE, TARGET :: Zout(:)

    REAL(SP), POINTER :: ZinP(:,:),ZoutP(:)

    IF (allocated(Zin)) THEN
       ZinP => Zin
    ELSE
       CALL FATAL_ERROR &
          &("INTERP: ZIN IS NOT ALLOCATED")
    END IF

    IF (allocated(Zout)) THEN
       ZoutP => Zout   !(lbound(zout):ubound(zout))
    ELSE
       CALL FATAL_ERROR &
          &("INTERP: ZOUT IS NOT ALLOCATED")
    END IF
 
    CALL INTERP_P(zinp,Weights,zoutp)

  END SUBROUTINE INTERP_A
  
  SUBROUTINE INTERP_P(zin,Weights,zout)
     REAL(SP), POINTER :: Zin(:,:)
     TYPE(INTERP_WEIGHTS), INTENT(in) :: WEIGHTS
     REAL(SP), POINTER :: Zout(:)

     REAL(SP), ALLOCATABLE :: Zvec(:)
     REAL(SP), ALLOCATABLE :: ZTEMP(:)
     
     integer :: I,status, msze,nsze, lb, ub

     IF (.NOT. ASSOCIATED(ZIN))CALL FATAL_ERROR &
          &("INTERP: ZIN IS NOT ALLOCATED")

     ! CAN'T CHECK ALLOCATION STATUS OF INTENT OUT VARIABLE!

     msze= size(Zin,1)
     nsze= size(Zin,2)

     IF (Weights%Nin /= nsze*msze) CALL FATAL_ERROR &
          &("INTERP: THE DATA INPUT SIZE DOES NOT MATCH THE WEIGHT!")

     ! MAKE VECTOR IN RIGHT ORDER...
     ALLOCATE(Zvec(0:Weights%Nin),stat=status)
     if(status /= 0) CALL FATAL_ERROR("INPTERP: COULD NOT ALLOCATE SPACE")

     ALLOCATE(ZTEMP(0:Weights%Nin),stat=status)
     if(status /= 0) CALL FATAL_ERROR("INPTERP: COULD NOT ALLOCATE SPACE")

     ! THERE MUST BE A BETTER WAY TO DO THIS?
     ! THE ORDER HERE MUST MATCH THAT IN SETUP_INTERP
     DO i = 1,nsze
        lb = (i-1)*msze+1
        ub = i*msze
        ZTEMP(lb:ub)=Zin(1:msze,i)
     END DO

     DO I = 1,WEIGHTS%Nin
        Zvec(I) = ZTEMP(WEIGHTS%INDEX(I))
     END DO

     DO I = 1,WEIGHTS%Nout

        CALL INTERP_WEIGH(ZVEC,WEIGHTS%PTW(I),Zout(I))

     END DO

   end subroutine INTERP_P

! VECTORIZE THIS - THERE IS NO REASON TO USE AN ARRAY TO STORE THE WEIGHTS!
   subroutine interp_weigh(Zvec,PTW,zval)
     implicit none
     real(SP), INTENT(OUT) :: ZVAL
     TYPE(R2PTS), INTENT(IN) :: PTW
     real(SP), ALLOCATABLE , INTENT(IN):: zvec(:)
     integer  :: i,j,n
     real(SP) :: tt

     ! WEIGHTS ARE NORMALIZED TO ONE ALREADY: DO NOT DIVIDE BY THE SUM!

     tt=0.0_SP
     do i=1,quad
        do j=1,nppq
           n=PTW%quad_num(j,i)
           tt=Zvec(n)*PTW%weights(j,i)+tt
        end do
     end do
     zval=tt

   end subroutine interp_weigh
   

 END MODULE MOD_INTERP
