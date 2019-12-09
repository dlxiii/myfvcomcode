










!************************************************************************
!                                                                       *
   SUBROUTINE SWBOUN ( )         
   END SUBROUTINE SWBOUN
 
!*********************************************************************
!                                                                    *
   SUBROUTINE BCFILE (FBCNAM, BCTYPE)                     

!  (This subroutine has not been used and tested yet)
!                                                                    *
!*********************************************************************
!
!     Reads file data for boundary condition
!
!*********************************************************************
!
   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE M_BNDSPEC                                                       
!
   IMPLICIT NONE

   CHARACTER FBCNAM *(*), BCTYPE *(*)

   INTEGER :: ISTATF, NDSL, NDSD, IOSTAT, IERR, NBOUNC, NANG, NFRE
   INTEGER :: IBOUNC, DORDER
   INTEGER :: IENT,IOPTT
   INTEGER :: NHEDF, NHEDT, NHEDS, IFRE , IANG
   INTEGER :: NQUANT, IQUANT, IBC, II, NBGRPT_PREV,IIPT2
   REAL    :: XP, YP, XP2, YP2
   REAL    :: FREQHZ, DIRDEG, DIRRD1,DIRRAD, EXCV
   CHARACTER BTYPE *4, HEDLIN *80
   LOGICAL         CCOORD                                            

!
   NDSL = 0
   IIPT2 = 0                                                           
!  open data file
   NDSD = 0
   IOSTAT = 0
   CALL FOR (NDSD, FILENM, 'OF', IOSTAT)
!
!     --- initialize array BFILED of BSPFIL                               
!      BSPFIL%BFILED = 0                                                   
!
!     start reading from the data file
      READ (NDSD, '(A)') HEDLIN
!      IF (EQCSTR(HEDLIN,'TPAR')) THEN                                     
!        BTYPE  = 'TPAR'
!        ISTATF = 1
!        IOPTT  = 1
!        NBOUNC = 1
!        NANG   = 0
!        NFRE   = 0
!        NHEDF  = 0
!        NHEDT  = 0
!        NHEDS  = 0
!        DORDER = 0
!        ALLOCATE(BSPFIL%BSPFRQ(NFRE))                                     
!        ALLOCATE(BSPFIL%BSPDIR(NANG))                                     
!        IF (NSTATM.EQ.0) CALL MSGERR (3,'time information not allowed in stationary mode')
!        NSTATM = 1
!      ELSE IF (EQCSTR(HEDLIN,'SWAN')) THEN                                
!      ELSE
!        CALL MSGERR (3, 'unsupported boundary data file')
!      ENDIF
!
!      ALLOCATE(BSPFIL%BSPLOC(NBOUNC))                                     
!      DO IBC = 1, NBOUNC
!         BSPFIL%BSPLOC(IBC) = NBSPEC + IBC                                
!      ENDDO
      NBSPEC = NBSPEC + NBOUNC
!
!     store file reading parameters in array BFILED
!
!      IF (ITEST.GE.80) WRITE(PRINTF,81) NBFILS, NBSPEC,(BSPFIL%BFILED(II), II=1,16)                                  
  81  FORMAT (' array BFILED: ', 2I4, 2(/,8I10))
!
      RETURN
      END SUBROUTINE BCFILE
