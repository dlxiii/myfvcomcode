! $Id: mod_prec.F,v 1.2.2.3 2008/04/03 16:53:50 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.2.2.3 $
!===============================================================================!
! DEFINE FLOATING POINT PRECISION USING KIND                                    !
!===============================================================================!
MODULE MOD_PREC
#  if defined (MULTIPROCESSOR)
   use mpi
#  endif

   IMPLICIT NONE

!--Single Precision Coding------------------------------------------------------!
#  if !defined (DOUBLE_PRECISION)
   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,30)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_F = MPI_REAL
# endif
#  else 
!--Double Precision Coding------------------------------------------------------!
   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(12,300)
#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_F = MPI_DOUBLE_PRECISION
# endif
#  endif

  
   INTEGER, PARAMETER :: SPA = SELECTED_REAL_KIND(6,30)
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(12,300)

#  if defined (MULTIPROCESSOR)
   INTEGER, PARAMETER :: MPI_DP = MPI_DOUBLE_PRECISION
# endif


END MODULE MOD_PREC

