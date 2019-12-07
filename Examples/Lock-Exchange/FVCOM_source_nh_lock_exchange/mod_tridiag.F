! $Id: mod_tridiag.F,v 1.2.2.2 2008/04/03 17:03:48 dstuebe Exp $
! $Name: New_Input $
! $Revision: 1.2.2.2 $
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mtridiagonal --- solver for tri-diagonal matrices \label{sec:tridiagonal}
!
! !INTERFACE:
   MODULE mtridiagonal
   use mod_prec
!
! !DESCRIPTION: 
!
!  Solves a linear system of equations with a tridiagonal matrix
!  using Gaussian elimination. 
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_tridiagonal,tridiagonal
!
! !PUBLIC DATA MEMBERS:
   real(sp), dimension(:), allocatable     :: au,bu,cu,du
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: mod_tridiag.F,v $
!  Revision 1.2.2.2  2008/04/03 17:03:48  dstuebe
!  Changed Header
!
!  Revision 1.2.2.1  2007/12/26 18:02:14  dstuebe
!  added cvs keyword tags
!
!  Revision 1.2  2007/02/28 21:40:05  dstuebe
!  Added comments about ICE
!
!  Revision 1.1.1.1  2007/02/28 19:48:18  dstuebe
!  Import of Gao's Ice model
!
!  Revision 1.1  2006/06/20 00:14:52  gcowles
!  *** empty log message ***
!
!  Revision 1.4  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:06:33  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
!  private data members
   real(sp), private, dimension(:),allocatable  ::  ru,qu
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory
!
! !INTERFACE:
   subroutine init_tridiagonal(N)
!
! !DESCRIPTION:
!  This routines allocates memory necessary to perform the Gaussian 
! elimination.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
!
!-----------------------------------------------------------------------
!BOC
   if(allocated(au))then
     deallocate(au)
     deallocate(bu)
     deallocate(cu)
     deallocate(du)
     deallocate(ru)
     deallocate(qu)
   endif

   allocate(au(0:N)) ; au = 0.
   allocate(bu(0:N)) ; bu = 0.
   allocate(cu(0:N)) ; cu = 0.
   allocate(du(0:N)) ; du = 0.
   allocate(ru(0:N)) ; ru = 0.
   allocate(qu(0:N)) ; qu = 0.

   return
   end subroutine init_tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Simplified Gaussian elimination 
!
! !INTERFACE:
   subroutine tridiagonal(N,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}. 
! The method used here is the simplified Gauss elimination, also called 
! \emph{Thomas algorithm}.  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N,fi,lt
!
! !OUTPUT PARAMETERS:
   real(sp)                                    :: value(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  $Log: mod_tridiag.F,v $
!  Revision 1.2.2.2  2008/04/03 17:03:48  dstuebe
!  Changed Header
!
!  Revision 1.2.2.1  2007/12/26 18:02:14  dstuebe
!  added cvs keyword tags
!
!  Revision 1.2  2007/02/28 21:40:05  dstuebe
!  Added comments about ICE
!
!  Revision 1.1.1.1  2007/02/28 19:48:18  dstuebe
!  Import of Gao's Ice model
!
!  Revision 1.1  2006/06/20 00:14:52  gcowles
!  *** empty log message ***
!
!  Revision 1.4  2003/03/28 09:20:36  kbk
!  added new copyright to files
!
!  Revision 1.3  2003/03/28 08:06:33  kbk
!  removed tabs
!
!  Revision 1.2  2003/03/10 08:54:16  gotm
!  Improved documentation and cleaned up code
!
!  Revision 1.1.1.1  2001/02/12 15:55:58  gotm
!  initial import into CVS
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)
   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do


   return
   end subroutine tridiagonal
!EOC

!-----------------------------------------------------------------------

   end module mtridiagonal

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------- 
