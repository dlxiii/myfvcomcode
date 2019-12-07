! CVS: $Id: ice_model_size.F,v 1.1.1.1.2.2 2008/01/03 18:39:34 dstuebe Exp $
!=======================================================================
!BOP
!
! !MODULE: ice_model_size 
!
! !DESCRIPTION:
!
! Defines the global domain size and number of categories and layers. \\
! Code originally based on model\_size.F in POP
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:
!
      module ice_model_size
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!=======================================================================

      use lims, only : MGL

      implicit none
      save

      integer (kind=int_kind), parameter ::  &
       imt_global = 1,     &   ! i-axis size 
       jmt_global = 100,  &    ! MGL,     &   ! j-axis size
       ncat       =   5,     &   ! number of categories
       nilyr      =   4,     &   ! number of layers in a single category
       ntilay     =  ncat*nilyr ! sum of number of layers in all categories

!=======================================================================

      end module ice_model_size

!=======================================================================
