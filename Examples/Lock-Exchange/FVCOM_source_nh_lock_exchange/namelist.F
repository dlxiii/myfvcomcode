SUBROUTINE NAMELIST
  USE MOD_UTILS
  USE CONTROL

  USE MOD_INPUT

  USE MOD_NESTING

# if defined (SEMI_IMPLICIT)
  USE MOD_SEMI_IMPLICIT
# endif
# if defined (BALANCE_2D)
  USE MOD_BALANCE_2D
# endif  
# if defined (ONE_D_MODEL)
  USE MOD_ONEDTIDE
# endif
# if defined (WATER_QUALITY)
  USE MOD_WQM
# endif    
# if defined (BioGen)
  USE MOD_BIO_3D
# endif    
# if defined (NH)
  USE NON_HYDRO, ONLY: NAME_LIST_INITIALIZE_NH, NAME_LIST_PRINT_NH, NAME_LIST_READ_NH
# endif
  IMPLICIT NONE


  !==============================================================================!
  ! SET DEFAULT VALUES IN NAME LIST                                                   
  !==============================================================================!
  CALL NAME_LIST_INITIALIZE

  CALL NAME_LIST_INITIALIZE_NEST

# if defined (SEMI_IMPLICIT)    
  CALL NAME_LIST_INITIALIZE_SEMI
# endif
# if defined (DYE_RELEASE)    
  CALL NAME_LIST_INITIALIZE_DYE
# endif
# if defined (BALANCE_2D)    
  CALL NAME_LIST_INITIALIZE_BALANCE
# endif
# if defined (ONE_D_MODEL)    
  CALL NAME_LIST_INITIALIZE_ONED
# endif
# if defined (WATER_QUALITY)    
  CALL NAME_LIST_INITIALIZE_WQM
# endif
# if defined (BioGen)    
  CALL NAME_LIST_INITIALIZE_BIO
# endif
# if defined (NH)
  CALL NAME_LIST_INITIALIZE_NH
# endif

  ! IF FVCOM IS ONLY PRINTING A BLANK NAME LIST FOR A NEW CASE:
  if (BLANK_NAMELIST) then
     CALL NAME_LIST_PRINT

     CALL NAME_LIST_PRINT_NEST

#    if defined (SEMI_IMPLICIT)    
     CALL NAME_LIST_PRINT_SEMI
#    endif
#    if defined (DYE_RELEASE)    
     CALL NAME_LIST_PRINT_DYE
#    endif
#    if defined (BALANCE_2D)    
     CALL NAME_LIST_PRINT_BALANCE
#    endif
#    if defined (ONE_D_MODEL)    
     CALL NAME_LIST_PRINT_ONED
#    endif
#    if defined (WATER_QUALITY)    
     CALL NAME_LIST_PRINT_WQM
#    endif
#    if defined (BioGen)    
     CALL NAME_LIST_PRINT_BIO
#    endif
#    if defined (NH)
     CALL NAME_LIST_PRINT_NH
#    endif

     CALL PSHUTDOWN
  end if

  !==============================================================================!
  !   SETUP MODEL RUN PARAMETERS                                                 !
  !==============================================================================!

  !READ DATA IN THE NAME LIST FILE
  CALL NAME_LIST_READ ! ALL PROCS READ THIS

  CALL NAME_LIST_READ_NEST

# if defined (SEMI_IMPLICIT)    
  CALL NAME_LIST_READ_SEMI
# endif
# if defined (DYE_RELEASE)    
  CALL NAME_LIST_READ_DYE
# endif
# if defined (BALANCE_2D)    
  CALL NAME_LIST_READ_BALANCE
# endif
# if defined (ONE_D_MODEL)    
  CALL GET_UMAX
# endif
# if defined (WATER_QUALITY)    
  CALL NAME_LIST_READ_WQM
# endif
# if defined (NH)
  CALL NAME_LIST_READ_NH
# endif

  !PRINT THE NAME LIST DATA TO THE SCREEN FOR THE LOG
  IF(DBG_SET(DBG_LOG)) CALL NAME_LIST_PRINT 


END SUBROUTINE NAMELIST
