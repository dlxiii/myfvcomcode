# $Id: make.inc,v 1.1.2.19 2008/05/28 17:54:44 dstuebe Exp $
# $Name: New_Input $
# $Revision: 1.1.2.19 $
#==========================================================================
#
#    WELCOME TO FVCOM
#
#    TO BUILD FVCOM, PLEASE SET THE FOLLOWING IN THIS FILE:
#        TOPDIR - the directory this make.inc file is in
#        LIBDIR - the directroy where libraries are installed
#        INCIDR - the directroy where include files are installed 
#
#        CHOOSE YOUR MODEL OPTIONAS - The Def Flags!
#
#        CHOOSE A COMPILER FROM THE LIST OR CREATE YOUR OWN
#
#    You can also use the makefile to build and install some of the libraries
#    used by fvcom. Set the INSTALLDIR and choose the 'LOCAL INSTALL' option
#    below. Select an non-mpi compiler from the list and run 'make libs_ext'
#
#    Good Luck!
#
#========== TOPDIR ========================================================
# TOPDIR is the directory in which this make file and the fvcom source reside

  TOPDIR  = ../FVCOM_source_nh_lock_exchange

# TOPDIR must be set!
#========== INSTALLDIR =====================================================
# INSTALLDIR is the directory where you wish to install external libraries
# The default is in the $(TOPDIR)/libs/install, but you can choose...
            INSTALLDIR    = $(TOPDIR)/libs/install
#===========================================================================
# PREPROCESSOR OPTIONS FOR CPP
            DEF_FLAGS     = -P -C -traditional 
#===========================================================================
#----------------------------------------------
# A few common environments for building FVCOM at SMAST: uncomment one of them...
#----------------------------------------------
####################### New Hydra #########################################
# Use the environmental variables, LIBPATH and INCLUDEPATH, set by MODULE 
# to find the packages used by this build.
colon=:
empty=
dashI= $(empty) -I
dashL= $(empty) -L
# ### UNCOMMENT HERE!
             LIBDIR       =  -L$(subst $(colon),$(dashL),$(LIBPATH))
             INCDIR       =  -I$(subst $(colon),$(dashI),$(INCLUDEPATH))

###########################################################################

# MAO SERVER
            # For em64t computers at SMAST (guppy)
#            LIBDIR        = -L/hosts/mao/usr/medm/lib64
#            INCDIR        = -I/hosts/mao/usr/medm/include64
            # For i386 computers at SMAST (salmon and minke)
#            LIBDIR        = -L/hosts/mao/usr/medm/lib
#            INCDIR        = -I/hosts/mao/usr/medm/include

# SALMON01 SERVER
            # For em64t computers at SMAST (guppy)
#            LIBDIR        = -L/hosts/salmon01/data00/medm/em64t/lib64
#            INCDIR        = -I/hosts/salmon01/data00/medm/em64t/include64
            # For i386 computers at SMAST (salmon and minke)
#            LIBDIR	= -L/hosts/salmon01/data00/medm/lib
#            INCDIR	= -I/hosts/salmon01/data00/medm/include

# Typhoeus and Hydra
#             LIBDIR       =  -L/usr/local/lib
#             INCDIR       =  -I/usr/local/include




# LOCAL INSTAL
#             LIBDIR       =  -L$(INSTALLDIR)/lib
#             INCDIR       =  -I$(INSTALLDIR)/include


#--------------------------------------------------------------------------
#	STANDARD LIBRARIES FOR DATA AND TIME IN fVCOM:
#
	    DTLIBS	= -ljulian
	    DTINCS	=  
#
#--------------------------------------------------------------------------
#        NETCDF OUTPUT      NETCDF IS NOW REQUIRED TO COMPILE FVCOM
#                           DUMP OUTPUT INTO NETCDF FILES (yes/no)
#                           REQUIRES SYSTEM DEPENDENT NETCDF LIBRARIES
#                           COMPILED WITH SAME F90 COMPILER
#                           SET PATH TO LIBRARIES WITH IOLIBS      
#                           SET PATH TO INCLUDE FILES (netcdf.mod) WITH IOINCS
#--------------------------------------------------------------------------
             IOLIBS       = -lnetcdf
             IOINCS       =
#--------------------------------------------------------------------------


#==========================================================================
#  BEGIN USER DEFINITION SECTION
#==========================================================================
#        SELECT MODEL OPTIONS
#          SELECT FROM THE FOLLOWING OPTIONS BEFORE COMPILING CODE
#          SELECT/UNSELECT BY COMMENTING/UNCOMMENTING LINE (#)
#          CODE MUST BE CLEANED (with "make clean") AND
#          RECOMPILED IF NEW SET OF OPTIONS IS DESIRED
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------

#            FLAG_1 =  -DDOUBLE_PRECISION 

#--------------------------------------------------------------------------
#        SPHERICAL          SELECT SPHERICAL COORDINATES FOR INTEGRATION
#                           DEFAULT: CARTESIAN
#                           UNCOMMENT TO SELECT SPHERICAL COORDINATES
#--------------------------------------------------------------------------

#            FLAG_2 = -DSPHERICAL

#--------------------------------------------------------------------------
#        FLOODYING/DRYING   INCLUDE WET/DRY TREATMENT OF DOMAIN
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH
#                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE
#                           DEFAULT: NO FLOODYING/DRYING INCLUDED
#                           UNCOMMENT TO INCLUDE FLOODYING/DRYING
#--------------------------------------------------------------------------

#             FLAG_3 = -DWET_DRY

#--------------------------------------------------------------------------
#        MULTI_PROCESSOR    INCLUDES PARALLELIZATION WITH MPI
#                           REQUIRES LINKING MPI LIBRARIES OR COMPILING
#                           WITH A PRELINKED SCRIPT (mpif90/mpf90/etc)
#                           DEFAULT: NO PARALLEL CAPABILITY
#                           UNCOMMENT TO INCLUDE MPI PARALLEL CAPABILITY
#--------------------------------------------------------------------------

             FLAG_4 = -DMULTIPROCESSOR

             PARLIB       = -lmetis

#--------------------------------------------------------------------------
#        WATER_QUALITY      INCLUDE EPA WATER QUALITY MOD
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH
#                           VARIABLE WQM_ON IN INPUT FILE
#                           DEFAULT: NO WATER QUALITY MODEL
#                           UNCOMMENT TO INCLUDE WATER QUALITY MODEL
#--------------------------------------------------------------------------

#             FLAG_5 = -DWATER_QUALITY

#--------------------------------------------------------------------------
#        PROJECTION         A Fortran90 wrapper for the Cartographic projection
#                           Software, proj4. 
#                           Proj can be obtained from:
#                           http://www.remotesensing.org/proj/
#                           Thanks to: USGS
#
#                           The Proj fortran bindings can be obtained from:
#                           http://forge.nesc.ac.uk/projects/glimmer/
#                           Thanks to: Magnus Hagdorn (Magnus.Hagdorn@ed.ac.uk)
#
#                           !! NOTE THAT THE PROJ 4 LIBRARY MUST BE IN YOUR
#                           LD_LIBRARY_PATH FOR DYNAMIC LOADING!!
#
#--------------------------------------------------------------------------

#             FLAG_6 = -DPROJ

#             PROJLIBS       = -lfproj4 -lproj -lm
#             PROJINCS       =

#--------------------------------------------------------------------------
#        DATA_ASSIMILATION  INCLUDE NUDGING BASED DATA ASSIMILATION FOR
#                           CURRENT/TEMP/SALINITY/SST
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH 
#                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE 
#                           DEFAULT: NO DATA ASSIMILATION INCLUDED 
#                           UNCOMMENT TO INCLUDE DATA ASSIMILATION 
#--------------------------------------------------------------------------

#             FLAG_7 = -DDATA_ASSIM
#--------------------------------------------------------------------------
#     DATA_ASSIMILATION_OI  INCLUDE OPTIMAL INTERPOLATION(OI) BASED DATA 
#                           ASSIMILATION FOR CURRENT/TEMP/SALINITY/SST
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH 
#                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE 
#                           DEFAULT: NO OI DATA ASSIMILATION INCLUDED 
#                           UNCOMMENT TO INCLUDE OI DATA ASSIMILATION 
#--------------------------------------------------------------------------

#             FLAG_7 = -DDATA_ASSIM_OI

#             OILIB  = -lmkl_lapack -lmkl_em64t -lguide -lpthread

# For Build on em64t computer (Guppy) 
#             LIBDIR = $(LIBDIR) -L/usr/local/lib64
# For Build on Cluster (Typhoeus and Hydra)
#             LIBDIR = $(LIBDIR) -L/usr/local/lib/em64t
# For i386 computers at SMAST (salmon and minke)
#   NO NEED TO ADD ANYTHING LIBS ARE IN THE DEFAULT PATH

#--------------------------------------------------------------------------
#        IN UPWIND LEAST SQUARE SCHEME:
#        LIMITED_NO: NO LIMITATION
#        LIMITED_1 : FIRST ORDER LIMITATION
#        LIMITED_2 : SECOND ORDER LIMITATION(  )
#   !!!!!! ONLY ONE OF THE FLAGS BELOW CAN BE AND MUST BE CHOSEN
#--------------------------------------------------------------------------

            FLAG_8 =  -DLIMITED_NO
#            FLAG_8 =  -DLIMITED_1
#            FLAG_8 =  -DLIMITED_2

#--------------------------------------------------------------------------
#  Semi-Implicit time stepping method
#--------------------------------------------------------------------------

#             FLAG_9 = -DSEMI_IMPLICIT
#             include ${PETSC_DIR}/bmake/common/variables


#---------------------------------------------------------------------------
#       SOLID BOUNDARY     IF GCN, NO GHOST CELL
#                          IF GCY1, GHOST CELL IS SYMMETRIC RELATIVE TO BOUNDARY
#                          CELL EDGE
#                          IF GCY2, GHOST CELL IS SYMMETRIC RELATIVE TO MIDDLE
#                          POINT OF THE BOUNDARY CELL EDGE
#   !!!!!! ONLY ONE OF THE FLAGS BELOW CAN BE AND MUST BE CHOSEN
#---------------------------------------------------------------------------

#             FLAG_10  = -DGCN
             FLAG_10  = -DGCY1
#             FLAG_10 = -DGCY2

#--------------------------------------------------------------------------
#        TURBULENCE MODEL   USE GOTM TURBULENCE MODEL INSTEAD OF THE ORIGINAL
#                           FVCOM MELLOR-YAMADA 2.5 IMPLEMENTATION
#                           UNCOMMENT TO USE GOTM TURBULENCE MODEL
#--------------------------------------------------------------------------

#            FLAG_11 = -DGOTM
#            GOTMLIB       = -L../GOTM_source -lturbulence -lutil
#            GOTMINCS      = -I../GOTM_source  

#--------------------------------------------------------------------------
#       EQUILIBRIUM TIDE
#--------------------------------------------------------------------------

#            FLAG_12 = -DEQUI_TIDE

#--------------------------------------------------------------------------
#       ATMOSPHERIC TIDE
#--------------------------------------------------------------------------

#            FLAG_13 = -DATMO_TIDE

#--------------------------------------------------------------------------
#        Using A fully multidimensional positive definite advection
#        transport algorithm with small implicit diffusion. 
#        Based on Smolarkiewicz, P. K; Journal of Computational
#        Physics, 54, 325-362, 1984
#--------------------------------------------------------------------------

#            FLAG_15 = -DMPDATA

#--------------------------------------------------------------------------
#         Run Two-D Barotropic Mode Only
#--------------------------------------------------------------------------

#             FLAG_16 = -DTWO_D_MODEL

#--------------------------------------------------------------------------
#         Output 2-D Momentum Balance Checking
#--------------------------------------------------------------------------

#            FLAG_17 = -DBALANCE_2D  

#--------------------------------------------------------------------------
#           OPEN BOUNDARY FORCING TYPE 
#           DEFAULT: OPEN BOUNDARY NODE WATER ELEVATION FORCING
#           UNCOMMENT TO SELECT BOTH OPEN BOUNDARY NODE WATER ELEVATION
#           FORCING AND OPEN BOUNDARY VOLUME TRANSPORT FORCING 
#---------------------------------------------------------------------------

#            FLAG_19 = -DMEAN_FLOW

#--------------------------------------------------------------------------
#           OUTPUT TIDAL INFORMATION AT NTIDENODE and NTIDECELL
#           FOR MEANFLOW CALCULATION.
#---------------------------------------------------------------------------

#            FLAG_20 = -DTIDE_OUTPUT

#--------------------------------------------------------------------------
#           dye release      
#---------------------------------------------------------------------------

#            FLAG_21 = -DDYE_RELEASE  

#--------------------------------------------------------------------------
#        SUSPENDED SEDIMENT MODEL:   UNCOMMENT TO INCLUDE MODEL
#--------------------------------------------------------------------------

#            FLAG_22 =  -DSEDIMENT

#--------------------------------------------------------------------------
#       Kalman Filters     TWO STEPS INCLUDED FOR REDUCED KALMAN FILTERS DATA
#                          ASSIMILATION. FIRSTLY USERS SHOULD OPEN FLAG_22 IN
#                          ORDER TO LET PROGRAM PREPARE INPUT DATA. AND THEN,
#                          SECOND OPTION FLAG_12 NEED BY APPLIED FOR ASSIMILATION
#                          BE NOTED: ONLY OPEN ONE OPTION AT ONE TIME                  
#--------------------------------------------------------------------------

#             FLAG_23 = -DENKF_ASSIM
#             FLAG_23 = -DRRK_PRE
#             FLAG_23 = -DRRK_ASSIM
#              FLAG_23 = -DENKF_ASSIM_EVE

#             KFLIB  = -lmkl_lapack -lmkl_em64t -lguide -llapack -lblas 

# For Build on em64t computer (Guppy) 
#             LIBDIR = $(LIBDIR) -L/usr/local/lib64
# For Build on Cluster (Typhoeus and Hydra)
#             LIBDIR = $(LIBDIR) -L/usr/local/lib/em64t
# For i386 computers at SMAST (salmon and minke)
#   NO NEED TO ADD ANYTHING LIBS ARE IN THE DEFAULT PATH

#--------------------------------------------------------------------------
#         Run One-D Mode with Biological Model
#--------------------------------------------------------------------------

#             FLAG_24 = -DONE_D_MODEL

#--------------------------------------------------------------------------
#       GENERAL BIOLOGICAL MODEL:    UNCOMMENT TO INCLUDE MODEL
#-------------------------------------------------------------------------- 
#            FLAG_25 = -DBioGen
#            BIOLIB       = -L$(TOPDIR)/BIO_source -lbiomodule
#            BIOINCS      = -I$(TOPDIR)/BIO_source/

#--------------------------------------------------------------------------
#         Dynamic/Thermodynamic Ice                     
#--------------------------------------------------------------------------

# NOTE: Must use -DSPHERICAL  and -DHEAT_FLUX
#            FLAG_26 = -DICE  

#--------------------------------------------------------------------------
#         CALCULATE THE NET HEAT FLUX IN MODEL                           
#--------------------------------------------------------------------------

#            FLAG_27 = -DAIR_PRESSURE  

#--------------------------------------------------------------------------
# Visit online simulation mode
#--------------------------------------------------------------------------  

#            FLAG_28 = -DVISIT

#      VISITLIB     = -lm -ldl -lsimf -lsim -lpthread
#      VISITLIBPATH = $(LIBDIR)/visit
#      VISITINCPATH = $(INCDIR)/visit


# USE DEVELOPER INSTALL VISIT
#      VISITLIBPATH = -L/hosts/salmon01/data00/medm/visit/src_64/visit/archives -L/hosts/salmon01/data00/medm/visit/src_64/visit/lib
#      VISITLIB     =  -lm -ldl -lsimf -lsim -lpthread
#      VISITINC     = -I/hosts/salmon01/data00/medm/visit/src_64/visit/include/visit

#--------------------------------------------------------------------------
#        NON-HYDROSTATIC MODEL:
#--------------------------------------------------------------------------

             FLAG_29 =  -DNH
             include ${PETSC_DIR}/bmake/common/variables

#--------------------------------------------------------------------------
#        SELECT COMPILER/PLATFORM SPECIFIC DEFINITIONS
#          SELECT FROM THE FOLLOWING PLATFORMS OR USE "OTHER" TO DEFINE
#          THE FOLLOWING VARIABLES:
#          CPP:  PATH TO C PREPROCESSOR 
#           FC:  PATH TO FORTRAN COMPILER (OR MPI COMPILE SCRIPT)
#          OPT:  COMPILER OPTIONS
#       MPILIB:  PATH TO MPI LIBRARIES (IF NOT LINKED THROUGH COMPILE SCRIPT)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#        COMPAQ/ALPHA Definitions
#--------------------------------------------------------------------------
#        COMPILER   = -DCOMPAQ
#        CPP        = /bin/cpp
#        FC         = f90
#        DEBFLGS    =  # -check bounds -check overflow -g 
#        OPT        =  -fast -arch ev6 -fpe1
#--------------------------------------------------------------------------
#  CRAY Definitions         
#--------------------------------------------------------------------------
#         COMPILER = -DCRAY 
#         CPP      = /opt/ctl/bin/cpp
#         FC       = f90
#         DEBFLGS  =  
#         OPT      =
#--------------------------------------------------------------------------
#   Linux/Portland Group Definitions 
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         COMPILER =  
#         FC       = pgf90
#         DEBFLGS  = -Mbounds -g -Mprof=func
#         OPT      = #-fast  -Mvect=assoc,cachesize:512000,sse  
#--------------------------------------------------------------------------
#  Intel Compiler Definitions
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         COMPILER = -DIFORT
#         FC       = ifort
#         CC       = icc
#         CXX      = icc
#         CFLAGS   = -O3
#         DEBFLGS  = #-check all
# Use 'OPT = -O0 -g' for fast compile to test the make 
# Use 'OPT = -xP' for fast run on em64t (Hydra, Typhoeus and Guppy)
# Use 'OPT = -xN' for fast run on em64t (Hydra, Typhoeus and Guppy)
#         OPT      = -O0 -g  
#         OPT      = -xN 

# Do not set static for use with visit!
#         VISOPT   = -Wl,--export-dynamic
#         LDFLAGS  = $(VISITLIBPATH)
#--------------------------------------------------------------------------
#  Intel/MPI Compiler Definitions (SMAST)      
#--------------------------------------------------------------------------
         CPP      = /usr/bin/cpp
         COMPILER = -DIFORT  
         CC       = mpicc
         CXX      = mpicxx
         CFLAGS   = -O3 
         FC       = mpif90
         DEBFLGS  = #-check all
# Use 'OPT = -O0 -g'  for fast compile to test the make 
# Use 'OPT = -xP'      for fast run on em64t (Hydra, Typhoeus and Guppy)
# Use 'OPT = -xN'      for fast run on em64t (Salmon and Minke)
#         OPT      = -O0 -g
         OPT      = -xP

# Do not set static for use with visit!
#         VISOPT   = -Wl,--export-dynamic
#         LDFLAGS  =  $(VISITLIBPATH)
#--------------------------------------------------------------------------
#  gfortran defs 
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp 
#         COMPILER = -DGFORTRAN
#         FC       = gfortran  -O3 
#         DEBFLGS  = 
#         OPT      = 
#         CLIB     = 
#--------------------------------------------------------------------------
#  absoft / mac os x defs
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp 
#         COMPILER = -DABSOFT  
#         FC       = f90  -O3  -lU77
#         DEBFLGS  = 
#         OPT      = 
#         CLIB     = 
#--------------------------------------------------------------------------
#  IBM/AIX Definitions
#--------------------------------------------------------------------------
#         COMPILER = -DAIX 
#         CPP      = /usr/local/bin/cpp
#         FC       = mpxlf90 -qsuffix=f=f90
#         DEBFLGS  =  # -qcheck -C -g
#         OPT      = -O -qarch=pwr4 -qtune=pwr4  -bmaxdata:0x80000000 -qhot -qmaxmem=8096
#--------------------------------------------------------------------------
#  APPLE OS X/XLF Definitions (G5)
#--------------------------------------------------------------------------
#         COMPILER = -DAIX
#         CPP      = /usr/bin/cpp
#         FC       =  /opt/ibmcmp/xlf/8.1/bin/xlf90 -qsuffix=f=f90
#         DEBFLGS  =  # -qcheck -C -g
#         OPT      = -O5 -qarch=g5 -qtune=g5  -qhot -qmaxmem=8096 -qunroll=yes -Wl,-stack_size,10000000
#--------------------------------------------------------------------------
#  OTHER                   
#--------------------------------------------------------------------------
#         COMPILER =
#         CPP      = 
#         FC       = 
#         DEBFLGS  = 
#         OPT      = 
#==========================================================================
#  END USER DEFINITION SECTION
#==========================================================================
         CPPFLAGS = $(DEF_FLAGS) $(COMPILER)
         FFLAGS = $(DEBFLGS) $(OPT)
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib
         AR = ar rc
#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
           CPPARGS =  $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) \
   $(FLAG_3) $(FLAG_4) $(FLAG_5) $(FLAG_6)  \
   $(FLAG_7) $(FLAG_8) $(FLAG_9) $(FLAG_10) \
   $(FLAG_11) $(FLAG_12) $(FLAG_13) $(FLAG_14) \
   $(FLAG_15) $(FLAG_16) $(FLAG_17) $(FLAG_18) \
   $(FLAG_19) $(FLAG_20) $(FLAG_21) $(FLAG_22) \
   $(FLAG_23) $(FLAG_24) $(FLAG_25) $(FLAG_26) \
   $(FLAG_27) $(FLAG_28) $(FLAG_29) $(FLAG_LAM)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------

            LIBS  =	$(LIBDIR) $(CLIB)  $(PARLIB) $(IOLIBS)	$(DTLIBS)\
			$(MPILIB) $(GOTMLIB) $(KFLIB) $(BIOLIB)	\
			$(OILIB) $(VISITLIB) $(PROJLIBS) $(PETSC_LIB)

            INCS  = 	$(INCDIR) $(IOINCS) $(GOTMINCS) $(BIOINCS)\
			 $(VISITINCPATH) $(PROJINCS) $(DTINCS) \
			 $(PETSC_FC_INCLUDES)
