#
# Configuration file for sst2grd Code
#
            SHELL = /bin/sh


#--------------------------------------------------------------------------
#  NETCDF STUFF
#--------------------------------------------------------------------------

         EXEC     = xrestart_file
         MYOPT      = # -static
         FVCOMINCS = -I$(TOPDIR)
         FVCOMLIBS = -L$(TOPDIR) -lfvcom

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 .viz

.F.o:
	$(CPP) $(CPPARGS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(MYINCS) $*.f90
	\rm $*.f90

#--------------------------------------------------------------------------
#  Filter Source Code.
#--------------------------------------------------------------------------

MYLIBS = $(FVCOMLIBS) $(LIBS)
MYINCS = $(FVCOMINCS) $(INCS)


MOD1  = mod_restart.F
MAIN  = main.F

SRCS = $(MOD1) $(MAIN) 

OBJS = $(SRCS:.F=.o)

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) $(MYOPT) $(OBJS) $(MYLIBS) -o $(EXEC)

#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		mv makedepends makedepends.orig
		makedepf90  $(SRCS) >> makedepends


#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
		tar cvf restart_file.tar *.F makefile 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod xrestart_file

clobber:	
		make clean
		/bin/rm -f *.f90

#--------------------------------------------------------------------------
#  Common rules for all Makefiles - do not edit.
#--------------------------------------------------------------------------

emptyrule::

#--------------------------------------------------------------------------
#  Empty rules for directories that do not have SUBDIRS - do not edit.
#--------------------------------------------------------------------------

install::
	@echo "install in $(CURRENT_DIR) done"

install.man::
	@echo "install.man in $(CURRENT_DIR) done"

Makefiles::

includes::
include ../../make.inc
include ./makedepends
# DO NOT DELETE
