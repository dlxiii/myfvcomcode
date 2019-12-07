# $Id: makefile,v 1.1.1.1 2010/01/03 19:36:16 jqi Exp $
# $Name:  $
# $Revision: 1.1.1.1 $
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh

SUBDIRS =	tools		cnvrt_coords	edit_file	\
		interp_time	weather_data	\
		fvcom2gnome
all:
	for dir in $(SUBDIRS); do (cd $$dir && make) || exit 1; done

clean:
	for dir in $(SUBDIRS); do (cd $$dir && make clean) || exit 1; done

clobber:
	for dir in $(SUBDIRS); do (cd $$dir && make clobber) || exit 1; done

