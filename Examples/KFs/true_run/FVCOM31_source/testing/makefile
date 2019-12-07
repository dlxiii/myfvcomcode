# $Id: makefile,v 1.1.1.1 2010/01/03 19:36:02 jqi Exp $
# $Name:  $
# $Revision: 1.1.1.1 $
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh

SUBDIRS = nctools parallel interp


all:
	for dir in $(SUBDIRS); do (cd $$dir && make) || exit 1; done

clean:
	for dir in $(SUBDIRS); do (cd $$dir && make clean) || exit 1; done

clobber:
	for dir in $(SUBDIRS); do (cd $$dir && make clobber) || exit 1; done

