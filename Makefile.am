## Process this file with automake to produce Makefile.in
## $Id: Makefile.am,v 1.29 2006-01-23 23:07:49 flux Exp $

SUBDIRS = . img graphs doc src include scripts rcl testing shed

VERSION = @VERSION@

AM_CPPFLAGS = -DMAN_PATH=\"$(mandir)\"

configure.ac: configure.ac.in VERSION
	shed/setversion

EXTRA_DIST = \
   LICENSE VERSION README COPYING THANKS INSTALL AUTHORS \
   configure.ac.in

