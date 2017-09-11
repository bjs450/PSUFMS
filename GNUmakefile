# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:42 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := simutesting
G4TARGET := $(name)
G4EXLIB := true

#BS The following line is needed because by default the G4TMP directory is set to $G4WORDIR/tmp, which is will cause the program to be created in the origional geant directory. open up new bash session to make sure this doesnt cause issues after finishing if you want to work in the old directory
G4WORKDIR := /home/branden/testing/

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin


#BS, Without these next two lines I can't compoile, check on their purpose
#ROOTINC = -I/home/branden/builds/root/include
ROOTINC = -I/builds/root/include
CPPFLAGS += $(ROOTINC)

# Add ROOT options for compilation
CPPFLAGS += `root-config --cflags`
CPPFLAGS += -I$(HOME)/builds/CLHEP-install/include/CLHEP
LDFLAGS += `root-config --libs`
LDFLAGS += -L$(G4LIB)/$(G4SYSTEM)
LDFLAGS += -L$(HOME)/builds/CLHEP-install/lib

LDFLAGS += `root-config --glibs`

# Add FMSSRC headers
CPPFLAGS += -I$(FMSSRC)
EXTRALIBS := $(FMSSRC)/Fpdchan.so
G4DEBUG = 1

#BS. NOT SURE HOW TO HAVE THE MAKEFILE COMPILE THE CLASS WITHOUT THIS. NEED TO FIX LATER
CPPFLAGS += -I/home/branden/testing/fmsu/myRootClass
EXTRALIBS += /home/branden/testing/fmsu/myRootClass/libG4QT.so
#EXTRALIBS += /home/branden/builds/root/lib/libTree.so
EXTRALIBS += /builds/root/lib/libTree.so
include $(G4INSTALL)/config/binmake.gmk
