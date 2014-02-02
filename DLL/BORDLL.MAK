# BORLAND MAKEFILE for Sfit DLL's
# To make your own DLL using the Borland OS/2 C++ compiler. 
# Change the SOURCENAME=GAUSS5 in this file to the name of your source file.
# type: make -fbordll.mak on the command line.
# NOTE ON SOURCE FILE:  To write your source code, modify for example CUSTFUNC.C

STARTUP = \bcos2\lib\c02.obj
STARTUPDLL = \bcos2\lib\c02d.obj
LIB = \bcos2\lib
INCLUDE = \bcos2\include

#===================================================================
#
# Compile switchs  that are enabled
# -c      compile don't link
# -v      debug info
# -sm     multi threaded

SOURCENAME=gauss5

CFLAGS  = -c -v -sm

LFLAGS  = -c -x -v  
LINKDLL = tlink $(LFLAGS) -L$(LIB) $(STARTUPDLL)

DLLLIBS = c2mt + os2

OBJECTS = $(SOURCENAME).obj

$(SOURCENAME).dll : $(OBJECTS)  custfunc.def
	$(LINKDLL) $(OBJECTS),$(SOURCENAME),,$(DLLLIBS),custfunc.def;

$(SOURCENAME).obj : $(SOURCENAME).c
   bcc $(CFLAGS) $(SOURCENAME).c
