# IBM MAKEFILE for Sfit DLL's
# To make your own DLL using the IBM C set compiler. 
# Change the SOURCENAME=GAUSS5 in this file to the name of your source file.
# type: nmake /f ibmdll.mak on the command line.
# NOTE ON SOURCE FILE:  To write your source code, modify for example CUSTFUNC.C

SOURCENAME=GAUSS5

.SUFFIXES:
.SUFFIXES: .c 

$(SOURCENAME).DLL:  \
  $(SOURCENAME).OBJ \
  IBMDLL.MAK
   ICC.EXE @<<
 /Fe"$(SOURCENAME).DLL" custfunc.def 
$(SOURCENAME).OBJ
<<

{.}.c.obj:
   ICC.EXE /Ge- /G4 /C   .\$*.c

$(SOURCENAME).OBJ:  $(SOURCENAME).c IBMDLL.MAK
