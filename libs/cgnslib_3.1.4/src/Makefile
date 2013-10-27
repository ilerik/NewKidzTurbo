# nmake makefile for the CGNS Library under Windows

!include make.defs

.SUFFIXES :
.SUFFIXES : .c .$(O) $(EXE)

OBJDIR  = lib
CGNSLIB = $(LIBCGNS)
INSTLIB = cgns.$(A)

COPTS   = $(CFLAGS) $(CFGFLAGS) /I.  

#----------

CGNSOBJS=\
	$(OBJDIR)\cgns_error.$(O) \
	$(OBJDIR)\cgns_internals.$(O) \
	$(OBJDIR)\cgns_io.$(O) \
	$(OBJDIR)\cgnslib.$(O)

# ADF/ADFH routines

ADFOBJS=\
	$(OBJDIR)\ADF_interface.$(O) \
	$(OBJDIR)\ADF_internals.$(O)

F2COBJS= $(OBJDIR)\cg_ftoc.$(O) $(OBJDIR)\cgio_ftoc.$(O)

#----------

all     : cgns tools 
cgns    : $(CGNSLIB)
test    : runtests
install : install-cgns

#----------

$(CGNSLIB) : $(OBJDIR) $(CGNSOBJS) $(ADFOBJS) $(F2COBJS)
	-@$(RM) $@
	$(AR) /lib /nologo $(AROUT)$@ $(CGNSOBJS) $(ADFOBJS) $(F2COBJS)

#----------

$(OBJDIR) :
	-$(MKDIR) $(OBJDIR)

#----------

tools : cgns
	-cd tools && nmake /nologo

tests : cgns
	-cd tests && nmake /nologo

cgnstools : cgns
	-cd cgnstools && nmake /nologo

runtests : cgns
	-cd tests && nmake /nologo test

#----------

clean :
	-cd $(OBJDIR) && $(RM) *.$(O)
	-cd tools && nmake /nologo clean
	-cd tests && nmake /nologo clean
	-cd cgnstools && nmake /nologo clean

allclean : distclean

distclean : clean
	-cd tools && nmake /nologo allclean
	-cd tests && nmake /nologo allclean
	-cd tools && $(RM) Makefile
	-cd tests && $(RM) Makefile
	-cd cgnstools && nmake /nologo distclean
	-$(RM) $(CGNSLIB)
	-$(RMDIR) $(OBJDIR)
	-$(RM) cgnstypes.h cgnstypes_f.h cgnslib_f.h
	-$(RM) *.pdb
	-$(RM) make.defs cgnsBuild.defs Makefile

install-cgns : cgns $(INCLUDEDIR) $(LIBDIR) 
	$(INSTALL_DATA) cgnstypes.h $(INCLUDEDIR)\cgnstypes.h
	$(INSTALL_DATA) cgnstypes_f.h $(INCLUDEDIR)\cgnstypes_f.h
	$(INSTALL_DATA) cgnslib.h $(INCLUDEDIR)\cgnslib.h
	$(INSTALL_DATA) cgnslib_f.h $(INCLUDEDIR)\cgnslib_f.h
	$(INSTALL_DATA) cgnswin_f.h $(INCLUDEDIR)\cgnswin_f.h
	$(INSTALL_DATA) cgns_io.h $(INCLUDEDIR)\cgns_io.h
	$(INSTALL_DATA) cgnsBuild.defs $(INCLUDEDIR)\cgnsBuild.defs
	$(INSTALL_DATA) $(CGNSLIB) $(LIBDIR)\$(INSTLIB)

install-all : install

install-tools :
	-cd tools && nmake /nologo install

install-cgnstools :
	-cd cgnstools && nmake /nologo install

$(INCLUDEDIR) : $(INSTALLDIR)
	-$(MKDIR) $(INCLUDEDIR)

$(LIBDIR) : $(INSTALLDIR)
	-$(MKDIR) $(LIBDIR)

$(INSTALLDIR) :
	-$(MKDIR) $(INSTALLDIR)

install-adf : $(INCLUDEDIR)\adf
	$(INSTALL_DATA) adf\ADF.h $(INCLUDEDIR)\adf\ADF.h

$(INCLUDEDIR)\adf : $(INCLUDEDIR)
	-$(MKDIR) $(INCLUDEDIR)\adf

install-adfh : $(INCLUDEDIR)\adfh
	$(INSTALL_DATA) adfh\ADFH.h $(INCLUDEDIR)\adfh\ADFH.h

$(INCLUDEDIR)\adfh : $(INCLUDEDIR)
	-$(MKDIR) $(INCLUDEDIR)\adfh

#---------- mid-level library

$(OBJDIR)\cgns_error.$(O) : cgns_error.c cgnslib.h cgns_header.h cgns_io.h
	$(CC) $(COPTS) $(COOUT)$@ /c cgns_error.c

$(OBJDIR)\cgns_internals.$(O) : cgns_internals.c cgnslib.h cgns_header.h cgns_io.h
	$(CC) $(COPTS) $(COOUT)$@ /c cgns_internals.c

$(OBJDIR)\cgns_io.$(O) : cgns_io.c cgnslib.h cgns_io.h \
	adf\ADF.h
	$(CC) $(COPTS) $(COOUT)$@ /c cgns_io.c

$(OBJDIR)\cgnslib.$(O) : cgnslib.c cgnslib.h cgns_header.h cgns_io.h
	$(CC) $(COPTS) $(HDF5INC) $(MPIINC) $(COOUT)$@ /c cgnslib.c

$(OBJDIR)\cg_ftoc.$(O) : cg_ftoc.c fortran_macros.h cgnslib.h cgns_header.h cgns_io.h
	$(CC) $(COPTS) $(F2CFLAGS) $(COOUT)$@ /c cg_ftoc.c
$(OBJDIR)\cgio_ftoc.$(O) : cgio_ftoc.c fortran_macros.h cgns_io.h
	$(CC) $(COPTS) $(F2CFLAGS) $(COOUT)$@ /c cgio_ftoc.c

cgnslib.h : cgnstypes.h
cgns_header.h : cgnstypes.h
cgns_io.h : cgnstypes.h

#---------- ADF

$(OBJDIR)\ADF_interface.$(O) : adf\ADF_interface.c \
	adf\ADF.h adf\ADF_internals.h
	$(CC) $(COPTS) /Iadf $(COOUT)$@ /c adf\ADF_interface.c

$(OBJDIR)\ADF_internals.$(O) : adf\ADF_internals.c \
	adf\ADF.h adf\ADF_internals.h
	$(CC) $(COPTS) /Iadf $(COOUT)$@ /c adf\ADF_internals.c

adf\ADF.h : cgnstypes.h
adf\ADF_internals.h : cgnstypes.h

#---------- HDF5

$(OBJDIR)\ADFH.$(O) : adfh\ADFH.c adfh\ADFH.h
	$(CC) $(COPTS) /Iadfh $(HDF5INC) $(HDF5DLL) $(MPIINC) $(COOUT)$@ /c adfh\ADFH.c

adfh\ADFH.h : cgnstypes.h
