# nmake makefile for Windows

CGNSDIR = ..
!include $(CGNSDIR)\make.defs

CGNSLIB = $(CGNSDIR)\$(LIBCGNS)

COPTS  = $(CFLAGS) /I$(CGNSDIR)
LDLIBS = $(CGNSLIB) $(BUILDLIBS)

#----------

ALL =	cgnslist$(EXE) \
	cgnscheck$(EXE) \
	cgnsversion$(EXE) \
	cgnsconvert$(EXE) \
	cgnscompress$(EXE) \
	cgnsdiff$(EXE) \
	cgnsreplace$(EXE) \
	cgnsnames$(EXE)

all : $(ALL)

#----------

cgnslist$(EXE) : cgnslist.$(O) getargs.$(O) $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ cgnslist.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)
cgnslist.$(O) : cgnslist.c getargs.h
	$(CC) $(COPTS) /c cgnslist.c

#----------

cgnscheck$(EXE) : cgnscheck.$(O) getargs.$(O) hash.$(O) cgnames.$(O) $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ cgnscheck.$(O) getargs.$(O) hash.$(O) cgnames.$(O) $(LDLIBS) $(CLIBS)
cgnscheck.$(O) : cgnscheck.c getargs.h hash.h cgnames.h
	$(CC) $(COPTS) /c cgnscheck.c

#----------

cgnsversion$(EXE) : cgnsversion.$(O) getargs.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsversion.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)
cgnsversion.$(O) : cgnsversion.c getargs.h
	$(CC) $(COPTS) /c cgnsversion.c

#----------

cgnsconvert$(EXE) : cgnsconvert.$(O) getargs.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsconvert.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)
cgnsconvert.$(O) : cgnsconvert.c getargs.h
	$(CC) $(COPTS) /c cgnsconvert.c

#----------

cgnscompress$(EXE) : cgnscompress.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ cgnscompress.$(O) $(LDLIBS) $(CLIBS)
cgnscompress.$(O) : cgnscompress.c
	$(CC) $(COPTS) /c cgnscompress.c

#----------

cgnsdiff$(EXE) : cgnsdiff.$(O) getargs.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsdiff.$(O) getargs.$(O) $(LDLIBS) $(CLIBS)
cgnsdiff.$(O) : cgnsdiff.c getargs.h
	$(CC) $(COPTS) /c cgnsdiff.c

#----------

cgnsreplace$(EXE) : cgnsreplace.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsreplace.$(O) $(LDLIBS) $(CLIBS)
cgnsreplace.$(O) : cgnsreplace.c
	$(CC) $(COPTS) /c cgnsreplace.c

#----------

cgnsnames$(EXE) : cgnsnames.$(O) cgnames.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ cgnsnames.$(O) cgnames.$(O) $(LDLIBS) $(CLIBS)
cgnsnames.$(O) : cgnsnames.c cgnames.h
	$(CC) $(COPTS) /c cgnsnames.c

#----------

getargs.$(O) : getargs.c getargs.h
	$(CC) $(COPTS) /c getargs.c

hash.$(O) : hash.c hash.h
	$(CC) $(COPTS) /c hash.c

cgnames.$(O) : cgnames.c cgnames.h
	$(CC) $(COPTS) /c cgnames.c

install : all $(BINDIR)
	$(INSTALL_PROG) cgnslist$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnscheck$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnsversion$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnsconvert$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnscompress$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnsdiff$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnsreplace$(EXE) $(BINDIR)
	$(INSTALL_PROG) cgnsnames$(EXE) $(BINDIR)
	$(INSTALL_PROG) adf2hdf.bat $(BINDIR)
	$(INSTALL_PROG) hdf2adf.bat $(BINDIR)
	$(INSTALL_PROG) cgnsupdate.bat $(BINDIR)

$(BINDIR) : $(INSTALLDIR)
	-$(MKDIR) $(BINDIR)

$(INSTALLDIR) :
	-$(MKDIR) $(INSTALLDIR)

clean :
	-$(RM) *.$(O)

allclean : clean
	-$(RM) *.exe
	-$(RM) *.pdb *.ilk
