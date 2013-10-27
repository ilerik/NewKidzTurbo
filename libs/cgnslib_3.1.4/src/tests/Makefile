# nmake makefile for Windows

CGNSDIR = ..
!include $(CGNSDIR)\make.defs

CGNSLIB = $(CGNSDIR)\$(LIBCGNS)

COPTS  = $(CFLAGS) /I$(CGNSDIR)
FOPTS  = $(FFLAGS) /I$(CGNSDIR)
LDLIBS = $(CGNSLIB) $(BUILDLIBS)

#----------

CALL = \
	elemtest$(EXE) \
	test_exts$(EXE) \
	test_partial$(EXE) \
	test_goto$(EXE) \
	test_ver31$(EXE) \
	write_array$(EXE) \
	write_links$(EXE) \
	write_bcdata$(EXE) \
	write_test$(EXE) \
	write_zones$(EXE) \
	write_rind$(EXE)
FALL =	cgwrite$(EXE) \
	cgread$(EXE) \
	cgzconn$(EXE) \
	cgsubreg$(EXE)
CALL64 = test64c$(EXE)
FALL64 = test64f$(EXE)

TESTS = $(CALL) $(FALL)
ALL64 = $(CALL64) $(FALL64)
ALL   = dbtest$(EXE) open_cgns$(EXE) $(TESTS)

#----------

all : $(ALL)
fortran : $(FALL)
test64 : $(ALL64)

#----------

test : $(TESTS)
	@echo === running tests ===
	-@runtest.bat elemtest$(EXE)
	-@runtest.bat test_exts$(EXE)
	-@runtest.bat test_partial$(EXE)
	-@runtest.bat test_goto$(EXE)
	-@runtest.bat test_ver31$(EXE)
	-@runtest.bat write_array$(EXE)
	-@runtest.bat write_links$(EXE)
	-@runtest.bat write_bcdata$(EXE)
	-@runtest.bat write_test$(EXE)
	-@runtest.bat write_zones$(EXE)
	-@runtest.bat write_rind$(EXE)
	-@runtest.bat cgwrite$(EXE)
	-@runtest.bat cgread$(EXE)
	-@runtest.bat cgzconn$(EXE)
	-@runtest.bat cgsubreg$(EXE)
	@echo === finished ===

#----------

dbtest$(EXE) : dbtest.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ dbtest.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
dbtest.$(O) : dbtest.c utils.h
	$(CC) $(COPTS) /c dbtest.c

#----------

elemtest$(EXE) : elemtest.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ elemtest.c $(LDLIBS) $(CLIBS)

#----------

open_cgns$(EXE) : open_cgns.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ open_cgns.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
open_cgns.$(O) : open_cgns.c utils.h
	$(CC) $(COPTS) /c open_cgns.c

#----------

test_exts$(EXE) : test_exts.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ test_exts.c $(LDLIBS) $(CLIBS)

#----------

test_partial$(EXE) : test_partial.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ test_partial.c $(LDLIBS) $(CLIBS)

#----------

test_goto$(EXE) : test_goto.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ test_goto.c $(LDLIBS) $(CLIBS)

#----------

test_ver31$(EXE) : test_ver31.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ test_ver31.c $(LDLIBS) $(CLIBS)

#----------

write_array$(EXE) : write_array.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ write_array.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
write_array.$(O) : write_array.c utils.h
	$(CC) $(COPTS) /c write_array.c

#----------

write_links$(EXE) : write_links.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ write_links.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
write_links.$(O) : write_links.c utils.h
	$(CC) $(COPTS) /c write_links.c

#----------

write_bcdata$(EXE) : write_bcdata.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ write_bcdata.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
write_bcdata.$(O) : write_bcdata.c utils.h
	$(CC) $(COPTS) /c write_bcdata.c

#----------

write_test$(EXE) : write_test.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ write_test.c $(LDLIBS) $(CLIBS)

#----------

write_zones$(EXE) : write_zones.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ write_zones.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
write_zones.$(O) : write_zones.c utils.h
	$(CC) $(COPTS) /c write_zones.c

#----------

write_rind$(EXE) : write_rind.c $(CGNSLIB)
	$(CC) $(COPTS) $(CEOUT)$@ write_rind.c $(LDLIBS) $(CLIBS)

#----------

cgwrite$(EXE) : cgwrite.F $(CGNSLIB)
	$(F77) $(FOPTS) $(FEOUT)$@ cgwrite.F $(LDLIBS) $(FLIBS)

#----------

cgread$(EXE) : cgread.F $(CGNSLIB)
	$(F77) $(FOPTS) $(FEOUT)$@ cgread.F $(LDLIBS) $(FLIBS)

#----------

cgzconn$(EXE) : cgzconn.F $(CGNSLIB)
	$(F77) $(FOPTS) $(FEOUT)$@ cgzconn.F $(LDLIBS) $(FLIBS)

#----------

cgsubreg$(EXE) : cgsubreg.F $(CGNSLIB)
	$(F77) $(FOPTS) $(FEOUT)$@ cgsubreg.F $(LDLIBS) $(FLIBS)

#----------

test64c$(EXE) : test64c.$(O) utils.$(O) $(CGNSLIB)
	$(CC) $(CFLAGS) $(CEOUT)$@ test64c.$(O) utils.$(O) $(LDLIBS) $(CLIBS)
	$(STRIP) $@
test64c.$(O) : test64c.c utils.h
	$(CC) $(COPTS) /c test64c.c

#----------

test64f$(EXE) : test64f.F $(CGNSLIB)
	$(F77) $(FOPTS) $(FEOUT)$@ test64f.F $(LDLIBS) $(FLIBS)
	$(STRIP) $@

#----------

utils.$(O) : utils.c utils.h
	$(CC) $(COPTS) /c utils.c

clean :
	-$(RM) *.$(O)

allclean : clean
	-$(RM) *.exe
	-$(RM) *.pdb *.ilk
	-$(RM) *.cgns *.cgio
