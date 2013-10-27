# nmake makefile for the CGNS tools under Windows

DOMAKE = nmake /nologo /f Makefile.win
!include ..\make.defs

defaults : cgnsview cgnscalc cgnsplot utilities
install : install-config install-cgnsview install-cgnscalc \
	install-cgnsplot install-utilities

all : cgnsview cgnscalc cgnsplot utilities
install-all : install-config install-cgnsview install-cgnscalc \
	install-cgnsplot install-utilities

clean :
	cd cgnsview && $(DOMAKE) clean
	cd cgnscalc && $(DOMAKE) clean
	cd cgnsplot && $(DOMAKE) clean
	cd utilities && $(DOMAKE) clean
	cd calclib && $(DOMAKE) clean
	cd tkogl && $(DOMAKE) clean

distclean : clean
	-del make.win Makefile cgconfig.bat

cgnsview  : prog-cgnsview
cgnscalc  : prog-cgnscalc
cgnsplot  : prog-cgnsplot
utilities : prog-utilities

prog-cgnsview :
	cd cgnsview && $(DOMAKE)

prog-cgnscalc : lib-calclib
	cd cgnscalc && $(DOMAKE)

prog-cgnsplot : lib-tkogl
	cd cgnsplot && $(DOMAKE)

prog-utilities : lib-calclib
	cd utilities && $(DOMAKE)

lib-calclib :
	cd calclib && $(DOMAKE)

lib-tkogl :
	cd tkogl && $(DOMAKE)

install-config : $(BINDIR)
	$(INSTALL) cgconfig.bat $(BINDIR)

$(BINDIR) : $(INSTALLDIR)
	-mkdir $(BINDIR)

$(INSTALLDIR) :
	-mkdir $(INSTALLDIR)

install-cgnsview :
	cd cgnsview && $(DOMAKE) install

install-cgnscalc : lib-calclib
	cd cgnscalc && $(DOMAKE) install

install-cgnsplot : lib-tkogl
	cd cgnsplot && $(DOMAKE) install

install-utilities : lib-calclib
	cd utilities && $(DOMAKE) install
