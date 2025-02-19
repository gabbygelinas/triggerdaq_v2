#
# Makefile for the DarkLight online analyzer
#

CXXFLAGS += -g -O2 -Wall -Wuninitialized -I. -std=c++11 $(USERFLAGS)

# build with MIDAS

ifdef MIDASSYS
CXXFLAGS += -I$(MIDASSYS)/include -I$(MIDASSYS)/manalyzer -I$(MIDASSYS)/midasio -I$(MIDASSYS)/mvodb -I$(MIDASSYS)/mjson
LIBS += -L$(MIDASSYS)/lib -lmanalyzer_main -lmanalyzer -lmidas -lrt -lutil
else ifdef ROOTANASYS
CXXFLAGS += -I$(ROOTANASYS)/include
LIBS += -L$(ROOTANASYS)/lib -lmanalyzer_main -lmanalyzer
else
norootanasys:
	@echo Error: MIDASSYS or ROOTANASYS should be defined
endif

# add ROOT

ifdef ROOTSYS
CXXFLAGS += -DHAVE_ROOT $(shell root-config --cflags)
RLIBS    += -L$(ROOTSYS)/lib -lCore -lHist -lRIO -lGraf -lGui -lGpad -lRHTTP -lMathCore -lImt -lMatrix -lThread -lMultiProc -lNet
endif

MODULES += unpack_cb_module.o ncfm.o unpack_cb.o cbko_module.o dltdc_module.o dltdc.o dltdc4_module.o dltdc8_module.o DlTdcEvent.o coinc_module.o

ALL     += dlana.exe
#ALL     += ncfm.exe

all:: $(MODULES)
all:: $(ALL)

CMAKE := cmake3

cmake: $(GIT_SUBMODULES)
	mkdir -p build
	cd build; $(CMAKE) ..; $(MAKE); $(MAKE) install

cclean:
	-rm -rf build bin lib include


%.exe: $(MODULES)
	$(CXX) -o $@ $(MODULES) $(CXXFLAGS) $(LIBS) $(RLIBS) -lm -lz -lpthread -Wl,-rpath,$(ROOTSYS)/lib

ncfm.exe: %.exe: %.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(LIBS) -lm -lz -lpthread

%.o: %.cxx
	$(CXX) -o $@ $(CXXFLAGS) -c $<

dltdc.o: dltdc.cxx dltdc.h
DlTdcEvent.o: dltdc.h DlTdcEvent.h
dltdc_module.o: dltdc.h DlTdcEvent.h
dltdc4_module.o: dltdc.h DlTdcEvent.h
dltdc8_module.o: dltdc.h DlTdcEvent.h
coinc_module.o: coinc_module.cxx

#%.o: ../chronobox_software/%.cxx
#	$(CXX) -o $@ $(CXXFLAGS) -c $<

html/index.html:
	-mkdir html
	-make -k dox
	touch html/index.html

dox:
	doxygen

clean::
	-rm -f *.o *.a *.exe

clean::
	-rm -f $(ALL)

clean::
	-rm -f $(EXAMPLE_ALL)

clean::
	-rm -rf *.exe.dSYM

clean::
	-rm -rf html

clean:: cclean

# end
