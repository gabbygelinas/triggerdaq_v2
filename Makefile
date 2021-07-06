#
# Makefile for the ALPHA-g online analyzer
#

CXXFLAGS += -g -O2 -Wall -Wuninitialized -I. -std=c++11

# build with MIDAS

ifdef MIDASSYS
CXXFLAGS += -I$(MIDASSYS)/include -I$(MIDASSYS)/manalyzer -I$(MIDASSYS)/midasio -I$(MIDASSYS)/mvodb
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
CXXFLAGS += -DHAVE_ROOT -I$(ROOTSYS)/include
RLIBS    += -L$(ROOTSYS)/lib -lCore -lHist -lRIO -lGraf -lGui -lGpad -lRHTTP -lMathCore -lImt -lMatrix -lThread -ltbb -lMultiProc -lNet
endif

MODULES += ncfm.o unpack_module.o adc_module.o bsc_module.o pwb_module.o Alpha16.o feam_module.o TsSync.o Feam.o Tdc.o FeamEVB.o FeamAsm.o PwbAsm.o AgEvent.o AgEVB.o TrgAsm.o Unpack.o AgAsm.o wfsuppress.o wfsuppress2.o wfsuppress_pwb.o wfsuppress_adc.o wfexport_module.o pulser_module.o final_module.o coinc_module.o display_module.o

ALL     += agana.exe
#ALL     += ncfm.exe

all:: $(MODULES)
all:: $(ALL)

%.exe: $(MODULES)
	$(CXX) -o $@ $(MODULES) $(CXXFLAGS) $(LIBS) $(RLIBS) -lm -lz -lpthread -Wl,-rpath,$(ROOTSYS)/lib

ncfm.exe: %.exe: %.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(LIBS) -lm -lz -lpthread

%.o: %.cxx
	$(CXX) -o $@ $(CXXFLAGS) -c $<

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

# end
