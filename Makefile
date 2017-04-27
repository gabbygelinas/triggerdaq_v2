#
# Makefile for the MIDAS analyzer component of ROOTANA
#
# Options:
# make NO_ROOT=1 # build without ROOT support
#

CXXFLAGS += -g -O2 -Wall -Wuninitialized -I.

ifndef ROOTANASYS
norootanasys:
	@echo Error: ROOTANASYS in not defined, please source thisrootana.{sh,csh}
endif

# check for ROOT

ifdef NO_ROOT
ROOTSYS:=
endif

ifneq ($(ROOTSYS),)
HAVE_ROOT:=1
endif

# get the rootana Makefile settings

CXXFLAGS += -I$(ROOTANASYS)/include $(shell cat $(ROOTANASYS)/include/rootana_cflags.txt)
LIBS     += -L$(ROOTANASYS)/lib -lrootana $(shell cat $(ROOTANASYS)/include/rootana_libs.txt) -L$(ANALYSIS_TPC) -lAGTPC -lAGDAQ -lAGUTILS -lGeom -lRGL

# select the main program - local custom main()
# or standard main() from rootana

#MAIN := manalyzer_main.o
MAIN := $(ROOTANASYS)/obj/manalyzer_main.o

# uncomment and define analyzer modules here

MODULES += ncfm.o unpack_module.o a16module.o Alpha16.o feam_module.o TsSync.o Feam.o FeamEVB.o AgEvent.o AgEVB.o Unpack.o reco_module.o final_module.o

ALL     += agana.exe

#ALL     += ncfm.exe

# examples

EXAMPLE_ALL += manalyzer.exe
EXAMPLE_ALL += manalyzer_example_cxx.exe
EXAMPLE_ALL += manalyzer_example_flow.exe
ifdef HAVE_ROOT
EXAMPLE_ALL += manalyzer_example_root.exe
EXAMPLE_ALL += manalyzer_example_root_graphics.exe
endif

all:: $(EXAMPLE_ALL)
all:: $(MODULES)
all:: $(ALL)

%.exe: %.o manalyzer_main.o
	$(CXX) -o $@ manalyzer_main.o $< $(CXXFLAGS) $(LIBS) -lm -lz -lpthread

$(EXAMPLE_ALL): %.exe:
	$(CXX) -o $@ $(ROOTANASYS)/obj/manalyzer_main.o $< $(CXXFLAGS) $(LIBS) -lm -lz -lpthread

%.exe: $(MAIN) $(MODULES)
	$(CXX) -o $@ $(MAIN) $(MODULES) $(CXXFLAGS) $(LIBS) -lm -lz -lpthread -lMathMore -lMinuit -lPhysics

ncfm.exe: %.exe: %.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(LIBS) -lm -lz -lpthread

reco_module.o: reco_module.cxx
	$(CXX) -o $@ -I$(AGTPC_ANALYSIS) -I$(ANALYSIS_TPC)/include -I$(GARFIELDPP) $(CXXFLAGS) -c $<

# Signals.o: HAVE_ROOT = ""
# Signals.o: $(AGTPC_ANALYSIS)/Signals.cc
# 	echo HAVE_ROOT = $(HAVE_ROOT)
# 	$(CXX) -o $@ -I$(AGTPC_ANALYSIS) -I$(GARFIELDPP) $(CXXFLAGS) -c $<

# SpacePoints.o: HAVE_ROOT = ""
# SpacePoints.o: $(AGTPC_ANALYSIS)/SpacePoints.cc
# 	echo HAVE_ROOT = $(HAVE_ROOT)
# 	$(CXX) -o $@ -I$(AGTPC_ANALYSIS) -I$(ANALYSIS_TPC)/include -I$(GARFIELDPP) $(CXXFLAGS) -c $<

# TLookUpTable.o: $(ANALYSIS_TPC)/src/TLookUpTable.cc
# 	echo HAVE_ROOT = $(HAVE_ROOT)
# 	$(CXX) -o $@ -I$(AGTPC_ANALYSIS) -I$(ANALYSIS_TPC)/include -I$(GARFIELDPP) $(CXXFLAGS) -c $<

# TPCBase.o: $(GARFIELDPP)/TPCBase.cc
# 	$(CXX) -o $@ -I$(GARFIELDPP) $(CXXFLAGS) -c $<

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
