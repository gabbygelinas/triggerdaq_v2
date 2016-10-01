# Makefile

CXXFLAGS = -g -O2 -Wall -Wuninitialized

# optional ZLIB library

CXXFLAGS += -DHAVE_ZLIB

# ROOT libraries

ifdef ROOTSYS
ROOTGLIBS = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs) -lXMLParser -lXMLIO -lThread
HAVE_RHTTP = $(shell $(ROOTSYS)/bin/root-config --has-http)
ifeq ($(HAVE_RHTTP),yes)
ROOTGLIBS += -lRHTTP
endif
CXXFLAGS += -DHAVE_ROOT $(shell $(ROOTSYS)/bin/root-config --cflags)
endif

ifdef MIDASSYS
MIDASLIBS = $(MIDASSYS)/linux/lib/libmidas.a -lutil -lrt
CXXFLAGS += -DHAVE_MIDAS -DOS_LINUX -Dextname -I$(MIDASSYS)/include

UNAME=$(shell uname)
RTLIB=
ifeq ($(UNAME),Darwin)
#CXXFLAGS += -DOS_LINUX -DOS_DARWIN
MIDASLIBS = $(MIDASSYS)/darwin/lib/libmidas.a
#RPATH=
else
RTLIB=-lrt
endif

endif

ifdef ROOTSYS
CXXFLAGS += -DHAVE_LIBNETDIRECTORY -IlibNetDirectory
endif


# Check for ROOTANASYS env variable.
ifdef ROOTANASYS
ROOTANAINC = -I$(ROOTANASYS)/include
ROOTANALIBS = $(ROOTANASYS)/lib/librootana.a
else
ROOTANAINC = -I../include
ROOTANALIBS = ../lib/librootana.a
endif

OBJS:=
#OBJS += TV792Histogram.o TV1190Histogram.o
#OBJS += TL2249Histogram.o TAgilentHistogram.o
#OBJS += TV1720Waveform.o TDT724Waveform.o
#OBJS += TV1730DppWaveform.o TV1730RawWaveform.o
#OBJS += TV1720Correlations.o
#OBJS += TAnaManager.o

all:: $(OBJS)
#all:: ana.exe
all:: analyzer.exe
all:: anaDisplay.exe
#all:: midas2root.exe
all:: agana.exe

Alpha16.o: Alpha16.h Alpha16.cxx
Unpack.o: Unpack.h Alpha16.h
anaDisplay.o: Alpha16.o

anaDisplay.o: anaCommon.cxx
a16module.o: anaCommon.cxx

OBJS += Alpha16.o
OBJS += Unpack.o

ana.exe: ana.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl $(RTLIB) -lutil

analyzer.exe: analyzer.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl $(RTLIB) -lutil

anaDisplay.exe: anaDisplay.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl $(RTLIB) -lutil

agana.exe: $(OBJS) a16module.o
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread -lssl -lutil

midas2root.exe: midas2root.cxx $(OBJS) 
	$(CXX) -o $@ $(CXXFLAGS) $(ROOTANAINC) $^ $(ROOTANALIBS) $(MIDASLIBS) $(ROOTGLIBS) -lm -lz -lpthread $(RPATH) -lssl $(RTLIB) -lutil

%Dict.o: %Dict.cxx
	$(CXX) $(CXXFLAGS) $(ROOTANAINC) -c $<

%Dict.cxx: %.h %_LinkDef.h
	LD_LIBRARY_PATH=$(ROOTSYS)/lib $(ROOTSYS)/bin/rootcint -f $@ -c -p $(CXXFLAGS) $(ROOTANAINC) $^


%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(ROOTANAINC) -c $<

dox:
	doxygen

clean::
	-rm -f *.o *.a *.exe 

# end

