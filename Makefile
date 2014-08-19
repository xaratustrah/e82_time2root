CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)
CINT = rootcint

OS_NAME:=$(shell uname -s | tr A-Z a-z)
ifeq ($(OS_NAME),darwin)
STDINCDIR := -I/opt/local/include
STDLIBDIR := -L/opt/local/lib
else
STDINCDIR := 
STDLIBDIR := 
endif

CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR)
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lgsl -lgslcblas -lfftw3 -lm

CPPFLAGS += -g

TARGET = time2root

SRC = FritzDPSS.cxx iqtdata.cxx multitaper.cxx header.cxx header_dict.cxx setinfo.c

OBJ = $(SRC:.cxx=.o)

all : $(TARGET) visualizer compare read_iqt

header_dict.cxx : header.h header_linkdef.h
	$(CINT) -f header_dict.cxx -c header.h header_linkdef.h

%.o : %.cxx
	$(CXX) $(CPPFLAGS) -o $@ -c $<

$(TARGET) : main.cxx $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) main.cxx $(OBJ) $(LDFLAGS)

compare: compare.cxx $(OBJ)
	$(CXX) $(CPPFLAGS) -o compare compare.cxx $(OBJ) $(LDFLAGS)

read_iqt: read_iqt.cxx $(OBJ)
	$(CXX) $(CPPFLAGS) -o read_iqt read_iqt.cxx $(OBJ) $(LDFLAGS)

visualizer : $(OBJ)
	$(CXX) $(CPPFLAGS) -o visualizer visualizer.C $(OBJ) $(LDFLAGS)

clean :
	rm -f *.o *_dict.* $(TARGET) visualizer compare read_iqt *~
