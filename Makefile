CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

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

SRC = FritzDPSS.cxx iqtdata.cxx fft.c multitaper.cxx setinfo.c

OBJ = $(SRC:.cxx=.o)

all : $(TARGET) vis

$(TARGET) : main.cxx $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) main.cxx $(OBJ) $(LDFLAGS)

%.o : %.cxx
	$(CXX) $(CPPFLAGS) -o $@ -c $<

compare: compare.cxx $(OBJ)
	$(CXX) $(CPPFLAGS) -o compare compare.cxx $(OBJ) $(LDFLAGS)

clean :
	rm -f *.o $(TARGET) *~
vis :
	$(CXX) $(CPPFLAGS) -o visualizer visualizer.C $(LDFLAGS)

