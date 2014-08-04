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
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lgsl -lgslcblas -lm

CPPFLAGS += -g

TARGET = time2root

SRC = FritzDPSS.cxx iqtdata.cxx main.cxx fft.c multitaper.c setinfo.c

OBJ = $(SRC:.cpp=.o)

all : $(TARGET) vis

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CPPFLAGS) -o $@ -c $<

clean :
	rm -f *.o $(TARGET) *~
vis :
	$(CXX) $(CPPFLAGS) -o visualizer visualizer.C $(LDFLAGS)

