TARGET= mandelbrot_gpu
OBJS	= main.o mandel.o writepng.o 
#mandel_offload.o

OPT	= -g -fast -Msafeptr -Minfo -acc -mp=gpu -gpu=mem:separate:pinnedalloc -gpu=cc90 -gpu=lineinfo -cuda -mp=noautopar
ISA	= 
PARA	= -fopenmp

PNGWRITERPATH = pngwriter
ARCH	      = $(shell uname -p)
PNGWRTLPATH   = $(PNGWRITERPATH)/lib/$(ARCH)
PNGWRTIPATH   = $(PNGWRITERPATH)/include
PNGWRITERLIB  = $(PNGWRTLPATH)/libpngwriter.a

CC	= nvc

CXX	= nvc++
CFLAGS= -I $(PNGWRTIPATH) $(OPT) $(ISA) $(PARA) $(XOPT)
CXXFLAGS= -I $(PNGWRTIPATH) $(OPT) $(ISA) $(PARA) $(XOPT)

F90C  	= gfortran
LIBS	= -L $(PNGWRTLPATH) -lpngwriter -lpng 

all: $(PNGWRITERLIB) $(TARGET)

$(TARGET): $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

$(PNGWRITERLIB):
	@cd pngwriter/src && $(MAKE) -f $(MAKEFILE_LIST)

clean:
	@/bin/rm -f *.o $(TARGET)

realclean: clean
	@cd pngwriter/src && $(MAKE) -f $(MAKEFILE_LIST) clean
	@rm -f $(PNGWRITERLIB)
	@rm -f $(TARGET)
	@rm -f mandelbrot.png

# dependencies
#
main.o  : main.cpp mandel.h 
#mandel_offload.h
#mandel.o: mandel.cpp
mandel.o: mandel.cpp mandel.h 
	$(CXX) $(CXXFLAGS) -c mandel.cpp
#mandel_offload.o: mandel_offload.cpp
writepng.o: writepng.cc writepng.h
	$(CXX) $(CXXFLAGS) -c writepng.cc
