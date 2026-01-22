TARGET= libmatmult.so
OBJS	= matmult.o

CC	= nvc
CXX	= nvc++

OPT	= -g -fast -Msafeptr -Minfo -mp -mp=gpu -gpu=mem:separate:pinnedalloc -gpu=lineinfo -gpu=cc90 -cuda -cudalib=cublas
PIC   = -fpic -shared
ISA	=
PARA	=
INC   =
LIBS	= -lblas -lcblas

CXXFLAGS= $(OPT) $(PIC) $(INC) $(ISA) $(PARA) $(XOPT)

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	@/bin/rm -f $(TARGET) $(OBJS)
