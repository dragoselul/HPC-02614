CC   = gcc
OPT  = -g -O3 -std=c99 -funroll-all-loops
CHIP = -march=native -mtune=native
ISA  =
DEFS = -DDATA_ANALYSIS

PIC    = -fPIC
CFLAGS = $(OPT) $(PIC) $(CHIP) $(ISA) $(DEFS) $(XOPTS)

LIBS = -lm

BLASFLAGS = -L/usr/lib64/atlas -lsatlas
TARGET  = libmatmult.so
SOFLAGS = -shared

# IMPORTANT: library sources must NOT include main.c
LIBSRCS = matmult.c
LIBOBJS = $(LIBSRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(LIBOBJS)
	$(CC) -o $@ $(SOFLAGS) $(LIBOBJS) $(BLASFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(LIBOBJS)

clean-all: clean
	rm -f $(TARGET)
