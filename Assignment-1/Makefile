CC	= gcc
OPT	= -g -O3 -std=c99
CHIP	= -march=native -mtune=native
ISA	=
DEFS	= -DDATA_ANALYSIS
LIBS	= -lm
CFLAGS  = $(OPT) $(CHIP) $(ISA) $(DEFS) $(XOPT)
BLAS ?= openblas

ifeq ($(BLAS),atlas)
    BLASFLAGS = -L/usr/lib64/atlas -lsatlas
else ifeq ($(BLAS),openblas)
    BLASFLAGS = -lopenblas
else
    $(error Unknown BLAS=$(BLAS) (use BLAS=atlas or BLAS=openblas))
endif

SOURCES = datatools.c main.c matmul.c
OBJECTS = $(SOURCES:.c=.o)

all: main

main: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(BLASFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS)

clean-all:
	rm -f main $(OBJECTS)
# DO NOT DELETE
