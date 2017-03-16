INC_PATH  = ../../include
LIB_PATH  = ../../../lib

include ../../../Makefile.in

SRC = cohesionLaw.cpp weightedAdhesion.cpp nocohesion.cpp simpleAdhesion.cpp \

OBJS = $(subst .cpp,.o,$(SRC))

.SUFFIXES: .cpp
.cpp.o:
	$(CC) -c $(CFLAGS) $<
	cp -f $@ $(LIB_PATH)

all: $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o

depend:
	makedepend -- $(CFLAGS) -- *.cpp

# DO NOT DELETE
