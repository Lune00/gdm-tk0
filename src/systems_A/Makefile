INC_PATH  = ../include
LIB_PATH  = ../../lib
AN_PATH   = ../analyse

CC      = g++
CFLAGS  = -O3 -Wall -Wno-write-strings -I$(INC_PATH)

SRC = system_A.cpp shearP_CD_A.cpp biaxial_A.cpp brazilian_A.cpp

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
