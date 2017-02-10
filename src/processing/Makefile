INC_PATH  = ../include
LIB_PATH  = ../../lib

CC      = g++
CFLAGS  = -O3 -Wall -Werror -I$(INC_PATH)

SRC = alteration.cpp grid.cpp

OBJS = $(subst .cpp,.o,$(SRC))

.SUFFIXES: .cpp
.cpp.o:
	$(CC) -c $(CFLAGS) $<
	cp -f $@ $(LIB_PATH)

all: $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o $(DATA_PACK).a

depend:
	makedepend -- $(CFLAGS) -- *.cpp
# DO NOT DELETE
