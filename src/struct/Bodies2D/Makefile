INC_PATH  = ../../include
LIB_PATH  = ../../../lib

CC      = g++ -arch x86_64
CFLAGS  = -O3 -Wall -Werror -I$(INC_PATH)

SRC = body2d.cpp \
      disk.cpp polyg.cpp rline.cpp massPoint.cpp 

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

