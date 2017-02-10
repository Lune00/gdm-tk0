INC_PATH  = ../../include
LIB_PATH  = ../../../lib

CC      = g++ 
CFLAGS  = -O3 -Wall -Werror -I$(INC_PATH)

SRC = inter2d.cpp \
      dkdk.cpp dkdkP.cpp \
	  pgpg.cpp dkrl.cpp pgrl.cpp dkpg.cpp 
	  
	  #mpmp.cpp

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

