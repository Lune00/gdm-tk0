INC_PATH  = ../include
LIB_PATH  = ../../lib

include ../../Makefile.in

SRC = sample.cpp network.cpp \
      simulation.cpp \
      control.cpp near.cpp \
      io.cpp talk.cpp \
      groupData.cpp groupRelationData.cpp\
	  dof.cpp

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

