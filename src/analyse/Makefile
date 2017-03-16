DATA_PACK = analyze
INC_PATH  = ../include
LIB_PATH  = ../../lib

include ../../Makefile.in

SRC = tensor.cpp stress.cpp\
      pointSet.cpp dataSet.cpp\
      anisotropy.cpp NRsource.cpp solidfraction.cpp\
      connectivity_A.cpp speedProfile.cpp \
	  grading.cpp cluster_A.cpp interdof.cpp

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
