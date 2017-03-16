INC_PATH  = ../include
LIB_PATH  = ../../lib

include ../../Makefile.in

SRC = system.cpp \
      generalCD.cpp \
      biaxial.cpp shear.cpp\
	  shearP_CD.cpp biaxial_dof.cpp\
	  check_friction.cpp\
	  brazilian.cpp


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
