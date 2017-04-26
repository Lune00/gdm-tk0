APP = continu
INC_PATH  = ../../src/include
SYS_PATH  = ../../src/systems
SYS_A_PATH = ../../src/systems_A
LIB_PATH = ../../lib

include ../../Makefile.in
SRC = grid.cpp config.cpp point_bis.cpp


OBJS = $(subst .cpp,.o,$(SRC)) 
OBJS2 = $(LIB_PATH)/*.o

.SUFFIXES: .cpp
.cpp.o:
	$(CC) -c $(CFLAGS) $<
	cp -f $@ $(LIB_PATH)
#Il faut rajouter compulation de continu.cpp ici
all: $(OBJS)
	$(CC) -o $(APP) $(APP).cpp $(CFLAGS) $(OBJS2) 


clean:
	rm -f *~ \#*\#
	rm -f *.o

depend:
	makedepend -- $(CFLAGS) -- *.cpp

#
