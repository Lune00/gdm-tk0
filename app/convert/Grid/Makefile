APP = Grid
CC = g++
INC_PATH = ../../../src/include
SYS_PATH = ../../../src/systems
LIB_PATH = ../../../lib

CFLAGS = -I/sw/include -I$(INC_PATH) -I$(SYS_PATH) -Wall -O3 -Werror

OBJS = $(APP).o  $(LIB_PATH)/*.o   

all: $(APP)


$(APP).o : $(APP).cpp 
	$(CC) $(CFLAGS) -c $(APP).cpp

$(APP) : $(OBJS)
	$(CC) -o $@ $(CFLAGS)  $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o $(APP) 

depend:
	makedepend -- $(CFLAGS) -- *.cpp


# DO NOT DELETE

Grid.o: ../../../src/include/io.hpp ../../../src/include/algo.hpp
Grid.o: ../../../src/include/sample.hpp ../../../src/include/body2d.hpp
Grid.o: ../../../src/include/dof.hpp /usr/include/math.h
Grid.o: /usr/include/architecture/i386/math.h /usr/include/sys/cdefs.h
Grid.o: ../../../src/include/talk.hpp ../../../src/include/inter2d.hpp
Grid.o: ../../../src/include/vertex.hpp ../../../src/include/disk.hpp
Grid.o: ../../../src/include/polyg.hpp ../../../src/include/rline.hpp
Grid.o: ../../../src/include/purge.h ../../../src/include/network.hpp
Grid.o: ../../../src/include/groupRelationData.hpp
Grid.o: ../../../src/include/fsafe.hpp ../../../src/include/near.hpp
Grid.o: ../../../src/include/system.hpp ../../../src/include/control.hpp
Grid.o: ../../../src/include/probe.hpp ../../../src/include/groupData.hpp
Grid.o: ../../../src/include/simulation.hpp /usr/include/assert.h
Grid.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
Grid.o: /usr/include/machine/types.h /usr/include/i386/types.h
Grid.o: /usr/include/i386/_types.h /usr/include/sys/_types.h
Grid.o: /usr/include/machine/_types.h /usr/include/machine/endian.h
Grid.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
Grid.o: /usr/include/stdint.h /usr/include/libkern/OSByteOrder.h
Grid.o: /usr/include/libkern/i386/OSByteOrder.h /usr/include/signal.h
Grid.o: /usr/include/_types.h /usr/include/sys/signal.h
Grid.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Grid.o: ../../../src/include/system_A.hpp /usr/include/sys/stat.h
Grid.o: /usr/include/stdio.h ../../../src/include/grid.hpp
Grid.o: ../../../src/include/NRsource.hpp
