APP = disk2cluster
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

disk2cluster.o: ../../../src/include/io.hpp ../../../src/include/algo.hpp
disk2cluster.o: ../../../src/include/sample.hpp
disk2cluster.o: ../../../src/include/body2d.hpp ../../../src/include/dof.hpp
disk2cluster.o: /usr/include/math.h /usr/include/architecture/i386/math.h
disk2cluster.o: /usr/include/sys/cdefs.h ../../../src/include/talk.hpp
disk2cluster.o: ../../../src/include/inter2d.hpp
disk2cluster.o: ../../../src/include/vertex.hpp ../../../src/include/disk.hpp
disk2cluster.o: ../../../src/include/polyg.hpp ../../../src/include/rline.hpp
disk2cluster.o: ../../../src/include/purge.h ../../../src/include/network.hpp
disk2cluster.o: ../../../src/include/groupRelationData.hpp
disk2cluster.o: ../../../src/include/cohesionLaw.hpp
disk2cluster.o: ../../../src/include/fsafe.hpp ../../../src/include/near.hpp
disk2cluster.o: ../../../src/include/system.hpp
disk2cluster.o: ../../../src/include/control.hpp
disk2cluster.o: ../../../src/include/groupData.hpp
disk2cluster.o: ../../../src/include/simulation.hpp /usr/include/assert.h
disk2cluster.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
disk2cluster.o: /usr/include/machine/types.h /usr/include/i386/types.h
disk2cluster.o: /usr/include/i386/_types.h /usr/include/sys/_types.h
disk2cluster.o: /usr/include/machine/_types.h /usr/include/machine/endian.h
disk2cluster.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
disk2cluster.o: /usr/include/libkern/_OSByteOrder.h
disk2cluster.o: /usr/include/libkern/i386/_OSByteOrder.h
disk2cluster.o: /usr/include/sys/_structs.h /usr/include/signal.h
disk2cluster.o: /usr/include/_types.h /usr/include/sys/signal.h
disk2cluster.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
disk2cluster.o: /usr/include/i386/_structs.h
disk2cluster.o: ../../../src/include/system_A.hpp /usr/include/sys/stat.h
disk2cluster.o: /usr/include/Availability.h
disk2cluster.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
disk2cluster.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
disk2cluster.o: ../../../src/include/alteration.hpp
disk2cluster.o: ../../../src/include/NRsource.hpp
