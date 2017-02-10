APP = disk2polyg
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

disk2polyg.o: ../../../src/include/io.hpp ../../../src/include/algo.hpp
disk2polyg.o: ../../../src/include/sample.hpp ../../../src/include/body2d.hpp
disk2polyg.o: ../../../src/include/dof.hpp /usr/include/math.h
disk2polyg.o: /usr/include/architecture/i386/math.h /usr/include/sys/cdefs.h
disk2polyg.o: ../../../src/include/talk.hpp ../../../src/include/inter2d.hpp
disk2polyg.o: ../../../src/include/vertex.hpp ../../../src/include/disk.hpp
disk2polyg.o: ../../../src/include/polyg.hpp ../../../src/include/rline.hpp
disk2polyg.o: ../../../src/include/purge.h ../../../src/include/network.hpp
disk2polyg.o: ../../../src/include/groupRelationData.hpp
disk2polyg.o: ../../../src/include/cohesionLaw.hpp
disk2polyg.o: ../../../src/include/fsafe.hpp ../../../src/include/near.hpp
disk2polyg.o: ../../../src/include/system.hpp
disk2polyg.o: ../../../src/include/control.hpp
disk2polyg.o: ../../../src/include/groupData.hpp
disk2polyg.o: ../../../src/include/simulation.hpp /usr/include/assert.h
disk2polyg.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
disk2polyg.o: /usr/include/machine/types.h /usr/include/i386/types.h
disk2polyg.o: /usr/include/i386/_types.h /usr/include/sys/_types.h
disk2polyg.o: /usr/include/machine/_types.h /usr/include/machine/endian.h
disk2polyg.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
disk2polyg.o: /usr/include/libkern/_OSByteOrder.h
disk2polyg.o: /usr/include/libkern/i386/_OSByteOrder.h
disk2polyg.o: /usr/include/sys/_structs.h /usr/include/signal.h
disk2polyg.o: /usr/include/_types.h /usr/include/sys/signal.h
disk2polyg.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
disk2polyg.o: /usr/include/i386/_structs.h ../../../src/include/system_A.hpp
disk2polyg.o: /usr/include/sys/stat.h /usr/include/Availability.h
disk2polyg.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
disk2polyg.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
disk2polyg.o: ../../../src/include/alteration.hpp
disk2polyg.o: ../../../src/include/NRsource.hpp
