APP = Fusion
INC_PATH  = ../../src/include
SYS_PATH  = ../../src/systems
SYS_A_PATH = ../../src/systems_analyze
LIB_PATH = ../../lib

CC      = g++
CFLAGS  = -I/sw/include -Wall -O3 -I$(SYS_A_PATH) -I$(INC_PATH) -I$(SYS_PATH) 


OBJS = $(APP).o  $(LIB_PATH)/*.o


all: $(APP)

$(APP).o : $(APP).cpp 
	$(CC) $(CFLAGS) -c $(APP).cpp

$(APP) : $(OBJS)
	$(CC) -o $@ $(CFLAGS)  $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o $(DATA_PACK).a

depend:
	makedepend -- $(CFLAGS) -- *.cpp
# DO NOT DELETE

Fusion.o: ../../src/include/simulation.hpp /usr/include/assert.h
Fusion.o: /usr/include/sys/cdefs.h /usr/include/sys/types.h
Fusion.o: /usr/include/sys/appleapiopts.h /usr/include/machine/types.h
Fusion.o: /usr/include/i386/types.h /usr/include/i386/_types.h
Fusion.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
Fusion.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Fusion.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Fusion.o: /usr/include/libkern/i386/_OSByteOrder.h
Fusion.o: /usr/include/sys/_structs.h /usr/include/signal.h
Fusion.o: /usr/include/_types.h /usr/include/sys/signal.h
Fusion.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Fusion.o: /usr/include/i386/_structs.h ../../src/include/algo.hpp
Fusion.o: ../../src/include/sample.hpp ../../src/include/body2d.hpp
Fusion.o: ../../src/include/dof.hpp /usr/include/math.h
Fusion.o: /usr/include/architecture/i386/math.h ../../src/include/talk.hpp
Fusion.o: ../../src/include/inter2d.hpp ../../src/include/vertex.hpp
Fusion.o: ../../src/include/disk.hpp ../../src/include/polyg.hpp
Fusion.o: ../../src/include/rline.hpp ../../src/include/purge.h
Fusion.o: ../../src/include/network.hpp
Fusion.o: ../../src/include/groupRelationData.hpp
Fusion.o: ../../src/include/cohesionLaw.hpp ../../src/include/fsafe.hpp
Fusion.o: ../../src/include/near.hpp ../../src/include/system.hpp
Fusion.o: ../../src/include/control.hpp ../../src/include/groupData.hpp
Fusion.o: ../../src/include/io.hpp ../../src/include/system_A.hpp
Fusion.o: /usr/include/sys/stat.h /usr/include/Availability.h
Fusion.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
Fusion.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
