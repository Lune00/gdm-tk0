APP = spl2cin
CC = g++
INC_PATH = ../../../src/include
LIB_PATH = ../../../lib

CFLAGS = -I/sw/include -I$(INC_PATH) -Wall -O3 -Werror

OBJS = $(APP).o $(LIB_PATH)/*.o   

all: $(APP)


$(APP).o : $(APP).cpp 
	$(CC) $(CFLAGS) -c $(APP).cpp

$(APP) : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o $(APP) 

depend:
	makedepend -- $(CFLAGS) -- *.cpp


# DO NOT DELETE

spl2cin.o: ../../../src/include/io.hpp ../../../src/include/algo.hpp
spl2cin.o: ../../../src/include/sample.hpp ../../../src/include/body2d.hpp
spl2cin.o: ../../../src/include/dof.hpp /usr/include/math.h
spl2cin.o: /usr/include/architecture/i386/math.h /usr/include/sys/cdefs.h
spl2cin.o: ../../../src/include/talk.hpp ../../../src/include/inter2d.hpp
spl2cin.o: ../../../src/include/vertex.hpp ../../../src/include/disk.hpp
spl2cin.o: ../../../src/include/polyg.hpp ../../../src/include/rline.hpp
spl2cin.o: ../../../src/include/purge.h ../../../src/include/network.hpp
spl2cin.o: ../../../src/include/groupRelationData.hpp
spl2cin.o: ../../../src/include/cohesionLaw.hpp
spl2cin.o: ../../../src/include/fsafe.hpp ../../../src/include/near.hpp
spl2cin.o: ../../../src/include/system.hpp ../../../src/include/control.hpp
spl2cin.o: ../../../src/include/groupData.hpp
spl2cin.o: ../../../src/include/simulation.hpp /usr/include/assert.h
spl2cin.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
spl2cin.o: /usr/include/machine/types.h /usr/include/i386/types.h
spl2cin.o: /usr/include/i386/_types.h /usr/include/sys/_types.h
spl2cin.o: /usr/include/machine/_types.h /usr/include/machine/endian.h
spl2cin.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
spl2cin.o: /usr/include/libkern/_OSByteOrder.h
spl2cin.o: /usr/include/libkern/i386/_OSByteOrder.h
spl2cin.o: /usr/include/sys/_structs.h /usr/include/signal.h
spl2cin.o: /usr/include/_types.h /usr/include/sys/signal.h
spl2cin.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
spl2cin.o: /usr/include/i386/_structs.h ../../../src/include/system_A.hpp
spl2cin.o: /usr/include/sys/stat.h /usr/include/Availability.h
spl2cin.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
spl2cin.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
