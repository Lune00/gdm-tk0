APP = dof2clust
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

dof2clust.o: ../../../src/include/io.hpp ../../../src/include/algo.hpp
dof2clust.o: ../../../src/include/sample.hpp ../../../src/include/body2d.hpp
dof2clust.o: ../../../src/include/dof.hpp /usr/include/math.h
dof2clust.o: /usr/include/architecture/i386/math.h /usr/include/sys/cdefs.h
dof2clust.o: ../../../src/include/talk.hpp ../../../src/include/inter2d.hpp
dof2clust.o: ../../../src/include/vertex.hpp ../../../src/include/disk.hpp
dof2clust.o: ../../../src/include/polyg.hpp ../../../src/include/rline.hpp
dof2clust.o: ../../../src/include/purge.h ../../../src/include/network.hpp
dof2clust.o: ../../../src/include/groupRelationData.hpp
dof2clust.o: ../../../src/include/cohesionLaw.hpp
dof2clust.o: ../../../src/include/fsafe.hpp ../../../src/include/near.hpp
dof2clust.o: ../../../src/include/system.hpp ../../../src/include/control.hpp
dof2clust.o: ../../../src/include/groupData.hpp
dof2clust.o: ../../../src/include/simulation.hpp /usr/include/assert.h
dof2clust.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
dof2clust.o: /usr/include/machine/types.h /usr/include/i386/types.h
dof2clust.o: /usr/include/i386/_types.h /usr/include/sys/_types.h
dof2clust.o: /usr/include/machine/_types.h /usr/include/machine/endian.h
dof2clust.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
dof2clust.o: /usr/include/libkern/_OSByteOrder.h
dof2clust.o: /usr/include/libkern/i386/_OSByteOrder.h
dof2clust.o: /usr/include/sys/_structs.h /usr/include/signal.h
dof2clust.o: /usr/include/_types.h /usr/include/sys/signal.h
dof2clust.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
dof2clust.o: /usr/include/i386/_structs.h ../../../src/include/system_A.hpp
dof2clust.o: /usr/include/sys/stat.h /usr/include/Availability.h
dof2clust.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
dof2clust.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
