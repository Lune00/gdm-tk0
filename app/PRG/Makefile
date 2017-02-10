APP = PRG
INC_PATH  = ../../src/include
SYS_PATH  = ../../src/systems
SYS_A_PATH = ../../src/systems_A
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

PRG.o: ../../src/include/simulation.hpp /usr/include/assert.h
PRG.o: /usr/include/sys/cdefs.h /usr/include/sys/types.h
PRG.o: /usr/include/sys/appleapiopts.h /usr/include/machine/types.h
PRG.o: /usr/include/i386/types.h /usr/include/i386/_types.h
PRG.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
PRG.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
PRG.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
PRG.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/sys/_structs.h
PRG.o: /usr/include/signal.h /usr/include/_types.h /usr/include/sys/signal.h
PRG.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
PRG.o: /usr/include/i386/_structs.h ../../src/include/algo.hpp
PRG.o: ../../src/include/sample.hpp ../../src/include/body2d.hpp
PRG.o: ../../src/include/dof.hpp /usr/include/math.h
PRG.o: /usr/include/architecture/i386/math.h ../../src/include/talk.hpp
PRG.o: ../../src/include/inter2d.hpp ../../src/include/vertex.hpp
PRG.o: ../../src/include/disk.hpp ../../src/include/polyg.hpp
PRG.o: ../../src/include/rline.hpp ../../src/include/purge.h
PRG.o: ../../src/include/network.hpp ../../src/include/groupRelationData.hpp
PRG.o: ../../src/include/cohesionLaw.hpp ../../src/include/fsafe.hpp
PRG.o: ../../src/include/near.hpp ../../src/include/system.hpp
PRG.o: ../../src/include/control.hpp ../../src/include/groupData.hpp
PRG.o: ../../src/include/io.hpp ../../src/include/system_A.hpp
PRG.o: /usr/include/sys/stat.h /usr/include/Availability.h
PRG.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
PRG.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
PRG.o: ../../src/include/biaxial_A.hpp ../../src/include/biaxial.hpp
PRG.o: ../../src/include/circularProbe.hpp ../../src/include/probe.hpp
PRG.o: ../../src/include/stress.hpp ../../src/include/tensor.hpp
PRG.o: ../../src/include/heightProbe.hpp ../../src/include/dkdkP.hpp
PRG.o: ../../src/include/dkdk.hpp ../../src/include/interdof.hpp
PRG.o: ../../src/include/rectangularProbe.hpp
PRG.o: ../../src/include/solidfraction.hpp ../../src/include/anisotropy.hpp
PRG.o: ../../src/include/pointSet.hpp ../../src/include/point.hpp
PRG.o: ../../src/include/dataSet.hpp ../../src/include/speedProfile.hpp
PRG.o: ../../src/include/NRsource.hpp ../../src/include/cluster_A.hpp
