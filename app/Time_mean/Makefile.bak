APP = Time_mean
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

Time_mean.o: ../../src/include/simulation.hpp /usr/include/assert.h
Time_mean.o: /usr/include/sys/cdefs.h /usr/include/sys/types.h
Time_mean.o: /usr/include/sys/appleapiopts.h /usr/include/machine/types.h
Time_mean.o: /usr/include/i386/types.h /usr/include/i386/_types.h
Time_mean.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
Time_mean.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Time_mean.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Time_mean.o: /usr/include/libkern/i386/_OSByteOrder.h
Time_mean.o: /usr/include/sys/_structs.h /usr/include/signal.h
Time_mean.o: /usr/include/_types.h /usr/include/sys/signal.h
Time_mean.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Time_mean.o: /usr/include/i386/_structs.h ../../src/include/algo.hpp
Time_mean.o: ../../src/include/sample.hpp ../../src/include/body2d.hpp
Time_mean.o: ../../src/include/dof.hpp /usr/include/math.h
Time_mean.o: /usr/include/architecture/i386/math.h ../../src/include/talk.hpp
Time_mean.o: ../../src/include/inter2d.hpp ../../src/include/vertex.hpp
Time_mean.o: ../../src/include/disk.hpp ../../src/include/polyg.hpp
Time_mean.o: ../../src/include/rline.hpp ../../src/include/purge.h
Time_mean.o: ../../src/include/network.hpp
Time_mean.o: ../../src/include/groupRelationData.hpp
Time_mean.o: ../../src/include/cohesionLaw.hpp ../../src/include/fsafe.hpp
Time_mean.o: ../../src/include/near.hpp ../../src/include/system.hpp
Time_mean.o: ../../src/include/control.hpp ../../src/include/groupData.hpp
Time_mean.o: ../../src/include/io.hpp ../../src/include/system_A.hpp
Time_mean.o: /usr/include/sys/stat.h /usr/include/Availability.h
Time_mean.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
Time_mean.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
Time_mean.o: ../../src/include/dataSet.hpp ../../src/include/pointSet.hpp
Time_mean.o: ../../src/include/point.hpp
