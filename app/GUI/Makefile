APP = GUI
CC = g++
INC_PATH = ../../src/include
SYS_PATH = ../../src/systems
FLTK_PATH = ../../../fltk-1.1.8rc1.1

CFLAGS = -I/sw/include -I$(FLTK_PATH) -I$(INC_PATH) -Wall -O3 

LIBS =  -L$(FLTK_PATH)/lib -lfltk -lfltk_gl \
        -framework AGL -framework OpenGL \
        -framework Carbon -framework ApplicationServices \
        -lm \
        -lsupc++


OBJS = $(APP).o display.o ../../lib/*.o 


all: $(APP)
	
	

display.o : display.cpp display.h
	$(CC) $(CFLAGS) -c display.cpp

$(APP).o : $(APP).cpp $(APP).h
	$(CC) $(CFLAGS) -c $(APP).cpp

$(APP).a : $(OBJS)
	ar -r $@ $(OBJS)

$(APP): $(APP).a
	$(CC) -o $@ $(APP).a $(LIBS)
	/Developer/Tools/Rez -t APPL -o $(APP) ./mac.r

#	/Developer/Tools/Rez -t APPL -o $(APP) /sw/include/FL/mac.r

clean:
	rm -f *~ \#*\#
	rm -f *.o $(APP) $(APP).a

depend:
	makedepend -- $(CFLAGS) -- *.cpp


# DO NOT DELETE

GUI.o: GUI.h glob.h ../../src/include/simulation.hpp /usr/include/assert.h
GUI.o: /usr/include/sys/cdefs.h /usr/include/sys/types.h
GUI.o: /usr/include/sys/appleapiopts.h /usr/include/machine/types.h
GUI.o: /usr/include/i386/types.h /usr/include/i386/_types.h
GUI.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
GUI.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
GUI.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
GUI.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/sys/_structs.h
GUI.o: /usr/include/signal.h /usr/include/_types.h /usr/include/sys/signal.h
GUI.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
GUI.o: /usr/include/i386/_structs.h ../../src/include/algo.hpp
GUI.o: ../../src/include/sample.hpp ../../src/include/body2d.hpp
GUI.o: ../../src/include/dof.hpp /usr/include/math.h
GUI.o: /usr/include/architecture/i386/math.h ../../src/include/talk.hpp
GUI.o: ../../src/include/inter2d.hpp ../../src/include/vertex.hpp
GUI.o: ../../src/include/disk.hpp ../../src/include/polyg.hpp
GUI.o: ../../src/include/rline.hpp ../../src/include/purge.h
GUI.o: ../../src/include/network.hpp ../../src/include/groupRelationData.hpp
GUI.o: ../../src/include/cohesionLaw.hpp ../../src/include/fsafe.hpp
GUI.o: ../../src/include/near.hpp ../../src/include/system.hpp
GUI.o: ../../src/include/control.hpp ../../src/include/groupData.hpp
GUI.o: ../../src/include/io.hpp ../../src/include/system_A.hpp
GUI.o: /usr/include/sys/stat.h /usr/include/Availability.h
GUI.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
GUI.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
GUI.o: ../../src/include/CDalgo.hpp display.h ../../src/include/shearP_CD.hpp
GUI.o: ../../src/include/shearP_CD_A.hpp ../../src/include/connectivity_A.hpp
GUI.o: ../../src/include/probe.hpp ../../src/include/heightProbe.hpp
GUI.o: ../../src/include/rectangularProbe.hpp
GUI.o: ../../src/include/solidfraction.hpp
GUI.o: ../../src/include/circularProbe.hpp ../../src/include/anisotropy.hpp
GUI.o: ../../src/include/tensor.hpp ../../src/include/dkdkP.hpp
GUI.o: ../../src/include/pointSet.hpp ../../src/include/point.hpp
GUI.o: ../../src/include/dataSet.hpp ../../src/include/stress.hpp
GUI.o: ../../src/include/dkdk.hpp ../../src/include/interdof.hpp
GUI.o: ../../src/include/speedProfile.hpp ../../src/include/NRsource.hpp
GUI.o: ../../src/include/biaxial.hpp ../../src/include/biaxial_A.hpp
GUI.o: ../../src/include/cluster_A.hpp ../../src/include/pgpg.hpp
GUI.o: ../../src/include/dkrl.hpp ../../src/include/pgrl.hpp
GUI.o: ../../src/include/dkpg.hpp
display.o: display.h ../../src/include/simulation.hpp /usr/include/assert.h
display.o: /usr/include/sys/cdefs.h /usr/include/sys/types.h
display.o: /usr/include/sys/appleapiopts.h /usr/include/machine/types.h
display.o: /usr/include/i386/types.h /usr/include/i386/_types.h
display.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
display.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
display.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
display.o: /usr/include/libkern/i386/_OSByteOrder.h
display.o: /usr/include/sys/_structs.h /usr/include/signal.h
display.o: /usr/include/_types.h /usr/include/sys/signal.h
display.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
display.o: /usr/include/i386/_structs.h ../../src/include/algo.hpp
display.o: ../../src/include/sample.hpp ../../src/include/body2d.hpp
display.o: ../../src/include/dof.hpp /usr/include/math.h
display.o: /usr/include/architecture/i386/math.h ../../src/include/talk.hpp
display.o: ../../src/include/inter2d.hpp ../../src/include/vertex.hpp
display.o: ../../src/include/disk.hpp ../../src/include/polyg.hpp
display.o: ../../src/include/rline.hpp ../../src/include/purge.h
display.o: ../../src/include/network.hpp
display.o: ../../src/include/groupRelationData.hpp
display.o: ../../src/include/cohesionLaw.hpp ../../src/include/fsafe.hpp
display.o: ../../src/include/near.hpp ../../src/include/system.hpp
display.o: ../../src/include/control.hpp ../../src/include/groupData.hpp
display.o: ../../src/include/io.hpp ../../src/include/system_A.hpp
display.o: /usr/include/sys/stat.h /usr/include/Availability.h
display.o: /usr/include/AvailabilityInternal.h /usr/include/stdio.h
display.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
display.o: ../../src/include/shearP_CD.hpp ../../src/include/shearP_CD_A.hpp
display.o: ../../src/include/connectivity_A.hpp ../../src/include/probe.hpp
display.o: ../../src/include/heightProbe.hpp
display.o: ../../src/include/rectangularProbe.hpp
display.o: ../../src/include/solidfraction.hpp
display.o: ../../src/include/circularProbe.hpp
display.o: ../../src/include/anisotropy.hpp ../../src/include/tensor.hpp
display.o: ../../src/include/dkdkP.hpp ../../src/include/pointSet.hpp
display.o: ../../src/include/point.hpp ../../src/include/dataSet.hpp
display.o: ../../src/include/stress.hpp ../../src/include/dkdk.hpp
display.o: ../../src/include/interdof.hpp ../../src/include/speedProfile.hpp
display.o: ../../src/include/NRsource.hpp ../../src/include/biaxial.hpp
display.o: ../../src/include/biaxial_A.hpp ../../src/include/cluster_A.hpp
display.o: ../../src/include/pgpg.hpp ../../src/include/dkrl.hpp
display.o: ../../src/include/pgrl.hpp ../../src/include/dkpg.hpp
