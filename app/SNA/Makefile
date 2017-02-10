APP = SNA
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
