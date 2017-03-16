APP = gdm_run
INC_PATH = ../../src/include
LIB_PATH = ../../lib

include ../../Makefile.in

OBJS = $(APP).o  $(LIB_PATH)/*.o

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
