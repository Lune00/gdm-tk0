APP = gdm_run
CC = g++ -arch x86_64
INC_PATH = ../../src/include
LIB_PATH = ../../lib

CFLAGS = -I/sw/include -I$(INC_PATH) -Wall -O3 -Werror -g

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
