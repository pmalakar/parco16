CC=mpixlc
CXX=mpixlcxx

CFLAGS=-O3 -g
DEBUG=-DCETUS #-DDEBUG
LIBMPI=#-L/soft/perftools/hpctw/lib -lmpihpm -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm -lrt -lstdc++ 
LIBUTILS=#-L/home/preeti/lib -lmem
LIBMPITRACE=-L/soft/perftools/hpctw/lib -lmpitrace
LIBS=$(LIBMPITRACE) $(LIBUTILS) $(LIBMPI) 
INC=#-I/home/preeti/inc

SRCS = 	route.cxx \
		personality.cxx \
		neighbour.cxx \
		iotree.cxx

OBJS = 	$(SRCS:.cxx=.o)

TARGET = iotree 

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cxx
		$(CXX) $(CFLAGS) -c $< -o $@ -DBGQ $(INC) $(DEBUG) $(LIBS)

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) -DBGQ $(INC) $(DEBUG) $(LIBS) $(OBJS) 

clean:
		$(RM) *.o *~

