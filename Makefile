ifdef mpi3
CC = /home/pkcoff/mpich-opensrc/install-xl.legacy.ndebug/bin/mpicc
CXX = /home/pkcoff/mpich-opensrc/install-xl.legacy.ndebug/bin/mpicxx
LIBS = -L/home/pkcoff/mpich-opensrc/install-xl.legacy.ndebug/lib
INC = -I/home/pkcoff/mpich-opensrc/install-xl.legacy.ndebug/include -I/bgsys/drivers/ppcfloor -I/bgsys/drivers/ppcfloor/spi/include/kernel/cnk -I/bgsys/drivers/ppcfloor/hwi/include/bqc
DEFINES = -DMPI3
else
CC=mpixlc
CXX=mpixlcxx
endif

#common
CFLAGS=-O3 -g
DEBUG=-DCETUS #-DDEBUG
LIBMPI = #-L/soft/perftools/hpctw/lib -lmpihpm -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm -lrt -lstdc++ 
LIBUTILS =#-L/home/preeti/lib -lmem
LIBMPITRACE =-L/soft/perftools/hpctw/lib -lmpitrace
LIBS += $(LIBMPITRACE) $(LIBUTILS) $(LIBMPI) 
INC += #-I/home/preeti/inc

SRCS = 	route.cxx \
		personality.cxx \
		neighbour.cxx \
		iotree.cxx

OBJS = 	$(SRCS:.cxx=.o)

ifdef mpi3
TARGET = iotree.mpi3
else
TARGET = iotree 
endif

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cxx
		$(CXX) $(CFLAGS) -c $< -o $@ -DBGQ $(INC) $(DEBUG) $(LIBS) $(DEFINES)

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) -DBGQ $(INC) $(DEBUG) $(LIBS) $(OBJS) $(DEFINES)   

clean:
		$(RM) *.o *~

