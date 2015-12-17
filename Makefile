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
CFLAGS=-O2 -g #-S #-O3

ifdef mpi3
DEFINES += -DBGQ -DMPI3 #-DDEBUG
else
DEFINES += -DBGQ #-DSTATS #-DDEBUG
endif

ifdef VESTA
DEFINES += -DVESTA
else
DEFINES += -DCETUS
endif

LIBHPM = -L/soft/perftools/hpctw/lib -lmpihpm 
LIBBGPM = -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm -lrt -lstdc++ 
LIBUTILS =-L/projects/Performance/preeti/utils -lbgqutils
LIBMPITRACE =-L/soft/perftools/hpctw/lib -lmpitrace
ifdef mpi3
LIBS += $(LIBUTILS) $(LIBBGPM) #$(LIBMPITRACE)
else
LIBS += $(LIBUTILS) $(LIBBGPM) #$(LIBMPITRACE)
endif
INC += -I/projects/Performance/preeti/utils	

SRCS = 	route.cxx \
		personality.cxx \
		neighbour.cxx \
		iotree.cxx

OBJS = 	$(SRCS:.cxx=.o)

ifdef mpi3
TARGET = iotree.mpi3.test
else
TARGET = iotree 
#TARGET = iotreeallnw #mu
#TARGET = iotreemod
endif

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cxx
		$(CXX) $(CFLAGS) -c $< -o $@ $(INC) $(LIBS) $(DEFINES)

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) $(OBJS) $(INC) $(LIBS) $(DEFINES)   

clean:
		$(RM) *.o *~

