CC=mpixlc
CXX=mpixlcxx

CFLAGS=-O3 -g
DEBUG=-DDEBUG

SRCS = 	route.cxx \
		personality.cxx \
		neighbour.cxx \
		iotree.cxx

OBJS = 	$(SRCS:.cxx=.o)

TARGET = iotree 

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cxx
		$(CXX) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) $(DEBUG) $(OBJS) $(LIBS)

clean:
		$(RM) *.o *~

