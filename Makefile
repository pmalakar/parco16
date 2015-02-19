CC=mpixlc
CXX=mpixlcxx

CFLAGS=-O3 -g

SRCS = 	route.cxx \
		personality.cxx \
		neighbour.cxx \
		iotree.cxx

OBJS = 	$(SRCS:.cxx=.o)

TARGET = iotree2 

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cxx
		$(CXX) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
		$(RM) *.o *~

