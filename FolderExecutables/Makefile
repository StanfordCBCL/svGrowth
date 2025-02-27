CC = g++ -std=c++11
CFLAGS =
LDFLAGS= 
LDLIBS = -lgsl -lgslcblas -lm -lboost_program_options -D_GLIBCXX_USE_CXX11_ABI=1
SOURCES= vessel.cpp matfun.cpp functions.cpp main.cpp layer.cpp constituent.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=gnr

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -c -o $@ $(LDLIBS)

clean:
	rm -f *.o *.mod *~ *_out_* $(EXECUTABLE)
