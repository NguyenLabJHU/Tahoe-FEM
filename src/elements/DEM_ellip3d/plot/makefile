SOURCES = $(wildcard *.cpp)
EXECUTABLE = $(SOURCES:.cpp=)
CXX = g++
CXXFLAGS =

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
clean:
	-rm -f *~ $(EXECUTABLE)


