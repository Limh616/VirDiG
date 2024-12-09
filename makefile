CXX = g++
CXXFLAGS = -Wall -g
LDFLAGS = -lpthread
SRCS = assemble2.cpp load_reads2.cpp transcript.cpp utility.cpp

TARGET = virdig

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)
