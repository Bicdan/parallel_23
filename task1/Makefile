CXX = g++
CXXFLAGS = -std=c++11 -Xpreprocessor -fopenmp
LDFLAGS = -Xpreprocessor -fopenmp -lomp -L /opt/homebrew/opt/libomp/lib

SRCS = main.cpp Handler.cpp Printer.cpp Solver.cpp
OBJS = $(SRCS:.cpp=.o)

TARGET = parallel_task

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)
