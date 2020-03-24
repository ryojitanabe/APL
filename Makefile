TARGET = de
SRCDIR    = ./src
INCLUDE = $(wildcard $(SRCDIR)/*.h)
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)
CXX = g++
CXXFLAGS = -O3 -lm -c -std=c++11 -MMD -MP # -lboost_thread

all: $(TARGET)

-include $(DEPS)

$(TARGET): $(OBJS)
	$(CXX) -o  $(TARGET) $^

%.o: %.cc
	$(CXX) $(CXXFLAGS)  $<

clean:
	rm -rf $(OBJS) $(DEPS) $(TARGET)
